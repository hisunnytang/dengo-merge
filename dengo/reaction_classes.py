import numpy as np
from .chemistry_constants import tevk, tiny, mh
import types
import os
import sympy
import h5py
import docutils.utils.roman as roman
from .periodic_table import \
    periodic_table_by_name, \
    periodic_table_by_number
from .mixin import ComparableMixin
import re

try:
    import ChiantiPy.core as ch
    chu = 1 # temporary to avoid the check
    #import ChiantiPy.util as chu
except (ImportError, KeyError):
    ch = None
    #chu = None

reaction_registry = {}
cooling_registry = {}
species_registry = {}

def registry_setup(func):
    """a decorator function register unseen species,
    cooling and reactions into the registry
    """
    def _wfunc(*args, **kwargs):
        old_names = [set(d.keys()) for d in (species_registry,
                                             cooling_registry,
                                             reaction_registry)]
        func(*args, **kwargs)
        nn = []
        for on, r in zip(old_names, (species_registry,
                                     cooling_registry,
                                     reaction_registry)):
            nn.append(set(r.keys()).difference(on))
        return nn
    return _wfunc


def ensure_reaction(r):
    if isinstance(r, Reaction): return r
    return reaction_registry[r]

def ensure_cooling(c):
    if isinstance(c, CoolingAction): return c
    return cooling_registry[c]

def ensure_species(s):
    if isinstance(s, Species): return s
    return species_registry[s]

count_m = sympy.Symbol('m', integer=True)
index_i = sympy.Idx('i', count_m)

class ReactionCoefficient(sympy.Symbol):
    #TODO: the functional derivatives
    # sympy.diff() of the sum/ mul of the combination of `ReacionCoefficient` class
    # gives zeros....
    energy = sympy.simplify("ge")

    def _eval_derivative(self, s):
        if s == self.energy:
            return sympy.Symbol("r%s" % self)
        else:
            return super(ReactionCoefficient, self)._eval_derivative(s)

    def _eval_derivative_n_times(self, s, n):
        if n == 1:
            return self._eval_derivative(s)
        else:
            return 0


    def diff(self, s):
        if s == self.energy:
            return sympy.Symbol("r%s" %self)
        else:
            return super(ReactionCoefficient, self).diff(s)

    @property
    def free_symbols(self):
        return super().free_symbols.union(set([self.energy]))


class Reaction(ComparableMixin):
    def __init__(self, name, coeff_fn, left_side, right_side):
        self.name = name
        self.coeff_fn = coeff_fn
        #self.coeff_sym = sympy.IndexedBase(name, (count_m,))
        self.coeff_sym = ReactionCoefficient("%s[i]" % name)
        self.left_side = left_side
        self.right_side = right_side
        self.considered = set( (s.name for n, s in left_side + right_side) )
        reaction_registry[name] = self # Register myself

    def __contains__(self, c):
        c = ensure_species(c)
        return c in self.considered

    @property
    def down_species(self):
        return [s for n, s in self.left_side]

    @property
    def up_species(self):
        return [s for n, s in self.right_side]

    @property
    def species(self):
        return set(self.down_species + self.up_species)

    def net_change(self, sname):
        up = sum( n for n, s in self.right_side if s.name == sname)
        down = sum( n for n, s in self.left_side if s.name == sname)
        if up == down:
            return up - down
        else:
            return up - down

    @property
    def nbody_reaction(self):
        # this is handy function for calculating
        # number of species involved in the reaction
        # this is used primarily in counting the
        # number of scale factors a^3 needed to
        # `calibrate`the reaction rates
        return sum(n for n, s in self.left_side)

    def __call__(self, quantities, up_derivatives, down_derivatives):
        # We just calculate our net derivatives and stick them in the right
        # place
        r = self.rate(quantities)
        for n, s in self.left_side:
            r *= s.number_density(quantities)**n
        for n, s in self.left_side:
            down_derivatives[s.name] += r * n * s.weight
        for n, s in self.right_side:
            up_derivatives[s.name] += r * n * s.weight
        return r

    def _cmpkey(self):
        return repr(self)

    def __repr__(self):
        a = "%s : " % self.name \
          + " + ".join( ["%s*%s" % (i, s.name) for i, s in self.left_side] ) \
          + " => " \
          + " + ".join( ["%s*%s" % (i, s.name) for i, s in self.right_side] )
        return a

    @classmethod
    def create_reaction(cls, name, left_side, right_side):
        """Initialize `Reaction`,
        the chemical reactions between species

        Parameters
        ----------
        name: str
            name of the reaction
        left_side:
            left hand side of the chemical reaction
        right_side:
            right hand side of the chemical reaction
        f:
            a function that takes state,
            and returns the reaction coefficient
        """
        def _w(f):
            rxn = cls(name, f, left_side, right_side)
        return _w

    def species_equation(self, species):
        if isinstance(species, str):
            species = species_registry[species]
        elif isinstance(species, (sympy.IndexedBase, sympy.Symbol)):
            species = species_registry[str(species)]
        if species not in self.species: return 0
        nr = self.net_change(species.name)
        return nr * self.lhs_equation
    @property
    def lhs_equation(self):
        #eq = self.coeff_sym[index_i]
        eq = self.coeff_sym
        for i, s in self.left_side:
            for ii in range(i):
                #eq *= s.symbol[index_i]
                eq *= s.symbol

        # a little hacky here to unfold
        # the algebraic power
        # i.e. H_1**3 -> H_1*H_1*H_1
        #eq = eq.replace(
        #        lambda x: x.is_Pow and x.exp > 0,
        #        lambda x: sympy.Symbol('*'.join([x.base.name]*x.exp)) )

        return eq

reaction = Reaction.create_reaction
def chianti_rate(atom_name, sm1, s, sp1):
    if ch is None: raise ImportError
    ion_name = s.name.lower()
    if "_" not in ion_name:
        print ("Name must be in ChiantiPy format.")
        raise RuntimeError
    de = species_registry['de']
    new_rates = []
    def ion_rate(network):
        ion = ch.ion(ion_name, temperature = network.T)
        try:
            ion.ionizRate()
        except AttributeError:
            print(f"{ion_name} is not defined in ChiantiPy master list")
            print(f"manually adding temperature to {ion_name}_ion")
            ion.Temperature = network.T
            ion.NTempDens   = len(network.T)
            ion.ionizRate()
        vals = ion.IonizRate['rate']
        return vals
    if sp1 is not None:
        Reaction("%s_i" % s.name, ion_rate,
                 [(1, s), (1, de)], # left side
                 [(1, sp1), (2, de)]) # right side
        new_rates.append("%s_i" % s.name)

    def rec_rate(network):
        ion = ch.ion(ion_name, temperature = network.T)
        # for some reason, the latest chiantipy
        # failed to update the tempeature for fully ionized
        ion.Temperature = network.T
        ion.recombRate()
        vals = ion.RecombRate['rate']
        return vals
    if sm1 is not None:
        Reaction("%s_r" % s.name, rec_rate,
                 [(1, s), (1, de)], # left side
                 [(1, sm1), ]) # right side
        new_rates.append("%s_r" % s.name)
    return new_rates

def ion_photoionization_rate(species, photo_background='HM12'):
    if chu is None: raise ImportError
    ion_name = species.name.lower()
    #ion_name = chu.zion2name(np.int(species.number),
    #                         np.int(species.free_electrons + 1))
    if "_" not in ion_name:
        print ("Name must be in 'Ion Species' format.")
        raise RuntimeError
    element_name = ion_name.split("_")[0]
    ion_state = int(ion_name.split("_")[1])
    species_pi = "%s_%s" % (element_name.capitalize(), ion_state + 1)
    de = species_registry['de']
    new_rates = []

    def photoionization_rate(network):
        # Read in photoheating rates generated from Ben Oppenheimer's data
        # (from: http://noneq.strw.leidenuniv.nl/)
        # and do linear interpolation and then recompute
        # the ends with either an extrapolation or falloff
        # NOTE: these rates do the interpolation as a function fo redshift
        fn = os.path.join(os.path.dirname(__file__),
                '..', 'input', 'photoionization',
                '%s_ion_by_ion_photoionization_%s.h5' %(element_name, photo_background))
        f = h5py.File(fn,'r')
        #f = h5py.File('../input/photoionization/%s_ion_by_ion_photoionization_%s.h5'
        #              %(element_name, photo_background))
        #print('../input/photoionization/%s_ion_by_ion_photoionization_%s.h5'
        #              %(element_name, photo_background))
        ### Intepolate values within table values ###
        vals = np.interp(network.z, f['z'], f['%s' %(ion_name)])

        end_method = 0 # 0 = extrapolation, 1 = gaussian falloff

        if end_method == 0:
            ### Extrapolation in logspace ###
            # convert to log space
            vals = np.log10(vals)
            logz = np.log10(network.z)
            logdataz = np.log10(f['z'])
            logdataS = np.log10(f['%s' %(ion_name)])

            # extrapolate
            extrapdown = logdataS[0] + \
                (logz - logdataz[0]) * (logdataS[0] - logdataS[1]) \
                / (logdataz[0] - logdataz[1])
            vals[logz < logdataz[0]] = extrapdown[logz < logdataz[0]]
            extrapup = logdataS[-1] + \
                (logz - logdataz[-1]) * (logdataS[-1] - logdataS[-2]) \
                / (logdataz[-1] - logdataz[-2])
            vals[logz > logdataz[-1]] = extrapup[logz > logdataz[-1]]

            # convert back to linear
            vals = 10.0**vals
        if end_method == 1:
            ### Gaussian falloff when values extend beyond table values ###
            # rename some variables to symplify code
            z = network.z
            dataz = f['z']
            dataS = f['%s' %(ion_name)]

            # compute gaussian tails
            gaussdown = dataS[0] * (tiny/dataS[0])**(((z - dataz[0])/(z[0] - dataz[0])))**2
            vals[z < dataz[0]] = gaussdown[z < dataz[0]]
            gaussup = dataS[-1] * (tiny/dataS[-1])**(((z - dataz[-1])/(z[-1] - dataz[-1])))**2
            vals[z > dataz[-1]] = gaussup[z > dataz[-1]]

        f.close()
        return vals

    if species_pi in species_registry:
        species_pi = species_registry[species_pi]
        Reaction("%s_pi" % species.name, photoionization_rate,
                 [(1, species), ], # left side
                 [(1, species_pi), (1, de)]) # right side
        new_rates.append("%s_pi" % species.name)
    return new_rates

class Species(ComparableMixin):
    def __init__(self, name, weight, pretty_name = None):
        self.pretty_name = pretty_name or name
        self.weight = weight
        self.name = name
        self.symbol = sympy.Symbol("%s" % name)
        species_registry[name] = self

    def _cmpkey(self):
        return repr(self)

    def __repr__(self):
        return "Species: %s" % (self.name)

class ChemicalSpecies(Species):
    """initialize chemical species object

    Parameters
    ----------
    name: str
        name of the chemical species
    weight: float
        the weight of the species in terms of atomic mass unit
    free_electrons: int
        the number of free electrons of the species
    """
    def __init__(self, name, weight, free_electrons = 0.0,
                 pretty_name = None):
        self.weight = weight
        self.free_electrons = free_electrons
        #self.elements = {}
        super(ChemicalSpecies, self).__init__(name, weight, pretty_name)

    def number_density(self, quantities):
        if self.weight == 0:
            return quantities[self.name]
        return quantities[self.name]/self.weight

    def add_to_dict(self, e, v):
        if  e in self._edict:
            self._edict[e] += v
        else:
            self._edict[e] = v

    def find_constituent(self):
        """find elements based on its name
        return: total weight
        """
        self._edict = {}
        name = self.original_name
        split = re.findall("\d+|\D+[a-z]|\D", name)
        for s in split:
            if s.isalpha():
                _ele = s
            if s.isnumeric():
                self.add_to_dict(_ele, int(s)-1)
            else:
                self.add_to_dict(_ele, 1)
        self.elements = self._edict
    def get_weight(self):
        w = 0
        for s, n in self.elements.items():
            w += periodic_table_by_name[s][1]*n
        self._weight = w

class AtomicSpecies(ChemicalSpecies):
    """Initialize the atomic species

    Parameters
    ----------
    atom_name: str
        name of the atom
    free_electrons: int
        number of free electron

    Note
    ----
    since it is an atom, the weight is taken directly from the periodic table.
    """
    def __init__(self, atom_name, free_electrons):
        num, weight, pn = periodic_table_by_name[atom_name]
        if free_electrons < 0:
            name = "%s_m%i" % (atom_name, np.abs(free_electrons + 1))
        else:
            name = "%s_%01i" % (atom_name, free_electrons + 1)
        pretty_name = "%s with %s free electrons" % (
            pn, free_electrons)
        # update the self.elements based on the input name
        self.original_name = atom_name
        self.find_constituent ()

        super(AtomicSpecies, self).__init__(name, weight,
            free_electrons, pretty_name)

class MolecularSpecies(ChemicalSpecies):
    """Initialize molecular species

    Parameters
    ----------
    molecule_name: str
        the name of the molecules
    weight: float
        weight of the molecule in amu
    free_electrons: int
        free electron in molecules

    """
    def __init__(self, molecule_name, weight, free_electrons,
                 original_name = None):
        name = "%s_%i" % (molecule_name, free_electrons + 1)
        pretty_name = "%s with %s free electrons" % (
            name, free_electrons)
        self.original_name = original_name or molecule_name
        # update the self.elements based on the input name
        self.find_constituent()
        super(MolecularSpecies, self).__init__(name, weight,
            free_electrons, pretty_name)


class CoolingAction(object):
    _eq = None
    def __init__(self, name, equation):
        self.name = name
        self._equation = equation
        self.tables = {}
        self.temporaries = {}
        self.table_symbols = {}
        self.temp_symbols = {}
        cooling_registry[name] = self # Register myself

    # so that we can sort the coolingaction
    def __gt__(self, other):
        return self.name > other.name

    @property
    def equation(self):
        if self._eq is not None: return self._eq
        symbols = dict((n, s.symbol) for n, s in species_registry.items())
        #ta_sym = dict((n, sympy.IndexedBase(n, (count_m,))) for n in self.tables))

        # instead of using sympy symbols, we declare

        ta_sym = dict((n, ReactionCoefficient("%s_%s[i]" % (self.name, n))) for n in self.tables)
        self.table_symbols.update(ta_sym)
        #tp_sym = dict((n, sympy.IndexedBase(n, (count_m,))) for n in self.temporaries))

        # temporaries are in fact function of of the tables...
        tp_sym = dict((n, sympy.Symbol("%s" % (n))) for n in self.temporaries)
        self.temp_symbols.update(tp_sym)

        for n, s in species_registry.items():
            try:
                name, Ilevel = n.split('_')
                sp_name = name + eval(Ilevel)*'I'
                temp_dict = {sp_name: s.symbol}
                symbols.update(temp_dict)
            except:
                pass

        symbols.update(self.table_symbols)
        symbols.update(self.temp_symbols)
        self._eq = eval(self._equation, symbols)
        for n, e in self.temporaries.items():
            e = eval(e, symbols)
            e = sympy.sympify(e)
            for n2, e2 in ta_sym.items():
                e = e.subs(n2, e2)
            self._eq = self._eq.subs(n, e)
        return self._eq

    @property
    def species(self):
        # self.equation
        bad = set(self.temp_symbols.values()).update( set(self.table_symbols.values()) )
        species = set([])
        #for s in self.equation.atoms(sympy.IndexedBase):
        for s in self.equation.atoms(sympy.Symbol):
            bad = set(self.temp_symbols.values()).update( set(self.table_symbols.values()) )
            if bad != None:
                if s not in bad:
                    species.add(species_registry[str(s).replace("[i]","")])
        return species

    def table(self, func):
        self.tables[func.__name__] = func

    def temporary(self, name, eq):
        self.temporaries[name] = eq
    @classmethod
    def create_cooling_action(cls, name, equation):
        obj = cls(name, equation)
        def _W(f):
            f(obj)
        return _W

cooling_action = CoolingAction.create_cooling_action

def ion_cooling_rate(species, atom_name):

    species_c = species.name
    ion_name = species.name.lower()
    de = species_registry['de']
    new_rates = []

    def cooling_rate(network):
        # Read in cooling rates from Gnat & Ferland 2012
        # and do linear interpolation and then recompute
        # the ends with either an extrapolation or falloff
        fn = os.path.join(os.path.dirname(__file__),
                '..', 'input', 'cooling',
                '%s_ion_by_ion_cooling.h5' % atom_name.lower())
        data = h5py.File(fn, 'r')

        ### Intepolate values within table values ###
        vals = np.interp(network.T, data['T'], data['%s' %(ion_name)])

        end_method = 1 # 0 = extrapolation, 1 = gaussian falloff

        if end_method == 0:
            ### Extrapolation in logspace ###
            # convert to log space
            vals = np.log10(vals)
            logT = np.log10(network.T)
            logdataT = np.log10(data['T'])
            logdataS = np.log10(data['%s' %(ion_name)])

            # extrapolate
            extrapdown = logdataS[0] + \
                (logT - logdataT[0]) * (logdataS[0] - logdataS[1]) \
                / (logdataT[0] - logdataT[1])
            vals[logT < logdataT[0]] = extrapdown[logT < logdataT[0]]
            extrapup = logdataS[-1] + \
                (logT - logdataT[-1]) * (logdataS[-1] - logdataS[-2]) \
                / (logdataT[-1] - logdataT[-2])
            vals[logT > logdataT[-1]] = extrapup[logT > logdataT[-1]]

            # convert back to linear
            vals = 10.0**vals
        if end_method == 1:
            ### Gaussian falloff when values extend beyond table values ###
            # rename some variables to symplify code
            T = network.T
            dataT = data['T']
            dataS = data['%s' %(ion_name)]

            # compute gaussian tails
            gaussdown = dataS[0] * (tiny/dataS[0])**(((T - dataT[0])/(T[0] - dataT[0])))**2
            vals[T < dataT[0]] = gaussdown[T < dataT[0]]
            gaussup = dataS[-1] * (tiny/dataS[-1])**(((T - dataT[-1])/(T[-1] - dataT[-1])))**2
            vals[T > dataT[-1]] = gaussup[T > dataT[-1]]

        data.close()
        return vals

    ion_cooling_action = CoolingAction("%s_c" % species.name, #name
             "-%s_c * %s * de" %(species.name, species.name)) #equation
    ion_cooling_action.tables["%s_c" % species.name] = cooling_rate
    new_rates.append("%s_c" % species.name)
    return new_rates

def ion_photoheating_rate(species, photo_background='HM12'):
    if chu is None: raise ImportError
    #ion_name = chu.zion2name(np.int(species.number),
    #                         np.int(species.free_electrons + 1))
    ion_name = species.name.lower()
    if "_" not in ion_name:
        print ("Name must be in 'Ion Species' format.")
        raise RuntimeError
    element_name = ion_name.split("_")[0]
    ion_state = int(ion_name.split("_")[1])
    species_ph = species.name
    de = species_registry['de']
    new_rates = []

    def photoheating_rate(network):
        # Read in photoheating rates generated from Ben Oppenheimer's data
        # (from: http://noneq.strw.leidenuniv.nl/)
        # and do linear interpolation and then recompute
        # the ends with either an extrapolation or falloff
        # NOTE: these rates do the interpolation as a function fo redshift
        fn = os.path.join(os.path.dirname(__file__),
                '..', 'input', 'photoheating',
                '%s_ion_by_ion_photoheating_%s.h5' %(element_name, photo_background))
        f = h5py.File(fn,'r')
        #f = h5py.File('../input/photoheating/%s_ion_by_ion_photoheating_%s.h5' %(element_name,
        #                                                         photo_background))

        ### Intepolate values within table values ###
        vals = np.interp(network.z, f['z'], f['%s' %(ion_name)])

        end_method = 0 # 0 = extrapolation, 1 = gaussian falloff

        if end_method == 0:
            ### Extrapolation in logspace ###
            # convert to log space
            vals = np.log10(vals)
            logz = np.log10(network.z)
            logdataz = np.log10(f['z'])
            logdataS = np.log10(f['%s' %(ion_name)])

            # extrapolate
            extrapdown = logdataS[0] + \
                (logz - logdataz[0]) * (logdataS[0] - logdataS[1]) \
                / (logdataz[0] - logdataz[1])
            vals[logz < logdataz[0]] = extrapdown[logz < logdataz[0]]
            extrapup = logdataS[-1] + \
                (logz - logdataz[-1]) * (logdataS[-1] - logdataS[-2]) \
                / (logdataz[-1] - logdataz[-2])
            vals[logz > logdataz[-1]] = extrapup[logz > logdataz[-1]]

            # convert back to linear
            vals = 10.0**vals
        if end_method == 1:
            ### Gaussian falloff when values extend beyond table values ###
            # rename some variables to symplify code
            z = network.z
            dataz = f['z']
            dataS = f['%s' %(ion_name)]

            # compute gaussian tails
            gaussdown = dataS[0] * (tiny/dataS[0])**(((z - dataz[0])/(z[0] - dataz[0])))**2
            vals[z < dataz[0]] = gaussdown[z < dataz[0]]
            gaussup = dataS[-1] * (tiny/dataS[-1])**(((z - dataz[-1])/(z[-1] - dataz[-1])))**2
            vals[z > dataz[-1]] = gaussup[z > dataz[-1]]

        f.close()
        return vals

    if species_ph in species_registry:
        species_ph = species_registry[species_ph]
        ion_cooling_action = CoolingAction("%s_ph" % species.name, #name
                                           "%s_ph * %s" %(species.name, species.name)) #equation
        ion_cooling_action.tables["%s_ph" % species.name] = photoheating_rate
        new_rates.append("%s_ph" % species.name)
    return new_rates

