"""
Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: UCSD
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2010 Matthew Turk.  All Rights Reserved.

  This file is part of the dengo package.

  This file is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import numpy as np
from chemistry_constants import tevk, tiny, mh
import types
import sympy
import h5py
import docutils.utils.roman as roman

try:
    import chianti.core as ch
    import chianti.util as chu
except (ImportError, KeyError):
    ch = None
    chu = None

reaction_registry = {}
cooling_registry = {}
species_registry = {}

count_m = sympy.Symbol('m', integer=True)
index_i = sympy.Idx('i', count_m)

class ReactionCoefficient(sympy.Symbol):
    def _eval_derivative(self, s):
        if s == self.energy:
            return sympy.Symbol("r%s" % self)
        else:
            return super(ReactionCoefficient, self)._eval_derivative(s)

    energy = None

    @property
    def free_symbols(self):
        if self.energy is not None:
            return set([self, self.energy])
        else:
            return set([sefl])

class Reaction(object):
    def __init__(self, name, coeff_fn, left_side, right_side):
        self.name = name
        self.coeff_fn = coeff_fn
        #self.coeff_sym = sympy.IndexedBase(name, (count_m,))
        self.coeff_sym = ReactionCoefficient("%s[i]" % name)
        self.left_side = left_side
        self.right_side = right_side
        self.considered = set( (s.name for n, s in left_side + right_side) )
        reaction_registry[name] = self # Register myself

    updated = False
    def update(self, e):
        if self.updated: return
        self.coeff_sym.energy = e
        self.updated = True

    def __contains__(self, c):
        if isinstance(c, types.StringTypes):
            return c in self.down_species + self.up_species
        return c in (s for n, s in self.left_side + self.right_side)

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
        return up - down

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

    def __repr__(self):
        a = "%s : " % self.name \
          + " + ".join( ["%s*%s" % (i, s.name) for i, s in self.left_side] ) \
          + " => " \
          + " + ".join( ["%s*%s" % (i, s.name) for i, s in self.right_side] )
        return a

    @classmethod
    def create_reaction(cls, name, left_side, right_side):
        def _w(f):
            rxn = cls(name, f, left_side, right_side)
        return _w

    def species_equation(self, species):
        if isinstance(species, types.StringTypes):
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
        return eq

reaction = Reaction.create_reaction
def chianti_rate(species):
    if ch is None: raise ImportError
    if chu is None: raise ImportError
    ion_name = chu.zion2name(np.int(species.number),
                             np.int(species.free_electrons + 1))
    if "_" not in ion_name:
        print "Name must be in ChiantiPy format."
        raise RuntimeError
    element_name = ion_name.split("_")[0]
    ion_state = int(ion_name.split("_")[1])
    species_i = "%s%s" % (element_name.capitalize(), roman.toRoman(ion_state + 1))
    if ion_state != 1:
        species_r = "%s%s" % (element_name.capitalize(), roman.toRoman(ion_state - 1))
    else:
        species_r = None
    de = species_registry['de']
    new_rates = []
    
    def ion_rate(network):
        ion = ch.ion(ion_name, temperature = network.T)
        ion.ionizRate()
        vals = ion.IonizRate['rate']
        return vals
    if species_i in species_registry:
        species_i = species_registry[species_i]
        Reaction("%s_i" % species.name, ion_rate,
                 [(1, species), (1, de)], # left side
                 [(1, species_i), (2, de)]) # right side
        new_rates.append("%s_i" % species.name)

    def rec_rate(network):
        ion = ch.ion(ion_name, temperature = network.T)
        ion.recombRate()
        vals = ion.RecombRate['rate']
        return vals
    if species_r in species_registry:
        species_r = species_registry[species_r]
        Reaction("%s_r" % species.name, rec_rate,
                 [(1, species), (1, de)], # left side
                 [(1, species_r), ]) # right side
        new_rates.append("%s_r" % species.name)
    return new_rates

def ion_photoionization_rate(species, photo_background='HM12'):
    if chu is None: raise ImportError
    ion_name = chu.zion2name(np.int(species.number),
                             np.int(species.free_electrons + 1))
    if "_" not in ion_name:
        print "Name must be in 'Ion Species' format."
        raise RuntimeError
    element_name = ion_name.split("_")[0]
    ion_state = int(ion_name.split("_")[1])
    species_pi = "%s%s" % (element_name.capitalize(), roman.toRoman(ion_state + 1))
    de = species_registry['de']
    new_rates = []

    def photoionization_rate(network):
        # Read in photoheating rates generated from Ben Oppenheimer's data
        # (from: http://noneq.strw.leidenuniv.nl/)
        # and do linear interpolation and then recompute
        # the ends with either an extrapolation or falloff
        # NOTE: these rates do the interpolation as a function fo redshift
        f = h5py.File('input/photoionization/%s_ion_by_ion_photoionization_%s.h5'
                      %(element_name, photo_background))
        
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

class Species(object):
    def __init__(self, name, number, weight, free_electrons = 0.0, equilibrium = False,
                 computed = False):
        self.name = name
        self.number = number
        self.weight = weight
        self.free_electrons = free_electrons
        self.equilibrium = equilibrium
        self.computed = computed
        #self.symbol = sympy.IndexedBase(name, (count_m,))
        self.symbol = sympy.Symbol("%s" % name)
        if equilibrium and computed: raise RuntimeError
        if equilibrium: raise RuntimeError
        species_registry[name] = self

    def number_density(self, quantities):
        return quantities[self.name]/self.weight

    def __repr__(self):
        return "Species: %s" % (self.name)

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

    @property
    def equation(self):
        if self._eq is not None: return self._eq
        symbols = dict((n, s.symbol) for n, s in species_registry.items())
        #ta_sym = dict((n, sympy.IndexedBase(n, (count_m,))) for n in self.tables))
        ta_sym = dict((n, sympy.Symbol("%s_%s[i]" % (self.name, n))) for n in self.tables)
        self.table_symbols.update(ta_sym)
        #tp_sym = dict((n, sympy.IndexedBase(n, (count_m,))) for n in self.temporaries))
        tp_sym = dict((n, sympy.Symbol("%s" % (n))) for n in self.temporaries)
        self.temp_symbols.update(tp_sym)
        symbols.update(self.table_symbols)
        symbols.update(self.temp_symbols)
        self._eq = eval(self._equation, symbols)
        for n, e in self.temporaries.items():
            e = sympy.sympify(e)
            for n2, e2 in ta_sym.items():
                e = e.subs(n2, e2)
            self._eq = self._eq.subs(n, e)
        return self._eq

    @property
    def species(self):
        self.equation
        bad = set(self.temp_symbols.values() + self.table_symbols.values())
        species = set([])
        #for s in self.equation.atoms(sympy.IndexedBase):
        for s in self.equation.atoms(sympy.Symbol):
            if s not in bad:
                species.add(species_registry[str(s).replace("[i]","")])
        return species

    def table(self, func):
        self.tables[func.func_name] = func

    def temporary(self, name, eq):
        self.temporaries[name] = eq
        
    @classmethod
    def create_cooling_action(cls, name, equation):
        obj = cls(name, equation)
        def _W(f):
            f(obj)
        return _W

cooling_action = CoolingAction.create_cooling_action

def ion_cooling_rate(species):
    if chu is None: raise ImportError
    ion_name = chu.zion2name(np.int(species.number),
                             np.int(species.free_electrons + 1))
    if "_" not in ion_name:
        print "Name must be in 'Ion Species' format."
        raise RuntimeError
    element_name = ion_name.split("_")[0]
    ion_state = int(ion_name.split("_")[1])
    species_c = species.name
    de = species_registry['de']
    new_rates = []

    def cooling_rate(network):
        # Read in cooling rates from Gnat & Ferland 2012
        # and do linear interpolation and then recompute
        # the ends with either an extrapolation or falloff
        data = h5py.File('input/cooling/%s_ion_by_ion_cooling.h5'
                         %(element_name))
        
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

    if species_c in species_registry:
        species_c = species_registry[species_c]
        ion_cooling_action = CoolingAction("%s_c" % species.name, #name
                                           "-%s_c * %s * de" %(species.name, species.name)) #equation
        ion_cooling_action.tables["%s_c" % species.name] = cooling_rate
        new_rates.append("%s_c" % species.name)
    return new_rates

def ion_photoheating_rate(species, photo_background='HM12'):
    if chu is None: raise ImportError
    ion_name = chu.zion2name(np.int(species.number),
                             np.int(species.free_electrons + 1))
    if "_" not in ion_name:
        print "Name must be in 'Ion Species' format."
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
        f = h5py.File('input/photoheating/%s_ion_by_ion_photoheating_%s.h5' %(element_name,
                                                                 photo_background))
        
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

