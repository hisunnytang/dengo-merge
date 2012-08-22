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

import numpy as na
from chemistry_constants import tevk, tiny, mh
import types
import sympy
import h5py

try:
    import chianti.core as ch
except ImportError:
    ch = None

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
    if "_" not in species.name:
        print "Name must be in ChiantiPy format."
        raise RuntimeError
    ion_name = species.name
    element_name = ion_name.split("_")[0]
    ion_state = int(ion_name.split("_")[1])
    species_i = "%s_%s" % (element_name, ion_state + 1)
    species_r = "%s_%s" % (element_name, ion_state - 1)
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

class Species(object):
    def __init__(self, name, weight, free_electrons = 0.0, equilibrium = False,
                 computed = False):
        self.name = name
        self.weight = weight
        self.free_electrons = free_electrons
        self.equilibrium = equilibrium
        self.computed = computed
        #self.symbol = sympy.IndexedBase(name, (count_m,))
        self.symbol = sympy.Symbol("%s[i]" % name)
        if equilibrium and computed: raise RuntimeError
        if equilibrium: raise RuntimeError
        species_registry[name] = self

    def number_density(self, quantities):
        return quantities[self.name]/self.weight

    def __repr__(self):
        return "Species: %s" % (self.name)

class Constraint(object):
    pass

class ChargeConservation(Constraint):
    def __call__(self, quantities, up_derivatives, down_derivatives, dt):
        quantities["de"] = (quantities["HII"]
            + quantities["HeII"] / 4.0
            + quantities["HeIII"] / 2.0
            + quantities["H2II"] / 2.0
            - quantities["HM"])
        return
        quantities["de"] = 0.0
        for q in quantities.species_list:
            quantities["de"] += q.free_electrons * q.number_density(quantities)

class Floor(Constraint):
    def __call__(self, quantities, up_derivatives, down_derivatives, dt):
        for s in quantities.species_list:
            quantities[s.name] = max(quantities[s.name], 1e-30)

class ChemicalHeating(Constraint):
    def __call__(self, quantities, up_derivatives, down_derivatives, dt):
        # Get the total mass
        rho = sum(quantities[i] for i in
                ["HI","HII","HM","H2I","H2II","HeI","HeII","HeIII"])
        dH2 = (up_derivatives["H2I"] - down_derivatives["H2I"])*dt
        quantities["T"] += dH2 * 51998.0/rho

constraints = [ChemicalHeating(), ChargeConservation(), Floor()]

class QuantitiesTable(object):
    def __init__(self, species_list, initial_values = None):
        self._names = {}
        for i,s in enumerate(species_list):
            self._names[s.name] = i
        self.species_list = species_list
        self.values = na.zeros(len(species_list), dtype='float64')
        if initial_values is not None:
            for s, v in initial_values.items():
                self.values[self._names[s]] = v

    def __getitem__(self, name):
        return self.values[self._names[name]]

    def __setitem__(self, name, value):
        self.values[self._names[name]] = value

    def __iter__(self):
        for i, s in enumerate(self.species_list):
            yield self.values[self._names[s.name]]

    def get_by_name(self, name):
        return self.species_list[self._names[name]]

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
        ta_sym = dict((n, sympy.Symbol("%s[i]" % n)) for n in self.tables)
        self.table_symbols.update(ta_sym)
        #tp_sym = dict((n, sympy.IndexedBase(n, (count_m,))) for n in self.temporaries))
        tp_sym = dict((n, sympy.Symbol("%s[i]" % n)) for n in self.temporaries)
        self.temp_symbols.update(tp_sym)
        symbols.update(self.table_symbols)
        symbols.update(self.temp_symbols)
        self._eq = eval(self._equation, symbols)
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
    if "_" not in species.name:
        print "Name must be in 'Ion Species' format."
        raise RuntimeError
    ion_name = species.name
    element_name = ion_name.split("_")[0]
    ion_state = int(ion_name.split("_")[1])
    species_c = ion_name
    de = species_registry['de']
    new_rates = []

    def cooling_rate(network):
        # Read in cooling rates from Gnat & Ferland 2012
        # and do linear interpolation
        f = h5py.File('dengo/%s_ion_by_ion_cooling.h5' %(element_name))
        data = f['Table']
        vals = na.interp(network.T, data['T'], data['%s' %(ion_name)])
        f.close()
        return vals

    if species_c in species_registry:
        species_c = species_registry[species_c]
        ion_cooling_action = CoolingAction("%s_c" % ion_name, #name
                                           "-%s_c * %s * de" %(ion_name, ion_name)) #equation
        ion_cooling_action.tables["%s_c" % ion_name] = cooling_rate
        new_rates.append("%s_c" % ion_name)
    return new_rates

class CVODEPrinter(sympy.printing.str.StrPrinter):
    def __init__(self, template_vars, *args, **kwargs):
        sympy.printing.str.StrPrinter.__init__(self, *args, **kwargs)
        self.template_vars = template_vars

    def _print_Symbol(self, symbol):
        s = str(symbol)
        if s in species_symbols:
            vname = self.template_vars["species_varnames"][s]
            return "%s" % vname
        elif s in cooling_rates_table:
            cid = self.template_vars["cooling_rate_ids"][s]
            vname = "data->cooling_storage[%s][cell]" % cid
            return "%s" % vname
        elif s in known_variables:
            return s
        elif s in user_data_symbols:
            return "data->%s" % s
        elif s in user_data_cell_vars:
            return "data->%s[cell]" % s
        else:
            raise RuntimeError

    def _print_Derivative(self, expr):
        # Using the sympy printing docs as a hint here!
        rate = str(expr.args[0].func)
        if rate not in cooling_rates_table:
            raise RuntimeError
        cid = self.template_vars["cooling_rate_ids"][rate]
        vname = "cooling_slopes[%s]" % cid
        return "%s" % vname

    def _print_Pow(self, expr):
        PREC = sympy.printing.precedence.precedence(expr)
        return 'pow(%s,%s)'%(self.parenthesize(expr.base, PREC),
                             self.parenthesize(expr.exp, PREC))

