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

try:
    import chianti.core as ch
except ImportError:
    ch = None

reaction_registry = {}
species_registry = {}

class Reaction(object):

    def __init__(self, name, coeff_fn, left_side, right_side):
        self.name = name
        self.coeff_fn = coeff_fn
        self.left_side = left_side
        self.right_side = right_side
        self.considered = set( (s.name for n, s in left_side + right_side) )
        reaction_registry[name] = self # Register myself

    def __contains__(self, c):
        if isinstance(c, types.StringTypes):
            return c in self.down_species + self.up_species
        return c in (s for n, s in self.left_side + self.right_side)

    @property
    def down_species(self):
        return [s.name for n, s in self.left_side]

    @property
    def up_species(self):
        return [s.name for n, s in self.right_side]

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

    def print_c_calculation(self):
        # This function assumes we're inside a loop
        st = ""
        st += "for (i = 0; i < data->nvals ; i++)\n"
        st += "{\n"
        st += "  tv = 0.0\n"
        for n, s in self.left_side:
            #r *= s.number_density(quantities)**n
            st += "  tv -= pow(y[%s], %s);\n" % (s.name, n)
        for n, s in self.right_side:
            #up_derivatives[s.name] += r * n * s.weight
            st += "  tv += y[%s]*%s;\n" % (s.name, n)
        st += "  ydot[i] * data->stride + ri;\n"
        st += "}\n"
        return st

    def print_c_reaction_value(self, species_varnames = None):
        st = []
        for n, s in self.left_side:
            st += [" * ".join([
                s.print_c_number_density(species_varnames[s.name])
                for i in xrange(n)])]
        return " * ".join(st)

    @classmethod
    def create_reaction(cls, name, left_side, right_side):
        def _w(f):
            rxn = cls(name, f, left_side, right_side)
        return _w

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
        if equilibrium and computed: raise RuntimeError
        if equilibrium: raise RuntimeError
        species_registry[name] = self

    def number_density(self, quantities):
        return quantities[self.name]/self.weight

    def print_c_convert_number_density(self, input, output):
        return "%s = %s / %s;" % (output, input, self.weight)

    def print_c_number_density(self, input):
        return "(%s / %s)" % (input, self.weight)

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

