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

reaction_rates_table = dict()

class ReactionRate(object):
    def __init__(self, name, values):
        self.values = values
        self.name = name

    def __call__(self, quantities):
        T = quantities['T']
        return na.interp(T, self.T, self.values)

    @classmethod
    def init_temperature(cls, T_bounds, n_bins=1024):
        cls.n_bins = 1024
        cls.T = na.logspace(
            na.log10(T_bounds[0]), na.log10(T_bounds[1]), n_bins)
        cls.logT = na.log(cls.T)
        cls.tev = cls.T / tevk
        cls.logtev = na.log(cls.tev)
        cls.T_bounds = T_bounds

class Reaction(object):
    def __init__(self, rate, left_side, right_side):
        self.name = rate.replace("k","r")
        self.rate = reaction_rates_table[rate]
        self.left_side = left_side
        self.right_side = right_side
        self.considered = set( (s.name for n, s in left_side + right_side) )

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

    def print_c_reaction_value(self, species_varnames):
        st = []
        for n, s in self.left_side:
            st += [" * ".join([
                s.print_c_number_density(species_varnames[s.name])
                for i in xrange(n)])]
        return " * ".join(st)

class Species(object):
    def __init__(self, name, weight, free_electrons = 0.0, equilibrium = False,
                 computed = False):
        self.name = name
        self.weight = weight
        self.free_electrons = free_electrons
        self.equilibrium = equilibrium
        self.computed = computed
        if equilibrium and computed: raise RuntimeError

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
