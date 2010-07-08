"""
Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: UCSD
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2010 Matthew Turk.  All Rights Reserved.

  This file is part of the primordial_chemistry package.

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
        print self.considered

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

class Species(object):
    name = None
    weight = 1.0
    free_electrons = 0.0

    def __init__(self, name, weight, free_electrons = 0.0):
        self.name = name
        self.weight = weight
        self.free_electrons = free_electrons

    def number_density(self, quantities):
        return quantities[self.name]/self.weight

    def __repr__(self):
        return "Species: %s" % (self.name)

class Constraint(object):
    pass

class ChargeConservation(Constraint):
    def __call__(self, quantities):
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
    def __call__(self, quantities):
        for s in quantities.species_list:
            quantities[s.name] = max(quantities[s.name], 1e-30)

constraints = [ChargeConservation(), Floor()]

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

