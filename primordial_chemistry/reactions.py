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
from initialize_rate_tables import reaction_rates_table

class Reaction(object):
    def __init__(self, rate, left_side, right_side):
        self.name = rate.replace("k","r")
        self.rate = reaction_rates_table[rate]
        self.left_side = left_side
        self.right_side = right_side

    def __call__(self, quantities, derivatives):
        # We just calculate our net derivatives and stick them in the right
        # place
        r = 0.0
        for n, s in self.left_side:
            r += s.number_density(quantities)**n
        k = self.rate(quantities)
        r *= k
        for n, s in self.left_side:
            derivatives[s.name] -= r * n * s.weight
        for n, s in self.right_side:
            derivatives[s.name] += r * n * s.weight
        return r

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

class Constraint(object):
    pass

class ChargeConservation(Constraint):
    def __call__(self, quantities):
        pass

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

species_table = dict(
    HI = Species("HI", 1.0),
    HII = Species("HII", 1.0, 1.0),
    HeI = Species("HeI", 4.0),
    HeII = Species("HeII", 4.0, 1.0),
    HeIII = Species("HeIII", 4.0, 2.0),
    de = Species("de", 1.0),
    HM = Species("HM", 1.0, -1.0),
    H2I = Species("H2I", 2.0),
    H2II = Species("H2II", 2.0, 1.0),
    T = Species("T", 0.0, 0.0),
)
for s in species_table:
    exec("%s = species_table['%s']" % (s,s))

reaction_table = dict(
    r01 = Reaction('k01', [   (1,HI),   (1,de)], [  (1,HII),   (2,de)]),
    r02 = Reaction('k02', [  (1,HII),   (1,de)], [   (1,HI),         ]),
    r03 = Reaction('k03', [  (1,HeI),   (1,de)], [ (1,HeII),   (2,de)]),
    r04 = Reaction('k04', [ (1,HeII),   (1,de)], [  (1,HeI),         ]),
    r05 = Reaction('k05', [ (1,HeII),   (1,de)], [(1,HeIII),   (2,de)]),
    r06 = Reaction('k06', [(1,HeIII),   (1,de)], [ (1,HeII),         ]),
    r07 = Reaction('k07', [   (1,HI),   (1,de)], [   (1,HM),         ]),
    r08 = Reaction('k08', [   (1,HM),   (1,HI)], [  (1,H2I),   (1,de)]),
    r09 = Reaction('k09', [   (1,HI),  (1,HII)], [ (1,H2II),         ]),
    r10 = Reaction('k10', [ (1,H2II),   (1,HI)], [  (1,H2I),  (1,HII)]),

    r11 = Reaction('k11', [  (1,H2I),  (1,HII)], [  (1,H2II),  (1,HI)]),
    r12 = Reaction('k12', [  (1,H2I),   (1,de)], [  (2,HII),   (1,de)]),
    r13 = Reaction('k13', [  (1,H2I),   (1,HI)], [   (3,HI),         ]), #3b
    r14 = Reaction('k14', [   (1,HM),   (1,de)], [   (1,HI),   (2,de)]),
    r15 = Reaction('k15', [   (1,HM),   (1,HI)], [   (2,HI),   (1,de)]),
    r16 = Reaction('k16', [   (1,HM),  (1,HII)], [   (2,HI),         ]),
    r17 = Reaction('k17', [   (1,HM),  (1,HII)], [ (1,H2II),   (1,de)]),
    r18 = Reaction('k18', [ (1,H2II),   (1,de)], [   (2,HI),         ]),
    r19 = Reaction('k19', [ (1,H2II),   (1,HM)], [   (1,HI),  (1,H2I)]),
    r21 = Reaction('k21', [   (2,HI),  (1,H2I)], [  (2,H2I),         ]), #3b
    r22 = Reaction('k22', [   (2,HI),   (1,HI)], [  (1,H2I),   (1,HI)]), #3b
    r23 = Reaction('k23', [  (1,H2I),  (1,H2I)], [   (2,HI),  (1,H2I)]), #3b
)

for s in reaction_table:
    exec("%s = reaction_table['%s']" % (s,s))
