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
import sympy
from primordial_rates import species_table

species_symbols = dict( [
    (sname, sympy.Symbol(sname)) for
     sname in species_table ])

# We set up all the symbols we are aware of here.

#class CVODEPrinter(sympy.Printer):
#    pass

# Note that we actually re-implement here rather than subclass.
# That's because we want our cooling rates to be able to have different min/max
# temperature, and they may also behave differently in other ways from normal
# reaction rates.

class CoolingRate(object):
    def __init__(self, name, values):
        self.values = values
        self.name = name
        self.s = sympy.Symbol(self.name)

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

        ReactionRate.__init__(self, name, values)

cooling_rates_table = dict()

class CoolingAction(object):
    def __init__(self, rates, species, equation):
        # Equation should be a string
        self.rates = rates
        self.species = species
        symbols = dict( [(a, cooling_rates_table[a].s) for a in rates])
        symbols.update(dict([(a, species_symbols[a]) for a in species]))
        self.equation = eval(equation, {}, symbols)

cooling_action_table = dict()

# Sample implementation of collisional excitation

cooling_rates_table["ceHI"] = CoolingRate("ceHI", na.zeros(1024, dtype='float64'))

# prepend with 'a' for 'active'
cooling_action_table["aceHI"] = CoolingAction(
    ["ceHI"], ["HI","de"], "ceHI * HI * de")


