"""
Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: Columbia University
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2012 Matthew Turk.  All Rights Reserved.

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
from .reaction_table import ReactionTable
import numpy as na
from chemistry_constants import tevk, tiny, mh

class ChemicalNetwork(object):

    def __init__(self):
        self.reactions = ReactionTable()
        self.rates = rates
        self.cooling = cooling
        self.required_species = set([])

    def add_reaction(self, reaction):
        if reaction.name in self.reactions:
            raise RuntimeError
        self.reactions[reaction.name] = reaction
        for n, s in reaction.left_side:
            self.required_species.add(s)
        for n, s in reaction.right_side:
            self.required_species.add(s)
    
    def init_temperature(self, T_bounds = (1, 1e8), n_bins=1024):
        self.n_bins = 1024
        self.T = na.logspace(na.log(T_bounds[0]),
                             na.log(T_bounds[1]),
                             n_bins, base = na.e)
        self.logT = na.log(self.T)
        self.tev = self.T / tevk
        self.logtev = na.log(self.tev)
        self.T_bounds = T_bounds

    def __iter__(self):
        for rxn in sorted(self.reactions, key=lambda a: a.name):
            yield rxn
