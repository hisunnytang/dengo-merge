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
import numpy as na
from chemistry_constants import tevk, tiny, mh
from .reaction_classes import reaction_registry
import types

class ChemicalNetwork(object):

    def __init__(self):
        self.reactions = {}
        self.required_species = set([])

    def add_reaction(self, reaction):
        if isinstance(reaction, types.StringTypes):
            reaction = reaction_registry[reaction]
        if reaction.name in self.reactions:
            raise RuntimeError
        self.reactions[reaction.name] = reaction
        
        for n, s in reaction.left_side:
            self.required_species.add(s)
        for n, s in reaction.right_side:
            self.required_species.add(s)
        print "Adding reaction: %s" % reaction
    
    def init_temperature(self, T_bounds = (1, 1e8), n_bins=1024):
        self.n_bins = 1024
        self.T = na.logspace(na.log(T_bounds[0]),
                             na.log(T_bounds[1]),
                             n_bins, base = na.e)
        self.logT = na.log(self.T)
        self.tev = self.T / tevk
        self.logtev = na.log(self.tev)
        self.T_bounds = T_bounds

    def species_reactions(self, species):
        tr = []
        for rn, rxn in sorted(self.reactions.items()):
            if species in rxn: tr.append(rxn)
        return tr

    def __iter__(self):
        for rname, rxn in sorted(self.reactions.items()):
            yield rxn
