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

# This is a very simple solver creator.

from .utils import ensure_species, ensure_reaction

def write_reaction(rxn, snames = None, rnames = None):
    if snames is None: snames = {}
    if rnames is None: rnames = {}
    rxn = ensure_reaction(rxn)
    sc = []
    for n, species in rxn.left_side:
        sn = snames.get(species.name, species.name)
        if species.weight != 1:
            sn = "( %s / %0.1f )" % (sn, species.weight)
        if n == 1:
            sc.append(sn)
        else:
            sc.append(" * ".join(n*[sn]))
    rn = rnames.get(rxn.name, rxn.name)
    c = "%s_v = %s * ( %s );" % (rn, rn, " * ".join(sc))
    print c

def write_species_deriv(species, reactions, snames = None, rnames = None):
    species = ensure_species(species)
    if snames is None: snames = {}
    if rnames is None: rnames = {}
    f, d = [], []
    for rxn in reactions:
        rxn = ensure_reaction(rxn)
        rn = "%s_v" % (rnames.get(rxn.name, rxn.name))
        net = rxn.net_change(species.name)
        if abs(net) == 1:
            nn = ""
        elif net != 0:
            nn = "%0.1f * " % net
        if net > 0:
            f.append("%s%s" % (nn, rn))
        elif net < 0:
            d.append("%s%s" % (nn, rn))
    if species.weight != 1.0:
        snn = "%0.1f  * " % (species.weight)
    else:
        snn = ""
    c = "%s_p = %s( %s ) - ( %s )" % (
            snames.get(species.name, species.name), snn,
            " + ".join(f), " + ".join(d))
    print c
