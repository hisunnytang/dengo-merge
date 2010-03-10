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
from scipy.integrate import ode

from reactions import QuantitiesTable, reaction_table
import reactions as rxn
from initialize_rate_tables import reaction_rates_table, tiny
species_list = [rxn.species_table[s] for s in
            ["HI","HII","HeI","HeII","HeIII","HM","de","H2I","H2II","T"]]

rho = 100.0
fracs = dict(HI    = 1.0 - 1e-4 - 1e-6,
             HII   = 1e-4,
             HeI   = tiny,
             HeII  = tiny,
             HeIII = tiny,
             HM    = tiny,
             de    = 1e-4,
             H2I   = 1e-6,
             H2II  = tiny)
derivs = dict( ( (i, tiny) for i in fracs ) )
iv = []
for f,v in fracs.items(): 
    fracs[f] *= rho
fracs["T"] = 400.0

s = QuantitiesTable(species_list, fracs)
d = QuantitiesTable(species_list, derivs)

def update(quantities, derivatives, dt):
    for species in quantities.species_list:
        quantities[species.name] += dt * derivatives[species.name]

def calc_derivs(quantities, derivatives, reactions):
    for n, reaction in sorted(reactions.items()):
        reaction(quantities, derivatives)

def plot_vals(vals, prefix="values"):
    import matplotlib;matplotlib.use("Agg");import pylab
    for v in vals:
        if v == 't': continue
        pylab.clf()
        pylab.loglog(vals['t'], vals[v], '-x')
        pylab.savefig("%s_%s.png" % (prefix, v))
        print "Saving: %s_%s.png" % (prefix, v)

tf = 1e9
ti = 0.0
dt = 1e5
vals = dict(((i.name, []) for i in species_list))
vals['t'] = []

def append_vals(all_vals, new_vals, t):
    all_vals['t'].append(t)
    for species in new_vals._names:
        if species not in all_vals: continue
        all_vals[species].append(new_vals[species])

while ti < tf:
    print ti/tf
    append_vals(vals, s, ti)
    calc_derivs(s, d, reaction_table)
    update(s, d, dt)
    d.values *= 0.0
    ti += dt
plot_vals(vals)
