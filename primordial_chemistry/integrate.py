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

from reactions import QuantitiesTable, reaction_table, constraints
import reactions as rxn
from initialize_rate_tables import reaction_rates_table, tiny

# Get our list of species, but it's important to note this is a list of Species
# instances, not just names.
species_list = [rxn.species_table[s] for s in
            ["HI","HII","HeI","HeII","HeIII","HM","de","H2I","H2II","T"]]

rho = 100.0

# This is the initial fraction of every species
fracs = dict(HI    = 1.0 - 1e-2,
             HII   = 1e-2,
             HeI   = tiny,
             HeII  = tiny,
             HeIII = tiny,
             HM    = tiny,
             de    = 1e-2,
             H2I   = tiny,
             H2II  = tiny)

# Initialize our derivatives as tiny
derivs = dict( ( (i, tiny) for i in fracs ) )

# Convert our fractions => densities
for f,v in fracs.items(): 
    fracs[f] *= rho

# Initialize our temperature
fracs["T"] = 100.0

s = QuantitiesTable(species_list, fracs)
ud = QuantitiesTable(species_list, derivs)
dd = QuantitiesTable(species_list, derivs)

def euler_update(quantities, up_derivatives, down_derivatives, dt):
    for species in quantities.species_list:
        net = (up_derivatives[species.name] - down_derivatives[species.name])
        quantities[species.name] += dt * net

def bdf_update(quantities, up_derivatives, down_derivatives, dt):
    for species in quantities.species_list:
        n = species.name
        q = quantities[n]
        quantities[species.name] = ((q + dt*up_derivatives[n]) / 
                                    (1.0 + dt*down_derivatives[n]/q))

def calc_derivs(quantities, up_derivatives, down_derivatives, reactions):
    for n, reaction in sorted(reactions.items()):
        reaction(quantities, up_derivatives, down_derivatives)

def plot_vals(vals, prefix="values", norm = 1.0):
    import matplotlib;matplotlib.use("Agg");import pylab
    for v in vals:
        if v == 't': continue
        pylab.clf()
        pylab.loglog(vals['t']/(365*24*3600), vals[v]/norm, '-x')
        pylab.savefig("%s_%s.png" % (prefix, v))
        print "Saving: %s_%s.png" % (prefix, v)

tf = 1e9 * (365*24*3600)
ti = 0.0 * (365*24*3600)
dt = 1e5 * (365*24*3600)
vals = dict(((i.name, []) for i in species_list))
vals['t'] = []

def append_vals(all_vals, new_vals, t):
    all_vals['t'].append(t)
    for species in new_vals._names:
        if species not in all_vals: continue
        all_vals[species].append(new_vals[species])

nsteps = 0
while ti < tf:
    if (nsteps % 1000 == 0): print ti/tf
    nsteps += 1
    append_vals(vals, s, ti)
    calc_derivs(s, ud, dd, reaction_table)
    bdf_update(s, ud, dd, dt)
    for c in constraints: c(s)
    ud.values *= 0.0 + tiny
    dd.values *= 0.0 + tiny
    ti += dt
for v in vals: vals[v] = na.array(vals[v])
plot_vals(vals, norm=rho)
