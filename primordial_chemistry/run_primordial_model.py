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
import primordial_rates # This initializes reaction_rates_table
from primordial_rates import reaction_table, species_table
from reaction_classes import QuantitiesTable, constraints, \
                reaction_rates_table
from chemistry_constants import tiny
from integration_helpers import euler_update, bdf_update, calc_derivs, \
    all_euler_update, all_bdf_update, all_calc_derivs, all_plot_vals, \
    update_vals, calculate_timestep

# Get our list of species, but it's important to note this is a list of Species
# instances, not just names.
species_list = [species_table[s] for s in
            ["HI","HII","HeI","HeII","HeIII","HM","de","H2I","H2II","T"]]

Temperature = 300
rho = 1.0e10 # total rho in amu/cc
X   = 0.01 # ionization fraction
fH2 = 0.50 # ionization fraction

# This is the initial fraction of every species
fracs = dict(HI    = 1.0 - X - fH2,
             HII   = X,
             HeI   = tiny,
             HeII  = tiny,
             HeIII = tiny,
             HM    = tiny,
             de    = X,
             H2I   = fH2,
             H2II  = tiny)

# Initialize our derivatives as tiny
derivs = dict( ( (i, tiny) for i in fracs ) )

# Convert our fractions => densities
for f,v in fracs.items(): 
    fracs[f] *= rho

# Initialize our temperature
fracs["T"] = Temperature

s = QuantitiesTable(species_list, fracs)
ud = QuantitiesTable(species_list, derivs)
dd = QuantitiesTable(species_list, derivs)

tf = 1e11 * (365*24*3600)
ti = 0.0 * (365*24*3600)
dt = 1e5 * (365*24*3600)
vals = dict(((i.name, []) for i in species_list))
vals['t'] = []

order_of_update = ["HII","HI","de","HM","H2II","H2I","HeI","HeII","HeIII"]
for rname, r in []: #reaction_table.items():
    if rname not in ["r02"]:
        print "Zeroing", rname
        r.rate.values *= 0.0

nsteps = 0
try:
    while ti < tf:
        nsteps += 1
        # Now we calculate our timestep, using ALL the derivatives
        # Note that this implicitly clears it
        all_calc_derivs(s, ud, dd, reaction_table)
        dt = min( tf-ti, calculate_timestep(s, ud, dd, rho))
        update_vals(vals, s, ti)
        for sname in order_of_update:
            ud.values *= 0.0 + tiny
            dd.values *= 0.0 + tiny
            calc_derivs(s, ud, dd, reaction_table, sname)
            bdf_update(s, ud, dd, dt, sname)
        for c in constraints: c(s)
        ti += dt
        if nsteps % 100 == 0: print ti/tf, vals["H2I"][-1]/rho
except KeyboardInterrupt:
    pass
print "Took %0.3e steps to get to %0.3e seconds" % (nsteps, ti)
update_vals(vals, s, ti)
for v in vals: vals[v] = na.array(vals[v])
#plot_vals(vals, norm=rho)

# Neutral fraction:
c1 = 1.0 / X
fk02 = reaction_rates_table['k02'](vals)
foft = 1.0 - 1.0/(c1 + fk02 * vals['t'] * rho)
import matplotlib;matplotlib.use("Agg");import pylab
pylab.clf()
pylab.loglog(vals['t']/(365*24*3600), vals['HI']/rho, '-k')
pylab.loglog(vals['t']/(365*24*3600), foft, '--b')
pylab.savefig("recombination.png")
#print foft/(vals['HI']/rho)

all_plot_vals(vals, 'values', rho)
