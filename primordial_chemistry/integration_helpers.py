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
import matplotlib;matplotlib.use("Agg");import pylab
from chemistry_constants import tiny

def euler_update(quantities, up_derivatives, down_derivatives, dt, name):
    # We only update the named species
    species = quantities.get_by_name(name)
    q = quantities[name]
    if species.equilibrium:
        quantities[name] = q * up_derivatives[name]/down_derivatives[name]
        return
    net = (up_derivatives[name] - down_derivatives[name])
    quantities[name] += dt * net

def bdf_update(quantities, up_derivatives, down_derivatives, dt, name):
    species = quantities.get_by_name(name)
    q = quantities[name]
    if species.equilibrium:
        quantities[name] = q * ((max(up_derivatives[name], tiny))
                              / (max(down_derivatives[name], tiny)))
        return
    quantities[name] = ((quantities[name] + dt*up_derivatives[name]) / 
                        (1.0 + dt*down_derivatives[name]/q))

def calc_derivs(quantities, up_derivatives, down_derivatives, reactions,
                name):
    for n, reaction in sorted(reactions.items()):
        if name not in reaction.considered: continue
        reaction(quantities, up_derivatives, down_derivatives)

def all_euler_update(quantities, up_derivatives, down_derivatives, dt):
    for species in quantities.species_list:
        n = species.name
        q = quantities[n]
        if species.equilibrium:
            ud = max(up_derivatives[n], tiny)
            dd = max(down_derivatives[n], tiny)
            quantities[n] = q * ud/dd
            print species.name, quantities[n], q
            continue
        net = (up_derivatives[species.name] - down_derivatives[species.name])
        quantities[species.name] += dt * net

def all_bdf_update(quantities, up_derivatives, down_derivatives, dt):
    for species in quantities.species_list:
        n = species.name
        q = quantities[n]
        if species.equilibrium:
            ud = max(up_derivatives[n], tiny)
            dd = max(down_derivatives[n], tiny)
            quantities[n] = q * ud/dd
            continue
        quantities[n] = ((q + dt*up_derivatives[n]) / 
                         (1.0 + dt*down_derivatives[n]/q))

def all_calc_derivs(quantities, up_derivatives, down_derivatives, reactions):
    up_derivatives.values *= 0.0 + tiny
    down_derivatives.values *= 0.0 + tiny
    for n, reaction in sorted(reactions.items()):
        reaction(quantities, up_derivatives, down_derivatives)

def all_plot_vals(vals, prefix="values", norm = 1.0):
    for v in vals:
        if v == 't': continue
        pylab.clf()
        pylab.loglog(vals['t']/(365*24*3600), vals[v]/norm, '-x')
        pylab.savefig("%s_%s.png" % (prefix, v))
        print "Saving: %s_%s.png" % (prefix, v)

def update_vals(all_vals, new_vals, t):
    all_vals['t'].append(t)
    for species in new_vals._names:
        if species not in all_vals: continue
        all_vals[species].append(new_vals[species])

def calculate_timestep(quantities, up_derivatives, down_derivatives, rho):
    dt = 1e30
    for species in quantities.species_list:
        if quantities[species.name]/rho < 1e-10: continue
        net = (up_derivatives[species.name] - down_derivatives[species.name])
        dt = min(dt, 0.1*abs(quantities[species.name] / net))
    return dt
