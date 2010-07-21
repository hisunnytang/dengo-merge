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

from chemistry_constants import tiny
import jinja2
import h5py

def create_tables(rate_list, solver_name):
    f = h5py.File("%s_rate_tables.h5" % solver_name, "w")
    for name, rate in rate_list.items():
        f.create_dataset("/%s" % name, data = rate.values.astype("float64"))
    f.close()

def create_cvode_solver(rate_list, reaction_table, species, solver_name, ics = ""):
    # Generate inverse tables
    ireaction_table = dict([(r.name, r.reaction_id) for r in reaction_table.values()])
    irate_list = dict([(r.name, r.rate_id) for r in rate_list.values()])
    ispecies = dict([(s.name, s.species_id) for s in species_table.values()])

    env = jinja2.Environment(extensions=['jinja2.ext.loopcontrols'])
    s = open("simple_cvode_solver/cvode_solver.c.template").read()
    template = env.template_class(s)
    out_s = template.render(rate_table = rate_list, 
                            irate_table = irate_list, 
                            reactions = reaction_table,
                            ireactions = ireaction_table,
                            solver_name = solver_name,
                            species = species,
                            ispecies = ispecies,
                            initial_conditions = ics)
    f = open("simple_cvode_solver/%s_cvode_solver.c" % solver_name, "w")
    f.write(out_s)

if __name__ == "__main__":
    from primordial_rates import reaction_rates_table, reaction_table, \
        species_table
    s = ""

    Temperature = 350
    rho = 1.0e10 # total rho in amu/cc
    X   = 0.01 # ionization fraction
    fH2 = 0.10 # ionization fraction

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


    s += "for (i = 0; i < data.ncells ; i++) {\n"
    for sname, val in sorted(fracs.items()):
        # Get the species ID
        i = species_table[sname].species_id
        s += "NV_Ith_S(y, i*offset + %s) = %0.5e; // %s\n" % (
            i, val * rho, sname)
    s += "NV_Ith_S(y, i*offset + %s) = %0.5e; // T\n" % (
            species_table["T"].species_id, Temperature)
    s += "data.rho[i] = %0.5e;\n" % (rho)
    s += "}\n"

    create_tables(reaction_rates_table, "primordial")
    create_cvode_solver(reaction_rates_table, reaction_table, species_table,
                        "primordial", ics = s)
