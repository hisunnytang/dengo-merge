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

from jinja2 import Template
import h5py

def create_tables(rate_list, solver_name):
    f = h5py.File("%s_rate_tables.h5" % solver_name, "w")
    for name, rate in rate_list.items():
        f.create_dataset("/%s" % name, data = rate.values.astype("float64"))
    f.close()

def create_table_reader(rate_list, solver_name):
    s = open("simple_cvode_solver/cvode_solver.c.template").read()
    template = Template(s)
    out_s = template.render(rate_table = rate_list, solver_name = solver_name)
    f = open("simple_cvode_solver/%s_cvode_solver.c" % solver_name, "w")
    f.write(out_s)

if __name__ == "__main__":
    from primordial_rates import reaction_rates_table
    create_tables(reaction_rates_table, "primordial")
    create_table_reader(reaction_rates_table, "primordial")
