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

from dengo.write_cvode_solver import create_tables, create_cvode_solver
from dengo.primordial_rates import reaction_rates_table, reaction_table, \
    species_table
from dengo.primordial_cooling import cooling_action_table, CVODEPrinter, \
    cooling_rates_table

create_tables(reaction_rates_table, cooling_rates_table, "enzo_be_primordial")
create_cvode_solver(reaction_rates_table, reaction_table, species_table,
                    "enzo_be_primordial",
                    cooling_rates_table,
                    cooling_action_table, CVODEPrinter)

