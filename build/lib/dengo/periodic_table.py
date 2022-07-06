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

import csv, os

periodic_table_by_name = {}
periodic_table_by_number = {}

fn = os.path.join(os.path.dirname(__file__), "periodictabledump.csv")

with open(fn, "r") as csvfile:
    pt_reader = csv.reader(csvfile)
    for row in pt_reader:
        num = int(row[0])
        weight = float(row[1])
        full_name = row[2]
        symbol = row[3]
        periodic_table_by_name[symbol] = (num, weight, full_name)
        periodic_table_by_number[num] = (weight, symbol, full_name)
