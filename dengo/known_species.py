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

from reaction_classes import Species

HI = Species("HI", 1.0, 1.0)
HII = Species("HII", 1.0, 1.0)
HeI = Species("HeI", 2.0, 4.0)
HeII = Species("HeII", 2.0, 4.0, 1.0)
HeIII = Species("HeIII", 2.0, 4.0, 2.0)
de = Species("de", 1.0, 1.0)
HM = Species("HM", 1.0, 1.0, -1.0, equilibrium = False)
H2I = Species("H2I", 1.0, 2.0)
H2II = Species("H2II", 1.0, 2.0, 1.0, equilibrium = False)
ge = Species("ge", 1.0, 1.0, 0.0)
