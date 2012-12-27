"""
Author: Devin Silvia <devin.silvia@gmail.com>
Affiliation: UC Boulder
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

from reaction_classes import Species, ion_cooling_rate, species_registry
import docutils.utils.roman as roman

# Note: the atomic/species properties have to be hard-coded for this
# we may want to come up with a better solution here...
atomicSymbol = 'C'
atomicNumber = 6
atomicWeight = 12
nIons = atomicNumber + 1

for i in range(nIons):
    ion_state = i + 1
    speciesName = "%s%s" %(atomicSymbol, roman.toRoman(ion_state))
    # Check if the species already exists
    # in the species registry, if it does
    # we don't want to create it again
    if (speciesName in species_registry) == False:
        s = Species(speciesName, atomicNumber, atomicWeight, i)
    else:
        s = species_registry[speciesName]
    ion_cooling_rate(s)
