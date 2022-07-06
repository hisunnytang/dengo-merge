"""
Author: Devin Silvia <devin.silvia@gmail.com>
Affiliation: Michigan State University
Homepage: https://bitbucket.org/MatthewTurk/dengo
License:
  Copyright (C) 2013 Devin Silvia.  All Rights Reserved.

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

from .periodic_table import periodic_table_by_name, periodic_table_by_number
from .reaction_classes import (
    AtomicSpecies,
    chianti_rate,
    ion_cooling_rate,
    ion_photoheating_rate,
    ion_photoionization_rate,
    registry_setup,
)


@registry_setup
def setup_ionization(atom_name, photo_background=None, cooling=True):
    num, weight, pn = periodic_table_by_name[atom_name]

    ions = {}
    for i in range(num + 1):
        ion_state = i + 1
        ions[i] = AtomicSpecies(atom_name, i)

    for i in range(num + 1):
        ion_state = i + 1
        sm1 = ions.get(i - 1, None)
        s = ions[i]
        sp1 = ions.get(i + 1, None)
        chianti_rate(atom_name, sm1, s, sp1)
        if cooling:
            ion_cooling_rate(s, atom_name)

        if ion_state != num + 1 and photo_background is not None:
            ion_photoionization_rate(s, photo_background=photo_background)
            if cooling:
                ion_photoheating_rate(s, photo_background=photo_background)


# Generate all the ion-by-ion rates
# Note: all the ones that that "None" for the background
# don't yet have the appropriate tables to allow for
# photo-terms
# ion_by_ion_rates('H', photo_background='HM12')
# ion_by_ion_rates('He', photo_background='HM12')
# ion_by_ion_rates('C', photo_background=None)
# ion_by_ion_rates('N', photo_background=None)
# ion_by_ion_rates('O', photo_background='HM12')
# ion_by_ion_rates('Ne', photo_background=None)
# ion_by_ion_rates('Mg', photo_background=None)
# ion_by_ion_rates('Si', photo_background=None)
# ion_by_ion_rates('S', photo_background=None)

# Generate all the ion-by-ion cooling rates
# Note: all the ones that that "None" for the background
# don't yet have the appropriate tables to allow for
# photo-terms
# ion_by_ion_cooling('H', 1, 1.00794, photo_background='HM12')
# ion_by_ion_cooling('He', 2, 4.002602, photo_background='HM12')
# ion_by_ion_cooling('C', 6, 12.0107, photo_background=None)
# ion_by_ion_cooling('N', 7, 14.0067, photo_background=None)
# ion_by_ion_cooling('O', 8, 15.9994, photo_background='HM12')
# ion_by_ion_cooling('Ne', 10, 20.1797, photo_background=None)
# ion_by_ion_cooling('Mg', 12, 24.3050, photo_background=None)
# ion_by_ion_cooling('Si', 14, 28.0855, photo_background=None)
# ion_by_ion_cooling('S', 16, 32.065, photo_background=None)
