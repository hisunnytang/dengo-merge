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

from reaction_classes import Species, Reaction, reaction_rates_table

from roman_numeral_convert import *

#
# Rates for oxygen collision ionization and recombination
#

# Define atomic number
atomNum = 8

# Define element identifying strings
capitalElement = 'O'
element = 'o'

# Get the last ion number as roman numeral
lastIonNumber = int_to_roman(atomNum+1)

# Set the first rate index for this species
# based off the one that was last assigned
rateIndex = len(reaction_rates_table.keys()) + 2

# store the starting index for the oxygen
# ion rates
oxygenFirstRate = rateIndex

# We'll use ChiantiPy to grab the rates
import chianti.core as ch

# Populate rates table with coll. ionization
# and recombination rates for all oxygen ions
for i in range(atomNum+1):

    # Grab the ion object from Chianti
    ion = ch.ion('%s_%i' %(element,i+1), temperature=T)

    # Get Chianti to compute the collisional ionization rate
    # and then dump the rates into the table
    if ion.Spectroscopic != '%s %s' %(capitalElement,lastIonNumber): # last ion can't be ionized
        ion.ionizRate()
        vals = ion.IonizRate['rate']
        reaction_rates_table['k%02i' %(rateIndex)] = ReactionRate('k%02i' %(rateIndex), vals)
        rateIndex += 1

    # Get Chianti to compute the recombination rate
    # and then dump the rates into the table    
    if ion.Spectroscopic != '%s I' %(capitalElement): # neutral state can't recombine
        ion.recombRate()
        vals = ion.RecombRate['rate']
        reaction_rates_table['k%02i' %(rateIndex)] = ReactionRate('k%02i' %(rateIndex), vals)
        rateIndex += 1

#    reaction_table['r%02i' %(oxygenFirstRate)] = Reaction('k%02i' %(oxygenFirstRate), [   (1,OI),   (1,de)], [  (1,OII),   (2,de)])
#    reaction_table['r%02i' %(oxygenFirstRate+1)] = Reaction('k%02i' %(oxygenFirstRate+1), [  (1,OII),   (1,de)], [ (1,OIII),   (2,de)])
#    reaction_table['r%02i' %(oxygenFirstRate+2)] = Reaction('k%02i' %(oxygenFirstRate+2), [  (1,OII),   (1,de)], [   (1,OI),         ])
#    reaction_table['r%02i' %(oxygenFirstRate+3)] = Reaction('k%02i' %(oxygenFirstRate+3), [ (1,OIII),   (1,de)], [  (1,OIV),   (2,de)])
#    reaction_table['r%02i' %(oxygenFirstRate+4)] = Reaction('k%02i' %(oxygenFirstRate+4), [ (1,OIII),   (1,de)], [  (1,OII),         ])
#    reaction_table['r%02i' %(oxygenFirstRate+5)] = Reaction('k%02i' %(oxygenFirstRate+5), [  (1,OIV),   (1,de)], [   (1,OV),   (2,de)])
#    reaction_table['r%02i' %(oxygenFirstRate+6)] = Reaction('k%02i' %(oxygenFirstRate+6), [  (1,OIV),   (1,de)], [ (1,OIII),         ])
#    reaction_table['r%02i' %(oxygenFirstRate+7)] = Reaction('k%02i' %(oxygenFirstRate+7), [   (1,OV),   (1,de)], [  (1,OVI),   (2,de)])
#    reaction_table['r%02i' %(oxygenFirstRate+8)] = Reaction('k%02i' %(oxygenFirstRate+8), [   (1,OV),   (1,de)], [  (1,OIV),         ])
#    reaction_table['r%02i' %(oxygenFirstRate+9)] = Reaction('k%02i' %(oxygenFirstRate+9), [  (1,OVI),   (1,de)], [ (1,OVII),   (2,de)])
#    reaction_table['r%02i' %(oxygenFirstRate+10)] = Reaction('k%02i' %(oxygenFirstRate+10), [  (1,OVI),   (1,de)], [   (1,OV),         ])
#    reaction_table['r%02i' %(oxygenFirstRate+11)] = Reaction('k%02i' %(oxygenFirstRate+11), [ (1,OVII),   (1,de)], [ (1,OVII),   (2,de)])
#    reaction_table['r%02i' %(oxygenFirstRate+12)] = Reaction('k%02i' %(oxygenFirstRate+12), [ (1,OVII),   (1,de)], [  (1,OVI),         ])
#    reaction_table['r%02i' %(oxygenFirstRate+13)] = Reaction('k%02i' %(oxygenFirstRate+13), [(1,OVIII),   (1,de)], [  (1,OIX),   (2,de)])
#    reaction_table['r%02i' %(oxygenFirstRate+14)] = Reaction('k%02i' %(oxygenFirstRate+14), [(1,OVIII),   (1,de)], [ (1,OVII),         ])
#    reaction_table['r%02i' %(oxygenFirstRate+15)] = Reaction('k%02i' %(oxygenFirstRate+15), [  (1,OIX),   (1,de)], [(1,OVIII),         ])
#
