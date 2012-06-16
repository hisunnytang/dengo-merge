"""
Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: UCSD
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2010 Matthew Turk.  All Rights Reserved.

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

import numpy as na

from chemistry_constants import tevk, tiny, mh
from reaction_classes import ReactionRate, Species, Reaction, reaction_rates_table

# First we set up our constants
ReactionRate.init_temperature((1.0, 1e7))
T = ReactionRate.T
logT = ReactionRate.logT
tev= ReactionRate.tev
logtev = ReactionRate.logtev

# Now we actually set up our values...

# -- k01 --
vals = na.exp(-32.71396786375 
             + 13.53655609057*logtev
             - 5.739328757388*logtev**2 
             + 1.563154982022*logtev**3
             - 0.2877056004391*logtev**4
             + 0.03482559773736999*logtev**5
             - 0.00263197617559*logtev**6
             + 0.0001119543953861*logtev**7
             - 2.039149852002e-6*logtev**8)
reaction_rates_table['k01'] = ReactionRate('k01', vals)

# -- k02 --

_i1 = (T > 5500)
_i2 = ~_i1
vals = na.exp(-28.61303380689232
              - 0.7241125657826851*logtev
              - 0.02026044731984691*logtev**2
              - 0.002380861877349834*logtev**3
              - 0.0003212605213188796*logtev**4
              - 0.00001421502914054107*logtev**5
              + 4.989108920299513e-6*logtev**6
              + 5.755614137575758e-7*logtev**7
              - 1.856767039775261e-8*logtev**8
              - 3.071135243196595e-9*logtev**9)
_i1 = (tev < 0.8)
vals[_i2] = (1.54e-9*(1.+0.3/na.exp(8.099328789667/tev))
              / (na.exp(40.49664394833662/tev)*tev**1.5)
              + 3.92e-13/tev**0.6353) 
vals[_i1] = 3.92e-13/tev**0.6353
reaction_rates_table['k02'] = ReactionRate('k02', vals)

# -- k03 --

_i1 = (tev > 0.8)
_i2 = ~_i1
vals = na.exp(-44.09864886561001
              + 23.91596563469*logtev
              - 10.75323019821*logtev**2
              + 3.058038757198*logtev**3
              - 0.5685118909884001*logtev**4
              + 0.06795391233790001*logtev**5
              - 0.005009056101857001*logtev**6
              + 0.0002067236157507*logtev**7
              - 3.649161410833e-6*logtev**8)
vals[_i2] = tiny
reaction_rates_table['k03'] = ReactionRate('k03', vals)

# -- k04 --

_i1 = (tev > 0.8)
_i2 = ~_i1
vals = (1.54e-9*(1.+0.3/na.exp(8.099328789667/tev))
               / (na.exp(40.49664394833662/tev)*tev**1.5)
               + 3.92e-13/tev**0.6353) 
vals[_i2] = tiny
reaction_rates_table['k04'] = ReactionRate('k04', vals)

# -- k05 --
_i1 = (tev > 0.8)
_i2 = ~_i1
vals = na.exp(-68.71040990212001
              + 43.93347632635*logtev
              - 18.48066993568*logtev**2
              + 4.701626486759002*logtev**3
              - 0.7692466334492*logtev**4
              + 0.08113042097303*logtev**5
              - 0.005324020628287001*logtev**6
              + 0.0001975705312221*logtev**7
              - 3.165581065665e-6*logtev**8)
vals[_i2] = tiny
reaction_rates_table['k05'] = ReactionRate(
    'k05', vals)

# -- k06 --
vals = 3.36e-10/na.sqrt(T)/(T/1.e3)**0.2/(1+(T/1.e6)**0.7)
reaction_rates_table['k06'] = ReactionRate('k06', vals)

# -- k07 --
vals = 6.77e-15*tev**0.8779
reaction_rates_table['k07'] = ReactionRate('k07', vals)

# -- k08 --
_i1 = (tev > 0.1)
_i2 = ~_i1
vals = na.exp(-20.06913897587003
              + 0.2289800603272916*logtev
              + 0.03599837721023835*logtev**2
              - 0.004555120027032095*logtev**3
              - 0.0003105115447124016*logtev**4
              + 0.0001073294010367247*logtev**5
              - 8.36671960467864e-6*logtev**6
              + 2.238306228891639e-7*logtev**7)
vals[_i2] = 1.43e-9
reaction_rates_table['k08'] = ReactionRate('k08', vals)

# -- k09 --

_i1 = (T > 6.7e3)
vals = 1.85e-23*T**1.8
vals[_i1] = 5.81e-16*(T/56200)**(-0.6657*na.log10(T/56200))
reaction_rates_table['k09'] = ReactionRate('k09', vals)

# -- k10 --
vals = T * 0.0 + 6.0e-10
reaction_rates_table['k10'] = ReactionRate('k10', vals)

# -- k11 --
_i1 = (tev > 0.3)
_i2 = ~_i1
vals = na.exp(-24.24914687731536
              + 3.400824447095291*logtev
              - 3.898003964650152*logtev**2
              + 2.045587822403071*logtev**3
              - 0.5416182856220388*logtev**4
              + 0.0841077503763412*logtev**5
              - 0.007879026154483455*logtev**6
              + 0.0004138398421504563*logtev**7
              - 9.36345888928611e-6*logtev**8)
vals[_i2] = tiny
reaction_rates_table['k11'] = ReactionRate('k11', vals)

# -- k12 --
_i1 = (tev > 0.3)
_i2 = ~_i1
vals = 5.6e-11*na.exp(-102124/T)*T**0.5
reaction_rates_table['k12'] = ReactionRate('k12', vals)

# -- k13 --
# NOTE: This is the Glover 2008 rate
vals = 10.0**(-178.4239 - 68.42243 * na.log10(T) 
                        + 43.20243 * na.log10(T)**2
                        - 4.633167 * na.log10(T)**3 
                        + 69.70086 * na.log10(1 + 40870.38 / T)
                        - (23705.7 / T))
reaction_rates_table['k13'] = ReactionRate('k13', vals)

# -- k14 --
_i1 = (tev > 0.04)
_i2 = ~_i1
vals = na.exp(-18.01849334273
              + 2.360852208681*logtev
              - 0.2827443061704*logtev**2
              + 0.01623316639567*logtev**3
              - 0.03365012031362999*logtev**4
              + 0.01178329782711*logtev**5
              - 0.001656194699504*logtev**6
              + 0.0001068275202678*logtev**7
              - 2.631285809207e-6*logtev**8)
vals[_i2] = tiny
reaction_rates_table['k14'] = ReactionRate('k14', vals)

# -- k15 --
_i1 = (tev > 0.1)
_i2 = ~_i1
vals = na.exp(-20.37260896533324
              + 1.139449335841631*logtev
              - 0.1421013521554148*logtev**2
              + 0.00846445538663*logtev**3
              - 0.0014327641212992*logtev**4
              + 0.0002012250284791*logtev**5
              + 0.0000866396324309*logtev**6
              - 0.00002585009680264*logtev**7
              + 2.4555011970392e-6*logtev**8
              - 8.06838246118e-8*logtev**9) 
vals[_i2] = 2.56e-9*tev**1.78186
reaction_rates_table['k15'] = ReactionRate('k15', vals)

# -- k16 --
k16 = 6.5e-9/na.sqrt(tev)
reaction_rates_table['k16'] = ReactionRate('k16', vals)

# -- k17 --
_i1 = (T < 1e4)
_i2 = ~_i1
vals = 1.0e-8*T**(-0.4)
vals[_i2] = 4.0e-4*T**(-1.4)*na.exp(-15100.0/T)
reaction_rates_table['k17'] = ReactionRate('k17', vals)

# -- k18 --
_i1 = (T > 617)
_i2 = ~_i1
vals = 1.32e-6 * T**(-0.76)
vals[_i2] = 1.e-8 
reaction_rates_table['k18'] = ReactionRate('k18', vals)
      
# -- k19 --
vals = 5.e-7*na.sqrt(100./T)
reaction_rates_table['k19'] = ReactionRate('k19', vals)

# -- k21 --
vals = 2.8e-31 * (T**(-0.6))
reaction_rates_table['k21'] = ReactionRate('k21', vals)

# -- k22 --
# NOTE: This is the Glover 2008 rate
vals = 7.7e-31 / T**0.464
reaction_rates_table['k22'] = ReactionRate('k22', vals)

# -- k23 --
vals = ((8.125e-8 / na.sqrt(T))
      * na.exp(-52000/T)
      * (1.0 - na.exp(-6000/T)))
reaction_rates_table['k23'] = ReactionRate('k23', vals)

del vals, _i1, _i2

#
# Rates for oxygen collision ionization and recombination
#

# We'll try using ChiantiPy to grab the rates
import chianti.core as ch

# -- k30 -- OV collisional ionization
o5 = ch.ion('o_5', temperature=T)
o5.ionizRate()
vals = o5.IonizRate['rate']
reaction_rates_table['k30'] = ReactionRate('k30', vals)


#
# Now we create a number of species tables
#

# This is the table of Species we know about
species_table = dict(
    HI = Species("HI", 1.0),
    HII = Species("HII", 1.0, 1.0),
    HeI = Species("HeI", 4.0),
    HeII = Species("HeII", 4.0, 1.0),
    HeIII = Species("HeIII", 4.0, 2.0),
    de = Species("de", 1.0),
    HM = Species("HM", 1.0, -1.0, equilibrium = False),
    H2I = Species("H2I", 2.0),
    H2II = Species("H2II", 2.0, 1.0, equilibrium = False),
    OV = Species("OV", 16.0, 4.0),
    OVI = Species("OVI", 16.0, 5.0),
    ge = Species("ge", 0.0, 0.0),
)

locals().update(species_table)

# This is the full set of Reactions we know about.
reaction_table = dict(
    r01 = Reaction('k01', [   (1,HI),   (1,de)], [  (1,HII),   (2,de)]),
    r02 = Reaction('k02', [  (1,HII),   (1,de)], [   (1,HI),         ]),
    r03 = Reaction('k03', [  (1,HeI),   (1,de)], [ (1,HeII),   (2,de)]),
    r04 = Reaction('k04', [ (1,HeII),   (1,de)], [  (1,HeI),         ]),
    r05 = Reaction('k05', [ (1,HeII),   (1,de)], [(1,HeIII),   (2,de)]),
    r06 = Reaction('k06', [(1,HeIII),   (1,de)], [ (1,HeII),         ]),
    r07 = Reaction('k07', [   (1,HI),   (1,de)], [   (1,HM),         ]),
    r08 = Reaction('k08', [   (1,HM),   (1,HI)], [  (1,H2I),   (1,de)]),
    r09 = Reaction('k09', [   (1,HI),  (1,HII)], [ (1,H2II),         ]),
    r10 = Reaction('k10', [ (1,H2II),   (1,HI)], [  (1,H2I),  (1,HII)]),

    r11 = Reaction('k11', [  (1,H2I),  (1,HII)], [  (1,H2II),  (1,HI)]),
    r12 = Reaction('k12', [  (1,H2I),   (1,de)], [  (2,HII),   (1,de)]),
    r13 = Reaction('k13', [  (1,H2I),   (1,HI)], [   (3,HI),         ]), #3b
    r14 = Reaction('k14', [   (1,HM),   (1,de)], [   (1,HI),   (2,de)]),
    r15 = Reaction('k15', [   (1,HM),   (1,HI)], [   (2,HI),   (1,de)]),
    r16 = Reaction('k16', [   (1,HM),  (1,HII)], [   (2,HI),         ]),
    r17 = Reaction('k17', [   (1,HM),  (1,HII)], [ (1,H2II),   (1,de)]),
    r18 = Reaction('k18', [ (1,H2II),   (1,de)], [   (2,HI),         ]),
    r19 = Reaction('k19', [ (1,H2II),   (1,HM)], [   (1,HI),  (1,H2I)]),
    r21 = Reaction('k21', [   (2,HI),  (1,H2I)], [  (2,H2I),         ]), #3b
    r22 = Reaction('k22', [   (2,HI),   (1,HI)], [  (1,H2I),   (1,HI)]), #3b
    r23 = Reaction('k23', [  (1,H2I),  (1,H2I)], [   (2,HI),  (1,H2I)]), #3b

    r30 = Reaction('k30', [   (1,OV),   (1,de)], [  (1,OVI),   (2,de)]),
)

locals().update(reaction_table)
