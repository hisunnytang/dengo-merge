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

import numpy as np
from chemistry_constants import tevk, tiny, mh
from reaction_classes import reaction
from .known_species import *

# -- k01 --
@reaction('k01', [   (1,HI),   (1,de)], [  (1,HII),   (2,de)])
def rxn(state):
    vals = np.exp(-32.71396786375 
                 + 13.53655609057*state.logtev
                 - 5.739328757388*state.logtev**2 
                 + 1.563154982022*state.logtev**3
                 - 0.2877056004391*state.logtev**4
                 + 0.03482559773736999*state.logtev**5
                 - 0.00263197617559*state.logtev**6
                 + 0.0001119543953861*state.logtev**7
                 - 2.039149852002e-6*state.logtev**8)
    return vals

# -- k02 --

@reaction('k02', [  (1,HII),   (1,de)], [   (1,HI),         ])
def rxn(state):
    _i1 = (state.T > 5500)
    _i2 = ~_i1
    vals = np.exp(-28.61303380689232
                  - 0.7241125657826851*state.logtev
                  - 0.02026044731984691*state.logtev**2
                  - 0.002380861877349834*state.logtev**3
                  - 0.0003212605213188796*state.logtev**4
                  - 0.00001421502914054107*state.logtev**5
                  + 4.989108920299513e-6*state.logtev**6
                  + 5.755614137575758e-7*state.logtev**7
                  - 1.856767039775261e-8*state.logtev**8
                  - 3.071135243196595e-9*state.logtev**9)
    _i1 = (state.tev < 0.8)
    vals[_i2] = (1.54e-9*(1.+0.3/np.exp(8.099328789667/state.tev[_i2]))
                  / (np.exp(40.49664394833662/state.tev[_i2])*state.tev[_i2]**1.5)
                  + 3.92e-13/state.tev[_i2]**0.6353) 
    vals[_i1] = 3.92e-13/state.tev[_i1]**0.6353
    return vals

# -- k03 --
@reaction('k03', [  (1,HeI),   (1,de)], [ (1,HeII),   (2,de)])
def rxn(state):
    _i1 = (state.tev > 0.8)
    _i2 = ~_i1
    vals = np.exp(-44.09864886561001
                  + 23.91596563469*state.logtev
                  - 10.75323019821*state.logtev**2
                  + 3.058038757198*state.logtev**3
                  - 0.5685118909884001*state.logtev**4
                  + 0.06795391233790001*state.logtev**5
                  - 0.005009056101857001*state.logtev**6
                  + 0.0002067236157507*state.logtev**7
                  - 3.649161410833e-6*state.logtev**8)
    vals[_i2] = tiny
    return vals

# -- k04 --

@reaction('k04', [ (1,HeII),   (1,de)], [  (1,HeI),         ])
def rxn(state):
    _i1 = (state.tev > 0.8)
    _i2 = ~_i1
    vals = (1.54e-9*(1.+0.3/np.exp(8.099328789667/state.tev))
                   / (np.exp(40.49664394833662/state.tev)*state.tev**1.5)
                   + 3.92e-13/state.tev**0.6353) 
    vals[_i2] = tiny
    return vals

# -- k05 --
@reaction('k05', [ (1,HeII),   (1,de)], [(1,HeIII),   (2,de)])
def rxn(state):
    _i1 = (state.tev > 0.8)
    _i2 = ~_i1
    vals = np.exp(-68.71040990212001
                  + 43.93347632635*state.logtev
                  - 18.48066993568*state.logtev**2
                  + 4.701626486759002*state.logtev**3
                  - 0.7692466334492*state.logtev**4
                  + 0.08113042097303*state.logtev**5
                  - 0.005324020628287001*state.logtev**6
                  + 0.0001975705312221*state.logtev**7
                  - 3.165581065665e-6*state.logtev**8)
    vals[_i2] = tiny
    return vals

# -- k06 --
@reaction('k06', [(1,HeIII),   (1,de)], [ (1,HeII),         ])
def rxn(state):
    vals = (3.36e-10/np.sqrt(state.T)/(state.T/1.e3)**0.2 /
            (1+(state.T/1.e6)**0.7))
    return vals

# -- k07 --
@reaction('k07', [   (1,HI),   (1,de)], [   (1,HM),         ])
def rxn(state):
    vals = 6.77e-15*state.tev**0.8779
    return vals

# -- k08 --
@reaction('k08', [   (1,HM),   (1,HI)], [  (1,H2I),   (1,de)])
def rxn(state):
    _i1 = (state.tev > 0.1)
    _i2 = ~_i1
    vals = np.exp(-20.06913897587003
                  + 0.2289800603272916*state.logtev
                  + 0.03599837721023835*state.logtev**2
                  - 0.004555120027032095*state.logtev**3
                  - 0.0003105115447124016*state.logtev**4
                  + 0.0001073294010367247*state.logtev**5
                  - 8.36671960467864e-6*state.logtev**6
                  + 2.238306228891639e-7*state.logtev**7)
    vals[_i2] = 1.43e-9
    return vals

# -- k09 --

@reaction('k09', [   (1,HI),  (1,HII)], [ (1,H2II),         ])
def rxn(state):
    _i1 = (state.T > 6.7e3)
    vals = 1.85e-23*state.T**1.8
    vals[_i1] = 5.81e-16*(state.T[_i1]/56200)**(-0.6657*np.log10(state.T[_i1]/56200))
    return vals

# -- k10 --
@reaction('k10', [ (1,H2II),   (1,HI)], [  (1,H2I),  (1,HII)])
def rxn(state):
    vals = state.T * 0.0 + 6.0e-10
    return vals

# -- k11 --
@reaction('k11', [  (1,H2I),  (1,HII)], [  (1,H2II),  (1,HI)])
def rxn(state):
    _i1 = (state.tev > 0.3)
    _i2 = ~_i1
    vals = np.exp(-24.24914687731536
                  + 3.400824447095291*state.logtev
                  - 3.898003964650152*state.logtev**2
                  + 2.045587822403071*state.logtev**3
                  - 0.5416182856220388*state.logtev**4
                  + 0.0841077503763412*state.logtev**5
                  - 0.007879026154483455*state.logtev**6
                  + 0.0004138398421504563*state.logtev**7
                  - 9.36345888928611e-6*state.logtev**8)
    vals[_i2] = tiny
    return vals

# -- k12 --
@reaction('k12', [  (1,H2I),   (1,de)], [  (2,HII),   (1,de)])
def rxn(state):
    _i1 = (state.tev > 0.3)
    _i2 = ~_i1
    vals = 5.6e-11*np.exp(-102124/state.T)*state.T**0.5
    return vals

# -- k13 --
# NOTE: This is the Glover 2008 rate
@reaction('k13', [  (1,H2I),   (1,HI)], [   (3,HI),         ])
def rxn(state):
    vals = 10.0**(-178.4239 - 68.42243 * np.log10(state.T) 
                            + 43.20243 * np.log10(state.T)**2
                            - 4.633167 * np.log10(state.T)**3 
                            + 69.70086 * np.log10(1 + 40870.38 / state.T)
                            - (23705.7 / state.T))
    return vals

# -- k14 --
@reaction('k14', [   (1,HM),   (1,de)], [   (1,HI),   (2,de)])
def rxn(state):
    _i1 = (state.tev > 0.04)
    _i2 = ~_i1
    vals = np.exp(-18.01849334273
                  + 2.360852208681*state.logtev
                  - 0.2827443061704*state.logtev**2
                  + 0.01623316639567*state.logtev**3
                  - 0.03365012031362999*state.logtev**4
                  + 0.01178329782711*state.logtev**5
                  - 0.001656194699504*state.logtev**6
                  + 0.0001068275202678*state.logtev**7
                  - 2.631285809207e-6*state.logtev**8)
    vals[_i2] = tiny
    return vals

# -- k15 --
@reaction('k15', [   (1,HM),   (1,HI)], [   (2,HI),   (1,de)])
def rxn(state):
    _i1 = (state.tev > 0.1)
    _i2 = ~_i1
    vals = np.exp(-20.37260896533324
                  + 1.139449335841631*state.logtev
                  - 0.1421013521554148*state.logtev**2
                  + 0.00846445538663*state.logtev**3
                  - 0.0014327641212992*state.logtev**4
                  + 0.0002012250284791*state.logtev**5
                  + 0.0000866396324309*state.logtev**6
                  - 0.00002585009680264*state.logtev**7
                  + 2.4555011970392e-6*state.logtev**8
                  - 8.06838246118e-8*state.logtev**9) 
    vals[_i2] = 2.56e-9*state.tev[_i2]**1.78186
    return vals

# -- k16 --
@reaction('k16', [   (1,HM),  (1,HII)], [   (2,HI),         ])
def rxn(state):
    k16 = 6.5e-9/np.sqrt(state.tev)
    return k16

# -- k17 --
@reaction('k17', [   (1,HM),  (1,HII)], [ (1,H2II),   (1,de)])
def rxn(state):
    _i1 = (state.T < 1e4)
    _i2 = ~_i1
    vals = 1.0e-8*state.T**(-0.4)
    vals[_i2] = 4.0e-4*state.T[_i2]**(-1.4)*np.exp(-15100.0/state.T[_i2])
    return vals

# -- k18 --
@reaction('k18', [ (1,H2II),   (1,de)], [   (2,HI),         ])
def rxn(state):
    _i1 = (state.T > 617)
    _i2 = ~_i1
    vals = 1.32e-6 * state.T**(-0.76)
    vals[_i2] = 1.e-8 
    return vals
      
# -- k19 --
@reaction('k19', [ (1,H2II),   (1,HM)], [   (1,HI),  (1,H2I)])
def rxn(state):
    vals = 5.e-7*np.sqrt(100./state.T)
    return vals

# -- k21 --
@reaction('k21', [   (2,HI),  (1,H2I)], [  (2,H2I),         ])
def rxn(state):
    vals = 2.8e-31 * (state.T**(-0.6))
    return vals

# -- k22 --
# NOTE: This is the Glover 2008 rate
@reaction('k22', [   (2,HI),   (1,HI)], [  (1,H2I),   (1,HI)])
def rxn(state):
    vals = 7.7e-31 / state.T**0.464
    return vals

# -- k23 --
@reaction('k23', [  (1,H2I),  (1,H2I)], [   (2,HI),  (1,H2I)])
def rxn(state):
    vals = ((8.125e-8 / np.sqrt(state.T))
          * np.exp(-52000/state.T)
          * (1.0 - np.exp(-6000/state.T)))
    return vals

