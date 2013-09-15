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

from chemistry_constants import tevk
import numpy as np
import sympy
from reaction_classes import cooling_action, reaction_registry

dhuge = 1.0e30

# Collisional excitations

# -- ceHI --
@cooling_action("ceHI", "-ceHI * HI * de")
def cool(eq):
    @eq.table
    def ceHI(state):
        vals = 7.5e-19*np.exp(-np.minimum(np.log(dhuge), 118348.0/state.T)) \
                    / (1.0 + np.sqrt(state.T/1.0e5))
        return vals

# -- ceHeI --
@cooling_action("ceHeI", "-ceHeI * HeII * de**2")
def cool(eq):
    @eq.table
    def ceHeI(state):
        vals = 9.1e-27*np.exp(-np.minimum(np.log(dhuge), 13179.0/state.T)) \
                      *state.T**(-0.1687)/(1.0 + np.sqrt(state.T/1.0e5))
        return vals

# -- ceHeIIa --
@cooling_action("ceHeII", "-ceHeII * HeII * de")
def cool(eq):
    @eq.table
    def ceHeII(state):
        vals = 5.54e-17*np.exp(-np.minimum(np.log(dhuge),473638.0/state.T)) \
                       *state.T**(-0.397)/(1.+np.sqrt(state.T/1.0e5))
        return vals

# Collisional ionizations

# -- ciHeIS --
@cooling_action("ciHeIS", "-ciHeIS*HeII*de**2")
def cool(eq):
    @eq.table
    def ciHeIS(state):
        vals = 5.01e-27*(state.T)**(-0.1687)/(1.+np.sqrt(state.T/1.0e5)) \
                * np.exp(-np.minimum(np.log(dhuge),55338.0/state.T))
        return vals

# -- ciHI --
@cooling_action("ciHI", "-ciHI*HI*de")
def cool(eq):
    @eq.table
    def ciHI(state):
        vals = 2.18e-11*reaction_registry['k01'].coeff_fn(state)
        return vals

# -- ciHeI --
@cooling_action("ciHeI", "-ciHeI*HeI*de")
def cool(eq):
    @eq.table
    def ciHeI(state):
        vals = 3.94e-11*reaction_registry['k03'].coeff_fn(state)
        return vals

# -- ciHeII --
@cooling_action("ciHeII", "-ciHeII*HeII*de")
def cool(eq):
    @eq.table
    def ciHeII(state):
        vals = 8.72e-11*reaction_registry['k05'].coeff_fn(state)
        return vals

# Recombinations

# -- reHII --
@cooling_action("reHII", "-reHII * HII * de")
def cool(eq):
    @eq.table
    def reHII(state):
        vals = 8.70e-27*np.sqrt(state.T)*(state.T/1000.0)**(-0.2) \
            / (1.0 + (state.T/1.0e6)**(0.7)) 
        return vals

# -- reHeII1 --
@cooling_action("reHeII1", "-reHeII1 * HeII * de")
def cool(eq):
    @eq.table
    def reHeII1(state):
        vals = 1.55e-26*state.T**0.3647
        return vals

# -- reHeII2 --
@cooling_action("reHeII2", "-reHeII2 * HeII * de")
def cool(eq):
    @eq.table
    def reHeII2(state):
        vals = 1.24e-13*state.T**(-1.5) \
            * np.exp(-np.minimum(np.log(dhuge),470000.0/state.T)) \
            * (1.+0.3*np.exp(-np.minimum(np.log(dhuge),94000.0/state.T))) 
        return vals

# -- reHeIII --
@cooling_action("reHeIII", "-reHeIII * HeIII * de")
def cool(eq):
    @eq.table
    def reHeIII(state):
        vals = 3.48e-26*np.sqrt(state.T)*(state.T/1000.0)**(-0.2) \
            / (1.0 + (state.T/1.0e6)**(0.7))
        return vals

# -- brema --
@cooling_action("brem", "-brem * (HII + HeII + 4.0*HeIII) * de")
def cool(eq):
    @eq.table
    def brem(state):
        vals = 1.43e-27*np.sqrt(state.T) \
             *(1.1+0.34*np.exp(-(5.5-np.log10(state.T))**2/3.0))
        return vals

# Galli & Palla 1999 cooling

# Glover and Abel 2008.  Note that this includes the unused Galli & Palla 1999
# low density table.
#@cooling_action("gloverabel08", "-(H2I*0.5)*gphdl/(1.0+gphdl1/galdl)")
@cooling_action("gloverabel08", "-(H2I)*gphdl/(1.0+gphdl/galdl)")
def cool(eq):
    @eq.table
    def gpldl(state):
        tm = np.minimum(np.maximum(state.T, 13.0), 1e5)
        lt = np.log10(tm)
        # Low density limit from Galli & Palla
        # -- gpldl --
        vals = 10.**(-103.0+97.59*lt-48.05*lt**2+10.80*lt*lt*lt
                       -0.9032*lt*lt*lt*lt)
        return vals

    @eq.table
    def gphdl(state):
        # high density limit from HM79 (typo corrected Aug 30/2007)
        # -- gphdl --
        tm  = np.maximum(state.T, 10.0e0)
        tm  = np.minimum(tm, 1.e4)
        t3 = tm/1000.
        # HDLR is from p31 of HM79.
        HDLR = ((9.5e-22*t3**3.76)/(1.+0.12*t3**2.1)*
                np.exp(-(0.13/t3)**3)+3.e-24*np.exp(-0.51/t3))
        HDLV = (6.7e-19*np.exp(-5.86/t3) + 1.6e-18*np.exp(-11.7/t3))
        vals  = (HDLR + HDLV) 
        return vals

    # Excitation by HI
    # -- gaHI --

    @eq.table
    def gaHI(state):
        # Low density rates from Glover & Abel 2008
        tm  = np.maximum(state.T, 10.0e0)
        tm  = np.minimum(tm, 1.e4)
        lt3 = np.log10(tm / 1.e3)  

        _i1 = (state.T < 100.0)
        _i2 = (state.T < 1000.0)
        # Default value
        vals = 10**(-24.311209e0
             + 4.6450521e0 * lt3
             - 3.7209846e0 * lt3**2
             + 5.9369081e0 * lt3**3
             - 5.5108047e0 * lt3**4
             + 1.5538288e0 * lt3**5)
        vals[_i2] = 10**(-24.311209e0
             + 3.5692468e0 * lt3[_i2]
             - 11.332860e0 * lt3[_i2]**2
             - 27.850082e0 * lt3[_i2]**3
             - 21.328264e0 * lt3[_i2]**4
             - 4.2519023e0 * lt3[_i2]**5)
        vals[_i1] = 10**(-16.818342e0
             + 37.383713e0 * lt3[_i1]
             + 58.145166e0 * lt3[_i1]**2
             + 48.656103e0 * lt3[_i1]**3
             + 20.159831e0 * lt3[_i1]**4
             + 3.8479610e0 * lt3[_i1]**5)
        return vals

    # Excitation by H2
    # -- gaH2 --
    @eq.table
    def gaH2(state):
        tm  = np.maximum(state.T, 10.0e0)
        tm  = np.minimum(tm, 1.e4)
        lt3 = np.log10(tm / 1.e3)  
        vals = 10**(-23.962112e0
             + 2.09433740e0  * lt3
             - 0.77151436e0  * lt3**2
             + 0.43693353e0  * lt3**3
             - 0.14913216e0  * lt3**4
             - 0.033638326e0 * lt3**5)
        return vals

    # Excitation by He
    # -- gaHe --
    @eq.table
    def gaHe(state):
        tm  = np.maximum(state.T, 10.0e0)
        tm  = np.minimum(tm, 1.e4)
        lt3 = np.log10(tm / 1.e3)  
        vals = 10**(-23.689237e0
             + 2.1892372e0  * lt3
             - 0.81520438e0 * lt3**2
             + 0.29036281e0 * lt3**3
             - 0.16596184e0 * lt3**4
             + 0.19191375e0 * lt3**5)
        return vals

    # Excitation by H+
    # -- gaHp --
    @eq.table
    def gaHp(state):
        tm  = np.maximum(state.T, 10.0e0)
        tm  = np.minimum(tm, 1.e4)
        lt3 = np.log10(tm / 1.e3)  
        vals = 10**(-21.716699e0
             + 1.3865783e0   * lt3
             - 0.37915285e0  * lt3**2
             + 0.11453688e0  * lt3**3
             - 0.23214154e0  * lt3**4
             + 0.058538864e0 * lt3**5) 
        return vals

    # Excitation by electrons
    # -- gael --
    @eq.table
    def gael(state):
        tm  = np.maximum(state.T, 10.0e0)
        tm  = np.minimum(tm, 1.e4)
        lt3 = np.log10(tm / 1.e3)  
        _i1 = (state.T < 200)
        vals = 10**(-22.190316
             + 1.5728955  * lt3
             - 0.21335100 * lt3**2
             + 0.96149759 * lt3**3
             - 0.91023195 * lt3**4
             + 0.13749749 * lt3**5)
        vals[_i1] = 10**(-34.286155e0
             - 48.537163e0  * lt3[_i1]
             - 77.121176e0  * lt3[_i1]**2
             - 51.352459e0  * lt3[_i1]**3
             - 15.169160e0  * lt3[_i1]**4
             - 0.98120322e0 * lt3[_i1]**5) 
        return vals

    eq.temporary("galdl", "gaHI*HI + gaH2*H2I + gaHe*HeI + gaHp*HII + gael*de")
    # do we have gphdl right?
    #eq.temporary("gphdl1", "gphdl")

# Compton cooling (Peebles 1971)
# -- comp --
# FIX THIS
#@cooling_action("compton", "-comp1*(T - comp2)*de + comp3")
def cool(eq):
    @eq.table
    def comp(state):
        vals = 5.65e-26 + state.T*0.0
        return vals
    eq.temporary("comp1", "comp * (1.0 * redshift)**4")
    eq.temporary("comp2", "2.73 * (1.0 + redshift)")
    eq.temporary("comp3", "-comp_xray*(T - comp_temp) * de")

#  Photoelectric heating by UV-irradiated dust (Wolfire 1995)
#  with epsilon=0.05, G_0=1.7 (rate in erg s^-1 cm^-3)
# -- gammah --
# FIX THIS
@cooling_action("gammah", "0.0*gammah * HI * HII")
def cool(eq):
    @eq.table
    def gammah(state):
        vals = 8.5e-26 + state.T
        return vals

@cooling_action("h2formation", "h2mheat*HI*HI*HI - h2mcool*H2I*HI")
def h2formation(eq):
    @eq.table
    def h2mheat(state):
        vals = 7.177e-12 * reaction_registry['k22'].coeff_fn(state)
        return vals

    @eq.table
    def h2mcool(state):
        vals = 7.177e-12 * reaction_registry['k13'].coeff_fn(state)
        return vals
