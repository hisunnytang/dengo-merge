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
import numpy as na
import sympy
from reaction_classes import reaction_registry, species_registry
import primordial_rates


species_symbols = dict( [
    (sname, sympy.Symbol(sname)) for
     sname in species_registry ])

known_variables = dict( 
    T = sympy.Symbol("T"),
)

user_data_symbols = dict(
    redshift = sympy.Symbol("redshift"),
)

user_data_cell_vars = dict(
)

other_known_symbols = {}
other_known_symbols.update(known_variables)
other_known_symbols.update(user_data_symbols)
other_known_symbols.update(user_data_cell_vars)

# We set up all the symbols we are aware of here.

class CVODEPrinter(sympy.printing.str.StrPrinter):
    def __init__(self, template_vars, *args, **kwargs):
        sympy.printing.str.StrPrinter.__init__(self, *args, **kwargs)
        self.template_vars = template_vars

    def _print_Symbol(self, symbol):
        s = str(symbol)
        if s in species_symbols:
            vname = self.template_vars["species_varnames"][s]
            return "%s" % vname
        elif s in cooling_rates_table:
            cid = self.template_vars["cooling_rate_ids"][s]
            vname = "data->cooling_storage[%s][cell]" % cid
            return "%s" % vname
        elif s in known_variables:
            return s
        elif s in user_data_symbols:
            return "data->%s" % s
        elif s in user_data_cell_vars:
            return "data->%s[cell]" % s
        else:
            raise RuntimeError

    def _print_Derivative(self, expr):
        # Using the sympy printing docs as a hint here!
        rate = str(expr.args[0].func)
        if rate not in cooling_rates_table:
            raise RuntimeError
        cid = self.template_vars["cooling_rate_ids"][rate]
        vname = "cooling_slopes[%s]" % cid
        return "%s" % vname

    def _print_Pow(self, expr):
        PREC = sympy.printing.precedence.precedence(expr)
        return 'pow(%s,%s)'%(self.parenthesize(expr.base, PREC),
                             self.parenthesize(expr.exp, PREC))

# Note that we actually re-implement here rather than subclass.
# That's because we want our cooling rates to be able to have different min/max
# temperature, and they may also behave differently in other ways from normal
# reaction rates.

class CoolingRate(object):
    def __init__(self, name, values):
        self.values = values
        self.name = name
        self.s = sympy.Symbol(self.name)
        self.f = sympy.Function(self.name)(known_variables['T'])

    def __call__(self, quantities):
        T = quantities['T']
        return na.interp(T, self.T, self.values)

    @classmethod
    def init_temperature(cls, T_bounds, n_bins=1024):
        cls.n_bins = 1024
        cls.T = na.logspace(
            na.log(T_bounds[0]), na.log(T_bounds[1]), n_bins,
                   base = na.e)
        cls.logT = na.log(cls.T)
        cls.tev = cls.T / tevk
        cls.logtev = na.log(cls.tev)
        cls.T_bounds = T_bounds

cooling_rates_table = dict()

class CoolingAction(object):
    def __init__(self, rates, species, equation,
                 temp_vars = None, other_symbols = None):
        if other_symbols is None: other_symbols = []
        # Equation should be a string
        self.rates = rates
        self.species = species
        self.other_symbols = other_symbols
        symbols = dict( [(a, cooling_rates_table[a].s) for a in rates]
                      + [(a, species_symbols[a]) for a in species]
                      + [(a, other_known_symbols[a]) for a in other_symbols])
        if temp_vars is not None:
            for t, eq in temp_vars:
                symbols[t] = eval(eq, {}, symbols)
        self.equation = eval(equation, {}, symbols)

        # But we also want one that has the functional form of the cooling
        # rates
        symbols = dict( [(a, cooling_rates_table[a].f) for a in rates]
                      + [(a, species_symbols[a]) for a in species]
                      + [(a, other_known_symbols[a]) for a in other_symbols])
        if temp_vars is not None:
            for t, eq in temp_vars:
                symbols[t] = eval(eq, {}, symbols)
        self.dequation = eval(equation, {}, symbols)

cooling_action_table = dict()

# First we set up our constants
CoolingRate.init_temperature((1.0, 1e7))
T = CoolingRate.T
logT = CoolingRate.logT
tev= CoolingRate.tev
logtev = CoolingRate.logtev

dhuge = 1.0e30

qtable = dict(T = T)

# -- ceHI --
vals = 7.5e-19*na.exp(-na.minimum(na.log(dhuge), 118348.0/T)) \
              / (1.0 + na.sqrt(T/1.0e5))
cooling_rates_table['ceHI'] = CoolingRate('ceHI', vals)

# -- ceHeI --
vals = 9.1e-27*na.exp(-na.minimum(na.log(dhuge), 13179.0/T)) \
              *T**(-0.1687)/(1.0 + na.sqrt(T/1.0e5))
cooling_rates_table['ceHeI'] = CoolingRate('ceHeI', vals)

# -- ceHeIIa --
vals = 5.54e-17*na.exp(-na.minimum(na.log(dhuge),473638.0/T)) \
               *T**(-0.397)/(1.+na.sqrt(T/1.0e5))
cooling_rates_table['ceHeII'] = CoolingRate('ceHeII', vals)

# -- ciHeIS --
val = 5.01e-27*(T)**(-0.1687)/(1.+na.sqrt(T/1.0e5)) \
    * na.exp(-na.minimum(na.log(dhuge),55338.0/T))
cooling_rates_table['ciHeIS'] = CoolingRate('ciHeIS', vals)

# -- ciHI --
vals = 2.18e-11*reaction_registry['k01'].coeff_fn(CoolingRate)
cooling_rates_table['ciHI'] = CoolingRate('ciHI', vals)

# -- ciHeI --
vals = 3.94e-11*reaction_registry['k03'].coeff_fn(CoolingRate)
cooling_rates_table['ciHeI'] = CoolingRate('ciHeI', vals)

# -- ciHeII --
vals = 8.72e-11*reaction_registry['k05'].coeff_fn(CoolingRate)
cooling_rates_table['ciHeII'] = CoolingRate('ciHeII', vals)

# -- reHII --
vals = 8.70e-27*na.sqrt(T)*(T/1000.0)**(-0.2) \
     / (1.0 + (T/1.0e6)**(0.7)) 
cooling_rates_table['reHII'] = CoolingRate('reHII', vals)

# -- reHeII1 --
vals = 1.55e-26*T**0.3647
cooling_rates_table['reHeII1'] = CoolingRate('reHeII1', vals)

# -- reHeII2 --
vals = 1.24e-13*T**(-1.5) \
     * na.exp(-na.minimum(na.log(dhuge),470000.0/T)) \
     * (1.+0.3*na.exp(-na.minimum(na.log(dhuge),94000.0/T))) 
cooling_rates_table['reHeII2'] = CoolingRate('reHeII2', vals)

# -- reHeIII --
vals = 3.48e-26*na.sqrt(T)*(T/1000.0)**(-0.2) \
     / (1.0 + (T/1.0e6)**(0.7))
cooling_rates_table['reHeIII'] = CoolingRate('reHeIII', vals)

# -- brema --
vals = 1.43e-27*na.sqrt(T) \
     *(1.1+0.34*na.exp(-(5.5-na.log10(T))**2/3.0))
cooling_rates_table['brem'] = CoolingRate('brem', vals)

# Galli & Palla 1999 cooling
tm = na.minimum(na.maximum(T, 13.0), 1e5)
lt = na.log10(tm)

# Low density limit from Galli & Palla
# -- gpldl --
vals = 10.**(-103.0+97.59*lt-48.05*lt**2+10.80*lt*lt*lt
               -0.9032*lt*lt*lt*lt)
cooling_rates_table['gpldl'] = CoolingRate('gpldl', vals)

# high density limit from HM79 (typo corrected Aug 30/2007)
# -- gphdl --
t3 = tm/1000.
# HDLR is from p31 of HM79.
HDLR = ((9.5e-22*t3**3.76)/(1.+0.12*t3**2.1)*
        na.exp(-(0.13/t3)**3)+3.e-24*na.exp(-0.51/t3))
HDLV = (6.7e-19*na.exp(-5.86/t3) + 1.6e-18*na.exp(-11.7/t3))
vals  = (HDLR + HDLV) 
cooling_rates_table['gphdl'] = CoolingRate('gphdl', vals)

# Low density rates from Glover & Abel 2008
tm  = na.maximum(T, 10.0e0)
tm  = na.minimum(tm, 1.e4)
lt3 = na.log10(tm / 1.e3)  

# Excitation by HI
# -- gaHI --
# Apply these in reverse order, so that one doesn't overwrite the other
_i1 = (T < 100.0)
_i2 = (T < 1000.0)
# Default value
vals = 10**(-24.311209e0
     + 4.6450521e0 * lt3
     - 3.7209846e0 * lt3**2
     + 5.9369081e0 * lt3**3
     - 5.5108047e0 * lt3**4
     + 1.5538288e0 * lt3**5)
vals[_i2] = 10**(-24.311209e0
     + 3.5692468e0 * lt3
     - 11.332860e0 * lt3**2
     - 27.850082e0 * lt3**3
     - 21.328264e0 * lt3**4
     - 4.2519023e0 * lt3**5)
vals[_i1] = 10**(-16.818342e0
     + 37.383713e0 * lt3
     + 58.145166e0 * lt3**2
     + 48.656103e0 * lt3**3
     + 20.159831e0 * lt3**4
     + 3.8479610e0 * lt3**5)
cooling_rates_table['gaHI'] = CoolingRate('gaHI', vals)

# Excitation by H2
# -- gaH2 --
vals = 10**(-23.962112e0
     + 2.09433740e0  * lt3
     - 0.77151436e0  * lt3**2
     + 0.43693353e0  * lt3**3
     - 0.14913216e0  * lt3**4
     - 0.033638326e0 * lt3**5)
cooling_rates_table['gaH2'] = CoolingRate('gaH2', vals)

# Excitation by He
# -- gaHe --
vals = 10**(-23.689237e0
     + 2.1892372e0  * lt3
     - 0.81520438e0 * lt3**2
     + 0.29036281e0 * lt3**3
     - 0.16596184e0 * lt3**4
     + 0.19191375e0 * lt3**5)
cooling_rates_table['gaHe'] = CoolingRate('gaHe', vals)

# Excitation by H+
# -- gaHp --
vals = 10**(-21.716699e0
     + 1.3865783e0   * lt3
     - 0.37915285e0  * lt3**2
     + 0.11453688e0  * lt3**3
     - 0.23214154e0  * lt3**4
     + 0.058538864e0 * lt3**5) 
cooling_rates_table['gaHp'] = CoolingRate('gaHp', vals)

# Excitation by electrons
# -- gael --
_i1 = (T < 200)
vals = 10**(-22.190316
     + 1.5728955  * lt3
     - 0.21335100 * lt3**2
     + 0.96149759 * lt3**3
     - 0.91023195 * lt3**4
     + 0.13749749 * lt3**5)
vals[_i1] = 10**(-34.286155e0
     - 48.537163e0  * lt3
     - 77.121176e0  * lt3**2
     - 51.352459e0  * lt3**3
     - 15.169160e0  * lt3**4
     - 0.98120322e0 * lt3**5) 
cooling_rates_table['gael'] = CoolingRate('gael', vals)

# Compton cooling (Peebles 1971)
# -- comp --
vals = 5.65e-26 + T*0.0
cooling_rates_table['comp'] = CoolingRate('comp', vals)

#  Photoelectric heating by UV-irradiated dust (Wolfire 1995)
#  with epsilon=0.05, G_0=1.7 (rate in erg s^-1 cm^-3)
# -- gammah --
vals = 8.5e-26 + T*0.0
cooling_rates_table['gammah'] = CoolingRate('gammah', vals)

# COOLING ACTION DEFINITIONS

# Collisional excitations
cooling_action_table["ceHI"] = CoolingAction(
    ["ceHI"], ["HI","de"], "-ceHI * HI * de")
cooling_action_table["ceHeI"] = CoolingAction(
    ["ceHeI"], ["HeII","de"], "-ceHeI * HeII * de**2/4.0")
cooling_action_table["ceHeII"] = CoolingAction(
    ["ceHeII"], ["HeII","de"], "-ceHeII * HeII * de/4.0")

# Collisional ionizations
cooling_action_table["ciHI"] = CoolingAction(
    ["ciHI"], ["HI","de"], "-ciHI*HI*de")
cooling_action_table["ciHeI"] = CoolingAction(
    ["ciHeI"], ["HeI","de"], "-ciHeI*HeI*de/4.0")
cooling_action_table["ciHeII"] = CoolingAction(
    ["ciHeII"], ["HeII","de"], "-ciHeII*HeII*de/4.0")
cooling_action_table["ciHeIS"] = CoolingAction(
    ["ciHeIS"], ["HeII","de"], "-ciHeIS*HeII*de**2/4.0")

# Recombinations
cooling_action_table["reHII"] = CoolingAction(
    ["reHII"], ["HII","de"], "-reHII * HII * de")
cooling_action_table["reHeII1"] = CoolingAction(
    ["reHeII1"], ["HeII","de"], "-reHeII1 * HeII * de/4.0")
cooling_action_table["reHeII2"] = CoolingAction(
    ["reHeII2"], ["HeII","de"], "-reHeII2 * HeII * de/4.0")
cooling_action_table["reHeIII"] = CoolingAction(
    ["reHeIII"], ["HeIII","de"], "-reHeIII * HeIII * de/4.0")

cooling_action_table["compton"] = CoolingAction(
    ["comp"], ["de",], "-comp1*(T - comp2)*de",
    temp_vars = [ ("comp1", "comp * (1.0 * redshift)**4"),
                   ("comp2", "2.73 * (1.0 + redshift)") ],
    other_symbols = ["redshift", "T"])

# SKIP FOR NOW
#cooling_action_table["compton_xray"] = CoolingAction(
#    ["comp_xray"]
"""
!                    X-ray compton heating

     &             - comp_xraya*(tgas(i)-comp_temp)*de(i,j,k)*dom_inv
"""

cooling_action_table["brem"] = CoolingAction(
    ["brem"], ["HII","HeII","HeIII","de"],
        "-brem*(HII * HeII/4.0 + HeIII * de)")

# Skipped
"""
!                    Photoelectric heating by UV-irradiated dust

     &             + float(igammah)*gammaha*(HI(i,j,k)+HII(i,j,k))
     &             *dom_inv)
"""

# Glover & Abel 2008 H2 cooling
# do we have gphdl right?
cooling_action_table["gloverabel08"] = CoolingAction(
    ["gaHI","gaH2","gaHe","gaHp","gael","gphdl"],
    ["HI","H2I","HeI","HII","de"],
    "-(H2I*0.5)*gphdl/(1.0+gphdl1/galdl)",
    temp_vars = [
        ("galdl", "gaHI*HI + gaH2*H2I + gaHe*HeI + gaHp*HII + gael*de"),
        ("gphdl1", "gphdl"),
        ])

if __name__ == "__main__":
    pp = CVODEPrinter()
    print pp.doprint(species_symbols["H2I"])
    print pp.doprint(cooling_action_table["ceHI"].equation)
    print pp.doprint(cooling_action_table["ceHeI"].equation)
    print pp.doprint(cooling_action_table["compton"].equation)
    #print pp.doprint(cooling_action_table["agloverabel08"].equation)
    #print pp.doprint(cooling_action_table["simple"].equation)
