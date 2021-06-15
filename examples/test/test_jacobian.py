import numpy as np
from dengo.chemical_network import \
    ChemicalNetwork, \
    reaction_registry, \
    cooling_registry, \
    species_registry
from dengo.chemistry_constants import tiny, kboltz, mh, G
from dengo.reaction_classes import ReactionCoefficient
import sympy
import dengo.primordial_rates, dengo.primordial_cooling
import pyximport
import os
import pylab

check_defined_envpath()

# this would register all the reaction defined in primodial_rates!
# such that these reaction are set up and can be founnd in reaction registry and species registry
dengo.primordial_rates.setup_primordial()

k01 = ReactionCoefficient("k01[i]")
k22 = ReactionCoefficient("k22[i]")
rk01 = k01.diff("ge")
rk22 = k22.diff("ge")

# h2formation rates
h2mheat = ReactionCoefficient("h2formation_h2mheat[i]")
h2mcool = ReactionCoefficient("h2formation_h2mcool[i]")
ncrn  = ReactionCoefficient("h2formation_ncrn[i]")
ncrd1 = ReactionCoefficient("h2formation_ncrd1[i]")
ncrd2 = ReactionCoefficient("h2formation_ncrd2[i]")
# drate dT
rh2mheat = h2mheat.diff("ge")
rh2mcool = h2mcool.diff("ge")
rncrn  = ncrn.diff("ge")
rncrd1 = ncrd1.diff("ge")
rncrd2 = ncrd2.diff('ge')

# gloverabel rate for h2
gpldl = ReactionCoefficient('gloverabel08_gpldl[i]')
gphdl = ReactionCoefficient('gloverabel08_gphdl[i]')
gaHI  = ReactionCoefficient('gloverabel08_gaHI[i]')
gaH2  = ReactionCoefficient('gloverabel08_gaH2[i]')
gaHe  = ReactionCoefficient('gloverabel08_gaHe[i]')
gaHp  = ReactionCoefficient('gloverabel08_gaHp[i]')
gael  = ReactionCoefficient('gloverabel08_gael[i]')
h2lte = ReactionCoefficient('gloverabel08_h2lte[i]')
# drate dT
rgpldl = gpldl.diff("ge")
rgphdl = gphdl.diff("ge")
rgaHI  = gaHI.diff("ge")
rgaH2  = gaH2.diff("ge")
rgaHe  = gaHe.diff("ge")
rgaHp  = gaHp.diff("ge")
rgael  = gael.diff("ge")
rh2lte = h2lte.diff("ge")
h2opticaldepth = sympy.sympify("h2_optical_depth_approx")

# define symbols for species
H1, H_2, H2, He_1, de, ge = sympy.symbols("H_1, H_2, H2_1, He_1, de, ge")

class TestChemicalNetworkSolution:
    # species_total
    def species_total_solution(self,sp):
        if sp.name == "H2_1":
            return sympy.simplify(k22*H1**3)
        if sp.name == "H_1":
            return sympy.simplify(-k01*H1*de - 2*k22*H1**3)
        if sp.name == "H_2":
            return sympy.simplify(k01*H1*de)
        if sp.name == "de":
            return sympy.simplify(k01*H1*de)
        return sympy.simplify("0")

    def cooling_term_solution(self, ca):
        if ca == "h2formation":
            eqH2formation = 0.5 * ( ncrn/ (ncrd1*H1 + ncrd2*H2) + 1. )**-1.0 *(-h2mcool*H2*H1 + h2mheat*H1**3)
            return sympy.simplify(eqH2formation)
        if ca == "gloverabel08":
            galdl = gaH2*H2 + gaHI*H1 + gaHe*He_1 + gaHp * H_2 + gael * de
            eqga08 = (-h2opticaldepth * h2lte  * H2 / ( 1.0 + h2lte/galdl ))
            return sympy.simplify(eqga08)

    def jacobian_solution(self, s1, s2):
        if s1.name == "ge" and s2.name == "ge":
            # h2 formation
            d1dge =  0.5*(ncrn/(ncrd1*H1 + ncrd2*H2) + 1.0)**(-2.0)* \
                    (-h2mcool*H2*H1 + h2mheat*H1**3)* \
                    (-1.0*ncrn*(-H2*rncrd2 - H1*rncrd1)/(ncrd1*H1 + ncrd2*H2)**2  \
                     - 1.0*rncrn/(ncrd1*H1 + ncrd2*H2)) + 0.5*(ncrn/(ncrd1*H1 + ncrd2*H2) + 1.0)**(-1.0) \
                    *(-H2*H1*rh2mcool + H1**3*rh2mheat)
            # gloverable08
            d2dge = -h2lte*H2*h2opticaldepth*(-h2lte*(-H2*rgaH2 - H1*rgaHI - H_2*rgaHp - He_1*rgaHe - de*rgael) \
                                                /(gaH2*H2 + gaHI*H1 + gaHe*He_1 + gaHp*H_2 + gael*de)**2 - rh2lte/ \
                                                (gaH2*H2 + gaHI*H1 + gaHe*He_1 + gaHp*H_2 + gael*de))/ \
                                                (h2lte/(gaH2*H2 + gaHI*H1 + gaHe*He_1 + gaHp*H_2 + gael*de) + 1.0)**2 \
                                                - H2*h2opticaldepth*rh2lte/(h2lte/(gaH2*H2 + gaHI*H1 + gaHe*He_1 + gaHp*H_2 + gael*de) + 1.0)
            return d1dge + d2dge

        if s1.name == "H_1" and s2.name == "H_1":
            return -k01*de - 6*k22*H1**2

        if s1.name == "H_1" and s2.name == "ge":
            return -rk01*de*H1 - 2*rk22*H1**3

class _ChemicalNetwork():
    def __init__(self):
        primordial_simplified = ChemicalNetwork()
        primordial_simplified.add_reaction("k01")
        primordial_simplified.add_reaction("k22")

        primordial_simplified.add_cooling("h2formation")
        primordial_simplified.add_cooling("gloverabel08")

        self.network = primordial_simplified
        self.solution = TestChemicalNetworkSolution()

    def test_species_total(self):
        # check the equations on the dydt / rhs function
        network = self.network
        solution = self.solution
        for sp in network.required_species:
            eq_sp = sympy.simplify(network.species_total(sp))
            eq_sp_sol = solution.species_total_solution(sp)
            assert eq_sp == eq_sp_sol

    def test_cooling_equation(self):
        network = self.network
        solution = self.solution
        for ca in network.cooling_actions:
            caeq = sympy.simplify(network.cooling_actions[ca].equation)
            caeq_sol = solution.cooling_term_solution(ca)
            assert caeq == caeq_sol

    def test_jacobian(self):
        network = self.network
        solution = self.solution

        for sp in network.required_species:
            if sp.name == "ge":
                ge = sp
                break
        for sp in network.required_species:
            if sp.name == "H_1":
                H_1 = sp
                break

        # test the jacobian term for H2 by H1
        dh1_dh1 = network.species_total("H_1").diff("H_1")
        dh1_dh1_sol = solution.jacobian_solution( H_1, H_1 )
        assert dh1_dh1 == dh1_dh1_sol

        # test the jacobian term for H2 by H1
        dh1_dge = network.species_total("H_1").diff("ge")
        dh1_dge_sol = solution.jacobian_solution( H_1, ge )
        assert dh1_dge == dh1_dge_sol

        # test the jacobian term fod dfge_dge
        dh2formation_dge     = network.cooling_actions["h2formation"].equation.diff("ge")
        dgloverabel08_dge     = network.cooling_actions["gloverabel08"].equation.diff("ge")
        dge_dge = dh2formation_dge  + dgloverabel08_dge
        dge_dge_sol =  solution.jacobian_solution(ge, ge)
        assert dge_dge == dge_dge_sol


def test_jac():
    network = _ChemicalNetwork()
    network.test_species_total()
    network.test_cooling_equation()
    network.test_jacobian()
