"""ChemicalNetwork Module"""
import os
import pkgutil
import warnings
from collections import defaultdict

import h5py
import jinja2
import numpy as np
import sympy
from sympy.printing import ccode
from sympy.utilities import lambdify

from .chemistry_constants import tevk
from .periodic_table import periodic_table_by_name
from .reaction_classes import (
    ChemicalSpecies,
    Species,
    cooling_registry,
    count_m,
    index_i,
    reaction_registry,
    species_registry,
)

ge = Species("ge", 1.0, "Gas Energy")
de = ChemicalSpecies("de", 1.0, pretty_name="Electrons")


class ChemicalNetwork:
    """A chemical network that gathers all the species and reactions.

    Parameters
    ----------

    Attributes
    ----------
    reactions: dict
        the dictionary that holds the reactions added through `self.add_reaction`
    cooling_actions: dict
        the dictionary that tracks the cooling tracked by the `ChemicalNetwork`
    required_species:
        the set of chemical species tracked by `ChemicalNetwork`

    """

    energy_term = None
    skip_weight = ("ge", "de")

    def __init__(self, write_intermediate=False, stop_time=3.1557e13):
        self.reactions = {}
        self.cooling_actions = {}
        self.required_species = set([])
        self.write_intermediate_solutions = write_intermediate
        self.energy_term = species_registry["ge"]
        self.required_species.add(self.energy_term)
        self.stop_time = stop_time
        self.z_bounds = (0.0, 10.0)

        # create a list of species where gamma has to be
        # integrate / calculate separately
        self.interpolate_gamma_species = set([])
        self.interpolate_gamma_species_name = set(["H2_1", "H2_2"])

        self.threebody = 4
        self.equilibrium_species = set()
        self.enforce_conservation = False

        # store paths locations
        self._hdf5_path = None
        self._dengo_install_path = None
        self._libtool_path = None
        self._cvode_path = None
        self._suitesparse_path = None

    def add_collection(self, species_names, cooling_names, reaction_names):
        for s in species_names:
            self.add_species(s)
        for c in cooling_names:
            self.add_cooling(c, False)
        for r in reaction_names:
            self.add_reaction(r, False)

    def set_equilibrium_species(self, species):
        sp = species_registry.get(species, species)
        self.equilibrium_species.add(sp)

    def update_ode_species(self):
        """this refers to the set of species that our ODE solver solves"""
        self.ode_species = self.required_species.copy()
        for s in self.equilibrium_species:
            self.ode_species.remove(s)

    @property
    def chemical_species(self):
        """helper function for code generation"""
        chemical_species = self.required_species.copy()
        chemical_species.remove(self.energy_term)
        return chemical_species

    def solve_equilibrium_abundance(self, species):
        """write the equilibrium abundance of species sp
        this is by setting dydt = 0

        Parameters
        ----------
            species: dengo.reaction_classes.ChemicalSpecies

        Return
        ------
            equilibrium_sol: str
                C-Code snippet that gives equilibrium species formula

        """
        from sympy.solvers import solve

        eq = self.species_total(species)
        equil_sol = solve(eq, species)
        return ccode(equil_sol[0], assign_to=species)

    def add_species(self, species):
        sp = species_registry.get(species, species)
        self.required_species.add(sp)

    def add_reaction(self, reaction, auto_add=True):
        """Add reaction to the `ChemicalNetwork`

        Parameters
        ----------
        reaction
            the string of the reactions that are previously registered in `reaction_registry`

        Returns
        -------
        None

        """
        reaction = reaction_registry.get(reaction, reaction)
        if reaction.name in self.reactions:
            raise RuntimeError
        if auto_add:
            for _, s in reaction.left_side:
                self.required_species.add(s)

                if s.name in self.interpolate_gamma_species_name:
                    self.interpolate_gamma_species.add(s)
            for _, s in reaction.right_side:
                self.required_species.add(s)

                if s.name in self.interpolate_gamma_species_name:
                    self.interpolate_gamma_species.add(s)
        else:
            for _, s in reaction.left_side:
                if s not in self.required_species:
                    raise RuntimeError
            for _, s in reaction.right_side:
                if s not in self.required_species:
                    raise RuntimeError
        self.reactions[reaction.name] = reaction
        print("Adding reaction: %s" % reaction)

    def add_cooling(self, cooling_term, auto_add=True):
        """Add cooling to the `ChemicalNetwork`

        Parameters
        ----------
        cooling term
            the string that corresponds to the cooling that are previously registered in `cooling_registry`

        Returns
        -------
        None

        """
        cooling_term = cooling_registry.get(cooling_term, cooling_term)
        if cooling_term.name in self.cooling_actions:
            raise RuntimeError
        if auto_add:
            self.required_species.update(cooling_term.species)
        else:
            for s in cooling_term.species:
                if s not in self.required_species:
                    print("SPECIES NOT FOUND", s)
                    raise RuntimeError
        self.cooling_actions[cooling_term.name] = cooling_term

    def init_temperature(self, T_bounds=(1, 1e8), n_bins=1024):
        """Initialize the range of temperature over which the rate tables are generated

        Parameters
        ----------
        T_bounds: List[Float, Float], optional (default=(1,1e8))
            the range over which the rates table is interpolated

        n_bins: int, optional (default=1024)

        """
        self.n_bins = n_bins
        self.T = np.logspace(
            np.log(T_bounds[0]), np.log(T_bounds[1]), n_bins, base=np.e
        )
        self.logT = np.log(self.T)
        self.tev = self.T / tevk
        self.logtev = np.log(self.tev)
        self.T_bounds = T_bounds

    def init_redshift(self, z_bounds=(0.0, 10.0), n_z_bins=100):
        """Initialize the range of redshift over which the rate tables are generated

        Parameters
        ----------
        z_bounds: List[Float, Float], optional (default=(0.0,10.0))
            the range over which the rates table is interpolated

        n_bins: int, optional (default=1024)

        """
        self.n_z_bins = n_z_bins
        self.z = np.logspace(
            np.log(z_bounds[0] + 1.0), np.log(z_bounds[1] + 1.0), n_z_bins, base=np.e
        )
        self.z -= 1.0
        self.logz = np.log(self.z)
        self.z_bounds = z_bounds

    def species_total(self, species):
        """Output the RHS function of the given species

        Parameters
        ----------
        species: str,

        Return
        ------
        eq: Sympy.symbols
            the dydt (RHS function) in sympy expressions
        """
        eq = sympy.sympify("0")
        for _, rxn in sorted(self.reactions.items()):
            eq += rxn.species_equation(species)
        return eq

    def species_reactions(self, species):
        tr = []
        for _, rxn in sorted(self.reactions.items()):
            if species in rxn:
                tr.append(rxn)
        return tr

    def species_list(self):
        species_list = []
        for s in sorted(self.required_species):
            species_list.append(s.name)
        return species_list

    def __iter__(self):
        for _, rxn in sorted(self.reactions.items()):
            yield rxn

    def print_ccode(self, species, assign_to=None):
        # assign_to = sympy.IndexedBase("d_%s" % species.name, (count_m,))
        if assign_to is None:
            assign_to = sympy.Symbol("d_%s[i]" % species.name)
        if species == self.energy_term:
            return self.print_cooling(assign_to)
        eq = self.species_total(species)
        # eq = eq.replace(
        #    lambda x: x.is_Pow and x.exp > 0,
        #    lambda x: sympy.Symbol("*".join([x.base.name] * x.exp)),
        # )

        return ccode(eq, assign_to=assign_to)

    def cie_optical_depth_approx(self):
        return sympy.Symbol("cie_optical_depth_approx")

    def print_cooling(self, assign_to):
        eq = sympy.sympify("0")
        for term in self.cooling_actions:
            # TODO: make it a more general check case?
            if term not in ["h2formation"] and "cie_cooling" in self.cooling_actions:
                # This is independent of the continuum optical depth
                # from the CIE
                eq += (
                    self.cooling_actions[term].equation
                    * self.cie_optical_depth_approx()
                )
            else:
                eq += self.cooling_actions[term].equation

        return ccode(eq, assign_to=assign_to)

    def get_conserved_dict(self):
        self.conserved_dict = defaultdict(lambda: [])
        species_considered = ["H", "He"]
        for s in self.species_list():
            sp = species_registry[s]
            if s in self.skip_weight:
                continue
            for e in sp.elements:
                if e in species_considered:
                    self.conserved_dict[e] += [sp]
        return self.conserved_dict

    def print_conserved_species(self, s, assign_to, is_mass_density=True):
        """Calculate the total {H, He, ...} atoms locked in species considered"""
        self.get_conserved_dict()
        v = self.conserved_dict[s]
        eq = sympy.sympify("0")
        if is_mass_density:
            for i in v:
                eq += i.symbol
        else:
            for i in v:
                eq += i.symbol * i.elements[s]
        return ccode(eq, assign_to=assign_to)

    def print_apply_conservation(self, s, assign_to):
        for ele, v in self.conserved_dict.items():
            _, w, _ = periodic_table_by_name[ele]
            vname = [i.name for i in v]
            if s in vname:
                return "{0} = {1}*f{2}*density/total_{2};".format(assign_to, s, ele)
        return ""

    def cie_optical_depth_correction(self):
        mdensity = sympy.Symbol("mdensity")
        tau = (mdensity / 3.3e-8) ** 2.8
        tau = sympy.Max(tau, 1e-5)
        ciefudge = sympy.Min((1.0 - sympy.exp(-tau)) / tau, 1.0)
        return ciefudge

    def print_jacobian_component(self, s1, s2, assign_to=None, print_zeros=True):
        """Prints the jacobian of d (ds1/dt)/ ds2"""

        if s1 == self.energy_term:
            st = sum(
                self.cooling_actions[ca].equation for ca in sorted(self.cooling_actions)
            )
        else:
            st = self.species_total(s1)

        if assign_to is None:
            assign_to = sympy.Symbol("d_%s_%s" % (s1.name, s2.name))
        if isinstance(st, (list, tuple)):
            codes = []
            for temp_name, temp_eq in st[0]:
                teq = sympy.sympify(temp_eq)
                codes.append(ccode(teq, assign_to=temp_name))
            codes.append(ccode(st[1], assign_to=assign_to))
            return "\n".join(codes)

        eq = sympy.diff(st, s2.symbol)
        # eq = eq.replace(
        #    lambda x: x.is_Pow and x.exp > 0 and x.exp == sympy.Integer,
        #    lambda x: sympy.Symbol("*".join([x.base.name] * x.exp)),
        # )

        if eq == sympy.sympify("0") and not print_zeros:
            return ""

        return ccode(eq, assign_to=assign_to)

    def get_sparse_matrix_component(
        self, sparse_type="CSR", return_type="component", assign_to="data"
    ):
        """Return a sparse CSR matrix for the jacobian, primarily used in templating"""

        k = 0
        colvals = []
        all_comp = []
        rowptrs = []
        s1_list = []
        s2_list = []
        i1_list = []
        i2_list = []

        # sp_list = list( sorted(self.required_species) )
        self.update_ode_species()
        sp_list = list(sorted(self.ode_species))
        s1_now = sp_list[0]

        # getting rowptrs is a little tricky....
        # its the counter
        for i1, s1 in enumerate(sp_list):
            rowptrs.append(k)
            for i2, s2 in enumerate(sp_list):
                jac_comp = self.print_jacobian_component(
                    s1, s2, print_zeros=False, assign_to=""
                )

                if jac_comp:
                    colvals.append(i2)
                    all_comp.append(jac_comp)

                    s1_list.append(s1)
                    s2_list.append(s2)

                    i1_list.append(i1)
                    i2_list.append(i2)

                    # if s1 == s1_now:
                    #    rowptrs.append(k)
                    #    if i1 < len(sp_list) - 1:
                    #        s1_now = sp_list[i1 + 1]
                    #   else:
                    #        s1_now = None
                    k += 1
        # rowptrs.append(k)
        if return_type == "component":
            return zip(colvals, all_comp, s1_list, s2_list, range(k))
        elif return_type == "indexptrs":
            return rowptrs
        elif return_type == "index":
            return zip(i1_list, i2_list, range(k))
        elif return_type == "nsparse":
            return k

    def print_JacTimesVec_component(self, s1, assign_to=None):
        """
        Compute the product of Jacobian * Vec for a given Vec
        Might be useful when we use the CVSpils solver
        """
        if s1 == self.energy_term:
            st = sum(
                self.cooling_actions[ca].equation for ca in sorted(self.cooling_actions)
            )

        else:
            st = self.species_total(s1)
        if assign_to is None:
            assign_to = sympy.Symbol("d_%s_dy_y" % (s1.name))
        if isinstance(st, (list, tuple)):
            codes = []
            for temp_name, temp_eq in st[0]:
                teq = sympy.sympify(temp_eq)
                codes.append(ccode(teq, assign_to=temp_name))
            codes.append(ccode(st[1], assign_to=assign_to))
            return "\n".join(codes)

        JtV_eq = sympy.sympify("0")
        mdensity = sympy.sympify("mdensity")
        T_energy = sympy.Symbol("T{0}[i]".format(self.energy_term.name))

        i = 0
        for s2 in sorted(self.required_species):

            vec = sympy.sympify("v{0}".format(i))
            newterm = sympy.diff(st, s2.symbol) * vec

            if s1.name == "ge":
                newterm /= mdensity
            if s2.name == "ge":
                newterm *= T_energy

            JtV_eq += newterm
            i += 1
        return ccode(JtV_eq, assign_to=assign_to)

    def print_mass_density(self):
        # Note: this assumes things are number density at this point
        eq = sympy.sympify("0")
        for s in sorted(self.required_species):
            if s.weight > 0:
                if s.name in self.skip_weight:
                    continue
                eq += s.symbol * s.weight
        return ccode(eq)

    def interpolate_species_gamma(self, sp, deriv=False):
        if sp.name in ("H2_1", "H2_2"):
            expr_gammaH2 = self.species_gamma(
                species_registry["H2_1"], temp=True, name=False
            )

            if deriv is True:
                expr_dgammaH2_dT = sympy.diff(expr_gammaH2, "T")
                f_dgammaH2_dT = lambdify("T", expr_dgammaH2_dT)
                dgammaH2_dT = f_dgammaH2_dT(self.T)
                _i1 = np.isnan(dgammaH2_dT)
                dgammaH2_dT[_i1] = 0.0

                return dgammaH2_dT

            f_gammaH2 = lambdify("T", expr_gammaH2)
            gammaH2_T = f_gammaH2(self.T)
            _i1 = np.isnan(gammaH2_T)
            gammaH2_T[_i1] = 7.0 / 5.0
            return gammaH2_T

        print("Provide your gamma function for {}".format(sp.name))
        raise RuntimeError

    def species_gamma(self, species, temp=False, name=True):
        """create an expression for the species gamma
        This should be depracated, as species gamma should be attached to each species!
        and move to ChemicalSpecies
        """
        if species in self.interpolate_gamma_species:
            sp_name = species.name
            T = sympy.Symbol("T")

            if temp and name:
                # so gamma enters as a function that depends on temperature
                # when gamma factor enters the temperature calculation which
                # involves derivatives of temperature (dgammaH2_dT) this term will not be dropped off
                # and be left as a function of T, which can later be supplied from the interpolated values
                # as gammaH2 does

                # this returns name of the gamma as a function of T
                # goes into the analytical differntiation for energy
                f_gammaH2 = sympy.Function("gamma%s" % sp_name)(T)
                return f_gammaH2
            elif temp and ~name:
                # x = 6100.0/T
                # expx = sympy.exp(x)
                # gammaH2_expr = 2.0 / (5.0 + 2.0*x*x*expx / (expx - 1 )**2.0 ) + 1

                T0 = T ** (1 / 6.5)
                a0 = 64.2416
                a1 = -9.36334
                a2 = -0.377266
                a3 = 69.8091
                a4 = 0.0493452
                a5 = 2.28021
                a6 = 0.115357
                a7 = 0.114188
                a8 = 2.90778
                a9 = 0.689708

                gammaH2_expr = (
                    sympy.exp(-a0 * T0**a1) * (a2 + T0**-a3)
                    + a4 * sympy.exp(-((T0 - a5) ** 2) / a6)
                    + a7 * sympy.exp(-((T0 - a8) ** 2) / a9)
                    + 5.0 / 3.0
                )
                # x = 6100.0/T
                # expx = sympy.exp(x)
                # gammaH2_expr = 2.0 / (5.0 + 2.0*x*x*expx / (expx - 1 )**2.0 ) + 1

                return gammaH2_expr

            gammaH2 = sympy.Symbol("gamma%s" % sp_name)
            return gammaH2

        gamma = sympy.Symbol("gamma")
        return gamma

    def gamma_factor(self, temp=False):
        """Return an sympy expression of sum ( 1 / (gamma - 1) )"""
        eq = sympy.sympify("0")
        for s in sorted(self.required_species):
            if s.name != "ge":
                eq += (sympy.sympify(s.name)) / (self.species_gamma(s, temp=temp) - 1.0)
        return eq

    def temperature_calculation(
        self, derivative=False, derivative_dge_dT=False, get_dge=False
    ):
        """Evaluate Temperature symbolically"""
        # If derivative=True, returns the derivative of
        # temperature with respect to ge.  Otherwise,
        # returns just the temperature function
        ge = sympy.Symbol("ge")
        function_eq = (sympy.Symbol("density") * ge * sympy.Symbol("mh")) / (
            sympy.Symbol("kb") * (self.gamma_factor())
        )
        # function_eq = (ge) / \
        #    (sympy.Symbol('kb') * (self.gamma_factor()))
        if derivative is True:
            deriv_eq = sympy.diff(function_eq, ge)
            return ccode(deriv_eq)
        elif derivative_dge_dT is True:
            # when H2 presents, the gamma is dependent on  temperature
            # therefore temperature must iterate to a convergence for a given ge
            # this part evaluates the derivatives of the function ge with respect to T
            T = sympy.Symbol("T")
            f = (
                self.gamma_factor(temp=True)
                * sympy.Symbol("kb")
                * sympy.Symbol("T")
                / sympy.Symbol("mh")
                / sympy.Symbol("density")
            )
            dge_dT = sympy.diff(f, T)
            tmp = sympy.Symbol("tmp")
            for sp in self.interpolate_gamma_species:
                # substitute the sympy function with sympy Symbols
                sym_fgamma = sympy.Function("gamma%s" % sp.name)(T)
                sym_dfgamma = sympy.diff(sym_fgamma, T)
                dgamma = sympy.Symbol("dgamma%s_dT" % sp.name)
                dge_dT = dge_dT.subs({sym_dfgamma: dgamma})

                fgamma = sympy.Symbol("gamma%s" % sp.name)
                dge_dT = dge_dT.subs({sym_fgamma: tmp})
                dge_dT = dge_dT.subs({tmp: fgamma})

                # substitute 1/ gamma - 1 with
                # the expression _gamma{{sp}}_m1
                _gamma_m1 = sympy.Symbol("_gamma%s_m1" % sp.name)
                dge_dT = dge_dT.subs({1 / (fgamma - 1): _gamma_m1})

            gamma = sympy.Symbol("gamma")
            _gamma_m1 = sympy.Symbol("_gamma_m1")
            dge_dT = dge_dT.subs({1 / (gamma - 1): _gamma_m1})

            # to unravel the algebric power
            dge_dT = dge_dT.replace(
                lambda x: x.is_Pow and x.exp > 0,
                lambda x: sympy.Symbol("*".join([x.base.name] * x.exp)),
            )

            return ccode(dge_dT)
        elif get_dge is True:
            T = sympy.Symbol("T")
            dge = self.gamma_factor(temp=True) * sympy.Symbol("kb") * T / sympy.Symbol(
                "mh"
            ) / sympy.Symbol("density") - sympy.Symbol("ge")

            tmp = sympy.Symbol("tmp")
            for sp in self.interpolate_gamma_species:
                sym_fgamma = sympy.Function("gamma%s" % sp.name)(T)
                fgamma = sympy.Symbol("gamma%s" % sp.name)
                dge = dge.subs({sym_fgamma: tmp})
                dge = dge.subs({tmp: fgamma})

                # substitute 1/ gamma - 1 with
                # the expression _gamma{{sp}}_m1
                _gamma_m1 = sympy.Symbol("_gamma%s_m1" % sp.name)
                dge = dge.subs({1 / (fgamma - 1): _gamma_m1})

            gamma = sympy.Symbol("gamma")
            _gamma_m1 = sympy.Symbol("_gamma_m1")
            dge = dge.subs({1 / (gamma - 1): _gamma_m1})

            # to unravel the algebric power
            dge = dge.replace(
                lambda x: x.is_Pow and x.exp > 0,
                lambda x: sympy.Symbol("*".join([x.base.name] * x.exp)),
            )

            return ccode(dge)

        gamma = sympy.Symbol("gamma")
        _gamma_m1 = sympy.Symbol("_gamma_m1")
        tmp = sympy.Symbol("tmp")
        T = sympy.Symbol("T")
        for sp in self.interpolate_gamma_species:
            # substitute the sympy function with sympy Symbols
            sym_fgamma = sympy.Function("gamma%s" % sp.name)(T)
            _fgamma_m1 = sympy.Symbol("_gamma%s_m1" % sp.name)
            function_eq = function_eq.subs({1 / (sym_fgamma - 1): tmp})
            function_eq = function_eq.subs({tmp: _fgamma_m1})

        function_eq = function_eq.subs({1 / (gamma - 1): _gamma_m1})
        gamma = sympy.Symbol("gamma")
        _gamma_m1 = sympy.Symbol("_gamma_m1")

        return ccode(function_eq)

    # This function computes the total number density
    def calculate_number_density(self, values):
        # values should be a dict with all of the required species in it
        # The values should be in *mass* density
        n = np.zeros_like(list(values.values())[0])
        for s in self.required_species:
            if s.name in self.skip_weight:
                continue
            n += values[s.name] / s.weight
        return n

    # This function counts up the total number of free electrons
    def calculate_free_electrons(self, values):
        """evaluate the number density of free electrons
        Parameters
        ----------
        values: dict[np.ndarray]
            the abundance vectors specified by the user, it is assumed to be in mass density

        Returns
        -------
        n: np.ndarray
            the numpy array that holds the number density of electrons
        """
        # values should be add dict with all of the required species in it
        # The values should be in *mass* density
        n = np.zeros_like(list(values.values())[0])
        for s in self.required_species:
            if s.name in self.skip_weight:
                continue
            n += (values[s.name] / s.weight) * s.free_electrons
        return n

    # This computes the total mass density from abundance fractions
    def calculate_mass_density(self, values):
        # values should be a dict with all of the required species in it
        # The values should be in *mass* density
        n = np.zeros_like(list(values.values())[0])
        for s in self.required_species:
            if s.name in self.skip_weight:
                continue
            n += values[s.name] * s.weight
        return n

    # This function sums the densities (mass or number depending on what
    # is fed in) of non-electron species
    def calculate_total_density(self, values):
        """evalue the total density by summing over species with mass
        Parameters
        ----------
        values: dict[values]
            the abundance vectors specified by the used

        Note
        ----
        Depending on what the

        Returns
        -------
        n: np.ndarray
            an numpy array of total density
        """
        # values should be a dict with all of the required species in it
        # The values should be in *mass* density
        n = np.zeros_like(list(values.values())[0])
        for s in self.required_species:
            if s.name in self.skip_weight:
                continue
            n += values[s.name]
        return n

    def convert_to_mass_density(self, values):
        """converts from number density to mass density"""
        for s in self.required_species:
            if s.name in self.skip_weight:
                continue
            values[s.name] = values[s.name] * s.weight
        return values

    def write_cuda_solver(
        self,
        solver_name,
        solver_template="cvode_cuda",
        output_dir=".",
        init_values=None,
        main_name="main",
        input_is_number=False,
    ):
        self.input_is_number = input_is_number
        self.update_ode_species()

        if not os.path.isdir(output_dir):
            os.makedirs(output_dir)

        root_path = os.path.join(os.path.dirname(__file__), "templates")
        template_path = os.path.join(root_path, "cvode_cuda")
        env = jinja2.Environment(
            extensions=["jinja2.ext.loopcontrols"],
            loader=jinja2.FileSystemLoader(template_path),
        )
        template_vars = dict(network=self, solver_name=solver_name)

        for suffix in (".cu", "_main.cu", ".h"):
            iname = "%s%s" % (solver_template, suffix)
            oname = os.path.join(output_dir, "%s_solver%s" % (solver_name, suffix))
            template_inst = env.get_template(iname + ".template")
            solver_out = template_inst.render(**template_vars)
            with open(oname, "w", encoding="UTF-8") as f:
                f.write(solver_out)

        # now we create a Makefile
        try:
            temp_dir, _ = os.path.split(solver_template)
            iname = os.path.join(temp_dir, "Makefile")
            oname = os.path.join(output_dir, "Makefile")
            template_inst = env.get_template(iname + ".template")
            solver_out = template_inst.render(**template_vars)
            with open(oname, "w", encoding="UTF-8") as f:
                f.write(solver_out)
        except:
            print("This does not have a Makefile.template")
            pass

        # This writes out the rates for the species in the
        # chemical network to HDF5 files which can later be
        # read by the C++ code that is output by the template
        ofn = os.path.join(output_dir, "%s_tables.h5" % solver_name)
        f = h5py.File(ofn, "w")

        for rxn in sorted(self.reactions.values()):
            f.create_dataset(
                "/%s" % rxn.name, data=rxn.coeff_fn(self).astype("float64")
            )
            if hasattr(rxn, "tables"):
                for tab in rxn.tables:
                    print(rxn.name, tab, rxn)
                    f.create_dataset(
                        "/%s_%s" % (rxn.name, tab),
                        data=rxn.tables[tab](self).astype("float64"),
                    )

        for action in sorted(self.cooling_actions.values()):
            for tab in action.tables:
                f.create_dataset(
                    "/%s_%s" % (action.name, tab),
                    data=action.tables[tab](self).astype("float64"),
                )

        for sp in sorted(self.interpolate_gamma_species):
            name = sp.name
            f.create_dataset(
                "/gamma%s" % name,
                data=self.interpolate_species_gamma(sp).astype("float64"),
            )
            f.create_dataset(
                "/dgamma%s_dT" % name,
                data=self.interpolate_species_gamma(sp, deriv=True).astype("float64"),
            )
        f.close()

    def set_library_paths(self):
        if "HDF5_PATH" in os.environ:
            self._hdf5_path = os.environ["HDF5_PATH"]
        else:
            raise ValueError("Need to supply HDF5_PATH")

        if "DENGO_INSTALL_PATH" in os.environ:
            self._dengo_install_path = os.environ["DENGO_INSTALL_PATH"]
        else:
            raise ValueError("Need to supply DENGO_INSTALL_PATH")

        if "LIBTOOL_PATH" in os.environ:
            self._libtool_path = os.environ["LIBTOOL_PATH"]

    def write_solver(
        self,
        solver_name,
        solver_option=None,
        solver_template="rates_and_rate_tables",
        ode_solver_source="BE_chem_solve.C",
        output_dir=".",
        init_values=None,
    ):
        """Write the chemistry solver, based on the specified solver templates

        Parameters
        ----------
        solver_name: str
            the name of the solver
        solver_template: str, default = rate_and_rate_tables
            the jinja2 template the `ChemicalNetwork` populates, read dengo/templates for examples
        ode_solver_source: str, default = BE_chem_solve.C
            the ODE solver used
        output_dir: str, default = .
            the directory at which these files would be generated
        init_values: Optional[dict, None], default=None
            a initial value hdf5 file would be generated that can be fed to the solver (for testing purpose)
        """
        self.update_ode_species()

        if solver_option == "CVODE":
            solver_template = "be_chem_solve/rates_and_rate_tables"
            ode_solver_source = "BE_chem_solve.C"
        elif solver_option == "BE_CHEM_SOLVE":
            solver_template = ("cv_omp/sundials_CVDls",)
            ode_solver_source = "initialize_cvode_solver.C"
        else:
            warnings.warn(
                "going to be replaced solver_template and ode_solver_source with solver_option",
                PendingDeprecationWarning,
            )

        # with Dan suggestions, I have included the path to write the solver!
        # please specify the environ path!

        # Only HDF5_PATH and DENGO_INSTALL_PATH is needed for
        # running dengo, but then it can only be coupled with the 1st order BDF solver
        self.set_library_paths()
        # CVODE is optional, but recommended for performance boost
        if "CVODE_PATH" in os.environ:
            self._cvode_path = os.environ["CVODE_PATH"]
            # SuiteSpare is optional as well
            if "SUITESPARSE_PATH" in os.environ:
                self._suitesparse_path = os.environ["SUITESPARSE_PATH"]
            else:
                print("Need to supply SUITESPARSE_PATH if CVODE is compiled with KLU")
        else:
            print("CVODE_PATH is optional")
            print(
                "OR: You can choose to use the first order BDF solver that comes with dengo"
            )
            print("In that case, please set ode_solver_source to 'BE_chem_solve.C'")
            print("and solver_template to 'rates_and_rate_tables'. ")
            solver_template = "be_chem_solve/rates_and_rate_tables"
            ode_solver_source = "BE_chem_solve.C"
            # raise ValueError()

        if not os.path.isdir(output_dir):
            os.makedirs(output_dir)

        self.write_templates(
            solver_name, ode_solver_source, solver_template, output_dir, init_values
        )
        self.write_rates_to_file(output_dir, solver_name)
        self.write_initial_condition_to_file(solver_name, init_values, output_dir)

    def write_initial_condition_to_file(self, solver_name, init_values, output_dir):
        """write initial value to a hdf5 file,this serves as an input to the C-solver

        Parameters
        ----------
        solver_name: str
            name of the solver, the output is of the form f"{solver_name}_initial_conditions.h5"
        init_values: dict[np.array]
            an dictionary that holds the initial abundance to feed into the solver

        Returns
        -------
        None

        """
        if init_values is None:
            return

        # This write outs a set of initial conditions for to be fed
        # into C++ code for testing purposes
        ofn = os.path.join(output_dir, f"{solver_name}_initial_conditions.h5")
        f = h5py.File(ofn, "w")
        for name, init_value in init_values.items():
            f.create_dataset("/%s" % name, data=init_value.astype("float64"))
        f.close()

    def write_templates(
        self,
        solver_name,
        ode_solver_source,
        solver_template,
        output_dir=".",
        init_values=None,
    ):
        """Write the ODE solver based on the specified ChemicalNetwork"""
        # What we are handed here is:
        #   * self, a python object which holds all of the species, reactions,
        #     rate, etc. that we're keeping track of and will be solving in Enzo
        #   * solver_name, an identifier to produce a unique template and to
        #     correctly grab the right HDF5 tables
        #   * ode_solver_source, the ODE solver we will be linking against
        #   * solver_template, the template for the solver we will write
        #   * output_dir, a place to stick everythign
        #   * init_values, a set of initial conditions
        #
        # To utilize these inside our template, we will generate convenience
        # handlers that will explicitly number them.
        root_path = os.path.join(os.path.dirname(__file__), "templates")
        env = jinja2.Environment(
            extensions=["jinja2.ext.loopcontrols"],
            loader=jinja2.FileSystemLoader(root_path),
        )
        template_vars = dict(
            network=self, solver_name=solver_name, init_values=init_values
        )

        for suffix in (
            ".C",
            "_main.C",
            ".h",
            "_run.pyx",
            "_run.pyxbld",
            "_run.pyxdep",
            "_run.pxd",
            "_main.py",
        ):
            iname = "%s%s" % (solver_template, suffix)
            oname = os.path.join(output_dir, "%s_solver%s" % (solver_name, suffix))
            template_inst = env.get_template(iname + ".template")
            solver_out = template_inst.render(**template_vars)
            with open(oname, "w", encoding="UTF-8") as f:
                f.write(solver_out)

        # now we create a Makefile
        try:
            temp_dir, _ = os.path.split(solver_template)
            iname = os.path.join(temp_dir, "Makefile")
            oname = os.path.join(output_dir, "Makefile")
            template_inst = env.get_template(iname + ".template")
            solver_out = template_inst.render(**template_vars)
            with open(oname, "w", encoding="UTF-8") as f:
                f.write(solver_out)
        except:
            print("This does not have a Makefile.template")
            pass

        env = jinja2.Environment(
            extensions=["jinja2.ext.loopcontrols"],
            loader=jinja2.PackageLoader("dengo", "solvers"),
        )

        # Now we copy over anything else we might need.
        if ode_solver_source is not None:

            # iname = os.path.join("solvers", ode_solver_source)
            # src = pkgutil.get_data("dengo", iname)
            # src = src.decode("utf-8")

            template_inst = env.get_template(ode_solver_source)
            solver_out = template_inst.render(**template_vars)

            with open(
                os.path.join(output_dir, ode_solver_source), "w", encoding="UTF-8"
            ) as f:
                f.write(solver_out)
        if solver_template.endswith(".c.template"):
            hfn = solver_template.rsplit(".", 2)[0] + ".h"
            src = pkgutil.get_data("dengo", os.path.join("templates", hfn))
            with open(os.path.join(output_dir, hfn), "w", encoding="UTF-8") as f:
                f.write(src)

    def write_rates_to_file(self, output_dir, solver_name):
        """This writes out the rates for the species in the
        chemical network to HDF5 files which can later be
        read by the C++ code that is output by the template

        f"{solver_name}_tables.h5" is the file that contains all the reaction rates and cooling rates


        """
        ofn = os.path.join(output_dir, f"{solver_name}_tables.h5")
        f = h5py.File(ofn, "w")

        for rxn in sorted(self.reactions.values()):
            f.create_dataset(
                "/%s" % rxn.name, data=rxn.coeff_fn(self).astype("float64")
            )
            if hasattr(rxn, "tables"):
                for tab in rxn.tables:
                    print(rxn.name, tab, rxn)
                    f.create_dataset(
                        "/%s_%s" % (rxn.name, tab),
                        data=rxn.tables[tab](self).astype("float64"),
                    )
        for action in sorted(self.cooling_actions.values()):
            for tab in action.tables:
                f.create_dataset(
                    "/%s_%s" % (action.name, tab),
                    data=action.tables[tab](self).astype("float64"),
                )

        for sp in sorted(self.interpolate_gamma_species):
            name = sp.name
            f.create_dataset(
                "/gamma%s" % name,
                data=self.interpolate_species_gamma(sp).astype("float64"),
            )
            f.create_dataset(
                "/dgamma%s_dT" % name,
                data=self.interpolate_species_gamma(sp, deriv=True).astype("float64"),
            )
        f.close()
