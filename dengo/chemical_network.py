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
import numpy as np
from .chemistry_constants import tevk, tiny, mh
from .reaction_classes import reaction_registry, cooling_registry, \
    count_m, index_i, species_registry, Species, ChemicalSpecies
import types
import sympy
import pkgutil
import os
import jinja2
import h5py
from sympy.printing import ccode
from sympy.utilities import lambdify

ge = Species("ge", 1.0, "Gas Energy")
de = ChemicalSpecies("de", 1.0, pretty_name = "Electrons")

class ChemicalNetwork(object):

    energy_term = None
    skip_weight = ("ge", "de")


    def __init__(self, write_intermediate = False, stop_time = 3.1557e13):
        self.reactions = {}
        self.cooling_actions = {}
        self.required_species = set([])
        self.write_intermediate_solutions = write_intermediate
        self.energy_term = species_registry["ge"]
        self.required_species.add(self.energy_term)
        self.stop_time = stop_time
        self.z_bounds = (0.0, 0.0)

        # create a list of species where gamma has to be
        # integrate / calculate separately
        self.interpolate_gamma_species = set([])
        self.interpolate_gamma_species_name = set(['H2_1', 'H2_2'])

        self.threebody = 4

    def add_collection(self, species_names, cooling_names, reaction_names):
        for s in species_names:
            self.add_species(s)
        for c in cooling_names:
            self.add_cooling(c, False)
        for r in reaction_names:
            self.add_reaction(r, False)

    def add_species(self, species):
        sp = species_registry.get(species, species)
        self.required_species.add(sp)

    def add_reaction(self, reaction, auto_add = True):
        reaction = reaction_registry.get(reaction, reaction)
        if reaction.name in self.reactions:
            raise RuntimeError
        if auto_add:
            for n, s in reaction.left_side:
                self.required_species.add(s)

                if s.name in self.interpolate_gamma_species_name:
                    self.interpolate_gamma_species.add(s)
            for n, s in reaction.right_side:
                self.required_species.add(s)

                if s.name in self.interpolate_gamma_species_name:
                    self.interpolate_gamma_species.add(s)
        else:
            for n, s in reaction.left_side:
                if s not in self.required_species:
                    raise RuntimeError
            for n, s in reaction.right_side:
                if s not in self.required_species:
                    raise RuntimeError
        self.reactions[reaction.name] = reaction
        reaction.coeff_sym.energy = self.energy_term
        print ("Adding reaction: %s" % reaction)

    def add_cooling(self, cooling_term, auto_add = True):
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

    def init_temperature(self, T_bounds = (1, 1e8), n_bins=1024):
        self.n_bins = n_bins
        self.T = np.logspace(np.log(T_bounds[0]),
                             np.log(T_bounds[1]),
                             n_bins, base = np.e)
        self.logT = np.log(self.T)
        self.tev = self.T / tevk
        self.logtev = np.log(self.tev)
        self.T_bounds = T_bounds

    def init_redshift(self, z_bounds = (0.0, 10.0), n_z_bins=100):
        self.n_z_bins = n_z_bins
        self.z = np.logspace(np.log(z_bounds[0] + 1.0),
                             np.log(z_bounds[1] + 1.0),
                             n_z_bins, base = np.e)
        self.z -= 1.0
        self.logz = np.log(self.z)
        self.z_bounds = z_bounds

    def species_total(self, species):
        eq = sympy.sympify("0")
        for rn, rxn in sorted(self.reactions.items()):
            eq += rxn.species_equation(species)
        return eq

    def species_reactions(self, species):
        tr = []
        for rn, rxn in sorted(self.reactions.items()):
            if species in rxn: tr.append(rxn)
        return tr

    def species_list(self):
        species_list = []
        for s in sorted(self.required_species):
            species_list.append(s.name)
        return species_list

    def __iter__(self):
        for rname, rxn in sorted(self.reactions.items()):
            yield rxn

    def print_ccode(self, species, assign_to = None):
        #assign_to = sympy.IndexedBase("d_%s" % species.name, (count_m,))
        if assign_to is None: assign_to = sympy.Symbol("d_%s[i]" % species.name)
        if species == self.energy_term:
            return self.print_cooling(assign_to)
        eq = self.species_total(species)
        eq = eq.replace(
                lambda x: x.is_Pow and x.exp > 0,
                lambda x: sympy.Symbol('*'.join([x.base.name]*x.exp)) )


        return ccode(eq , assign_to = assign_to)

    def cie_optical_depth_approx(self):
        return sympy.Symbol("cie_optical_depth_approx")

    def print_cooling(self, assign_to):
        eq = sympy.sympify("0")
        for term in self.cooling_actions:
            if term not in ['h2formation']:
                # This is independent of the continuum optical depth
                # from the CIE
                eq += self.cooling_actions[term].equation * self.cie_optical_depth_approx()
            else:
                eq += self.cooling_actions[term].equation

        return ccode(eq, assign_to = assign_to)

    def cie_optical_depth_correction(self):
        mdensity = sympy.Symbol('mdensity')
        tau = ( mdensity/ 3.3e-8 )**2.8
        tau = sympy.Max( tau, 1e-5 )
        ciefudge = sympy.Min((1.0 - sympy.exp(-tau))/ tau, 1.0)
        return ciefudge

    def print_jacobian_component(self, s1, s2, assign_to = None, print_zeros = True):

        if s2 == self.energy_term:
            st = self.print_ccode(s1, assign_to = assign_to)
            for k in self.reactions.keys():
                k = str(k)
                st = st.replace(k + "[i]", "r" + k + "[i]")
            return st


        if s1 == self.energy_term:
            st = sum(self.cooling_actions[ca].equation
                     for ca in sorted(self.cooling_actions))
        else:
            st = self.species_total(s1)


        if assign_to is None:
            assign_to = sympy.Symbol("d_%s_%s" % (s1.name, s2.name))
        if isinstance(st, (list, tuple)):
            codes = []
            for temp_name, temp_eq in st[0]:
                teq = sympy.sympify(temp_eq)
                codes.append(ccode(teq, assign_to = temp_name))
            codes.append(ccode(st[1], assign_to = assign_to))
            return "\n".join(codes)

        eq = sympy.diff(st, s2.symbol)
        eq = eq.replace(
                lambda x: x.is_Pow and x.exp > 0 and x.exp == sympy.Integer,
                lambda x: sympy.Symbol('*'.join([x.base.name]*x.exp)) )

        if eq == sympy.sympify('0') and not print_zeros:
            return

        return ccode(eq , assign_to = assign_to)

    def get_sparse_matrix_component(self, sparse_type = "CSR", return_type = "component", assign_to = "data"  ):

        k = 0
        colvals     = []
        all_comp    = []
        rowptrs     = []
        s1_list     = []
        s2_list     = []
        i1_list     = []
        i2_list     = []

        sp_list = list( sorted(self.required_species) )
        s1_now  = sp_list[0]

        for i1, s1 in enumerate(sp_list):
            for i2, s2 in enumerate(sp_list):
                jac_comp = self.print_jacobian_component(s1, s2, print_zeros = False, assign_to = "" )

                if jac_comp:
                    colvals.append( i2 )
                    all_comp.append(jac_comp)

                    s1_list.append(s1)
                    s2_list.append(s2)

                    i1_list.append(i1)
                    i2_list.append(i2)

                    if s1 == s1_now:
                        rowptrs.append(k)
                        if i1 < len(sp_list) - 1:
                            s1_now = sp_list[i1 + 1]
                        else:
                            s1_now = None
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



    def print_JacTimesVec_component(self, s1, assign_to = None):
        """
        Compute the product of Jacobian * Vec for a given Vec
        Might be useful when we use the CVSpils solver
        """
        if s1 == self.energy_term:
            st = sum(self.cooling_actions[ca].equation
                     for ca in sorted(self.cooling_actions))

        else:
            st = self.species_total(s1)
        if assign_to is None:
            assign_to = sympy.Symbol("d_%s_dy_y" % (s1.name))
        if isinstance(st, (list, tuple)):
            codes = []
            for temp_name, temp_eq in st[0]:
                teq = sympy.sympify(temp_eq)
                codes.append(ccode(teq, assign_to = temp_name))
            codes.append(ccode(st[1], assign_to = assign_to))
            return "\n".join(codes)

        JtV_eq = sympy.sympify("0")
        mdensity = sympy.sympify("mdensity")
        T_energy = sympy.Symbol("T{0}[i]".format(self.energy_term.name) )

        i = 0
        for s2 in sorted(self.required_species):

            vec    = sympy.sympify("v{0}".format(i))
            newterm = sympy.diff(st,s2.symbol) * vec

            if s1.name == "ge":
                newterm /= mdensity
            if s2.name == "ge":
                newterm *= T_energy

            JtV_eq += newterm
            i += 1
        return ccode( JtV_eq , assign_to = assign_to)



    def print_mass_density(self):
        # Note: this assumes things are number density at this point
        eq = sympy.sympify("0")
        for s in sorted(self.required_species):
            if (s.weight > 0):
                if s.name in self.skip_weight: continue
                eq += s.symbol * s.weight
        return ccode(eq)

    def interpolate_species_gamma(self, sp, deriv=False):
        if (sp.name == 'H2_1') or (sp.name == 'H2_2') :
            expr_gammaH2 = self.species_gamma( species_registry['H2_1'], temp=True, name=False )

            if deriv is True:
                expr_dgammaH2_dT = sympy.diff(expr_gammaH2, 'T')
                f_dgammaH2_dT = lambdify('T', expr_dgammaH2_dT)
                dgammaH2_dT = f_dgammaH2_dT(self.T)
                _i1 = np.isnan(dgammaH2_dT)
                dgammaH2_dT[_i1] = 0.0

                return dgammaH2_dT
            else:
                f_gammaH2 = lambdify('T',expr_gammaH2)
                gammaH2_T =  f_gammaH2(self.T)
                _i1 = np.isnan(gammaH2_T)
                gammaH2_T[_i1] = 7./5.
                return gammaH2_T

        else:
            print('Provide your gamma function for {}'.format(sp.name) )
            raise RuntimeError


    def species_gamma(self, species, temp=False, name=True):
        if species in self.interpolate_gamma_species:
            sp_name = species.name
            T = sympy.Symbol('T')

            if temp and name:
                # so gamma enters as a function that depends on temperature
                # when gamma factor enters the temperature calculation which
                # involves derivatives of temperature (dgammaH2_dT) this term will not be dropped off
                # and be left as a function of T, which can later be supplied from the interpolated values
                # as gammaH2 does

                # this returns name of the gamma as a function of T
                # goes into the analytical differntiation for energy
                f_gammaH2 = sympy.Function('gamma%s' %sp_name)(T)
                return f_gammaH2
            elif temp and ~name:
                # x = 6100.0/T
                # expx = sympy.exp(x)
                # gammaH2_expr = 2.0 / (5.0 + 2.0*x*x*expx / (expx - 1 )**2.0 ) + 1

                T0 = T**(1/6.5)
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

                gammaH2_expr = sympy.exp(-a0*T0**a1)*(a2+T0**-a3) \
                        + a4*sympy.exp( - (T0-a5)**2 / a6)\
                        + a7*sympy.exp( -(T0-a8)**2 / a9) \
                        + 5./3.
                #x = 6100.0/T
                #expx = sympy.exp(x)
                #gammaH2_expr = 2.0 / (5.0 + 2.0*x*x*expx / (expx - 1 )**2.0 ) + 1

                return gammaH2_expr
            else:
                gammaH2 = sympy.Symbol('gamma%s' %sp_name)
                return gammaH2

        else:
            gamma = sympy.Symbol('gamma')
        return gamma


    def gamma_factor(self, temp=False):
        eq = sympy.sympify("0")
        for s in sorted(self.required_species):
            if s.name != 'ge':
                eq += (sympy.sympify(s.name)) / \
                    (self.species_gamma(s, temp=temp) - 1.0)
        return eq

    def temperature_calculation(self, derivative=False, derivative_dge_dT=False, get_dge=False):
        # If derivative=True, returns the derivative of
        # temperature with respect to ge.  Otherwise,
        # returns just the temperature function
        ge = sympy.Symbol('ge')
        function_eq = (sympy.Symbol('density') * ge * sympy.Symbol('mh')) / \
            (sympy.Symbol('kb') * (self.gamma_factor()))
        #function_eq = (ge) / \
        #    (sympy.Symbol('kb') * (self.gamma_factor()))
        if derivative == True:
            deriv_eq = sympy.diff(function_eq, ge)
            return ccode(deriv_eq)
        elif derivative_dge_dT == True:
            # when H2 presents, the gamma is dependent on  temperature
            # therefore temperature must iterate to a convergence for a given ge
            # this part evaluates the derivatives of the function ge with respect to T
            T = sympy.Symbol('T')
            f = self.gamma_factor(temp=True) * sympy.Symbol('kb') * sympy.Symbol('T') \
                    / sympy.Symbol('mh') / sympy.Symbol('density')
            dge_dT = sympy.diff(f, T)
            tmp = sympy.Symbol('tmp')
            for sp in self.interpolate_gamma_species:
                # substitute the sympy function with sympy Symbols
                sym_fgamma = sympy.Function('gamma%s' %sp.name)(T)
                sym_dfgamma = sympy.diff(sym_fgamma, T)
                dgamma = sympy.Symbol('dgamma%s_dT' %sp.name)
                dge_dT = dge_dT.subs({sym_dfgamma: dgamma})

                fgamma = sympy.Symbol('gamma%s' %sp.name)
                dge_dT = dge_dT.subs({sym_fgamma: tmp})
                dge_dT = dge_dT.subs({tmp : fgamma})

                # substitute 1/ gamma - 1 with
                # the expression _gamma{{sp}}_m1
                _gamma_m1 = sympy.Symbol('_gamma%s_m1' %sp.name)
                dge_dT = dge_dT.subs( {1/(fgamma - 1) : _gamma_m1})

            gamma = sympy.Symbol('gamma')
            _gamma_m1 = sympy.Symbol('_gamma_m1')
            dge_dT = dge_dT.subs( { 1/ (gamma-1) : _gamma_m1 } )

            # to unravel the algebric power
            dge_dT = dge_dT.replace(
                lambda x: x.is_Pow and x.exp > 0,
                lambda x: sympy.Symbol('*'.join([x.base.name]*x.exp)) )


            return ccode(dge_dT)
        elif get_dge == True:
            T = sympy.Symbol('T')
            dge = self.gamma_factor(temp=True) * sympy.Symbol('kb') * T / sympy.Symbol('mh') / sympy.Symbol('density') - sympy.Symbol('ge')

            tmp = sympy.Symbol('tmp')
            for sp in self.interpolate_gamma_species:
                sym_fgamma = sympy.Function('gamma%s' %sp.name)(T)
                fgamma = sympy.Symbol('gamma%s' %sp.name)
                dge = dge.subs({sym_fgamma: tmp})
                dge = dge.subs({tmp: fgamma})

                # substitute 1/ gamma - 1 with
                # the expression _gamma{{sp}}_m1
                _gamma_m1 = sympy.Symbol('_gamma%s_m1' %sp.name)
                dge = dge.subs( {1/(fgamma - 1) : _gamma_m1})

            gamma = sympy.Symbol('gamma')
            _gamma_m1 = sympy.Symbol('_gamma_m1')
            dge = dge.subs( { 1/ (gamma-1) : _gamma_m1 } )

            # to unravel the algebric power
            dge = dge.replace(
                lambda x: x.is_Pow and x.exp > 0,
                lambda x: sympy.Symbol('*'.join([x.base.name]*x.exp)) )


            return ccode(dge)
        else:
            gamma = sympy.Symbol('gamma')
            _gamma_m1 = sympy.Symbol('_gamma_m1')
            tmp = sympy.Symbol('tmp')
            T = sympy.Symbol('T')
            for sp in self.interpolate_gamma_species:
                # substitute the sympy function with sympy Symbols
                sym_fgamma = sympy.Function('gamma%s' %sp.name)(T)
                _fgamma_m1 = sympy.Symbol('_gamma%s_m1' %sp.name)
                function_eq = function_eq.subs({ 1 / (sym_fgamma - 1) : tmp})
                function_eq = function_eq.subs({tmp : _fgamma_m1})

            function_eq = function_eq.subs( { 1/ (gamma-1) : _gamma_m1 } )
            gamma = sympy.Symbol('gamma')
            _gamma_m1 = sympy.Symbol('_gamma_m1')

            return ccode(function_eq)

    # This function computes the total number density
    def calculate_number_density(self, values, skip = ()):
        # values should be a dict with all of the required species in it
        # The values should be in *mass* density
        n = np.zeros_like(values.values()[0])
        for s in self.required_species:
            if s.name in self.skip_weight: continue
            n += values[s.name] / s.weight
        return n

    # This function counts up the total number of free electrons
    def calculate_free_electrons(self, values):
        # values should be a dict with all of the required species in it
        # The values should be in *mass* density
        n = np.zeros_like(values.values()[0])
        for s in self.required_species:
            if s.name in self.skip_weight: continue
            n += ( values[s.name] / s.weight ) * s.free_electrons
        return n

    # This computes the total mass density from abundance fractions
    def calculate_mass_density(self, values):
        # values should be a dict with all of the required species in it
        # The values should be in *mass* density
        n = np.zeros_like(values.values()[0])
        for s in self.required_species:
            if s.name in self.skip_weight: continue
            n += values[s.name] * s.weight
        return n

    # This function sums the densities (mass or number depending on what
    # is fed in) of non-electron species
    def calculate_total_density(self, values):
        # values should be a dict with all of the required species in it
        # The values should be in *mass* density
        n = np.zeros_like(values.values()[0])
        for s in self.required_species:
            if s.name in self.skip_weight: continue
            n += values[s.name]
        return n

    # This function converts from fractional abundance to mass density
    def convert_to_mass_density(self, values, skip = ()):
        for s in self.required_species:
            if s.name in self.skip_weight: continue
            values[s.name] = values[s.name] * s.weight
        return values

    def write_solver(self, solver_name,
                     solver_template = "rates_and_rate_tables",
                     ode_solver_source = "BE_chem_solve.C",
                     output_dir = ".", init_values = None,
                     main_name = "main",
                     input_is_number = False):
        self.input_is_number = input_is_number

        if not os.path.isdir(output_dir): os.makedirs(output_dir)
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
        env = jinja2.Environment(extensions=['jinja2.ext.loopcontrols'],
                loader = jinja2.PackageLoader("dengo", "templates"))
        template_vars = dict(network = self, solver_name = solver_name)


        for suffix in (".C", "_main.C", ".h", "_run.pyx", "_run.pyxbld",
                       "_run.pyxdep", "_main.py"):
            iname = "%s%s" % (solver_template, suffix)
            oname = os.path.join(output_dir,
                            "%s_solver%s" % (solver_name, suffix))
            template_inst = env.get_template(iname + ".template")
            solver_out = template_inst.render(**template_vars)
            with open(oname ,"w") as f:
                f.write(solver_out)

        # now we create a Makefile
        temp_dir, _ = os.path.split(solver_template)
        iname = os.path.join(temp_dir, "Makefile")
        oname = os.path.join(output_dir, "Makefile")
        template_inst = env.get_template(iname + ".template")
        solver_out = template_inst.render(**template_vars)
        with open(oname ,"w") as f:
            f.write(solver_out)


        env = jinja2.Environment(extensions=['jinja2.ext.loopcontrols'],
                loader = jinja2.PackageLoader("dengo", "solvers"))

        # Now we copy over anything else we might need.
        if ode_solver_source is not None:

            #iname = os.path.join("solvers", ode_solver_source)
            #src = pkgutil.get_data("dengo", iname)
            #src = src.decode("utf-8")

            template_inst = env.get_template(ode_solver_source)
            solver_out = template_inst.render(**template_vars)

            with open(os.path.join(output_dir, ode_solver_source), "w") as f:
                f.write(solver_out)
        if solver_template.endswith(".c.template"):
            hfn = solver_template.rsplit(".", 2)[0] + ".h"
            src = pkgutil.get_data("dengo", os.path.join("templates", hfn))
            with open(os.path.join(output_dir, hfn), "w") as f:
                f.write(src)
        # This writes out the rates for the species in the
        # chemical network to HDF5 files which can later be
        # read by the C++ code that is output by the template
        ofn = os.path.join(output_dir, "%s_tables.h5" % solver_name)
        f = h5py.File(ofn, "w")

        # construct a 2d array of reaction rates
        # reaction_rates[ T ][ k0i ]
        all_rates = []
        for rxn in sorted(self.reactions.values()):
            data = rxn.coeff_fn(self).astype("float64")
            all_rates.append(data)
        all_rates = np.array(all_rates).transpose().flatten()
        f.create_dataset("/all_reaction_rates" , data = all_rates)

        # construct a 2d array of cooling rates
        # cooling_rates [ T ][ cooling ]
        all_cooling = []
        for action in sorted(self.cooling_actions.values()):
            for tab in sorted(action.tables):
                data = action.tables[tab](self).astype("float64")
                all_cooling.append(data)
        all_cooling = np.array(all_cooling).transpose().flatten()
        f.create_dataset("/all_cooling_rates", data= all_cooling)

        for rxn in sorted(self.reactions.values()):
            f.create_dataset("/%s" % rxn.name, data = rxn.coeff_fn(self).astype("float64"))
        for action in sorted(self.cooling_actions.values()):
            for tab in action.tables:
                f.create_dataset("/%s_%s" % (action.name, tab),
                    data=action.tables[tab](self).astype("float64"))

        for sp in sorted(self.interpolate_gamma_species):
            name = sp.name
            f.create_dataset("/gamma%s" %name,
                    data = self.interpolate_species_gamma(sp).astype("float64"))
            f.create_dataset("/dgamma%s_dT" %name,
                    data = self.interpolate_species_gamma(sp,deriv=True).astype("float64"))
        f.close()


        if init_values is None: return

        # This write outs a set of initial conditions for to be fed
        # into C++ code for testing purposes
        ofn = os.path.join(output_dir, "%s_initial_conditions.h5" % solver_name)
        f = h5py.File(ofn, "w")
        for name, init_value in init_values.items():
            f.create_dataset("/%s" % name, data=init_value.astype('float64'))
        f.close()
