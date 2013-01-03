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
from chemistry_constants import tevk, tiny, mh
from .reaction_classes import reaction_registry, cooling_registry, \
    count_m, index_i, species_registry
import types
import sympy
import pkgutil
import os
import jinja2
import h5py
from sympy.printing import ccode

class ChemicalNetwork(object):

    energy_term = None

    def __init__(self):
        self.reactions = {}
        self.cooling_actions = {}
        self.required_species = set([])
        self.write_intermediate_solutions = False

    def add_reaction(self, reaction):
        if isinstance(reaction, types.StringTypes):
            reaction = reaction_registry[reaction]
        if reaction.name in self.reactions:
            raise RuntimeError
        self.reactions[reaction.name] = reaction
        
        for n, s in reaction.left_side:
            self.required_species.add(s)
        for n, s in reaction.right_side:
            self.required_species.add(s)
        print "Adding reaction: %s" % reaction
    def add_energy_term(self):
        if self.energy_term is None:
            self.energy_term = species_registry["ge"]
            self.required_species.add(self.energy_term)

    def add_cooling(self, cooling_term):
        self.add_energy_term()
        if cooling_term.name in self.cooling_actions:
            raise RuntimeError
        self.cooling_actions[cooling_term.name] = cooling_term
        self.required_species.update(cooling_term.species)
    
    def init_temperature(self, T_bounds = (1, 1e8), n_bins=1024):
        self.n_bins = 1024
        self.T = np.logspace(np.log(T_bounds[0]),
                             np.log(T_bounds[1]),
                             n_bins, base = np.e)
        self.logT = np.log(self.T)
        self.tev = self.T / tevk
        self.logtev = np.log(self.tev)
        self.T_bounds = T_bounds

    def species_total(self, species):
        eq = sympy.sympify("0")
        for rn, rxn in sorted(self.reactions.items()):
            rxn.update(self.energy_term.symbol)
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
        return ccode(self.species_total(species), assign_to = assign_to)

    def print_cooling(self, assign_to):
        eq = sympy.sympify("0")
        for term in self.cooling_actions:
            eq += self.cooling_actions[term].equation
        return ccode(eq, assign_to = assign_to)

    def print_jacobian_component(self, s1, s2, assign_to = None):
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
        return ccode(sympy.diff(st, s2.symbol), assign_to = assign_to)

    def print_number_density(self):
        eq = sympy.sympify("0")
        for s in sorted(self.required_species):
            if s.name != 'ge' and s.name != 'de': 
                eq += s.symbol
        return ccode(eq)

    def print_mass_density(self):
        eq = sympy.sympify("0")
        for s in sorted(self.required_species):
            if s.name not in ('ge', 'de'): 
                eq += s.symbol * s.weight
        return ccode(eq)

    def species_gamma(self, species):
        if species.name == 'H2I' or species.name == 'H2II':
            gamma = sympy.Symbol('gammaH2')
        else:
            gamma = sympy.Symbol('gamma')
        return gamma

    def gamma_factor(self):
        eq = sympy.sympify("0")
        for s in sorted(self.required_species):
            if s.name != 'ge':
                eq += (sympy.sympify(s.name)) / \
                    (self.species_gamma(s) - 1.0)
        return eq

    def temperature_calculation(self, derivative=False):
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
        else:
            return ccode(function_eq)
        return ccode(eq)

    def write_solver(self, solver_name,
                     solver_template = "rates_and_rate_tables.c.template",
                     ode_solver_source = "BE_chem_solve.C",
                     output_dir = ".", init_values = None):
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
        
        template_inst = env.get_template(solver_template)
        template_vars = dict(network = self, solver_name = solver_name)
        solver_out = template_inst.render(**template_vars)
        f = open(os.path.join(output_dir, "%s_solver.c" % solver_name), "w")
        f.write(solver_out)
        f.close()

        # Now we copy over anything else we might need.
        if ode_solver_source is not None:
            src = pkgutil.get_data("dengo", os.path.join("solvers", ode_solver_source))
            with open(os.path.join(output_dir, ode_solver_source), "w") as f:
                f.write(src)
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
        for rxn in sorted(self.reactions.values()):
            f.create_dataset("/%s" % rxn.name, data = rxn.coeff_fn(self).astype("float64"))
        for action in sorted(self.cooling_actions.values()):
            for tab in action.tables:
                f.create_dataset("/%s_%s" % (action.name, tab),
                    data=action.tables[tab](self).astype("float64"))
        f.close()

        if init_values is None: return

        # This write outs a set of initial conditions for to be fed
        # into C++ code for testing purposes
        ofn = os.path.join(output_dir, "%s_initial_conditions.h5" % solver_name)
        f = h5py.File(ofn, "w")
        for name, init_value in init_values.items():
            f.create_dataset("/%s" % name, data=init_value.astype('float64'))
        f.close()
