from chemistry_constants import tiny, kboltz, mh
import jinja2
import h5py
import numpy as na
import os

def create_rate_tables(network, solver_name):
    # This writes out the rates for the species in the
    # chemical network to HDF5 files which can later be
    # read by the C++ code that is output by the template
    if not os.path.isdir("output"): os.makedirs("output")
    f = h5py.File("output/%s_rate_tables.h5" % solver_name, "w")
    for rxn in sorted(network.reactions.values()):
        f.create_dataset("/%s" % rxn.name, data = rxn.coeff_fn(network).astype("float64"))
    f.close()
    f = h5py.File("output/%s_cooling_tables.h5" % solver_name, "w")
    for action in sorted(network.cooling_actions.values()):
        for tab in action.tables:
            f.create_dataset("/%s_%s" % (action.name, tab),
                data=action.tables[tab](network).astype("float64"))
    f.close()

def create_rate_reader(network, solver_name):
    # What we are handed here is:
    #   * network, a python object which holds all of the species, reactions,
    #     rate, etc. that we're keeping track of and will be solving in Enzo
    #   * solver_name, an identifier to produce a unique template and to
    #     correctly grab the right HDF5 tables
    #
    # To utilize these inside our template, we will generate convenience
    # handlers that will explicitly number them.
    env = jinja2.Environment(extensions=['jinja2.ext.loopcontrols'],
            loader = jinja2.FileSystemLoader(["templates/","."]))
    solver_template = env.get_template(
        "templates/rates_and_rate_tables.c.template")
    template_vars = dict(network = network, solver_name = solver_name)
    solver_out = solver_template.render(**template_vars)
    if not os.path.isdir("output"): os.makedirs("output")
    f = open("output/%s_solver.c" % solver_name, "w")
    f.write(solver_out)

def create_initial_conditions(init_values, solver_name):
    # This write outs a set of initial conditions for to be fed
    # into C++ code for testing purposes
    f = h5py.File("output/%s_initial_conditions.h5" % solver_name, "w")
    for name, init_value in init_values.items():
        f.create_dataset("/%s" % name, data=init_value.astype('float64'))
    f.close()
