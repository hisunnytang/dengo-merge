from chemistry_constants import tiny, kboltz, mh
import jinja2
import h5py
import numpy as na
import os

def create_rate_tables(network, solver_name):
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

def create_rate_reader(network, solver_name): #, pp_class):
    # This comment block is currently out of date -- Devin 7/27/2012
    # What we are handed here is:
    #   * rate_table, which is a dict of "kXX" to ReactionRate objects.  These
    #     will be used in the reading of rates from disk, but will not directly
    #     be used in the calculation of the derivatives.
    #   * reaction_table, which is a dict of "rXX" to Reaction objects.  These
    #     objects have left_side, right_side, and considered attributes that
    #     describe which species are members of each and how they contribute.
    #     Note that *left_side* and *right_side* are tuples of (n, S), where S
    #     is an instance of Species, but that *considered* is a set of strings.
    #   * species_table, which is a dict of "H2I" (etc) to Species objects.
    #     The Species objects are notable for having information about whether
    #     they are in *equilibrium*, as well as their atomic weight.
    #
    # To utilize these inside our template, we will generate convenience
    # handlers that will explicitly number them.
    env = jinja2.Environment(extensions=['jinja2.ext.loopcontrols'],
            loader = jinja2.FileSystemLoader(["templates/","."]))
    solver_template = env.get_template(
        "templates/rates_and_rate_tables.c.template")
    # template_vars = dict(num_solved_species = num_solved_species,
    #                      num_total_species = num_total_species,
    template_vars = dict(network = network, solver_name = solver_name)
    #template_vars['pp'] = pp_class(template_vars)
    solver_out = solver_template.render(**template_vars)
    if not os.path.isdir("output"): os.makedirs("output")
    f = open("output/%s_solver.c" % solver_name, "w")
    f.write(solver_out)
