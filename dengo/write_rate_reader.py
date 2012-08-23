from chemistry_constants import tiny, kboltz, mh
import jinja2
import h5py
import numpy as na

def create_rate_tables(network, solver_name):
    f = h5py.File("%s_rate_tables.h5" % solver_name, "w")
    for rxn in sorted(network.reactions.values()):
        f.create_dataset("/%s" % rxn.name, data = rxn.coeff_fn(network.T).astype("float64"))
    f.close()
    f = h5py.File("%s_cooling_tables.h5" % solver_name, "w")
    for action in sorted(network.cooling_actions.values()):
        f.create_dataset("/%s" % action.name, data = action.tables[action.name](network.T).astype("float64"))
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
    rate_table = network.reactions
    species_table = network.required_species
    cooling_rate_table = network.cooling_actions
    
    ireaction_table = dict([(rid, rname) 
            for rid, rname in enumerate(sorted(rate_table.values()))])
    reaction_ids = dict([(a, b) for b, a in ireaction_table.items()])
    reaction_varnames = dict([(a, "r_%02i" % i)
            for i, a in enumerate(sorted(rate_table.values()))])
    irate_table = dict([(rid, rname) 
            for rid, rname in enumerate(sorted(rate_table.values()))])
    rate_ids = dict([(a, b) for b, a in irate_table.items()])

    # Now the cooling stuff
    icooling_rate_table = dict([(cid, cname) 
            for cid, cname in enumerate(sorted(cooling_rate_table))])
    cooling_rate_ids = dict([(a, b) for b, a in icooling_rate_table.items()])

    env = jinja2.Environment(extensions=['jinja2.ext.loopcontrols'],
            loader = jinja2.FileSystemLoader(["enzo_templates/","."]))
    solver_template = env.get_template(
        "enzo_templates/%s_data.c.template" % (solver_name))
    template_vars = dict(network = network,
                         rate_ids = rate_ids,
                         solver_name = solver_name,
                         cooling_rate_ids = cooling_rate_ids)
    solver_out = solver_template.render(**template_vars)
    f = open("enzo_templates/%s_data.c" % solver_name, "w")
    f.write(solver_out)
