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

    # Note here that we have three mechanisms for numbering species.  We have
    # our equilibrium species, which are not calculated in the derivative
    # calculation (this should include charge conservation calculations, as
    # well) and we have our actual local variable names.  These variables are
    # local to the cell, so this may give us better cache performance, and it
    # allows us to address things more cleanly.
    #   * non_eq_species_table: maps species name to Species object
    #   * eq_species_table: maps species name to Species object
    #   * species_varnames: maps species name to local variable name
    # non_eq_species_table = dict([ (a, b)
    #         for a, b in species_table.items() if not b.equilibrium
    #                 and not b.computed])
    # non_eq_species_ids = dict([ (a, b)
    #         for b, a in enumerate(sorted(non_eq_species_table))])
    # eq_species_table = dict([ (a, b)
    #         for a, b in species_table.items() if b.equilibrium 
    #                 and not b.computed])
    # eq_species_ids = dict([ (a, b)
    #         for b, a in enumerate(sorted(eq_species_table))])
    # species_varnames = dict([(a, "s_%02i" % i)
    #         for i, a in enumerate(sorted(species_table))])
    # num_solved_species = len(non_eq_species_table)
    # num_total_species = len(species_varnames)

    # Now the cooling stuff
    icooling_rate_table = dict([(cid, cname) 
            for cid, cname in enumerate(sorted(cooling_rate_table))])
    cooling_rate_ids = dict([(a, b) for b, a in icooling_rate_table.items()])

    env = jinja2.Environment(extensions=['jinja2.ext.loopcontrols'],
            loader = jinja2.FileSystemLoader(["enzo_templates/","."]))
    solver_template = env.get_template(
        "enzo_templates/%s_data.c.template" % (solver_name))
    # template_vars = dict(num_solved_species = num_solved_species,
    #                      num_total_species = num_total_species,
    template_vars = dict(network = network,
                         # rate_table = rate_table,
                         rate_ids = rate_ids,
                         # irate_table = irate_table, 
                         # reaction_table = reaction_table,
                         # reaction_ids = reaction_ids,
                         # reaction_varnames = reaction_varnames,
                         # ireaction_table = ireaction_table,
                         solver_name = solver_name,
                         # species_table = species_table,
                         # non_eq_species_table = non_eq_species_table,
                         # non_eq_species_ids = non_eq_species_ids,
                         # eq_species_table = eq_species_table,
                         # eq_species_ids = eq_species_ids,
                         # species_varnames = species_varnames,
                         # cooling_rate_table = cooling_rate_table,
                         cooling_rate_ids = cooling_rate_ids)
                         # icooling_rate_table = icooling_rate_table)
                         # cooling_action_table = cooling_action_table)
    #template_vars['pp'] = pp_class(template_vars)
    solver_out = solver_template.render(**template_vars)
    f = open("enzo_templates/%s_data.c" % solver_name, "w")
    f.write(solver_out)
