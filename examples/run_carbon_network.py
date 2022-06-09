from dengo.chemical_network import \
    ChemicalNetwork, \
    reaction_registry, \
    cooling_registry
import dengo.primordial_rates, dengo.primordial_cooling
import dengo.carbon_rates, dengo.carbon_cooling
from dengo.write_rate_reader import \
    create_rate_tables, \
    create_rate_reader, \
    create_initial_conditions
from dengo.chemistry_constants import tiny, kboltz, mh
import os

# If only a subset of species are wanted put them here
# and change the commented lines below
want = ('CIII', 'CIV', 'CV', 'de', 'ge')

carbon = ChemicalNetwork()
carbon.add_energy_term()

# for ca in cooling_registry.values():
#     # The following line can be used to specify a subset of species
#     #if all(sp.name in want for sp in ca.species):
#     if ca.name.startswith("C"):
#        carbon.add_cooling(ca)

for s in reaction_registry.values():
    # The following line can be used to specify a subset of species
    #if all(sp.name in want for sp in s.species):
    if s.name.startswith("C"):
        carbon.add_reaction(s)

# This defines the temperature range for the rate tables
carbon.init_temperature((1e4, 1e8))

# Set to false if you don't want intermediate solution output
carbon.write_intermediate_solutions = True

# Write the rate tables and the corresponding C++ code
create_rate_tables(carbon, "carbon")
create_rate_reader(carbon, "carbon")

# Generate initial conditions (switch to False to disable this)
generate_initial_conditions = True
start_neutral = False
initial_state = 'CVII' # if using start_neutral = True, this should be the _1 state
                      # if starting in a close to CIE state, this should be the
                      # ion species that will be closest for one for the given Temp.
                      # note, this will get hairy when we're not testing a single T

if generate_initial_conditions:
    import numpy as np
    NCELLS = 1
    density = 1
    init_array = np.ones(NCELLS)
    init_values = dict()
    init_values['density'] = density * init_array
    init_values[initial_state] = init_array.copy() # use conservation to set this below

    # set up initial temperatures values used to define ge
    temperature = np.logspace(4, 6.7, NCELLS)
    temperature[:] = 1e7; # need to remove this line for the above one to matter
    init_values['T'] = temperature

    if start_neutral:
        X = 1e-6

        # populate initial fractional values for the other species
        for s in sorted(carbon.required_species):
            if s.name != 'ge' and s.name != initial_state:
                if s.name == 'de':
                    continue
                else:
                    init_values[s.name] = X * init_array
                init_values[initial_state] -= init_values[s.name]
        init_values['de'] = init_array * 0.0
        for s in sorted(carbon.required_species):
            if s.name in ("ge", "de"): continue
            init_values['de'] += init_values[s.name] * s.free_electrons
            print ("Adding %0.5e to electrons from %s" %
                (init_values[s.name] * s.free_electrons)[0], s.name))
        print("Total de: %0.5e" % (init_values['de'][0]))

    else:
        # Unless neutral is desired, we'll start in a perturbed CIE solution
        import chianti.core as ch
        import chianti.util as chu

        # populate initial fractional values for the species
        for s in sorted(carbon.required_species):
            if s.name != 'ge' and s.name != initial_state:
                if s.name == 'de':
                    continue
                else:
                    ion_name = chu.zion2name(s.number, s.free_electrons + 1)
                    ion = ch.ion(ion_name, temperature=init_values['T'])
                    ion.ioneqOne()
                    ion_frac = ion.IoneqOne
                    if ion_frac == 0.0:
                        init_values[s.name] = tiny * init_array
                    else:
                        ion_frac[ion_frac < tiny] = tiny
                        init_values[s.name] = ion_frac * init_array

                # add some random noise
                init_values[s.name] += 0.5 * init_values[s.name] * np.random.randn(NCELLS)
                # in case something went negative:
                init_values[s.name][init_values[s.name] < tiny] = tiny

                # conservation...
                init_values[initial_state] -= init_values[s.name]

        init_values[initial_state][init_values[initial_state] < tiny] = tiny
        init_values['de'] = init_array * 0.0
        for s in sorted(carbon.required_species):
            if s.name in ("ge", "de"): continue
            init_values['de'] += init_values[s.name] * s.free_electrons
            print ("Adding %0.5e to electrons from %s" % (
                (init_values[s.name] * s.free_electrons)[0], s.name))
        print ("Total de: %0.5e" % (init_values['de'][0]))

    # convert to masses to multiplying by the density factor and the species weight
    for s in carbon.required_species:
        if s.name == 'ge': continue
        if s.name == 'de':
            init_values[s.name] *= (density) # this is still just number density
        else:
            init_values[s.name] *= (density * s.weight)

    #compute new total density and number density
    density = 0.0
    number_density = 0.0
    for s in sorted(carbon.required_species):
        if s.name == 'ge': continue
        number_density += init_values[s.name][0]/s.weight
        if s.name == 'de': continue
        density += init_values[s.name][0]

    # calculate ge (very crudely)
    gamma = 5.e0/3.e0
    init_values['ge'] = ((temperature * number_density * kboltz)
                         / (density * mh * (gamma - 1)))

    # Write the initial conditions file
    create_initial_conditions(init_values, 'carbon')

    includepath = "/Users/devinsilvia/Research/code/yt-x86_64/include/"
    enzoincludepath = "/Users/devinsilvia/Research/code/enzo-dev_chemistry/src/enzo/"
    libpath = "/Users/devinsilvia/Research/code/yt-x86_64/lib/"
    os.system("mpic++ -O0 -ggdb -DLINUX -DH5_USE_16_API -DUSE_SQLITE -DSAB -D__max_subgrids=100000 -D__max_baryons=20 -D__max_cpu_per_node=8 -D__memory_pool_size=100000 -DINITS64 -DSMALL_INTS -DCONFIG_PINT_8 -DIO_32 -DUSE_MPI -DCONFIG_PFLOAT_8 -DCONFIG_BFLOAT_8 -DUSE_HDF5_GROUPS -DTRANSFER -DNEW_GRID_IO -DFAST_SIB -DBITWISE_IDENTICALITY -DFLUX_FIX -DSET_ACCELERATION_BOUNDARY -I%s -I%s -L%s -lhdf5_hl -lhdf5 -lm output/carbon_solver.c templates/BE_chem_solve.C -o output/carbon.out" %(includepath, enzoincludepath, libpath))
    os.chdir("output")
    os.system("ipython clean.py")
    os.system("./carbon.out >& estd.out")
    os.system("ipython ../plot_intermediate.py carbon")
