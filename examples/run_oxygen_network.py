from dengo.chemical_network import \
    ChemicalNetwork, \
    reaction_registry, \
    cooling_registry
import dengo.primordial_rates, dengo.primordial_cooling
import dengo.oxygen_rates, dengo.oxygen_cooling
from dengo.write_rate_reader import \
    create_rate_tables, \
    create_rate_reader, \
    create_initial_conditions
from dengo.chemistry_constants import tiny, kboltz, mh
import os

# If only a subset of species are wanted put them here
# and change the commented lines below
want = ('o_5', 'o_6', 'o_7', 'de', 'ge')

oxygen = ChemicalNetwork()
oxygen.add_energy_term()

# for ca in cooling_registry.values():
#     # The following line can be used to specify a subset of species
#     #if all(sp.name in want for sp in ca.species):
#     if ca.name.startswith("o_"):
#        oxygen.add_cooling(ca)

for s in reaction_registry.values():
    # The following line can be used to specify a subset of species
    #if all(sp.name in want for sp in s.species):
    if s.name.startswith("o_"):
        oxygen.add_reaction(s)

# This defines the temperature range for the rate tables
oxygen.init_temperature((1e4, 1e8))

# Set to false if you don't want intermediate solution output
oxygen.write_intermediate_solutions = True

# Write the rate tables and the corresponding C++ code
create_rate_tables(oxygen, "oxygen")
create_rate_reader(oxygen, "oxygen")

# Generate initial conditions (switch to False to disable this)
generate_initial_conditions = True

if generate_initial_conditions:
    import numpy as na
    NCELLS = 1
    density = 1.0
    init_array = na.ones(NCELLS) 
    X = 1e-6

    init_values = dict()
    init_values['density'] = density * init_array
    initial_state = 'o_1'
    init_values[initial_state] = init_array.copy() # use conservation to set this below
    # populate initial fractional values for the other species
    for s in sorted(oxygen.required_species):
        if s.name != 'ge' and s.name != initial_state:
            if s.name == 'de':
                continue
            else:
                init_values[s.name] = X * init_array
            init_values[initial_state] -= init_values[s.name]
    init_values['de'] = init_array * 0.0
    for s in sorted(oxygen.required_species):
        if s.name in ("ge", "de"): continue
        init_values['de'] += init_values[s.name] * s.free_electrons
        print "Adding %0.5e to electrons from %s" % (
            (init_values[s.name] * s.free_electrons)[0], s.name)
    print "Total de: %0.5e" % (init_values['de'][0])

    # convert to masses to multiplying by the density factor and the species weight
    for s in oxygen.required_species:
        if s.name == 'ge': continue
        if s.name == 'de':
            init_values[s.name] *= (density) # this is still just number density
        else:
            init_values[s.name] *= (density * s.weight)
    
    #compute new total density and number density
    density = 0.0
    number_density = 0.0
    for s in sorted(oxygen.required_species):
        if s.name == 'ge': continue
        number_density += init_values[s.name][0]/s.weight
        if s.name == 'de': continue
        density += init_values[s.name][0]

    # set up initial temperatures values used to define ge
    temperature = na.logspace(4, 6.7, NCELLS)
    temperature[:] = 1e7; # need to remove this line for the above one to matter
    init_values['T'] = temperature

    # calculate ge (very crudely)
    gamma = 5.e0/3.e0
    init_values['ge'] = ((temperature * number_density * kboltz)
                         / (density * mh * (gamma - 1)))
    
    # Write the initial conditions file
    create_initial_conditions(init_values, 'oxygen')

    includepath = "/Users/devinsilvia/Research/code/yt-x86_64/include/"
    enzoincludepath = "/Users/devinsilvia/Research/code/enzo-dev_chemistry/src/enzo/"
    libpath = "/Users/devinsilvia/Research/code/yt-x86_64/lib/"
    os.system("mpic++ -O0 -ggdb -DLINUX -DH5_USE_16_API -DUSE_SQLITE -DSAB -D__max_subgrids=100000 -D__max_baryons=20 -D__max_cpu_per_node=8 -D__memory_pool_size=100000 -DINITS64 -DSMALL_INTS -DCONFIG_PINT_8 -DIO_32 -DUSE_MPI -DCONFIG_PFLOAT_8 -DCONFIG_BFLOAT_8 -DUSE_HDF5_GROUPS -DTRANSFER -DNEW_GRID_IO -DFAST_SIB -DBITWISE_IDENTICALITY -DFLUX_FIX -DSET_ACCELERATION_BOUNDARY -I%s -I%s -L%s -lhdf5_hl -lhdf5 -lm output/oxygen_solver.c templates/BE_chem_solve.C -o output/oxygen.out" %(includepath, enzoincludepath, libpath))
    os.chdir("output")
    os.system("ipython clean.py")
    os.system("./oxygen.out >& estd.out")
    os.system("ipython ../plot_intermediate.py oxygen")
