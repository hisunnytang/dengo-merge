from dengo.chemical_network import \
    ChemicalNetwork, \
    reaction_registry, \
    cooling_registry
import dengo.primordial_rates, dengo.primordial_cooling
from dengo.write_rate_reader import \
    create_rate_tables, \
    create_rate_reader, \
    create_initial_conditions
from dengo.chemistry_constants import tiny, kboltz, mh
from dengo.known_species import *

# If only a subset of species are wanted put them here
# and change the commented lines below
want = ("HI", "HII", "de", "ge")

primordial = ChemicalNetwork()
primordial.add_energy_term()

# Set to false if intermediate solution output is not wanted
primordial.write_intermediate_solutions = True

for ca in cooling_registry.values():
    #if not all(sp.name in want for sp in ca.species): continue
    primordial.add_cooling(ca)


for i, rname in enumerate(sorted(reaction_registry)):
    s = reaction_registry[rname]
    #if not all(sp.name in want for sp in s.species): continue
    primordial.add_reaction(s)

# This defines the temperature range for the rate tables
primordial.init_temperature((1e0, 1e8))

# Write the rate tables and corresponding C++ code
create_rate_tables(primordial, "primordial")
create_rate_reader(primordial, "primordial")

# Generate initial conditions (switch to False to disable this)
generate_initial_conditions = True

if generate_initial_conditions:
    import numpy as na
    NCELLS = 1
    density = 1e7
    init_array = na.ones(NCELLS) 
    X = 0.5

    init_values = dict()
    init_values['density'] = density * init_array
    init_values['HI'] = init_array.copy() # use conservation to set this below
    # populate initial fractional values for the other species
    for s in sorted(primordial.required_species):
        if s.name != 'ge' and s.name != 'HI':
            if s.name == 'de':
                continue
            elif s.name == 'HII':
                init_values[s.name] = X * init_array
            else:
                init_values[s.name] = tiny * init_array
            init_values['HI'] -= init_values[s.name]
    init_values['de'] = init_array * 0.0
    for s in sorted(primordial.required_species):
        if s.name in ('ge', 'de'): continue
        init_values['de'] += init_values[s.name] * s.free_electrons
        print "Adding %0.5e to electrons from %s" % (
            (init_values[s.name] * s.free_electrons/s.weight)[0], s.name)
    print "Total de: %0.5e" % (init_values['de'][0])

    # convert to masses to multiplying by the density factor and the species weight
    for s in primordial.required_species:
        if s.name == 'ge': continue
        if s.name == 'de':
            init_values[s.name] *= (density) # this is still just number density
        else:
            init_values[s.name] *= (density * s.weight)

    #compute new total density and number density
    density = 0.0
    number_density = 0.0
    for s in sorted(primordial.required_species):
        if s.name == 'ge': continue
        number_density += init_values[s.name][0]/s.weight
        if s.name == 'de': continue
        density += init_values[s.name][0]
    
    # set up initial temperatures values used to define ge
    temperature = na.logspace(2, 4, NCELLS)
    temperature[:] = 2e2
    init_values['T'] = temperature

    # calculate ge (very crudely)
    gamma = 5.e0/3.e0
    init_values['ge'] = ((temperature * density * kboltz)
                         / (number_density * mh * (gamma - 1)))
    
    # Write the initial conditions file
    create_initial_conditions(init_values, 'primordial')
