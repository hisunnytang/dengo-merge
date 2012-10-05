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

oxygen = ChemicalNetwork()
for ca in cooling_registry.values():
    if ca.name.startswith("o_"):
        oxygen.add_cooling(ca)

for s in reaction_registry.values():
    if s.name.startswith("o_"):
        oxygen.add_reaction(s)

oxygen.init_temperature((1e0, 1e8))

create_rate_tables(oxygen, "oxygen")
create_rate_reader(oxygen, "oxygen")

# Generate initial conditions (switch to False to disable this)
generate_initial_conditions = True

if generate_initial_conditions:
    import numpy as na
    NCELLS = 64
    density = 1.0e4
    init_array = na.ones(NCELLS) 
    X = 1.0/9.0

    init_values = dict()
    init_values['density'] = density * init_array
    init_values['o_1'] = init_array.copy() # use conservation to set this below
    # populate initial fractional values for the other species
    for s in sorted(oxygen.required_species):
        if s.name != 'ge' and s.name != 'o_1':
            if s.name == 'de':
                continue
            if s.name == 'o_2':
                init_values[s.name] = X * init_array
            else:
                init_values[s.name] = X * init_array
            init_values['o_1'] -= init_values[s.name]
    init_values['de'] = init_array * 0.0
    for s in sorted(oxygen.required_species):
        if s.name == "ge": continue
        init_values['de'] += init_values[s.name] * s.free_electrons / s.weight
    #print init_values['de'][0] / (init_values['o_2'][0] / 16.0)

    # convert to masses to multiplying by the density factor and the species weight
    for s in oxygen.required_species:
        if s.name == 'ge': continue
        if s.name == 'de':
            init_values[s.name] *= (density)
        else:
            init_values[s.name] *= (density * s.weight)
    
    #compute new total density and number density
    density = 0.0
    number_density = 0.0
    for s in sorted(oxygen.required_species):
        if s.name == 'ge': continue
        density += init_values[s.name][0]
        number_density += init_values[s.name][0]/s.weight

    # set up initial temperatures values used to define ge
    temperature = na.logspace(4, 8, NCELLS)

    # calculate ge (very crudely)
    gamma = 5.e0/3.e0
    init_values['ge'] = ((temperature * number_density * kboltz)
                         / (density * mh * (gamma - 1)))
    
    create_initial_conditions(init_values, 'oxygen')
