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
oxygen.add_energy_term()
for ca in cooling_registry.values():
    #if ca.name in ('o_1', 'o_2', 'de'):
    if ca.name.startswith("o_"):
        oxygen.add_cooling(ca)
        print ca

for s in reaction_registry.values():
    #if all(sp.name in ('o_1', 'o_2', 'de', 'ge') for sp in s.species):
    if s.name.startswith("o_"):
        oxygen.add_reaction(s)
        print s

oxygen.init_temperature((1e0, 1e7))
oxygen.write_intermediate_solutions = True


create_rate_tables(oxygen, "oxygen")
create_rate_reader(oxygen, "oxygen")

# Generate initial conditions (switch to False to disable this)
generate_initial_conditions = True

if generate_initial_conditions:
    import numpy as na
    NCELLS = 4
    density = 1.0
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
        if s.name in ("ge", "de"): continue
        init_values['de'] += init_values[s.name] * s.free_electrons# / s.weight
        print "Adding %0.5e to electrons from %s" % (
            (init_values[s.name] * s.free_electrons/s.weight)[0],
            s.name,
        )
    print "Total de: %0.5e" % (init_values['de'][0])
    #print init_values['de'][0] / (init_values['o_2'][0] / 16.0)

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
        if s.name == 'de': continue #don't want this in either density
        density += init_values[s.name][0]

    # set up initial temperatures values used to define ge
    temperature = na.logspace(4, 6, NCELLS)
    init_values['T'] = temperature

    # calculate ge (very crudely)
    gamma = 5.e0/3.e0
    init_values['ge'] = ((temperature * number_density * kboltz)
                         / (density * mh * (gamma - 1)))
    #init_values['ge'] = ((temperature * number_density * kboltz)
    #                     / ((gamma - 1)))
    
    create_initial_conditions(init_values, 'oxygen')
