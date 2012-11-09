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

primordial = ChemicalNetwork()
for ca in cooling_registry.values():
    if any(sp.name in ("H2I", "H2II", "HeI", "HeII", "HeIII") for sp in ca.species): continue
    primordial.add_cooling(ca)
#primordial.add_cooling(cooling_registry["ceHI"])

for i, rname in enumerate(sorted(reaction_registry)):
    s = reaction_registry[rname]
    if any(sp.name in ("H2I", "H2II", "HeI", "HeII", "HeIII") for sp in s.species): continue
    primordial.add_reaction(s)
#s = reaction_registry["k01"]
#primordial.add_reaction(s)
#s = reaction_registry["k02"]
#primordial.add_reaction(s)

primordial.init_temperature((1e0, 1e8))
primordial.write_intermediate_solutions = True

create_rate_tables(primordial, "primordial")
create_rate_reader(primordial, "primordial")

# Generate initial conditions (switch to False to disable this)
generate_initial_conditions = True

if generate_initial_conditions:
    import numpy as na
    NCELLS = 8
    density = 1.0e3
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
                print s.name, init_values[s.name][0]
            else:
                init_values[s.name] = tiny * init_array
            init_values['HI'] -= init_values[s.name]
    init_values['de'] = init_array * 0.0
    for s in sorted(primordial.required_species):
        if s.name == "ge": continue
        init_values['de'] += init_values[s.name] * s.free_electrons / s.weight

    print init_values['de'][0] / init_values['HII'][0],
    print init_values['de'][0],  init_values['HII'][0]
    # convert to masses to multiplying by the density factor
    for iv in init_values:
        init_values[iv] *= density
    print init_values['de'][0] / init_values['HII'][0],
    print init_values['de'][0],  init_values['HII'][0]
    
    # set up initial temperatures values used to define ge
    temperature = na.logspace(2, 4, NCELLS)

    # calculate the number density
    number_density = 0.0
    for s in sorted(primordial.required_species):
        if s.name != 'ge':
            number_density += init_values[s.name]/s.weight

    # calculate ge (very crudely)
    gamma = 5.e0/3.e0
    init_values['ge'] = ((temperature * density * kboltz)
                         / (number_density * mh * (gamma - 1)))
    
    create_initial_conditions(init_values, 'primordial')
