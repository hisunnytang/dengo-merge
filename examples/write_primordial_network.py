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

primordial = ChemicalNetwork()
for ca in cooling_registry.values():
    primordial.add_cooling(ca)

for s in reaction_registry.values():
    primordial.add_reaction(s)

primordial.init_temperature((1e4, 1e8))

create_rate_tables(primordial, "primordial")
create_rate_reader(primordial, "primordial")

# Generate initial conditions (switch to False to disable this)
generate_initial_conditions = True

if generate_initial_conditions:
    import numpy as na
    NCELLS = 64
    density = 1.0e3
    init_array = na.ones(NCELLS) 
    X = 1.0e-4

    init_values = dict()
    init_values['density'] = density * init_array
    init_values['HI'] = init_array # use conservation to set this below
    # populate initial fractional values for the other species
    for s in sorted(primordial.required_species):
        if s.name != 'ge' and s.name != 'HI':
            if s.name == 'de' or s.name == 'HII':
                init_values[s.name] = X * init_array
            else:
                init_values[s.name] = tiny * init_array
            init_values['HI'] -= init_values[s.name]

    # convert to masses to multiplying by the density factor
    for iv in init_values:
        init_values[iv] *= density
    
    # set up initial temperatures values used to define ge
    temperature = na.logspace(2, 6, NCELLS)

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
