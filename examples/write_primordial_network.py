from dengo.chemical_network import \
    ChemicalNetwork, \
    reaction_registry, \
    cooling_registry
import dengo.primordial_rates, dengo.primordial_cooling
from dengo.chemistry_constants import tiny, kboltz, mh
# from dengo.known_species import *

# If only a subset of species are wanted put them here
# and change the commented lines below
want = ("HI", "HII", "de", "ge")

primordial = ChemicalNetwork()
#primordial.add_energy_term()

# Set to false if intermediate solution output is not wanted
primordial.write_intermediate_solutions = True

#for ca in cooling_registry.values():
    #if not all(sp.name in want for sp in ca.species): continue
#    primordial.add_cooling(ca)

for i, rname in enumerate(sorted(reaction_registry)):
    s = reaction_registry[rname]
    #if not all(sp.name in want for sp in s.species): continue
    primordial.add_reaction(s)

# This defines the temperature range for the rate tables
primordial.init_temperature((1e0, 1e8))

# Generate initial conditions (switch to False to disable this)
generate_initial_conditions = True

if generate_initial_conditions:
    import numpy as np
    NCELLS = 4
    density = 1e7
    X = 0.5
    tiny = 1e-10

    init_array = np.ones(NCELLS) * density
    init_values = dict()
    init_values['HII']     = X * init_array
    init_values['HM']      = init_array * tiny
    init_values['HeI']     = init_array * tiny
    init_values['HeII']    = init_array * tiny
    init_values['HeIII']   = init_array * tiny
    init_values['H2I']     = init_array * tiny
    init_values['H2II']    = init_array * tiny
    init_values['de'] = init_array * 0.0

    total_density = primordial.calculate_total_density(init_values)
    init_values["HI"] = init_array.copy() - total_density
    init_values = primordial.convert_to_mass_density(init_values)
    init_values['de'] = primordial.calculate_free_electrons(init_values)
    init_values['density'] = primordial.calculate_total_density(init_values)
    number_density = primordial.calculate_number_density(init_values)
    temperature = np.logspace(2, 4, NCELLS)
    temperature[:] = 1e3
    init_values['T'] = temperature

    # calculate ge (very crudely, no H2 help here)
    gamma = 5.0/3.0
    init_values['ge'] = ((temperature * number_density * kboltz)
                         / (init_values['density'] * mh * (gamma - 1)))

primordial.write_solver("primordial", output_dir = ".",
                        init_values=init_values)
