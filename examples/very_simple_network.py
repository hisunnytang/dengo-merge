import numpy as np
from dengo.chemical_network import \
    ChemicalNetwork, \
    reaction_registry, \
    cooling_registry
import dengo.primordial_rates, dengo.primordial_cooling
from dengo.chemistry_constants import tiny, kboltz, mh

NCELLS = 4
density = 1e7
temperature = np.logspace(2, 4, NCELLS)
temperature[:] = 2e3
X = 0.75

dengo.primordial_rates.setup_primordial()

primordial = ChemicalNetwork()

primordial.add_reaction("k01")
primordial.add_reaction("k02")
primordial.add_reaction("k03")
primordial.add_reaction("k04")
primordial.add_reaction("k05")
primordial.add_reaction("k06")

# This defines the temperature range for the rate tables
primordial.init_temperature((1e0, 1e8))

tiny = 1e-10

init_array = np.ones(NCELLS) * density
init_values = dict()
init_values['H_2']     = X * init_array
init_values["H_1"]     = X * init_array
init_values["He_1"] = ( 1.0 - X )*init_array
init_values["He_2"] = ( 1.0 - X )*init_array
init_values["He_3"] = ( 1.0 - X )*init_array


total_density = primordial.calculate_total_density(init_values)


init_values = primordial.convert_to_mass_density(init_values)
init_values['de'] = primordial.calculate_free_electrons(init_values)
init_values['density'] = primordial.calculate_total_density(init_values)
number_density = primordial.calculate_number_density(init_values)

# set up initial temperatures values used to define ge
init_values['T'] = temperature

# calculate ge (very crudely, no H2 help here)
gamma = 5.0/3.0
init_values['ge'] = ((temperature * number_density * kboltz)
                     / (init_values['density'] * mh * (gamma - 1)))

primordial.write_solver("sunny", output_dir = ".",
                        init_values=init_values,
                        input_is_number=False)
