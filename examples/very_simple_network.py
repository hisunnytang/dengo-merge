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
X = 0.5

dengo.primordial_rates.setup_primordial()

primordial = ChemicalNetwork()

primordial.add_reaction("k01")
primordial.add_reaction("k02")

# This defines the temperature range for the rate tables
primordial.init_temperature((1e0, 1e8))

tiny = 1e-10

init_array = np.ones(NCELLS) * density
init_values = dict()
init_values['HII']     = X * init_array
init_values['de'] = init_array * 0.0

#total_density = primordial.calculate_total_density(init_values, ("HI",))
#init_values["HI"] = init_array.copy() - total_density
#init_values = primordial.convert_to_mass_density(init_values)
#init_values['de'] = primordial.calculate_free_electrons(init_values)
#init_values['density'] = primordial.calculate_total_density(init_values)
#number_density = primordial.calculate_number_density(init_values)

# set up initial temperatures values used to define ge
init_values['T'] = temperature

# calculate ge (very crudely, no H2 help here)
gamma = 5.0/3.0
#init_values['ge'] = ((temperature * number_density * kboltz)
                     #/ (init_values['density'] * mh * (gamma - 1)))

primordial.write_solver("sunny", output_dir = ".",
                        init_values=init_values,
                        input_is_number=False)
