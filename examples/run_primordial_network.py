import numpy as np
from dengo.chemical_network import \
    ChemicalNetwork, \
    reaction_registry, \
    cooling_registry
import dengo.primordial_rates, dengo.primordial_cooling
from dengo.chemistry_constants import tiny, kboltz, mh
from dengo.known_species import *

NCELLS = 128
density = 1e7
X = 0.5

primordial = ChemicalNetwork()
primordial.add_energy_term()

for ca in cooling_registry.values():
    primordial.add_cooling(ca)

for r in reaction_registry.values():
    primordial.add_reaction(r)

# This defines the temperature range for the rate tables
primordial.init_temperature((1e0, 1e8))

init_array = np.ones(NCELLS) 
init_values = dict()
init_values['HII']     = X * init_array
init_values['HM']      = init_array * tiny
init_values['HeI']     = init_array * tiny
init_values['HeII']    = init_array * tiny
init_values['HeIII']   = init_array * tiny
init_values['H2I']     = init_array * tiny
init_values['H2II']    = init_array * tiny
init_values['de'] = init_array * 0.0

total_density = primordial.calculate_total_density(init_values, ("HI",))
init_values["HI"] = init_array.copy() - total_density
init_values = primordial.convert_to_mass_density(init_values)
init_values['de'] = primordial.calculate_free_electrons(init_values)
init_values['density'] = primordial.calculate_total_density(init_values)
number_density = primordial.calculate_number_density(init_values)

# set up initial temperatures values used to define ge
temperature = np.logspace(2, 4, NCELLS)
temperature[:] = 1e3
init_values['T'] = temperature

# calculate ge (very crudely, no H2 help here)
gamma = 5.0/3.0
init_values['ge'] = ((temperature * init_values['density'] * kboltz)
                     / (number_density * mh * (gamma - 1)))

# Write the initial conditions file
primordial.write_solver("primordial", output_dir = ".")

import pyximport
pyximport.install(setup_args={"include_dirs":np.get_include()})

primordial_solver_run = pyximport.load_module("primordial_solver_run",
                            "primordial_solver_run.pyx")
primordial_solver_run.run_primordial(init_values, 1e8)
