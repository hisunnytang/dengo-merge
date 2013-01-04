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

primordial.add_cooling("brem")
primordial.add_cooling("reHII")
primordial.add_cooling("reHeIII")
primordial.add_cooling("gloverabel08")
primordial.add_cooling("ceHI")
primordial.add_cooling("h2formation")
primordial.add_cooling("reHeII2")
primordial.add_cooling("reHeII1")
primordial.add_cooling("ciHeIS")
primordial.add_cooling("ceHeII")
primordial.add_cooling("ciHI")
primordial.add_cooling("ceHeI")
primordial.add_cooling("gammah")
primordial.add_cooling("ciHeI")
primordial.add_cooling("ciHeII")

primordial.add_reaction("k01")
primordial.add_reaction("k02")
primordial.add_reaction("k03")
primordial.add_reaction("k04")
primordial.add_reaction("k05")
primordial.add_reaction("k06")
primordial.add_reaction("k07")
primordial.add_reaction("k08")
primordial.add_reaction("k09")
primordial.add_reaction("k10")
primordial.add_reaction("k11")
primordial.add_reaction("k12")
primordial.add_reaction("k13")
primordial.add_reaction("k14")
primordial.add_reaction("k15")
primordial.add_reaction("k16")
primordial.add_reaction("k17")
primordial.add_reaction("k18")
primordial.add_reaction("k19")
primordial.add_reaction("k21")
primordial.add_reaction("k22")
primordial.add_reaction("k23")

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
pyximport.install(setup_args={"include_dirs":np.get_include()},
                  reload_support=True, inplace=True)

primordial_solver_run = pyximport.load_module("primordial_solver_run",
                            "primordial_solver_run.pyx",
                            build_inplace = True, pyxbuild_dir = "_dengo_temp")
primordial_solver_run.run_primordial(init_values, 1e8)
