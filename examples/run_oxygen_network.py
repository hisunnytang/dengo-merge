import numpy as np

import dengo.oxygen_cooling
import dengo.oxygen_rates
import dengo.primordial_cooling
import dengo.primordial_rates
from dengo.chemical_network import ChemicalNetwork, cooling_registry, reaction_registry
from dengo.chemistry_constants import kboltz, mh, tiny
from dengo.known_species import *

NCELLS = 4
density = 1.0
temperature = np.logspace(4, 6.7, NCELLS)
temperature[:] = 1e7
X = 1e-3

oxygen = ChemicalNetwork()
oxygen.add_energy_term()

for ca in cooling_registry.values():
    if ca.name.startswith("O"):
        oxygen.add_cooling(ca)

for s in reaction_registry.values():
    if s.name.startswith("O"):
        oxygen.add_reaction(s)

# This defines the temperature range for the rate tables
oxygen.init_temperature((1e0, 1e8))

tiny = 1e-10

init_array = np.ones(NCELLS) * density
init_values = dict()
init_values["OII"] = X * init_array
init_values["OIII"] = init_array * X
init_values["OIV"] = init_array * X
init_values["OV"] = init_array * X
init_values["OVI"] = init_array * X
init_values["OVII"] = init_array * X
init_values["OVIII"] = init_array * X
init_values["OIX"] = init_array * X
init_values["de"] = init_array * 0.0

total_density = oxygen.calculate_total_density(init_values, ("OI",))
init_values["OI"] = init_array.copy() - total_density
init_values = oxygen.convert_to_mass_density(init_values)
init_values["de"] = oxygen.calculate_free_electrons(init_values)
init_values["density"] = oxygen.calculate_total_density(init_values)
number_density = oxygen.calculate_number_density(init_values)

# set up initial temperatures values used to define ge
init_values["T"] = temperature

# calculate ge (very crudely, no H2 help here)
gamma = 5.0 / 3.0
init_values["ge"] = (temperature * number_density * kboltz) / (
    init_values["density"] * mh * (gamma - 1)
)

# Write the initial conditions file
oxygen.write_solver("oxygen", output_dir=".")

import pyximport

pyximport.install(
    setup_args={"include_dirs": np.get_include()}, reload_support=True, inplace=True
)

oxygen_solver_run = pyximport.load_module(
    "oxygen_solver_run",
    "oxygen_solver_run.pyx",
    build_inplace=True,
    pyxbuild_dir="_dengo_temp",
)
rv, rv_int = oxygen_solver_run.run_oxygen(init_values, 1e16)

import pylab

pylab.clf()

mask = rv_int["successful"]
for name in sorted(rv_int):
    if len(rv_int[name].shape) == 1:
        rv_int[name] = rv_int[name][mask]
    else:
        rv_int[name] = rv_int[name][0, mask]

skip = ("successful", "dt", "t", "ge")
for n, v in sorted(rv_int.items()):
    if n in skip:
        continue
    pylab.loglog(rv_int["t"], v, label=n)

pylab.ylim(density * 1e-20, density * 10)
pylab.xlabel("time [s]")
pylab.legend(loc="best", fontsize="xx-small")
pylab.savefig("plot.png")
