import os

import numpy as np

import dengo.umist_rates
from dengo.chemical_network import ChemicalNetwork, cooling_registry, reaction_registry
from dengo.chemistry_constants import kboltz, mh, tiny
from dengo.get_rates import setup_umist_reactions, setup_umist_species

NCELLS = 1
density = 1.0e0
temperature = np.logspace(1, 3, NCELLS)
temperature[:] = 1e4
X = 1e-2

umist = ChemicalNetwork()
umist.skip_weight += ("us_e_0",)

# This defines the temperature range for the rate tables
umist.init_temperature((1e0, 1e5))

# Get UMIST rates for a bunch of species for an example network

desired_species = [("H", 1), ("H2", 2), ("e-", 1)]
added_species = set([])

for name, weight in desired_species:
    s, c, r = setup_umist_species(name, weight)
    print(("ADDED", s, c, r))
    added_species.update(s)
    umist.add_collection(s, c, r)
# Add ionic species by hand, since we need correct atomic weights
print(added_species)
s, c, r = setup_umist_reactions(added_species)
umist.add_collection(s, c, r)
print(reaction_registry)

tiny = 1e-10

init_array = np.ones(NCELLS) * density
init_values = dict()
init_values["us_H_1"] = init_array * X
init_values["us_H2_1"] = init_array * X
init_values["us_e_0"] = init_array * 1e-10

print(init_values)
# print sorted(umist.reactions.values())

for species in umist.required_species:
    if species.name not in init_values:
        init_values[species.name] = init_array * 0.0

total_density = umist.calculate_total_density(init_values)
init_values = umist.convert_to_mass_density(init_values)
init_values["us_e_0"] = init_array * 1e-10
# umist.calculate_free_electrons(init_values)
init_values["density"] = umist.calculate_total_density(init_values)
number_density = umist.calculate_number_density(init_values)


# set up initial temperatures values used to define ge
init_values["T"] = temperature


# calculate ge (very crudely, no H2 help here)
gamma = 5.0 / 3.0
init_values["ge"] = (temperature * number_density * kboltz) / (
    init_values["density"] * mh * (gamma - 1)
)

print(init_values)


# import pdb; pdb.set_trace()

# Write the initial conditions file
umist.write_solver(
    "umist",
    output_dir=".",
    solver_template="cv_omp/sundials_CVDls",
    ode_solver_source="initialize_cvode_solver.C",
)


import pyximport

pyximport.install(
    setup_args={"include_dirs": np.get_include()}, reload_support=True, inplace=True
)

umist_solver_run = pyximport.load_module(
    "umist_solver_run",
    "umist_solver_run.pyx",
    build_inplace=True,
    pyxbuild_dir="_dengo_temp",
)
rv, rv_int = umist_solver_run.run_umist(init_values, 1e16, niter=1e4)

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
