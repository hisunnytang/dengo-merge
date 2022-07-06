import os

import numpy as np

from dengo.chemical_network import (
    ChemicalNetwork,
    cooling_registry,
    reaction_registry,
    species_registry,
)
from dengo.chemistry_constants import kboltz, mh, tiny
from dengo.ion_by_ion import setup_ionization
from dengo.reaction_classes import AtomicSpecies

os.environ["HDF5_DIR"] = "/home/kwoksun2/anaconda3"

NCELLS = 1
density = 1.0
# * 1.67e-24
temperature = np.logspace(2, 8, NCELLS)
temperature[:] = 5e6
X = 1e-8

ion_by_ion = ChemicalNetwork(write_intermediate=False, stop_time=3.1557e13)
HI = AtomicSpecies("H", 0.0)
HII = AtomicSpecies("H", 1.0)
HeI = AtomicSpecies("He", 0.0)
HeII = AtomicSpecies("He", 1.0)
de = species_registry["de"]

OI = AtomicSpecies("O", 0.0)
OII = AtomicSpecies("O", 1.0)
OIII = AtomicSpecies("O", 2.0)
OIV = AtomicSpecies("O", 3.0)
OV = AtomicSpecies("O", 4.0)
OVI = AtomicSpecies("O", 5.0)
ion_by_ion.add_species(de)
ion_by_ion.add_species(HI)
ion_by_ion.add_species(HII)
ion_by_ion.add_species(HeI)
ion_by_ion.add_species(HeII)
ion_by_ion.add_species(OI)
ion_by_ion.add_species(OII)
ion_by_ion.add_species(OIII)
ion_by_ion.add_species(OIV)
ion_by_ion.add_species(OV)
ion_by_ion.add_species(OVI)

for atom in ["H", "He", "O"]:  # "C", "N", "O", "Ne", "Si"]:
    s, c, r = setup_ionization(atom, photo_background="HM12")
    ion_by_ion.add_collection(s, c, r)

# This defines the temperature range for the rate tables
ion_by_ion.init_temperature((1e0, 1e12))

# This defines the redsfhit range for the rate tables
ion_by_ion.init_redshift((0.0, 9.0))

tiny = 1e-20

init_array = np.ones(NCELLS) * density
init_values = dict()

# set up initial temperatures values used to define ge
init_values["T"] = temperature

start_neutral = False

import ChiantiPy.core as ch

for s in ion_by_ion.required_species:
    print(s, type(s))

if start_neutral:
    for s in ion_by_ion.required_species:
        if getattr(s, "free_electrons", -1) == 0:
            init_values[s.name] = init_array.copy()
        else:
            init_values[s.name] = X * init_array
        # Scale to solar abundances
        if s.name not in ["de", "ge"]:
            ion_name = s.name.lower()
            ion = ch.ion(ion_name, temperature=init_values["T"])
            init_values[s.name] *= ion.Abundance

    init_values["de"][:] = 1e-30
    init_values = ion_by_ion.convert_to_mass_density(init_values)
else:
    # start CIE

    for s in sorted(ion_by_ion.required_species):
        if s.name != "ge":
            if s.name == "de":
                continue
            else:
                ion_name = s.name.lower()
                ion = ch.ion(ion_name, temperature=init_values["T"])
                print(ion_name)
                # ion.ioneqOne()
                # this calcuate the equilirbium abundance @ T
                # however this attr is not defined for fully ionized species...
                # fix that later, take this as 1 for now first
                # ion_frac = ion.IoneqOne
                ion_frac = 1.0
                # ion.Abundance is the elemental abundance relative to hydrogen
                init_values[s.name] = ion_frac * init_array * ion.Abundance

                # in case something is negative or super small:
                init_values[s.name][init_values[s.name] < tiny] = tiny

    init_values["de"] = init_array * 0.0
    init_values = ion_by_ion.convert_to_mass_density(init_values)

init_values["de"] = ion_by_ion.calculate_free_electrons(init_values)
init_values["density"] = ion_by_ion.calculate_total_density(init_values)
number_density = ion_by_ion.calculate_number_density(init_values)

# calculate ge (very crudely)
gamma = 5.0 / 3.0
init_values["ge"] = (temperature * number_density * kboltz) / (
    init_values["density"] * mh * (gamma - 1)
)


# Write the initial conditions file
# IF you need to use the Makefile, and c-library
# you will have to specified the library_path


# Write the initial conditions file
ion_by_ion.write_solver(
    "ion_by_ion",
    output_dir=".",
    init_values=init_values,
    input_is_number=False,
    solver_template="cv_omp/sundials_CVDls",
    ode_solver_source="initialize_cvode_solver.C",
)

import pyximport

pyximport.install(
    setup_args={"include_dirs": np.get_include()}, reload_support=True, inplace=True
)

ion_by_ion_solver_run = pyximport.load_module(
    "ion_by_ion_solver_run",
    "ion_by_ion_solver_run.pyx",
    build_inplace=True,
    pyxbuild_dir="_dengo_temp",
)
rv, rv_int = ion_by_ion_solver_run.run_ion_by_ion(init_values, 1e16, 100000, z=0.0)
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

pylab.ylim(density * 1e-30, density * 10)
pylab.xlabel("time [s]")
pylab.legend(loc="best", fontsize="xx-small")
pylab.savefig("plot.png")

pylab.clf()
pylab.loglog(rv_int["t"], rv_int["T"], label="T")
pylab.xlabel("time [s]")
pylab.savefig("plot_temp.png")
