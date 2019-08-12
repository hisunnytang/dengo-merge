import numpy as np
from dengo.chemical_network import \
    ChemicalNetwork, \
    reaction_registry, \
    cooling_registry
import dengo.primordial_rates, dengo.primordial_cooling
from dengo.chemistry_constants import tiny, kboltz, mh

import os
# this is required to compiled cython
os.environ["HDF5_DIR"] = "/home/kwoksun2/anaconda3"

NCELLS = 4
density = 1e12
temperature = np.logspace(2, 4, NCELLS)
temperature[:] = 2e3
X = 1.0e-4

dengo.primordial_rates.setup_primordial()
primordial = ChemicalNetwork()
#primordial.add_energy_term()

# add the neccessary cooling and reactions
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

# This defines the temperature range for the rate tables
primordial.init_temperature((1e0, 1e8))

# setting initial conditions
init_array = np.ones(NCELLS) * density
init_values = dict()
init_values["H_1"]      = (1 - X) * init_array
init_values['H_2']     = X * init_array
init_values['H_m0']      = init_array * tiny
init_values['He_1']     = init_array * tiny
init_values['He_2']    = init_array * tiny
init_values['He_3']   = init_array * tiny
init_values['H2_1']     = init_array * tiny
init_values['H2_2']    = init_array * tiny
init_values['de'] = init_array * 0.0

# update and calculate electron density and etc with the handy functions
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

# Write the initial conditions file
# IF you need to use the Makefile, and c-library
# you will have to specified the library_path
library_path = {}
library_path["CVODE_PATH"] = "/home/kwoksun2/cvode-3.1.0/instdir"
library_path["HDF5_PATH"] = "/home/kwoksun2/anaconda3"
library_path["SUITESPARSE_PATH"] = "/home/kwoksun2/SuiteSparse"
library_path["DENGO_INSTALL_PATH"] = "/home/kwoksun2/dengo_install"

primordial.write_solver("primordial", output_dir = ".", solver_template = "cv_omp/sundials_CVDls",ode_solver_source = "initialize_cvode_solver.C", library_path = library_path)

import pyximport
pyximport.install(setup_args={"include_dirs":np.get_include()},
                  reload_support=True, inplace=True)

primordial_solver_run = pyximport.load_module("primordial_solver_run",
                            "primordial_solver_run.pyx",
                            build_inplace = True, pyxbuild_dir = "_dengo_temp")
rv, rv_int = primordial_solver_run.run_primordial(init_values, 1e10, niter=100)

import pylab
pylab.clf()

mask = rv_int['successful']
for name in sorted(rv_int):
    if len(rv_int[name].shape) == 1:
        rv_int[name] = rv_int[name][mask]
    else:
        rv_int[name] = rv_int[name][0, mask]

skip = ('successful', 'dt', 't', 'ge')
for n, v in sorted(rv_int.items()):
    if n in skip: continue
    pylab.loglog(rv_int['t'], v, label = n)

pylab.ylim(density * 1e-20, density * 10)
pylab.xlabel("time [s]")
pylab.legend(loc='best', fontsize='xx-small')
pylab.savefig("plot.png")
