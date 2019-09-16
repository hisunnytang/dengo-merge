import numpy as np
from dengo.chemical_network import \
    ChemicalNetwork, \
    reaction_registry, \
    cooling_registry
import dengo.new_rates
from dengo.chemistry_constants import tiny, kboltz, mh
#from dengo.known_species import *

import os
os.environ["HDF5_DIR"] = "/home/kwoksun2/anaconda3"

NCELLS = 1
density = 1e0
temperature = np.logspace(2, 4, NCELLS)
temperature[:] = 2e3
X = 0.5

dengo.new_rates.setup_new()
new = ChemicalNetwork()
new.add_reaction("I01")
new.add_reaction("I01_re")
new.add_reaction("I02")
new.add_reaction("I03")
new.add_reaction("I04")
new.add_reaction("I05")
new.add_reaction("C05")
new.add_reaction("D01")


# This defines the temperature range for the rate tables
new.init_temperature((1e0, 1e8))

init_array = np.ones(NCELLS) * density
init_values = dict()
init_values["HOI_1"]      = 8.0e-11 * init_array * (200)
init_values['I_m0']     = 1.0e-10 * init_array * (126)
init_values['H_2']     = 0.056 * init_array
init_values['I2_1']     = 8.0e-8 * init_array * (126*2)
init_values['H2O_1']     = 1.0 * init_array * (18)

init_values['HOIO_1']     = 9.0e-11 * init_array * (126 + 32 + 1)
init_values['IO3m_0']     = 0.01 * init_array * ( 1+16*3 )
init_values['O2_1']     = 2.5e3 * init_array * ( 16*2 )
init_values['CH2_COOH2_1']     = 0.0015 * init_array * ( 12+2+(12+16*2+1)*2 )
init_values['CHI_COOH2_1']     = 1.0e-20 * init_array * (  12+126+(12+16*2+1)*2 )
init_values["H2O2_1"] = 0.33*init_array*(2+16)*2

print(new.required_species)

total_density = new.calculate_total_density(init_values)
init_values = new.convert_to_mass_density(init_values)
init_values['de'] = new.calculate_free_electrons(init_values)
init_values['density'] = init_array # new.calculate_total_density(init_values)
number_density = new.calculate_number_density(init_values)

# set up initial temperatures values used to define ge
init_values['T'] = temperature

# calculate ge (very crudely, no H2 help here)
gamma = 5.0/3.0
init_values['ge'] = ((temperature * number_density * kboltz) / (init_values['density'] * mh * (gamma - 1)))


# Write the initial conditions file
# IF you need to use the Makefile, and c-library
# you will have to specified the library_path
library_path = {}
library_path["CVODE_PATH"] = "/home/kwoksun2/cvode-3.1.0/instdir"
library_path["HDF5_PATH"] = "/home/kwoksun2/anaconda3"
library_path["SUITESPARSE_PATH"] = "/home/kwoksun2/SuiteSparse"
library_path["DENGO_INSTALL_PATH"] = "/home/kwoksun2/dengo_install"


# Write the initial conditions file
"""
new.write_solver("new", output_dir = ".",
                 solver_template = "cv_omp/sundials_CVDls",
                 ode_solver_source = "initialize_cvode_solver.C",
                 library_path = library_path)
"""

import pyximport
pyximport.install(setup_args={"include_dirs":np.get_include()},
                  reload_support=True, inplace=True, language_level=3)

new_solver_run = pyximport.load_module("new_solver_run",
                            "new_solver_run.pyx",
                            build_inplace = True, pyxbuild_dir = "_dengo_temp")
rv, rv_int = new_solver_run.run_new(init_values, 5e1, niter=1e2)

import pylab
pylab.clf()

mask = rv_int['successful']
for name in sorted(rv_int):
    if len(rv_int[name].shape) == 1:
        rv_int[name] = rv_int[name][mask]
    else:
        rv_int[name] = rv_int[name][0, mask]

skip = ('successful', 'dt', 't', 'ge', "T")
for n, v in sorted(rv_int.items()):
    if n in skip: continue
    if n in ["I2_1", "I_m0"]:
        pylab.semilogy(rv_int['t'], v, label = n, marker = '.')


#print(rv_int)
#pylab.ylim(1e-20*density, 1e1*density)
pylab.xlabel("time [s]")
pylab.legend(loc='best', fontsize='xx-small')
pylab.savefig("plot_new.png")
