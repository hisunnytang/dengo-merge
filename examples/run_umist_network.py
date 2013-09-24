import numpy as np
from dengo.chemical_network import \
    ChemicalNetwork, \
    reaction_registry, \
    cooling_registry
import dengo.umist_rates
from dengo.chemistry_constants import tiny, kboltz, mh
from dengo.known_species import *

NCELLS = 4
density = 1.0
temperature = np.logspace(4, 6.7, NCELLS)
temperature[:] = 1e7
X = 1e-3

umist = ChemicalNetwork()
umist.add_energy_term()

for r in reaction_registry.values():
        umist.add_reaction(s)

# This defines the temperature range for the rate tables
umist.init_temperature((1e0, 1e3))

tiny = 1e-10

init_array = np.ones(NCELLS) * density
init_values = dict()
init_values['O']     = X * init_array
init_values['C']    = init_array * X
init_values['H']     = init_array * X
init_values['H2O']      = init_array * X
init_values['H2']     = init_array * X
init_values['CO']    = init_array * X
init_values['OH']   = init_array * X
init_values['O2']    = init_array * X
init_values['de']      = init_array * 0.0

total_density = umist.calculate_total_density(init_values, ("OI",))
init_values["OI"] = init_array.copy() - total_density
init_values = umist.convert_to_mass_density(init_values)
init_values['de'] = umist.calculate_free_electrons(init_values)
init_values['density'] = umist.calculate_total_density(init_values)
number_density = umist.calculate_number_density(init_values)

# set up initial temperatures values used to define ge
init_values['T'] = temperature

# calculate ge (very crudely, no H2 help here)
gamma = 5.0/3.0
init_values['ge'] = ((temperature * number_density * kboltz)
                     / (init_values['density'] * mh * (gamma - 1)))

# Write the initial conditions file
umist.write_solver("umist", output_dir = ".")

import pyximport
pyximport.install(setup_args={"include_dirs":np.get_include()},
                  reload_support=True, inplace=True)

umist_solver_run = pyximport.load_module("umist_solver_run",
                            "umist_solver_run.pyx",
                            build_inplace = True, pyxbuild_dir = "_dengo_temp")
rv, rv_int = umist_solver_run.run_umist(init_values, 1e16)

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
