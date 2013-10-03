import numpy as np
from dengo.chemical_network import \
    ChemicalNetwork, \
    reaction_registry, \
    cooling_registry
import dengo.umist_rates
from dengo.get_rates import get_rates
from dengo.chemistry_constants import tiny, kboltz, mh
from dengo.known_species import *

NCELLS = 4
density = 1.0
temperature = np.logspace(4, 6.7, NCELLS)
temperature[:] = 1e7
X = 1e-3

umist = ChemicalNetwork()
umist.add_energy_term()

# This defines the temperature range for the rate tables
umist.init_temperature((1e1, 1e3))

get_rates('CO', 28, -1, umist)
get_rates('C', 12, -1, umist)
get_rates('O', 16, -1, umist)
get_rates('OH', 17, -1, umist)
get_rates('H', 1, -1, umist)
get_rates('H2', 2, -2, umist)
get_rates('H2O', 18, -1, umist)
get_rates('O2', 32, -1, umist)

# Define small subset of species for restricted calculation
sub = set(['us_CO', 'us_COm','us_em', 'us_C','us_Cm','us_Cp','us_Om','us_Op','us_O','us_OH','us_OHp','us_OHm','us_H','us_Hm','us_Hp','us_H2','us_H2p','us_H2m'])

# Add an if statement - want to restrict the species to the 8 above (and their ions and electrons)
for r in reaction_registry.values():
    s = r.considered
    if s.issubset(sub) == True:
        print 'Found one!!!' 
        umist.add_reaction(r)
print 'Finished looking through reaction registry.'

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
init_values['em']      = init_array * 0.0

print init_values
print sorted(umist.reactions.values())

for species in umist.required_species:
    if species.name not in init_values:
        init_values[species.name] = init_array * 0.0

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
