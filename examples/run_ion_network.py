import numpy as np
from dengo.chemical_network import \
    ChemicalNetwork, \
    reaction_registry, \
    cooling_registry
from dengo.ion_by_ion import setup_ionization
from dengo.chemistry_constants import tiny, kboltz, mh
import numpy as np

NCELLS = 1
density = 1e-3
#* 1.67e-24
temperature = np.logspace(2, 8, NCELLS)
temperature[:] = 5e6
X = 1e-8

ion_by_ion = ChemicalNetwork(write_intermediate = False,
                             stop_time = 3.1557e13)
ion_by_ion.add_species("de")

for atom in ["H", "He", "O"]:#"C", "N", "O"]:# , "Ne", "Mg", "Si", "S"]:
    s, c, r = setup_ionization(atom)
    ion_by_ion.add_collection(s, c, r)

# This defines the temperature range for the rate tables
ion_by_ion.init_temperature((1e0, 1e12))

# This defines the redsfhit range for the rate tables
ion_by_ion.init_redshift((0.0, 9.0))

tiny = 1e-10

init_array = np.ones(NCELLS) * density
init_values = dict()

# set up initial temperatures values used to define ge
init_values['T'] = temperature

start_neutral = True

if start_neutral:
    for s in ion_by_ion.required_species:
        if getattr(s, 'free_electrons', -1) == 0:
            init_values[s.name] = init_array.copy()
        else:
            init_values[s.name] = X * init_array
    init_values['de'][:] = 1e-30
    init_values = ion_by_ion.convert_to_mass_density(init_values)
else:
    # start CIE
    import chianti.core as ch
    import chianti.util as chu

    for s in sorted(ion_by_ion.required_species):
            if s.name != 'ge':
                if s.name == 'de':
                    continue
                else:
                    print s.name, s.free_electrons + 1
                    ion_name = chu.zion2name(np.int(s.number),
                                             np.int(s.free_electrons + 1))
                    ion = ch.ion(ion_name, temperature=init_values['T'])
                    ion.ioneqOne()
                    ion_frac = ion.IoneqOne
                    init_values[s.name] = ion_frac * init_array * ion.Abundance
                
                # in case something is negative or super small:
                init_values[s.name][init_values[s.name] < tiny] = tiny

    init_values['de']      = init_array * 0.0
    init_values = ion_by_ion.convert_to_mass_density(init_values)

init_values['de'] = ion_by_ion.calculate_free_electrons(init_values)
init_values['density'] = ion_by_ion.calculate_total_density(init_values)
number_density = ion_by_ion.calculate_number_density(init_values)

# calculate ge (very crudely)
gamma = 5.0/3.0
init_values['ge'] = ((temperature * number_density * kboltz)
                     / (init_values['density'] * mh * (gamma - 1)))

# Write the initial conditions file
ion_by_ion.write_solver("ion_by_ion", output_dir = ".",
                        init_values=init_values,
                        input_is_number=False)

import pyximport
pyximport.install(setup_args={"include_dirs":np.get_include()},
                  reload_support=True, inplace=True)

ion_by_ion_solver_run = pyximport.load_module("ion_by_ion_solver_run",
                            "ion_by_ion_solver_run.pyx",
                            build_inplace = True, pyxbuild_dir = "_dengo_temp")
rv, rv_int = ion_by_ion_solver_run.run_ion_by_ion(init_values, 1e16, 100000)
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

pylab.ylim(density * 1e-30, density * 10)
pylab.xlabel("time [s]")
pylab.legend(loc='best', fontsize='xx-small')
pylab.savefig("plot.png")
