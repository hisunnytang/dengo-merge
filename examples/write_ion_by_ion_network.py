from dengo.chemical_network import \
    ChemicalNetwork, \
    reaction_registry, \
    cooling_registry
from dengo.primordial_rates import setup_primordial
from dengo.ion_by_ion import setup_ionization
from dengo.chemistry_constants import tiny, kboltz, mh
import numpy as np

NCELLS = 1
density = 1e-3
temperature = np.logspace(4, 6.7, NCELLS)
temperature[:] = 5e6
X = 1e-3

ion_by_ion = ChemicalNetwork()

s, c, r = setup_ionization("O")

ion_by_ion.add_species("de")
ion_by_ion.add_collection(s, c, r)

#s, c, r = setup_primordial()
#ion_by_ion.add_collection(s, c, r)

# This defines the temperature range for the rate tables
ion_by_ion.init_temperature((1e0, 1e12))

# This defines the redsfhit range for the rate tables
ion_by_ion.init_redshift((0.0, 9.0))

# Want intermediate output?
#combined.write_intermediate_solutions = True

tiny = 1e-10

init_array = np.ones(NCELLS) * density
init_values = dict()

# set up initial temperatures values used to define ge
init_values['T'] = temperature

start_neutral = True

if start_neutral:
    init_values['O_2']     = X * init_array
    init_values['O_3']    = init_array * X
    init_values['O_4']     = init_array * X
    init_values['O_5']      = init_array * X
    init_values['O_6']     = init_array * X
    init_values['O_7']    = init_array * X
    init_values['O_8']   = init_array * X
    init_values['O_9']    = init_array * X
    init_values['de']      = init_array * 0.0

    total_density = ion_by_ion.calculate_total_density(init_values, ("O_1",))
    init_values["O_1"] = init_array.copy() - total_density
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
#    total_density = combined.calculate_total_density(init_values, ("OI",))
#    init_values["OI"] = init_array.copy() - total_density
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
                    init_values=init_values)
