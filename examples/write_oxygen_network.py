from dengo.chemical_network import \
    ChemicalNetwork, \
    reaction_registry, \
    cooling_registry
import dengo.primordial_rates, dengo.primordial_cooling
import dengo.oxygen_rates, dengo.oxygen_cooling
from dengo.chemistry_constants import tiny, kboltz, mh
import numpy as np

NCELLS = 1
density = 1.0
temperature = np.logspace(4, 6.7, NCELLS)
temperature[:] = 5e6
X = 1e-3

oxygen = ChemicalNetwork(write_intermediate = True)
oxygen.add_energy_term()

for ca in cooling_registry.values():
    if ca.name.startswith("O"):
       oxygen.add_cooling(ca)

for s in reaction_registry.values():
    if s.name.startswith("O"):
        oxygen.add_reaction(s)

# This defines the temperature range for the rate tables
oxygen.init_temperature((1e0, 1e8))

# Want intermediate output?

tiny = 1e-10

init_array = np.ones(NCELLS) * density
init_values = dict()

# set up initial temperatures values used to define ge
init_values['T'] = temperature

start_neutral = False

if start_neutral:
    init_values['OII']     = X * init_array
    init_values['OIII']    = init_array * X
    init_values['OIV']     = init_array * X
    init_values['OV']      = init_array * X
    init_values['OVI']     = init_array * X
    init_values['OVII']    = init_array * X
    init_values['OVIII']   = init_array * X
    init_values['OIX']    = init_array * X
    init_values['de']      = init_array * 0.0

    total_density = oxygen.calculate_total_density(init_values, ("OI",))
    init_values["OI"] = init_array.copy() - total_density
    init_values = oxygen.convert_to_mass_density(init_values)
else:
    # start CIE
    import chianti.core as ch
    import chianti.util as chu

    for s in sorted(oxygen.required_species):
            if s.name != 'ge':
                if s.name == 'de':
                    continue
                else:
                    ion_name = chu.zion2name(s.number, s.free_electrons + 1)
                    ion = ch.ion(ion_name, temperature=init_values['T'])
                    ion.ioneqOne()
                    ion_frac = ion.IoneqOne
                    init_values[s.name] = ion_frac * init_array * ion.Abundance
                
                # in case something is negative or super small:
                init_values[s.name][init_values[s.name] < tiny] = tiny

    init_values['de']      = init_array * 0.0
#    total_density = oxygen.calculate_total_density(init_values, ("OI",))
#    init_values["OI"] = init_array.copy() - total_density
    init_values = oxygen.convert_to_mass_density(init_values)

init_values['de'] = oxygen.calculate_free_electrons(init_values)
init_values['density'] = oxygen.calculate_total_density(init_values)
number_density = oxygen.calculate_number_density(init_values)

# calculate ge (very crudely)
gamma = 5.0/3.0
init_values['ge'] = ((temperature * number_density * kboltz)
                     / (init_values['density'] * mh * (gamma - 1)))


    # Write the initial conditions file
oxygen.write_solver("oxygen", output_dir = ".",
                    init_values=init_values)
