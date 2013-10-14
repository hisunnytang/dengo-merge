from dengo.chemical_network import \
    ChemicalNetwork, \
    reaction_registry, \
    cooling_registry
import dengo.primordial_rates, dengo.primordial_cooling
import dengo.ion_by_ion
from dengo.chemistry_constants import tiny, kboltz, mh
import numpy as np

NCELLS = 1
density = 1e-3
temperature = np.logspace(4, 6.7, NCELLS)
temperature[:] = 5e6
X = 1e-3

ion_by_ion = ChemicalNetwork()
ion_by_ion.add_energy_term()

for ca in cooling_registry.values():
   if ca.name.startswith("O"):
      ion_by_ion.add_cooling(ca)

ion_by_ion.add_cooling("brem")
ion_by_ion.add_cooling("reHII")
ion_by_ion.add_cooling("reHeIII")
ion_by_ion.add_cooling("ceHI")
ion_by_ion.add_cooling("reHeII2")
ion_by_ion.add_cooling("reHeII1")
ion_by_ion.add_cooling("ciHeIS")
ion_by_ion.add_cooling("ceHeII")
ion_by_ion.add_cooling("ciHI")
ion_by_ion.add_cooling("ceHeI")
ion_by_ion.add_cooling("ciHeI")
ion_by_ion.add_cooling("ciHeII")
ion_by_ion.add_cooling("compton")

for r in reaction_registry.values():
    if r.name.startswith("O"): 
        ion_by_ion.add_reaction(r)

ion_by_ion.add_reaction("k01")
ion_by_ion.add_reaction("k02")
ion_by_ion.add_reaction("k03")
ion_by_ion.add_reaction("k04")
ion_by_ion.add_reaction("k05")
ion_by_ion.add_reaction("k06")

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

    total_density = ion_by_ion.calculate_total_density(init_values, ("OI",))
    init_values["OI"] = init_array.copy() - total_density
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
                    print s.name, s.number, s.free_electrons + 1
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
