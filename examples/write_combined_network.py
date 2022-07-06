import numpy as np

import dengo.carbon_cooling
import dengo.carbon_rates
import dengo.magnesium_cooling
import dengo.magnesium_rates
import dengo.neon_cooling
import dengo.neon_rates
import dengo.nitrogen_cooling
import dengo.nitrogen_rates
import dengo.oxygen_cooling
import dengo.oxygen_rates
import dengo.primordial_cooling
import dengo.primordial_rates
import dengo.silicon_cooling
import dengo.silicon_rates
import dengo.sulfur_cooling
import dengo.sulfur_rates
from dengo.chemical_network import ChemicalNetwork, cooling_registry, reaction_registry
from dengo.chemistry_constants import kboltz, mh, tiny

NCELLS = 1
density = 1.0
temperature = np.logspace(4, 6.7, NCELLS)
temperature[:] = 5e6
X = 1e-3

combined = ChemicalNetwork()
combined.add_energy_term()

for ca in list(cooling_registry.values()):
    if (
        ca.name.startswith("C")
        or ca.name.startswith("N")
        or ca.name.startswith("O")
        or ca.name.startswith("Ne")
        or ca.name.startswith("Mg")
        or ca.name.startswith("Si")
        or ca.name.startswith("S")
    ):
        combined.add_cooling(ca)

combined.add_cooling("brem")
combined.add_cooling("reHII")
combined.add_cooling("reHeIII")
combined.add_cooling("ceHI")
combined.add_cooling("reHeII2")
combined.add_cooling("reHeII1")
combined.add_cooling("ciHeIS")
combined.add_cooling("ceHeII")
combined.add_cooling("ciHI")
combined.add_cooling("ceHeI")
combined.add_cooling("ciHeI")
combined.add_cooling("ciHeII")

for r in list(reaction_registry.values()):
    if (
        r.name.startswith("C")
        or r.name.startswith("N")
        or r.name.startswith("O")
        or r.name.startswith("Ne")
        or r.name.startswith("Mg")
        or r.name.startswith("Si")
        or r.name.startswith("S")
    ):
        combined.add_reaction(r)

combined.add_reaction("k01")
combined.add_reaction("k02")
combined.add_reaction("k03")
combined.add_reaction("k04")
combined.add_reaction("k05")
combined.add_reaction("k06")

# This defines the temperature range for the rate tables
combined.init_temperature((1e0, 1e8))

# Want intermediate output?
combined.write_intermediate_solutions = True

tiny = 1e-10

init_array = np.ones(NCELLS) * density
init_values = dict()

# set up initial temperatures values used to define ge
init_values["T"] = temperature

start_neutral = False

if start_neutral:
    init_values["OII"] = X * init_array
    init_values["OIII"] = init_array * X
    init_values["OIV"] = init_array * X
    init_values["OV"] = init_array * X
    init_values["OVI"] = init_array * X
    init_values["OVII"] = init_array * X
    init_values["OVIII"] = init_array * X
    init_values["OIX"] = init_array * X
    init_values["de"] = init_array * 0.0

    total_density = combined.calculate_total_density(init_values, ("OI",))
    init_values["OI"] = init_array.copy() - total_density
    init_values = combined.convert_to_mass_density(init_values)
else:
    # start CIE
    import chianti.core as ch
    import chianti.util as chu

    for s in sorted(combined.required_species):
        if s.name != "ge":
            if s.name == "de":
                continue
            else:
                print(s.name, s.number, s.free_electrons + 1)
                ion_name = chu.zion2name(np.int(s.number), np.int(s.free_electrons + 1))
                ion = ch.ion(ion_name, temperature=init_values["T"])
                ion.ioneqOne()
                ion_frac = ion.IoneqOne
                init_values[s.name] = ion_frac * init_array * ion.Abundance

            # in case something is negative or super small:
            init_values[s.name][init_values[s.name] < tiny] = tiny

    init_values["de"] = init_array * 0.0
    #    total_density = combined.calculate_total_density(init_values, ("OI",))
    #    init_values["OI"] = init_array.copy() - total_density
    init_values = combined.convert_to_mass_density(init_values)

init_values["de"] = combined.calculate_free_electrons(init_values)
init_values["density"] = combined.calculate_total_density(init_values)
number_density = combined.calculate_number_density(init_values)

# calculate ge (very crudely)
gamma = 5.0 / 3.0
init_values["ge"] = (temperature * number_density * kboltz) / (
    init_values["density"] * mh * (gamma - 1)
)


# Write the initial conditions file
combined.write_solver("combined", output_dir=".", init_values=init_values)
