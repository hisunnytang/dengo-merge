from dengo.chemical_network import \
    ChemicalNetwork, \
    reaction_registry, \
    cooling_registry
from dengo.reaction_classes import index_i, count_m
import dengo.primordial_rates, dengo.primordial_cooling
import dengo.oxygen_rates, dengo.oxygen_cooling
from dengo.write_rate_reader import create_rate_tables, create_rate_reader
import sympy
from sympy.printing import print_ccode

oxygen = ChemicalNetwork()
for ca in list(cooling_registry.values()):
    if ca.name.startswith("o_"):
        oxygen.add_cooling(ca)

print("These species are required for cooling alone:")
print(("\n".join([s.name for s in sorted(oxygen.required_species)])))

for s in list(reaction_registry.values()):
    if s.name.startswith("o_"):
        oxygen.add_reaction(s)

print("These species are required for chemistry and cooling:")
print(("\n".join([s.name for s in sorted(oxygen.required_species)])))

functions = []

cooling = sum(v.equation for n,v in sorted(oxygen.cooling_actions.items()))

for species in oxygen.required_species:
    print()
    print(("// HANDLING SPECIES", species.name))
    print()
    eq = oxygen.species_total(species)
    oxygen.print_ccode(species)
    oxygen.print_jacobian(species)
    #ds_dt = sympy.IndexedBase("d_%s[i]" % species.name, (count_m,))
    #print_ccode(eq, assign_to = ds_dt)

T_bounds = [1.0e4, 1.0e8]
oxygen.init_temperature(T_bounds)
print()
print("// WRITING RATE TABLES TO HDF5 FILES")
print()
create_rate_tables(oxygen, "oxygen")

print()
print("// WRITING C CODE TO READ HDF5 FILES")
print()
create_rate_reader(oxygen, "oxygen")
