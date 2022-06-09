from dengo.chemical_network import \
    ChemicalNetwork, \
    reaction_registry, \
    cooling_registry
from dengo.reaction_classes import index_i, count_m
import dengo.primordial_rates, dengo.primordial_cooling
from dengo.write_rate_reader import create_rate_tables, create_rate_reader
import sympy
from sympy.utilities.codegen import codegen
from sympy.printing import print_ccode

primordial = ChemicalNetwork()
for ca in list(cooling_registry.values()):
    primordial.add_cooling(ca)

print("These species are required for cooling alone:")
print("\n".join([s.name for s in sorted(primordial.required_species)]))

for s in list(reaction_registry.values()):
    primordial.add_reaction(s)

print("These species are required for chemistry and cooling:")
print("\n".join([s.name for s in sorted(primordial.required_species)]))

functions = []

cooling = sum(v.equation for n,v in sorted(primordial.cooling_actions.items()))

for species in primordial.required_species:
    print()
    print("// HANDLING SPECIES", species.name)
    print()
    eq = primordial.species_total(species)
    primordial.print_ccode(species)
    primordial.print_jacobian(species)

#primordial.print_ccode(primordial.energy_term)
#primordial.print_jacobian(primordial.energy_term)

    #ds_dt = sympy.IndexedBase("d_%s" % species.name, (count_m,))
    #expr = sympy.Equality(ds_dt[index_i], eq)
    #functions.append(("d_%s" % species.name, expr))

#print codegen(functions, "C", "codegen", to_files=True)

#from sympy.abc import x, y, z
#print codegen(("f", x+y*z), "C", "test")

T_bounds = [1.0e4, 1.0e8]
primordial.init_temperature(T_bounds)
print()
print("// WRITING RATE TABLES TO HDF5 FILES")
print()
create_rate_tables(primordial, "primordial")

print()
print("// WRITING C CODE TO READ HDF5 FILES")
print()
create_rate_reader(primordial, "primordial")
