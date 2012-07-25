from dengo.chemical_network import \
    ChemicalNetwork, \
    reaction_registry, \
    cooling_registry
import dengo.primordial_rates, dengo.primordial_cooling
from dengo.solver_writer import write_reaction, write_species_deriv
from sympy.printing import print_ccode

primordial = ChemicalNetwork()
for ca in cooling_registry.values():
    primordial.add_cooling(ca)

print "These species are required for cooling alone:"
print "\n".join([s.name for s in sorted(primordial.required_species)])

for s in reaction_registry.values():
    primordial.add_reaction(s)

print "These species are required for chemistry and cooling:"
print "\n".join([s.name for s in sorted(primordial.required_species)])

for species in primordial.required_species:
    eq = primordial.species_total(species)
    print_ccode(eq, assign_to = "d_%s" % species.name)
