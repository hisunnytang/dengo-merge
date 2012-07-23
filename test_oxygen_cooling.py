from dengo.chemical_network import \
    ChemicalNetwork, \
    reaction_registry, \
    cooling_registry
import dengo.primordial_rates, dengo.primordial_cooling
import dengo.oxygen_rates, dengo.oxygen_cooling
from dengo.solver_writer import write_reaction, write_species_deriv

oxygen = ChemicalNetwork()
for ca in cooling_registry.values():
    if ca.name.startswith("o_"):
        print "ca = %s" %ca
        oxygen.add_cooling(ca)

print "These species are required for cooling alone:"
print "\n".join([s.name for s in sorted(oxygen.required_species)])

for s in reaction_registry.values():
    if s.name.startswith("o_"):
        oxygen.add_reaction(s)

print "These species are required for chemistry and cooling:"
print "\n".join([s.name for s in sorted(oxygen.required_species)])
