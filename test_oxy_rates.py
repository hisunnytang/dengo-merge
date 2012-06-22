from dengo.chemical_network import ChemicalNetwork
from dengo.reaction_classes import Species, chianti_rate, reaction_registry
import dengo.primordial_rates
import dengo.oxygen_rates

oxygen_network = ChemicalNetwork()

for rate in sorted(reaction_registry):

    if rate.startswith("o_"):
        oxygen_network.add_reaction(rate)

print "These species are required:"
print "\n".join([s.name for s in sorted(oxygen_network.required_species)])

for species in oxygen_network.required_species:
    print "Querying: ", species
    for rxn in oxygen_network.species_reactions(species):
        print "  ", rxn
