from dengo.chemical_network import ChemicalNetwork, reaction_registry
import dengo.primordial_rates
from dengo.solver_writer import write_reaction, write_species_deriv

primordial = ChemicalNetwork()
for i in range(23):
    n = "k%02i" % (i+1)
    if n in reaction_registry:
        primordial.add_reaction(n)

print "These species are required:"
print "\n".join([s.name for s in sorted(primordial.required_species)])

for species in primordial.required_species:
    print "Querying: ", species
    for rxn in primordial.species_reactions(species):
        print "  ", rxn

for rxn in primordial:
    write_reaction(rxn)

for species in primordial.required_species:
    write_species_deriv(species, primordial.species_reactions(species))
