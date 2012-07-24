from dengo.chemical_network import \
    ChemicalNetwork, \
    reaction_registry, \
    cooling_registry
import dengo.primordial_rates, dengo.primordial_cooling
import dengo.oxygen_rates, dengo.oxygen_cooling
from dengo.solver_writer import write_reaction, write_species_deriv
from matplotlib import pylab

oxygen = ChemicalNetwork()
for ca in cooling_registry.values():
    if ca.name.startswith("o_"):
        oxygen.add_cooling(ca)

print "These species are required for cooling alone:"
print "\n".join([s.name for s in sorted(oxygen.required_species)])

for s in reaction_registry.values():
    if s.name.startswith("o_"):
        oxygen.add_reaction(s)

print "These species are required for chemistry and cooling:"
print "\n".join([s.name for s in sorted(oxygen.required_species)])

oxygen.init_temperature(T_bounds=(1.0e4, 1.0e8))

# Let's create a cooling function plot for all oxygen ion species
pylab.ion()
pylab.figure(1)
for cooling_action in oxygen.cooling_actions.values():
    pylab.loglog(oxygen.T, cooling_action.tables[cooling_action.name](oxygen), label=cooling_action.name)

pylab.legend()

