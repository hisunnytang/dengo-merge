import numpy as na
import h5py
from dengo.chemical_network import ChemicalNetwork
from dengo.reaction_classes import Species, chianti_rate, reaction_registry, species_registry
import dengo.primordial_rates
import dengo.oxygen_rates
from dengo.oxygen_cooling import cooling_rates_table, cooling_action_table

oxygen_network = ChemicalNetwork()

for rate in sorted(reaction_registry):
    if rate.startswith("o_"):
        oxygen_network.add_reaction(rate)

print("These species are required:")
print("\n".join([s.name for s in sorted(oxygen_network.required_species)]))

for species in oxygen_network.required_species:
    print("Querying: ", species)
    for rxn in oxygen_network.species_reactions(species):
        print("  ", rxn)

for rxn in oxygen_network:
    write_reaction(rxn)

for species in oxygen_network.required_species:
    print(species.name)
    write_species_deriv(species, oxygen_network.species_reactions(species))

# Define the temperature we want rate values for
Tlog_bounds = (1, 9) 
T = na.logspace(Tlog_bounds[0], Tlog_bounds[1], 1024)

# Write the ionization/recombination rate values to file
f = h5py.File('oxygen_rates_tables.h5', 'w')
for rxn_name in na.sort(list(oxygen_network.reactions.keys())):
    print('Reaction %s is being written to HDF5 File...' %rxn_name)
    dset = f.create_dataset('/%s' %rxn_name,
                            data=oxygen_network.reactions[rxn_name].coeff_fn(T).astype('float64'))
    # Store the reaction equation as an attribute of the dataset
    dset.attrs.create('reaction_equation', data=oxygen_network.reactions[rxn_name].__str__())
f.close()

# Based on the required species, create a list of cooling rates
cooling_rate_list = []
for s in oxygen_network.required_species:
    if s.name != 'de':
        cooling_rate_list.append('%s_c' %s.name)
cooling_rate_list = na.sort(cooling_rate_list)

# Write cooling rates to HDF5 file
f = h5py.File('oxygen_cooling_tables.h5', 'w')
for r in cooling_rate_list:
    print('Cooling rate for %s is being written to HDF5 File...' %r)
    dset = f.create_dataset('/%s' %r,
                            data=cooling_rates_table['%s' %r].values.astype('float64'))
f.close()
