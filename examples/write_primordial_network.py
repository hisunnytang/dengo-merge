from dengo.chemical_network import \
    ChemicalNetwork, \
    reaction_registry, \
    cooling_registry
import dengo.primordial_rates, dengo.primordial_cooling
from dengo.write_rate_reader import create_rate_tables, create_rate_reader
from dengo.write_rate_reader import create_rate_tables, create_rate_reader

primordial = ChemicalNetwork()
for ca in cooling_registry.values():
    primordial.add_cooling(ca)

for s in reaction_registry.values():
    primordial.add_reaction(s)

primordial.init_temperature((1e4, 1e8))

create_rate_tables(primordial, "primordial")
create_rate_reader(primordial, "primordial")
