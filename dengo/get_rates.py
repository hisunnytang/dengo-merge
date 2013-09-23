# Gets rates for all formation reactions for a given species in the UMIST database

from reaction_classes import Species, species_registry
from umist_rates import umist_rates

def get_rates(atomicSymbol, atomicWeight, atomicNumber):
    speciesName = "%s" %(atomicSymbol)

    if speciesName not in species_registry:
        s = Species(speciesName, atomicNumber, atomicWeight, 0)
    else:
        s = species_registry[speciesName]

    umist_rates(s)
    

get_rates("CO", 28, -1)
get_rates("C", 12, -1)
get_rates("O", 16, -1)
get_rates("OH", 17, -1)
get_rates("H", 1, -1)
get_rates("H2", 2, -2)
get_rates("H2O", 18, -1)
get_rates("O2", 32, -1)
