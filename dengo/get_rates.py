# Gets rates for all formation reactions for a given species in the UMIST database

from reaction_classes import Species, species_registry
from umist_rates import umist_rates

def get_rates(atomicSymbol, atomicWeight, atomicNumber, network):
    speciesName = "us_%s" %(atomicSymbol)
    free_electrons = speciesName.count('+')-speciesName.count('-')
    if speciesName not in species_registry:
        s = Species(speciesName, atomicNumber, atomicWeight, free_electrons)
    else:
        s = species_registry[speciesName]

    umist_rates(s, network)
