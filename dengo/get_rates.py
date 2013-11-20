# Gets rates for all formation reactions for a given species in the UMIST database

import numpy as np
import re

from reaction_classes import \
    reaction, \
    AtomicSpecies, \
    ChemicalSpecies, \
    MolecularSpecies, \
    registry_setup, \
    species_registry

from umist_rates import umist_rates

def _get_species(sn):
    species_name = "us_%s" % (sn)
    free_electrons = species_name.count('+')-species_name.count('-')
    return species_registry["us_%s_%s" % (sn, free_electrons + 1)]

def _create_reaction(rline):
    rtype = rline[1]
    rA = _get_species(rline[2])
    rB = _get_species(rline[3])
    pC = _get_species(rline[4])
    pD = _get_species(rline[5])
    a = float(rline[9])
    b = float(rline[10])
    g = float(rline[11])
    T_lower = float(rline[12])
    T_upper = float(rline[13])
    reactants = [(1, rA), (1, rB)]
    products = [(1, pC), (1, pD)]
    rname = "%s_plus_%s" % (rA.name, rB.name)

    if rline[6] != '':
        pE = _get_species(rline[6])
        products.append((1,pE))

    if rline[7] != '':
        pF = _get_species(rline[7])
        products.append((1,pF))

    if rtype == 'CP':
        @reaction(rname, reactants, products)
        def rxn(state):
            return a
    elif rtype == "CR":
        @reaction(rname, reactants, products)
        def rxn(state):
            return 0.0
    else:
        @reaction(rname, reactants, products)
        def rxn(state):
            # rate coefficient with units cm^3 / s
            _i1 = (state.T < T_upper) & (state.T > T_lower)
            _i2 = ~_i1
            vals = a*(state.T**b)*(np.exp(-g / state.T)) 
            vals[_i2] = 0.0
            return vals

@registry_setup
def setup_umist_species(species_symbol, atomic_weight):
    species_name = "us_%s" % (species_symbol)
    free_electrons = species_name.count('+')-species_name.count('-')
    my_sym = MolecularSpecies(species_name, atomic_weight, free_electrons)

@registry_setup
def setup_umist_reactions(allowed_species):
    umist_names = set([s[3:].rsplit("_", 1)[0] for s in allowed_species])
    #umist_names.add("PHOTON")
    for line in open("RATE12.txt", "r"):
        rline = re.split(":?", line)
        if not all(r in umist_names for r in rline[2:6]): continue
        _create_reaction(rline)

def get_rates(atomicSymbol, atomicNumber, atomicWeight, network):
    speciesName = "us_%s" %(atomicSymbol)
    free_electrons = speciesName.count('+')-speciesName.count('-')
    if speciesName not in species_registry:
        s = Species(speciesName, atomicNumber, atomicWeight, free_electrons)
    else:
        s = species_registry[speciesName]

    umist_rates(s, network)
