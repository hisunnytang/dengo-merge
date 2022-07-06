import numpy as np
from reaction_classes import Species, species_registry, Reaction

# Check whether the species exists and if not, add to species registry
def _ensure_species(sp):
    # "us" => "umist species"
    sp = "us_" + sp
    sp = sp.replace("+", "p").replace("-", "m")
    if sp not in species_registry:
        i = sp.count("p") - sp.count("m")
        sp = Species(sp, 1.0, 1.0, i)
    else:
        sp = species_registry[sp]
    return sp


def get_rate(reaction, network):
    type = reaction[1]
    rA = _ensure_species(reaction[2])
    rB = _ensure_species(reaction[3])
    pC = _ensure_species(reaction[4])
    pD = _ensure_species(reaction[5])
    a = float(reaction[9])
    b = float(reaction[10])
    g = float(reaction[11])
    T_lower = int(reaction[12])
    T_upper = int(reaction[13])
    reactants = [(1, rA), (1, rB)]
    products = [(1, pC), (1, pD)]

    temp = network.T

    if reaction[6] != "":
        pE = _ensure_species(reaction[6])
        products.append((1, pE))

    if reaction[7] != "":
        pF = _ensure_species(reaction[7])
        products.append((1, pF))

    if np.all(T_lower <= temp) and np.all(temp <= T_upper):
        if type == "CP":
            rate = lambda network: a  # rate coefficient with units 1/s
            units = "1/s"
            # return rate, units
        elif type == "CR":
            rate = 0.0
            units = ""
            # return rate, units
        else:
            rate = (
                lambda network: a * (network.T**b) * (np.exp(-g / network.T))
            )  # rate coefficient with units cm^3 / s
            units = "cm^3 / s"
            # return rate, units
        return Reaction("%s_plus_%s" % (rA.name, rB.name), rate, reactants, products)
    else:
        rate = 0
        units = ""
        return rate, units
