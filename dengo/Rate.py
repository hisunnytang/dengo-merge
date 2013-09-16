import math
from reaction_classes import Species, species_registry, Reaction

#Check whether the species exists and if not, add to species registry
def _ensure_species(sp):
    if sp not in species_registry:
        sp = Species(sp,-1,-1,-1)
    else:
        sp = speces_registry[sp]
    return sp

def get_rate(reaction, temp):
    type = reaction[1]
    reactantA = _ensure_species(reaction[2])
    reactantB = _ensure_species(reaction[3])
    prodC = _ensure_species(reaction[4])
    prodD = _ensure_species(reaction[5])
    prodE = _ensure_species(reaction[6])
    prodF = _ensure_species(reaction[7])
    a = float(reaction[9])
    b = float(reaction[10])
    g = float(reaction[11])
    T_lower = int(reaction[12])
    T_upper = int(reaction[13])
    
    if (temp >= T_lower and temp <= T_upper):
        
        if type == 'CP':
            rate = a # rate coefficient with units 1/s
            units = "1/s"
            return rate, units
        elif type == 'CR':
            rate = 0.0
            units = ''
            return rate, units
        else:
            t = float(temp) / float(300)
            rate = a*(math.pow(t,b))*(math.exp(-g / temp)) # rate coefficient with units cm^3 / s
            units = "cm^3 / s"
            return rate, units
        return Reaction("%s+%s" % (rA.name, rB.name), rate, [(1, rA.name), (1, rB.name)], [(1, pC.name), (1, pD.name), (1, pE.name), (1,pF.name)]) 
    else:
        rate = 0
        units = ''
        return rate, units