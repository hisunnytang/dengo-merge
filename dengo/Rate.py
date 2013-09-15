import math
from reaction_classes import Species, species_registry, Reaction

def get_rate(reaction, temp):
    type = reaction[1]
    reactantA = reaction[2]
    reactantB = reaction[3]
    prodC = reaction[4]
    prodD = reaction[5]
    prodE = reaction[6]
    prodF = reaction[7]
    a = float(reaction[9])
    b = float(reaction[10])
    g = float(reaction[11])
    T_lower = int(reaction[12])
    T_upper = int(reaction[13])
    
    if (temp >= T_lower and temp <= T_upper):
        if (reactantA in species_registry) == False:
            rA = Species(reactantA, -1, -1, -1)
        else:
            rA = species_registry[reactantA]
        if (reactantB in species_registry) == False:
            rB = Species(reactantB, -1, -1, -1)
        else:
            rB = species_registry[reactantB]
        if (prodC in species_registry) == False:
            pC = Species(prodC, -1, -1, -1)
        else:
            pC = species_registry[prodC]
    
        if (prodD in species_registry) == False:
            pD = Species(prodD, -1, -1, -1)
        else:
            pD = species_registry[prodD]
    
        if (prodE in species_registry) == False:
            pE = Species(prodE, -1, -1, -1)
        else:
            pE = species_registry[prodE]
   
        if (prodF in species_registry) == False:
            pF = Species(prodF, -1, -1, -1)
        else:
            pF = species_registry[prodF]
        
        
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