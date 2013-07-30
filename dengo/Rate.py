import math
from reaction_classes import Species, species_registry

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
    
    if (reactantA in species_registry) == False:
        rA = Species(speciesName, atomicNumber, atomicWeight, i)
    else:
        rA = species_registry[speciesName]
    if (reactantB in species_registry) == False:
        rB = Species(speciesName, atomicNumber, atomicWeight, i)
    else:
        rB = species_registry[speciesName]
    if (productC in species_registry) == False:
        pC = Species(speciesName, atomicNumber, atomicWeight, i)
    else:
        pC = species_registry[speciesName]
    
    if (productD in species_registry) == False:
        pD = Species(speciesName, atomicNumber, atomicWeight, i)
    else:
        pD = species_registry[speciesName]
    
    if (productE in species_registry) == False:
        pE = Species(speciesName, atomicNumber, atomicWeight, i)
    else:
        pE = species_registry[speciesName]
   
    if (productF in species_registry) == False:
        pF = Species(speciesName, atomicNumber, atomicWeight, i)
    else:
        pF = species_registry[speciesName]
   
    if (reactantB in species_registry) == False:
        rB = Species(speciesName, atomicNumber, atomicWeight, i)
    else:
        rB = species_registry[speciesName]
        
    
        
    if type == 'CP':
        rate = a # rate coefficient with units 1/s
        units = "1/s"
        return rate, units
    elif type == 'CR':
        pass
    else:
        t = float(temp) / float(300)
        rate = a*(math.pow(t,b))*(math.exp(-g / temp)) # rate coefficient with units cm^3 / s
        units = "cm^3 / s"
        return rate, units
    Reaction("%s+%s" % (rA.name, rB.name), rate, [(1, rA.name), (1, rB.name)], [(1, prodC), (1, prodD), (1, prodE), (1,prodF)]) 