#Parsing the UMIST database

import re
import math

class umist_reaction(object):

    # Holds information for reactions from the UMIST database.

    # type = code for type of reaction, e.g. IN (ion-neutral)
    # reactantA = reactant
    # reactantB = another reactant
    # prodC = product
    # prodD = another product
    # prodE = another product
    # prodF = another product
    # alpha, beta, gamma = coefficients for calculation of rate coefficient from temperature

    def __init__(self):
        self.type = None
        self.reactantA = None
        self.reactantB = None
        self.prodC = None
        self.prodD = None
        self.prodE = None
        self.prodF = None
        self.alpha = None
        self.beta = None
        self.gamma = None
 


def get_umistrates(network):
    
    with open("RATE12.txt", "r") as f:
        umist = f.read()
        reactions = []
        reactions = re.split("\n+", umist)

        for line in reactions:

            reaction = re.split(":?", line)
            umist_reaction = umist_reaction()

            umist_reaction.type = reaction[1]
            umist_reaction.reactantA = reaction[2]
            umist_reaction.reactantB = reaction[3]
            umist_reaction.prodC = reaction[4]
            umist_reaction.prodD = reaction[5]
            umist_reaction.prodE = reaction[6]
            umist_reaction.prodF = reaction[7]
            umist_reaction.alpha = float(reaction[9])
            umist_reaction.beta = float(reaction[10])
            umist_reaction.gamma = float(reaction[11])


            temp = network.T
            if umist_reaction.type == 'CP':
                rate = umist_rxn.alpha # rate coefficient with units 1/s
                return rate
            elif umist_reaction.type == 'CR':
                pass
            else:
                a = umist_reaction.alpha
                b = umist_reaction.beta
                g = umist_rxn.gamma
                t = float(temp) / float(300)
                rate = a*(math.pow(t,b))*(math.exp(-g / temp)) # rate coefficient with units cm^3 / s
                return rate
