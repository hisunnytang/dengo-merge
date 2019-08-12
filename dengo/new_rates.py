import numpy as np
from .chemistry_constants import tevk, tiny, mh
from .reaction_classes import \
    reaction, \
    AtomicSpecies, \
    ChemicalSpecies, \
    MolecularSpecies, \
    registry_setup, \
    species_registry

@registry_setup
def setup_new():

#    HI = AtomicSpecies("H", 0.0)
#    HII = AtomicSpecies("H", 1.0)
#    HeI = AtomicSpecies("He", 0.0)
#    HeII = AtomicSpecies("He", 1.0)
#    de = species_registry['de']
#

#   Let's implement the Briggs Rauscher Reaction!
#   http://www.math.udel.edu/~rossi/Math512/2003/br5.pdf

    HOI = MolecularSpecies("HOI", 1+16+126 ,0.0 )
    Im  = AtomicSpecies("I", -1.0)
    Hp  = AtomicSpecies("H", 1.0)

    I2  = MolecularSpecies("I2", 126*2, 0.0)
    H2O = MolecularSpecies("H2O", 2 + 16, 0.0)

    HOIO = MolecularSpecies("HOIO", 1 + 16 + 126 + 16, 0.0)
    IO3m = MolecularSpecies("IO3m", 1 + 16*3, -1)

    O2 = MolecularSpecies("O2", 16*2, 0.0)

    CH2_COOH2 = MolecularSpecies("CH2_COOH2", 12+2+(12+16*2+1)*2, 0)
    CHI_COOH2 = MolecularSpecies("CHI_COOH2", 12+126+(12+16*2+1)*2, 0)
    H2O2 = MolecularSpecies("H2O2", (2+16)*2, 0)

    # -- I01 --
    @reaction('I01', [   (1,HOI),   (1,Im), (1, Hp)], [  (1,I2),   (1, H2O)])
    def rxn(state):
       return 3.1e12* np.ones_like(state.T)


    # -- I01_re --
    @reaction('I01_re', [(1,I2),   (1, H2O)], [(1,HOI),   (1,Im), (1, Hp)] )
    def rxn(state):
        return 2.2* np.ones_like(state.T)

    # -- I02 --
    @reaction('I02', [(1, Hp),   (1, HOIO), (1, Im)], [(2, HOI)] )
    def rxn(state):
        return 5.0e9* np.ones_like(state.T)

    # -- I03 --
    @reaction('I03', [(2, Hp),   (1, IO3m ), (1, Im)], [(1, HOIO), (1, HOI)] )
    def rxn(state):
        return 1.4e3* np.ones_like(state.T)

    # -- I04 --
    @reaction('I04', [(2, HOIO) ], [(1, Hp), (1, IO3m), (1, HOI)] )
    def rxn(state):
        return 3.0e9* np.ones_like(state.T)

    # -- I05 --
    @reaction('I05', [(1, Hp), (1, IO3m), (1, HOIO) ], [(2, HOIO), (0.5, O2) ] )
    def rxn(state):
        return 2.6e5* np.ones_like(state.T)

    # -- C05 --
    @reaction('C05', [(1, CH2_COOH2), (1, I2) ], [(1, CHI_COOH2), (1, Im), (1, Hp)] )
    def rxn(state):
        return 3.494* np.ones_like(state.T)

    # -- D01 --
    @reaction('D01', [(1, HOI), (1, H2O2) ], [(1, Im), (1, O2), (1, Hp), (1, H2O)] )
    def rxn(state):
        return 2.0e3* np.ones_like(state.T)





