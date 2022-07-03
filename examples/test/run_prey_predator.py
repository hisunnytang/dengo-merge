import numpy as np
from dengo.reaction_classes import (
    reaction,
    ChemicalSpecies,
    registry_setup,
    species_registry,
)
from dengo.chemical_network import ChemicalNetwork
import os
import pyximport
import matplotlib
import matplotlib.pyplot as plt
import copy
import pytest
import h5py
from utilities import (
    setup_solver_options,
    write_network,
    run_solver,
    check_defined_envpath,
    run_c_solver,
    write_init_to_file,
)

plt.switch_backend("agg")


solver_name = "predator_prey"
output_dir = "temp_prey_predator"
pytest_dir = os.getcwd()


# parameters for the prey-predator model
alpha = 2.0 / 3.0
beta = 4.0 / 3.0
gamma = 1.0
delta = 1.0

# initial conditions
predator0 = 2.0
prey0 = 2.0


@registry_setup
def setup_predator_prey_rates():
    predator = ChemicalSpecies("predator", 1.0)
    dead_predator = ChemicalSpecies("dead_predator", 1.0)
    prey = ChemicalSpecies("prey", 1.0)
    dead_prey = ChemicalSpecies("dead_prey", 1.0)

    # predator-prey model
    @reaction(
        "exp_growth_prey",
        [
            (1, prey),
        ],
        [
            (2, prey),
        ],
    )
    def rxn(state):
        return alpha * np.ones_like(state.T)

    @reaction(
        "predation",
        [(1, predator), (1, prey)],
        [(1, dead_prey), (gamma / beta + 1, predator)],
    )
    def rxn(state):
        return beta * np.ones_like(state.T)

    @reaction(
        "natural_death_predator",
        [
            (1, predator),
        ],
        [
            (1, dead_predator),
        ],
    )
    def rxn(state):
        return gamma * np.ones_like(state.T)


def predator_prey_network():
    setup_predator_prey_rates()
    cN = ChemicalNetwork()
    cN.add_reaction("exp_growth_prey")
    cN.add_reaction("predation")
    cN.add_reaction("natural_death_predator")

    # this shouldnt be compulsory...
    cN.init_temperature((1e0, 1e8))
    return cN

from test_prey_predator import run_prey_predator
if __name__ == "__main__":
    solver_options = {"output_dir": "temp",
                      "solver_name": "primordial",
                      "use_omp": False,
                      "use_cvode": False,
                      "use_suitesparse": False,
                      "niters": 1e3,
                      "NCELLS": 128,
                      "reltol": 1.0e-6}
    run_prey_predator(
        predator_prey_network(),
        solver_options
    )
