import numpy as np
from dengo.chemistry_constants import tevk, tiny, mh
from dengo.reaction_classes import \
    reaction, \
    AtomicSpecies, \
    ChemicalSpecies, \
    MolecularSpecies, \
    registry_setup, \
    species_registry
from dengo.chemical_network import ChemicalNetwork
import os
import pyximport
import matplotlib
import matplotlib.pyplot as plt
import copy
import pytest

matplotlib.use("Agg")

def set_env_variables(var, path):
    if var not in os.environ:
        os.environ[var] = path


set_env_variables("HDF5_DIR", "/home/kwoksun2/anaconda3")
set_env_variables("CVODE_PATH", "/home/kwoksun2/cvode-3.1.0/instdir")
set_env_variables("HDF5_PATH", "/home/kwoksun2/anaconda3")
set_env_variables("SUITESPARSE_PATH", "/home/kwoksun2/SuiteSparse")
set_env_variables("DENGO_INSTALL_PATH", "/home/kwoksun2/dengo_install")

solver_name = "predator_prey"
output_dir = "temp_prey_predator"

alpha = 2./3.
beta = 4./3.
gamma = 1.0
delta = 1.0

predator0 = 2.0
prey0 = 2.0


@pytest.fixture
def init_values(request):
    print(request.param)
    setup_predator_prey(*request.param)
    cN = write_predator_prey_model(write_file=True)
    yield write_initial_conditions(cN)


@pytest.fixture
def setup_solver_options(update_options={}):
    solver_options = {"output_dir": "./{}".format(output_dir),
                      "solver_name": "primordial",
                      "use_omp": False,
                      "use_cvode": False,
                      "use_suitesparse": False,
                      "niters": 1e3,
                      "NCELLS": 128,
                      "reltol": 1.0e-6}
    solver_options.update(update_options)
    return solver_options



@registry_setup
def setup_predator_prey(alpha, beta, gamma, delta):

    predator = ChemicalSpecies("predator", 1.0)
    dead_predator = ChemicalSpecies("dead_predator", 1.0)
    prey = ChemicalSpecies("prey", 1.0)
    dead_prey = ChemicalSpecies("dead_prey", 1.0)

    # predator-prey model
    @reaction("exp_growth_prey", [(1, prey), ], [(2, prey), ])
    def rxn(state):
        return alpha * np.ones_like(state.T)

    @reaction("predation", [(1, predator), (1, prey)], [
              (1, dead_prey), (gamma/beta + 1, predator)])
    def rxn(state):
        return beta * np.ones_like(state.T)

    @reaction(
        "natural_death_predator", [(1, predator), ],
        [(1, dead_predator), ])
    def rxn(state):
        return gamma * np.ones_like(state.T)


def write_predator_prey_model(write_file=True):
    cN = ChemicalNetwork()
    cN.add_reaction("exp_growth_prey")
    cN.add_reaction("predation")
    cN.add_reaction("natural_death_predator")

    # this shouldnt be compulsory...
    cN.init_temperature((1e0, 1e8))
    init_values = write_initial_conditions(cN)
    if write_file:
        cN.write_solver(solver_name, output_dir=output_dir,
                        solver_template="cv_omp/sundials_CVDls",
                        ode_solver_source="initialize_cvode_solver.C",
                        init_values = init_values)
    return cN


def write_initial_conditions(cN):
    init_values = {}
    N = 1

    init_values["predator"] = np.ones((N)) * predator0
    init_values["prey"] = np.ones((N)) * prey0

    init_values["dead_predator"] = np.ones((N)) * 0.0
    init_values["dead_prey"] = np.ones((N)) * 0.0

    init_values['density'] = cN.calculate_total_density(init_values)
    init_values['ge'] = np.ones((N)) * 10.0

    return init_values



def run_model(init_values, reltol=1.0e-5, make_plot=True, dtf=5.0e1):
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    os.chdir(output_dir)
    print(os.getcwd())
    pyximport.install(setup_args={"include_dirs": np.get_include()},
                      reload_support=True, inplace=True, language_level=3)

    _solver_run = pyximport.load_module(
        "{}_solver_run".format(solver_name),
        "{}_solver_run.pyx".format(solver_name),
        build_inplace=True, pyxbuild_dir="_dengo_temp")
    print("finished loading!")
    rv, rv_int = _solver_run.run_predator_prey(
        init_values, dtf, niter=5e2, adaptive_step=False, reltol=reltol)
    mask = rv_int['successful']
    for name in sorted(rv_int):
        if len(rv_int[name].shape) == 1:
            rv_int[name] = rv_int[name][mask]
        else:
            rv_int[name] = rv_int[name][0, mask]

    if make_plot:
        plt.clf()
        skip = (
            'successful',
            'dt',
            't',
            'ge',
            "T",
            "dead_predator",
            "dead_prey")
        for n, v in sorted(rv_int.items()):
            if n in skip:
                continue
            plt.plot(rv_int['t'], v, label=n)
        plt.xlabel("time [s]")
        plt.legend(loc='best', fontsize='small')
        plt.savefig("prey_predator.png")
    os.chdir("../")
    return rv_int


def prey_predator_conserved_variables(prey, predator):
    V = delta*prey - gamma*np.log(prey) + beta * \
                                  predator - alpha*np.log(predator)
    return V


def phase_plot(results, filename=None):
    prey = results["prey"]
    predator = results["predator"]
    V0 = prey_predator_conserved_variables(prey0, predator0)

    def Vt(p0, p1): return prey_predator_conserved_variables(p0, p1) - V0

    x = np.linspace(0, 5, 1000)
    y = np.linspace(0, 5, 1000)
    X, Y = np.meshgrid(x, y)
    res2d = Vt(X, Y)
    tolerance = 1.0e-2

    flag = np.logical_and(res2d > -tolerance, res2d < tolerance)
    plt.clf()
    plt.scatter(X[flag], Y[flag])
    plt.scatter(prey, predator, color="C1", s=1)
    plt.xlabel("Prey")
    plt.ylabel("Predator")
    if filename is None:
        plt.savefig("prey_predator_phase.png")
    else:
        plt.savefig(filename)


def conserved_variables(results, filename=None, make_plot=True):
    prey = results["prey"]
    predator = results["predator"]
    V0 = prey_predator_conserved_variables(prey0, predator0)
    Vt = prey_predator_conserved_variables(prey, predator)
    if not make_plot:
        return V0, Vt
    plt.clf()
    plt.plot(np.abs((Vt - V0)/V0*100))
    plt.xlabel("Time")
    plt.ylabel("Percentage Error of the conserved variable")
    if filename is None:
        plt.savefig("prey_predator_conserved_variables.png")
    else:
        plt.savefig(filename)


def TestConvergence(init_values, rtol_array):
    f, ax = plt.subplots()
    init = copy.deepcopy(init_values)
    for reltol in rtol_array:
        results = run_model(init, reltol=reltol, make_plot=False)
        v0, vt = conserved_variables(results, make_plot=False)
        ax.semilogy(np.abs((vt - v0)/v0),
                    label="rtol = {0:0.1e}".format(reltol))
    ax.legend()
    f.savefig("prey_predator_convergence.png")


@pytest.mark.parametrize(
    "init_values",
    ([alpha,beta,gamma,delta],),
    indirect=True
)
def test_prey_predator(init_values):
    results = run_model(init_values)
    phase_plot(results)
    TestConvergence(init_values, rtol_array=np.logspace(-8, -3, 6))
