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
import h5py
from utilities import setup_solver_options, write_network,\
    run_solver, set_env_variables, run_c_solver, write_init_to_file

matplotlib.use("Agg")

set_env_variables("HDF5_DIR", "/home/kwoksun2/anaconda3")
set_env_variables("CVODE_PATH", "/home/kwoksun2/cvode-3.1.0/instdir")
set_env_variables("HDF5_PATH", "/home/kwoksun2/anaconda3")
set_env_variables("SUITESPARSE_PATH", "/home/kwoksun2/SuiteSparse")
set_env_variables("DENGO_INSTALL_PATH", "/home/kwoksun2/dengo_install")

solver_name = "predator_prey"
output_dir = "temp_prey_predator"

pytest_dir = os.getcwd()

alpha = 2./3.
beta = 4./3.
gamma = 1.0
delta = 1.0

predator0 = 2.0
prey0 = 2.0


@registry_setup
def setup_predator_prey_rates():
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


@pytest.fixture
def setup_predator_prey_network():
    setup_predator_prey_rates()
    cN = ChemicalNetwork()
    cN.add_reaction("exp_growth_prey")
    cN.add_reaction("predation")
    cN.add_reaction("natural_death_predator")

    # this shouldnt be compulsory...
    cN.init_temperature((1e0, 1e8))
    return cN


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


def TestConvergence(init_values, rtol_array, setup_solver_options):
    f, ax = plt.subplots()
    init = copy.deepcopy(init_values)
    for reltol in rtol_array:
        setup_solver_options["reltol"] = reltol
        results = run_solver(
            init,
            setup_solver_options,
            make_plot=False,
            dtf=10.0)
        v0, vt = conserved_variables(results, make_plot=False)
        ones = np.ones_like(vt)
        # this shows the tolerance level of the solver
        np.testing.assert_array_almost_equal(ones, vt/v0, decimal=2)
        ax.semilogy(np.abs((vt - v0)/v0),
                    label="rtol = {0:0.1e}".format(reltol))
    ax.legend()
    f.savefig("prey_predator_convergence.png")

def check_output_conservation(solver_options):
    solver_name = solver_options["solver_name"]
    solver_dir = solver_options["output_dir"]
    os.chdir(solver_dir)

    fn = "{}_sol.h5".format(solver_name)
    outfile = h5py.File(fn,'r')
    results = {}
    for k in outfile.keys():
        results[k] = outfile[k].value
    v0, vt = conserved_variables(results, make_plot=False)
    ones = np.ones_like(vt)
    np.testing.assert_array_almost_equal(ones, vt/v0, decimal=2)

    outfile.close()

    os.chdir("../")

@pytest.mark.parametrize('setup_solver_options',
                         ({"use_omp": False, "use_cvode": False,
                           "use_suitesparse": False,
                           "output_dir": "be_chem_solve"},
                          {"use_omp": False, "use_cvode": True,
                           "use_suitesparse": False,
                           "output_dir": "cvode_dls"},
                          {"use_omp": True, "use_cvode": True,
                           "use_suitesparse": True,
                           "output_dir": "cvode_klu"}),
                         indirect=True)
def test_prey_predator(setup_predator_prey_network, setup_solver_options):
    setup_solver_options["solver_name"] = solver_name
    write_network(setup_predator_prey_network,
                  setup_solver_options)
    init_values = write_initial_conditions(setup_predator_prey_network)

    os.chdir(pytest_dir)
    os.chdir(setup_solver_options["output_dir"])
    results = run_solver(
        init_values,
        setup_solver_options,
        make_plot=False,
        dtf=10.0)
    phase_plot(results)
    os.chdir("../")


#@pytest.mark.xfail(reason="relative level of the tolerance?")
@pytest.mark.parametrize(
    'setup_solver_options',
    (pytest.param(
         {"use_omp": False, "use_cvode": False, "use_suitesparse": False,
          "output_dir": "be_chem_solve"},
         marks=pytest.mark.xfail(
             reason="BE_chem_solve is unstable to stiff equations?")),
     {"use_omp": False, "use_cvode": True, "use_suitesparse": False,
      "output_dir": "cvode_dls"},
     {"use_omp": True, "use_cvode": True, "use_suitesparse": True,
      "output_dir": "cvode_klu"}),
    indirect=True)
def test_prey_predator_convergence(
        setup_predator_prey_network, setup_solver_options):
    init_values = write_initial_conditions(setup_predator_prey_network)
    setup_solver_options["solver_name"] = solver_name
    #
    os.chdir(pytest_dir)
    os.chdir(setup_solver_options["output_dir"])
    results = run_solver(
        init_values,
        setup_solver_options,
        make_plot=False,
        dtf=10.0)
    phase_plot(results)
    TestConvergence(
        init_values, rtol_array=np.logspace(-8, -4, 5),
        setup_solver_options=setup_solver_options)
    os.chdir("../")


@pytest.mark.last
@pytest.mark.parametrize('setup_solver_options',
                         ({"use_omp": False, "use_cvode": True,
                           "use_suitesparse": False,
                           "output_dir": "cvode_dls"},
                          {"use_omp": True, "use_cvode": True,
                           "use_suitesparse": True,
                           "output_dir": "cvode_klu"}),
                         indirect=True)
def test_c_solver(setup_predator_prey_network, setup_solver_options):
    setup_solver_options["solver_name"] = "prey_predator"
    init_values = write_initial_conditions(setup_predator_prey_network)

    os.chdir(pytest_dir)
    write_init_to_file(init_values, setup_solver_options)
    run_c_solver(setup_solver_options)
    check_output_conservation(setup_solver_options)

def write_initial_conditions(cN, N=1024):
    init_values = {}

    init_values["predator"] = np.ones((N)) * predator0
    init_values["prey"] = np.ones((N)) * prey0

    init_values["dead_predator"] = np.ones((N)) * 1.0e-4
    init_values["dead_prey"] = np.ones((N)) * 1.0e-4

    init_values['density'] = cN.calculate_total_density(init_values)
    init_values['ge'] = np.ones((N)) * 10.0

    return init_values
