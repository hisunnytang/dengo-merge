import os
import numpy as np
from dengo.chemical_network import \
    ChemicalNetwork, \
    reaction_registry, \
    cooling_registry, species_registry
import dengo.primordial_rates
import dengo.primordial_cooling
from dengo.chemistry_constants import tiny, kboltz, mh, G
import yt
import yt.units as u
import numpy
import pickle
import time
import pyximport
import sys
import matplotlib.pyplot as plt
import pytest
import h5py

def set_env_variables(var, path):
    if var not in os.environ:
        os.environ[var] = path


set_env_variables("HDF5_DIR", "/home/kwoksun2/anaconda3")
set_env_variables("CVODE_PATH", "/home/kwoksun2/cvode-3.1.0/instdir")
set_env_variables("HDF5_PATH", "/home/kwoksun2/anaconda3")
set_env_variables("SUITESPARSE_PATH", "/home/kwoksun2/SuiteSparse")
set_env_variables("DENGO_INSTALL_PATH", "/home/kwoksun2/dengo_install")


def freefall_time(density):
    return 1.0 / np.sqrt(G*mh*density)


@pytest.fixture
def setup_primordial_network():
    """Initial a ChemicalNetwork object
       for primordial network 9-species model
    Return:
        primordial: ChemicalNetwork with primordial reactions and cooling
    """
    # this register all the rates specified in `primordial_rates.py`
    dengo.primordial_rates.setup_primordial()

    # initialize the chmical network object
    primordial = ChemicalNetwork()

    # add all the reactions
    primordial.add_reaction("k01")
    primordial.add_reaction("k02")
    primordial.add_reaction("k03")
    primordial.add_reaction("k04")
    primordial.add_reaction("k05")
    primordial.add_reaction("k06")
    primordial.add_reaction("k07")
    primordial.add_reaction("k08")
    primordial.add_reaction("k09")
    primordial.add_reaction("k10")
    primordial.add_reaction("k11")
    primordial.add_reaction("k12")
    primordial.add_reaction("k13")
    primordial.add_reaction("k14")
    primordial.add_reaction("k15")
    primordial.add_reaction("k16")
    primordial.add_reaction("k17")
    primordial.add_reaction("k18")
    primordial.add_reaction("k19")
    primordial.add_reaction("k21")
    primordial.add_reaction("k22")
    primordial.add_reaction("k23")

    primordial.add_cooling("brem")
    primordial.add_cooling("reHII")
    primordial.add_cooling("reHeIII")
    primordial.add_cooling("gloverabel08")
    primordial.add_cooling("ceHI")
    primordial.add_cooling("h2formation")
    primordial.add_cooling("reHeII2")
    primordial.add_cooling("reHeII1")
    primordial.add_cooling("ciHeIS")
    primordial.add_cooling("ceHeII")
    primordial.add_cooling("ciHI")
    primordial.add_cooling("ceHeI")
    primordial.add_cooling("gammah")
    primordial.add_cooling("ciHeI")
    primordial.add_cooling("ciHeII")

    # This defines the temperature range for the rate tables
    primordial.init_temperature((1e0, 1e8))

    return primordial


@pytest.fixture
def setup_solver_options(request):
    solver_options = {"output_dir": "temp",
                      "solver_name": "primordial",
                      "use_omp": False,
                      "use_cvode": False,
                      "use_suitesparse": False,
                      "niters": 1e3,
                      "NCELLS": 128,
                      "reltol": 1.0e-6}
    solver_options.update(request.param)
    return solver_options


def setup_initial_conditions(network, density, temperature, h2frac, NCELLS):
    # setting initial conditions
    temperature = np.ones((NCELLS))*temperature
    init_array = np.ones(NCELLS) * density
    init_values = dict()
    init_values["H_1"] = init_array * 0.76
    init_values['H_2'] = init_array * tiny
    init_values['H_m0'] = init_array * tiny
    init_values['He_1'] = init_array * 0.24
    init_values['He_2'] = init_array * tiny
    init_values['He_3'] = init_array * tiny
    init_values['H2_1'] = init_array * h2frac
    init_values['H2_2'] = init_array * tiny
    init_values['de'] = init_array * tiny

    # update and calculate electron density and etc with the handy functions
    # init_values = primordial.convert_to_mass_density(init_values)
    init_values['de'] = network.calculate_free_electrons(init_values)
    # init_values['density'] = primordial.calculate_total_density(init_values)
    init_values['density'] = np.ones((NCELLS))*density
    number_density = network.calculate_number_density(init_values)

    # set up initial temperatures values used to define ge
    init_values['T'] = temperature

    # calculate ge (very crudely, no H2 help here)
    gamma = 5.0/3.0
    mH = 1.67e-24
    init_values["ge"] = 3.0 / 2.0 * temperature * kboltz / mH

    return init_values


def write_network(network, solver_options={"output_dir": "test_dir",
                                           "solver_name": "primordial",
                                           "use_omp": False,
                                           "use_cvode": False,
                                           "use_suitesparse": False}):
    """Write solver based on the ChemicalNetwork
    """
    # Write the initial conditions file
    # IF you need to use the Makefile, and c-library
    # you will have to specified the library_path

    output_dir = solver_options["output_dir"]
    solver_name = solver_options["solver_name"]
    use_omp = solver_options["use_omp"]
    use_cvode = solver_options["use_cvode"]
    use_suitesparse = solver_options["use_suitesparse"]

    if use_cvode:
        network.write_solver(solver_name, output_dir=output_dir,
                             solver_template="cv_omp/sundials_CVDls",
                             ode_solver_source="initialize_cvode_solver.C")
        return

    if not(use_omp and use_cvode and use_suitesparse):
        network.write_solver(
            solver_name, output_dir=output_dir,
            solver_template="be_chem_solve/rates_and_rate_tables",
            ode_solver_source="BE_chem_solve.C")
        return


def run_solver(init_values, solver_options, dtf=None,
               make_plot=True, intermediate=True):

    print(os.getcwd())
    solver_name = solver_options["solver_name"]
    solver_dir = solver_options["output_dir"]
    niters = solver_options["niters"]
    reltol = solver_options["reltol"]
    pyximport.install(setup_args={"include_dirs": np.get_include()},
                      reload_support=True, inplace=True)

    _solver_run = pyximport.load_module(
        "{}_solver_run".format(solver_name),
        "{}_solver_run.pyx".format(solver_name),
        build_inplace=True, pyxbuild_dir="_dengo_temp")

    if dtf is None:
        dtf = freefall_time(init_values["density"])
    rv, rv_int = eval(
        "_solver_run.run_{}(init_values, dtf, niter={}, reltol = {}, floor_value = 1.0e-20)".format(
            solver_name, niters, reltol))

    mask = rv_int['successful']
    for name in sorted(rv_int):
        if len(rv_int[name].shape) == 1:
            rv_int[name] = rv_int[name][mask]
        else:
            rv_int[name] = rv_int[name][0, mask]

    if make_plot:
        plt.clf()
        skip = ('successful', 'dt', 't', 'ge')
        for n, v in sorted(rv_int.items()):
            if n in skip:
                continue
            plt.loglog(rv_int['t'], v, label=n)
            print(n, "{0:0.5g}".format(v[-1]))
        # plt.ylim(density * 1e-20, density * 100)
        plt.xlabel("time [s]")
        plt.legend(loc='best', fontsize='xx-small')
        plt.savefig("plot_{}_network.png".format(solver_name))

    if intermediate:
        return rv_int
    else:
        try:
            for n, v in sorted(rv_int.items()):
                rv_int[n] = np.array([v[-1]])
            return rv_int
        except BaseException:
            return rv_int


def run_c_solver(solver_options, dt=100.0):
    solver_dir = solver_options["output_dir"]
    solver_name = solver_options["solver_name"]
    os.chdir(solver_dir)

    os.system("make clean")
    os.system("make")
    os.system("make test")
    out = os.system(
        "./run_dengo {}_ic.h5 {}_sol.h5 {}".format(solver_name, solver_name, dt))
    print(out)
    os.chdir("../")


def write_init_to_file(init_values, solver_options):
    solver_name = solver_options["solver_name"]
    os.chdir(solver_options["output_dir"])
    print(solver_options["output_dir"])
    print(os.getcwd())
    icname = "{}_ic.h5".format(solver_name)
    f = h5py.File(icname, "w")
    for k, v in init_values.items():
        print(k,v)
        f.create_dataset(k, data=v)
    f.close()
    os.chdir("../")
