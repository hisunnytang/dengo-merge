import numpy as np
from dengo.chemical_network import \
    ChemicalNetwork, \
    reaction_registry, \
    cooling_registry
import dengo.primordial_rates
import dengo.primordial_cooling
from dengo.chemistry_constants import tiny, kboltz, mh, G
import pyximport
import os
import pylab
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pytest
from utilities import setup_primordial_network, setup_solver_options, write_network, setup_initial_conditions, run_solver,freefall_time


output_dir = "test_primordial"
pytest_dir = os.getcwd()


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
def test_tolerance_convergence(setup_primordial_network, setup_solver_options):
    """Change number of iteration to check if the results converges
    """
    os.chdir(pytest_dir)
    write_network(
        setup_primordial_network,
        solver_options=setup_solver_options)
    output_dir = setup_solver_options["output_dir"]
    print(setup_solver_options)
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    os.chdir(output_dir)
    density = 1.0e0
    temperature = 2000.0
    h2frac = 1.0e-5
    NCELLS = 1

    rtol_array = np.logspace(-9, -4, 6)
    start = True
    skip = ('successful', 'dt', 't')
    for r in rtol_array:
        setup_solver_options["reltol"] = r
        init_values = setup_initial_conditions(
            setup_primordial_network, density, temperature, h2frac, NCELLS)
        r = run_solver(init_values, setup_solver_options, make_plot=False, intermediate=False)
        if start:
            for k, v in r.items():
                if k in skip:
                    continue
                r[k] = [float(r[k])]
            results = r
            print(results)
            start = False
        else:
            for k, v in r.items():
                if k in skip:
                    continue
                results[k].append(v[0])
    f, ax = plt.subplots()
    for k, v in results.items():
        if k in skip:
            continue
        # here we treat reltol = 1.0e-9 as ground truth
        # and compare it to the rest to see if it converges
        v = np.array(v)
        error = np.abs((v[1:] - v[0]) / v[0])
        ax.loglog(rtol_array[1:], error, label=k)
    ax.plot(rtol_array, rtol_array, ls='--', color='k')
    #ax.set_ylim(1.0e-9, 1e-2)
    #ax.set_xlim(1.0e-8, 1.0e-4)
    ax.set_xlabel("Relative Tolerance of the Solver")
    ax.set_ylabel("Fractional Differnce to rtol =1e-9")
    ax.legend(loc="best", fontsize="xx-small")
    f.savefig("primordial_tolerance_convergence.png")
    os.chdir("../")


def TestConservation(cN, results, density):
    """Check convergence of electrons, and mass, and species
    """
    # test charge conservation
    # de_ is the net charge of the rest of the species
    percentage_error = []
    de_ = cN.calculate_free_electrons(results)
    de_solver = results["de"]
    pdiff = (de_ - de_solver) / de_
    percentage_error.append(pdiff)

    # test mass conservation
    # conservation of Hydrogen mass
    true_hfrac = 0.76 * density
    hmass = results["H2_1"] + results["H2_2"] + \
        results["H_1"] + results["H_2"] + results["H_m0"]
    pdiff = (hmass - true_hfrac) / hmass
    percentage_error.append(pdiff)

    # conservation of helium mass
    true_hefrac = 0.24 * density
    hemass = results["He_1"] + results["He_2"] + results["He_3"]
    pdiff = (hemass - true_hefrac) / hemass
    percentage_error.append(pdiff)

    return percentage_error


def run_grid(network, solver_options, density, temperature, h2frac):
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    os.chdir(output_dir)
    NCELLS = solver_options["NCELLS"]
    perror = []
    start = True
    for d in density:
        for T in temperature:
            for f in h2frac:
                print("{:0.2E} {:0.2E} {:.2E}".format(d, T, f))
                init_values = setup_initial_conditions(
                    network, d, T, f, NCELLS)
                r = run_solver(
                    init_values,
                    solver_options=solver_options,
                    intermediate = False,
                    make_plot=False)
                perror.append(
                    TestConservation(
                        network,
                        r,
                        init_values["density"]))
                if start:
                    # there should be a smarter way to do this...
                    for k, v in r.items():
                        r[k] = [float(r[k])]
                    results = r
                    start = False
                else:
                    for k, v in results.items():
                        results[k].append(v[0])

    return results, perror


def main(density, temperature, h2frac,
         solver_options={"output_dir": ".",
                         "solver_name": "primordial",
                         "use_omp": False,
                         "use_cvode": False,
                         "use_suitesparse": False,
                         "niters": 1e3,
                         "NCELLS": 128,
                         "reltol": 1.0e-6}):
    NCELLS = solver_options["NCELLS"]

    network = setup_primordial_network()
    write_network(network, solver_options=solver_options)
    init_values = setup_initial_conditions(
        network, density, temperature, h2frac, NCELLS)
    print(init_values)
    results = run_solver(
        init_values,
        solver_options=solver_options,
        make_plot=False)
    TestConservation(network, results, density)


def runtime_ncells(setup_primordial_network, setup_solver_options):
    density = 1.0e14
    temperature = 2000.0
    h2frac = 1.0e-3
    solver_name = setup_solver_options["solver_name"]
    niters = setup_solver_options["niters"]
    reltol = setup_solver_options["reltol"]

    write_network(
        setup_primordial_network,
        solver_options=setup_solver_options)
    output_dir = setup_solver_options["output_dir"]
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    os.chdir(output_dir)
    pyximport.install(setup_args={"include_dirs": np.get_include()},
                      reload_support=True, inplace=True)

    _solver_run = pyximport.load_module(
        "{}_solver_run".format(solver_name),
        "{}_solver_run.pyx".format(solver_name),
        build_inplace=True, pyxbuild_dir="_dengo_temp")
    ncells = 1
    init_values = setup_initial_conditions(
           setup_primordial_network, density, temperature, h2frac, ncells)
    dtf = freefall_time(init_values["density"])

    from timeit import default_timer as timer
    for ncells in [2**(i*2) for i in range(6)]:
        init_values = setup_initial_conditions(
           setup_primordial_network, density, temperature, h2frac, ncells)

        start = timer()
        rv, rv_int = eval(
            "_solver_run.run_{}(init_values, dtf, niter={}, reltol = {}, floor_value = 1.0e-16)".format(
                solver_name, niters, reltol))
        end = timer()
        print("ncells = {:d}, with {:0.2e} s\n".format(
                  ncells, (end - start) / ncells))
        TestConservation(setup_primordial_network, rv_int, density)
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
def test_different_solvers(setup_primordial_network, setup_solver_options):
    runtime_ncells(setup_primordial_network, setup_solver_options)
    return


@pytest.mark.parametrize('setup_solver_options',
                         ({"use_omp": False, "use_cvode": True,
                           "use_suitesparse": False,
                           "output_dir": "cvode_dls"},
                          {"use_omp": True, "use_cvode": True,
                           "use_suitesparse": True,
                           "output_dir": "cvode_klu"}),
                         indirect=True)
@pytest.mark.parametrize(('nd', 'nT', 'nf'),
                         ([16, 16, 3],))
def test_run_grid(setup_primordial_network, setup_solver_options, nd, nT, nf):
    os.chdir(pytest_dir)
    # test on grid of initial conditions
    density_array = np.logspace(0, 15, nd)
    temp_array = np.logspace(2, 3.7,  nT)
    h2frac_array = np.logspace(-6, -4, nf)

    setup_solver_options["niters"] = 1
    setup_solver_options["NCELLS"] = 1
    write_network(setup_primordial_network,
                  solver_options=setup_solver_options)
    rgrid, error_grid = run_grid(
        setup_primordial_network,
        setup_solver_options, density_array, temp_array, h2frac_array)

    # visualize the percentage error plot
    error_grid = np.array(error_grid).reshape(nd, nT, nf, 3)
    fig, axes = plt.subplots(nf, 2, figsize=(nf*4, 2*4))
    error_name = ["$e^-$ ", "H ", "He "]
    X, Y = np.meshgrid(density_array, temp_array)
    for i in range(2):
        for j in range(nf):
            im = axes[j, i].pcolor(X, Y,
                                   np.log10(np.abs(error_grid[:, :, j, i])))
            axes[j, i].set_title(error_name[i] +
                                 " {0:0.2e}".format(h2frac_array[j]))
            axes[j, i].set_xscale("log")
            axes[j, i].set_yscale("log")
            div = make_axes_locatable(axes[j, i])
            cax = div.append_axes(position="right", size="20%", pad=0.2)
            cbar = plt.colorbar(im, cax=cax)
    fig.tight_layout()
    fig.savefig("primordial_error_rate.png")
