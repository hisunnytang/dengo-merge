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

# this is required to compiled cython
os.environ["HDF5_DIR"] = "/home/kwoksun2/anaconda3"
# this is to fill in relative path in the templates
os.environ["CVODE_PATH"] = "/home/kwoksun2/cvode-3.1.0/instdir"
os.environ["HDF5_PATH"] = "/home/kwoksun2/anaconda3"
os.environ["SUITESPARSE_PATH"] = "/home/kwoksun2/SuiteSparse"
os.environ["DENGO_INSTALL_PATH"] = "/home/kwoksun2/dengo_install"

output_dir = "./test_primordial"

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
    # init_values['ge'] = ((temperature * number_density * kboltz) /
    # (init_values['density'] * mh * (gamma - 1)))
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

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    if use_omp and use_cvode and use_suitesparse:
        network.write_solver(solver_name, output_dir=output_dir,
                             solver_template="cv_omp/sundials_CVDls",
                             ode_solver_source="initialize_cvode_solver.C")
        return

    if not (use_omp and use_cvode and use_suitesparse):
        network.write_solver(
            solver_name, output_dir=output_dir,
            solver_template="be_chem_solve/rates_and_rate_tables",
            ode_solver_source="BE_chem_solve.C")
        return


def freefall_time(density):
    return 1.0 / np.sqrt(density[0] * G * mh)


def run_solver(init_values, solver_options, make_plot=True):

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
        pylab.clf()
        skip = ('successful', 'dt', 't', 'ge')
        for n, v in sorted(rv_int.items()):
            if n in skip:
                continue
            pylab.loglog(rv_int['t'], v, label=n)
            print(n, "{0:0.5g}".format(v[-1]))
        # pylab.ylim(density * 1e-20, density * 100)
        pylab.xlabel("time [s]")
        pylab.legend(loc='best', fontsize='xx-small')
        pylab.savefig("plot_primordial_network.png")
    else:
        try:
            for n, v in sorted(rv_int.items()):
                rv_int[n] = np.array([v[-1]])
            return rv_int
        except BaseException:
            return rv_int


def test_tolerance_convergence(setup_primordial_network, setup_solver_options):
    """Change number of iteration to check if the results converges
    """
    density = 1.0e10
    temperature = 2000.0
    h2frac = 1.0e-3
    NCELLS = 1

    rtol_array = np.logspace(-9, -4, 6)
    start = True
    skip = ('successful', 'dt', 't')
    for r in rtol_array:
        setup_solver_options["reltol"] = r
        init_values = setup_initial_conditions(
            setup_primordial_network, density, temperature, h2frac, NCELLS)
        r = run_solver(init_values, setup_solver_options, make_plot=False)
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
    ax.set_ylim(1.0e-9, 1e-2)
    ax.set_xlim(1.0e-8, 1.0e-4)
    ax.set_xlabel("Relative Tolerance of the Solver")
    ax.set_ylabel("Fractional Differnce to rtol =1e-9")
    ax.legend(loc="best", fontsize="xx-small")
    f.savefig("primordial_tolerance_convergence.png")


def test_iteration_convergence(init_values, solver_options):
    """Check how the results changes with reltol
    """
    pass


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
    return
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

    network = setup_primorial_network()
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
        # TestConservation(network, results, density)

# },


@pytest.mark.parametrize(
    'setup_solver_options',
    ({"use_omp": False, "use_cvode": False, "use_suitesparse": False},
     {"use_omp": False, "use_cvode": True, "use_suitesparse": False},
     {"use_omp": True, "use_cvode": True, "use_suitesparse": True}),
    indirect=True)
def test_different_solvers(setup_primordial_network, setup_solver_options):
    os.chdir(output_dir)
    runtime_ncells(setup_primordial_network, setup_solver_options)
    os.chdir(output_dir)
    return


@pytest.mark.parametrize(('nd', 'nT', 'nf'),
                         ([16, 16, 3],))
def test_run_grid(setup_primordial_network, setup_solver_options, nd, nT, nf):
    # test on grid of initial conditions
    density_array = np.logspace(0, 15, nd)
    temp_array = np.logspace(2, 3.7,  nT)
    h2frac_array = np.logspace(-6, -4, nf)

    setup_solver_options["use_omp"] = True
    setup_solver_options["use_cvode"] = True
    setup_solver_options["use_suitesparse"] = True
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


if __name__ == "__main__":
    test_primordial_suite()
