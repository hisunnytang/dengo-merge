import os
import pickle
import sys
import time

import matplotlib.pyplot as plt
import numpy
import numpy as np
import pyximport
import sympy
from unyt import amu_cgs

import dengo.primordial_cooling
import dengo.primordial_rates
from dengo.chemical_network import (
    ChemicalNetwork,
    cooling_registry,
    reaction_registry,
    species_registry,
)
from dengo.chemistry_constants import kboltz, mh, tiny

pytest_dir = os.getcwd()


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
    #primordial.add_reaction("k21")
    primordial.add_reaction("k22")
    #primordial.add_reaction("k23")

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
    primordial.add_cooling("ciHeI")
    primordial.add_cooling("ciHeII")
    primordial.add_cooling("cie_cooling")
    primordial.add_cooling("compton")
    # This defines the temperature range for the rate tables
    primordial.init_temperature((1e0, 1e8))

    return primordial


def setup_solver_options(update_options={}):
    solver_options = {
        "output_dir": "./{}".format(output_dir),
        "solver_name": "primordial",
        "use_omp": False,
        "use_cvode": False,
        "use_suitesparse": False,
        "niters": 1e3,
        "NCELLS": 128,
        "reltol": 1.0e-6,
    }
    solver_options.update(update_options)
    return solver_options


def setup_initial_conditions(network, density, temperature, h2frac, NCELLS):
    # setting initial conditions
    temperature = np.ones((NCELLS)) * temperature
    init_array = np.ones(NCELLS) * density
    init_values = dict()
    init_values["H_1"] = init_array * 0.76
    init_values["H_2"] = init_array * tiny
    init_values["H_m0"] = init_array * tiny
    init_values["He_1"] = init_array * 0.24
    init_values["He_2"] = init_array * tiny
    init_values["He_3"] = init_array * tiny
    init_values["H2_1"] = init_array * h2frac
    init_values["H2_2"] = init_array * tiny
    init_values["de"] = init_array * tiny

    # update and calculate electron density and etc with the handy functions
    # init_values = primordial.convert_to_mass_density(init_values)
    init_values["de"] = network.calculate_free_electrons(init_values)
    # init_values['density'] = primordial.calculate_total_density(init_values)
    init_values["density"] = np.ones((NCELLS)) * density
    number_density = network.calculate_number_density(init_values)

    # set up initial temperatures values used to define ge
    init_values["T"] = temperature

    # calculate ge (very crudely, no H2 help here)
    gamma = 5.0 / 3.0
    # init_values['ge'] = ((temperature * number_density * kboltz) /
    # (init_values['density'] * mh * (gamma - 1)))
    mH = 1.67e-24
    init_values["ge"] = 3.0 / 2.0 * temperature * kboltz / mH

    return init_values


def write_network(
    network,
    solver_options={
        "output_dir": "test_dir",
        "solver_name": "primordial",
        "use_omp": False,
        "use_cvode": False,
        "use_suitesparse": False,
    },
):
    """Write solver based on the ChemicalNetwork"""
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
        network.write_solver(
            solver_name,
            output_dir=output_dir,
            solver_template="cv_omp/sundials_CVDls",
            ode_solver_source="initialize_cvode_solver.C",
        )
        return

    if not (use_omp and use_cvode and use_suitesparse):
        network.write_solver(
            solver_name,
            output_dir=output_dir,
            solver_template="be_chem_solve/rates_and_rate_tables",
            ode_solver_source="BE_chem_solve.C",
        )
        return


def save_obj(obj, name):
    with open(name + ".pkl", "wb") as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)


def load_obj(name):
    with open(name + ".pkl", "rb") as f:
        return pickle.load(f)


def calculate_pressure(init, primordial):
    P = numpy.zeros((1))
    T = init["T"]
    for sp in primordial.required_species:
        if sp.name != "ge":
            n_sp = init[sp.name] / sp.weight
            P += n_sp * kboltz * T
    return P


def calculate_collapse_factor(pressure, density):
    # Calculate the effective adiabatic index, dlog(p)/dlog(rho).
    if len(pressure) < 3:
        return 0.0

    # compute dlog(p) / dlog(rho) using last two timesteps
    gamma_eff = np.log10(pressure[-1] / pressure[-2]) / np.log10(
        density[-1] / density[-2]
    )

    # compute a higher order derivative if more than two points available
    if len(pressure) > 2:
        gamma_eff += 0.5 * (
            (
                np.log10(pressure[-2] / pressure[-3])
                / np.log10(density[-2] / density[-3])
            )
            - gamma_eff
        )

    gamma_eff = min(gamma_eff, 4.0 / 3.0)

    # Equation 9 of Omukai et al. (2005)
    if gamma_eff < 0.83:
        force_factor = 0.0
    elif gamma_eff < 1.0:
        force_factor = (
            0.6 + 2.5 * (gamma_eff - 1) - 6.0 * np.power((gamma_eff - 1.0), 2.0)
        )
    else:
        force_factor = (
            1.0
            + 0.2 * (gamma_eff - (4.0 / 3.0))
            - 2.9 * np.power((gamma_eff - (4.0 / 3.0)), 2.0)
        )

    force_factor = max(force_factor, 0.0)
    force_factor = min(force_factor, 0.95)
    return force_factor


def evaluate_gamma_factor(init, primordial, temperature):
    gamma = 5.0 / 3.0
    num_den = {}
    for sp in primordial.required_species:
        num_den[sp.name] = init[sp.name].item()

    x = 6100.0 / temperature
    gammaH2 = (
        2.0 / (5.0 + 2.0 * x * x * numpy.exp(x) / (numpy.exp(x) - 1) ** 2.0) + 1
    ).item()

    gamma_fac = primordial.gamma_factor()
    gamma_factor = gamma_fac.subs(num_den).subs(
        {"gamma": gamma, "gammaH2_1": gammaH2, "gammaH2_2": gammaH2}
    )
    return gamma_factor


def calculate_gamma(init, primordial, temperature):
    #    gamma = 5.0/3.0
    #    temperature = init["T"]
    #    for sp in primordial.required_species:
    #        if sp.name == 'H2_1':
    #            sp_H2 = sp
    #            break
    #    gammaH2 = primordial.species_gamma(
    #        sp, temp=True, name=False).subs(
    #        {'T': temperature})
    #
    #    gammaH2 = 7./5.
    #    gamma_fac = primordial.gamma_factor()
    #    gamma_factor = gamma_fac.subs(init).subs(
    #        {'gamma': gamma}).subs(
    #        {'gammaH2_1': gammaH2,
    #         "gammaH2_2": gammaH2})
    #    print(gamma_factor)

    gamma_factor = evaluate_gamma_factor(init, primordial, temperature)
    n_density = 0.0
    for sp in primordial.required_species:
        if sp.name != "ge":
            n_density += init[sp.name]

    gamma_ad = n_density / gamma_factor + 1
    gamma_ad = float(gamma_ad)
    return gamma_ad


def calculate_temperature(init, primordial):
    dT = 10.0
    temperature = init["T"]
    nden = {}
    for sp in primordial.required_species:
        nden[sp.name] = init[sp.name] / sp.weight
    while dT > 0.1:
        x = 6100.0 / temperature
        # update the gammaH2 which is dependent on temperature
        # gammaH2 = 7./5.
        # gamma_factor = primordial.gamma_factor().subs(nden).subs(
        #    {'gammaH2_1': gammaH2, "gammaH2_2": gammaH2, 'gamma': 5./3., 'T': temperature})

        gamma_factor = evaluate_gamma_factor(nden, primordial, temperature)
        # with ge updated from compressional heating
        ge = init["ge"]

        new_T = numpy.array([float(init["density"] * ge * mh / kboltz / gamma_factor)])
        dT = numpy.abs(new_T - temperature)
        temperature = new_T

    return new_T


def calculate_energy(init, primordial):
    """Calculate energy from the abundance and temperature"""
    num_den = {}
    for sp in primordial.required_species:
        num_den[sp.name] = init[sp.name] / sp.weight

    # set up initial temperatures values used to define ge
    temperature = init["T"]
    # calculate gammaH2
    x = 6100.0 / temperature
    gammaH2 = 2.0 / (5.0 + 2.0 * x * x * numpy.exp(x) / (numpy.exp(x) - 1) ** 2.0) + 1

    # gamma_factor = primordial.gamma_factor().subs(num_den).subs(
    #    {'gammaH2_1': gammaH2, "gammaH2_2": gammaH2, 'gamma': 5./3., 'T': temperature})
    gamma_factor = evaluate_gamma_factor(num_den, primordial, temperature)

    ge = (temperature * kboltz) * gamma_factor / (init["density"] * mh)

    T = init["density"] * ge * mh / kboltz / gamma_factor

    return numpy.array([numpy.float64(ge)])


def update_initial_condition(
    init, primordial, pressure_array, density_array, safety_factor=0.01
):

    # should be in cgs units
    # dyne / cm^-2
    current_pressure = calculate_pressure(init, primordial)
    pressure_array = numpy.append(pressure_array, current_pressure)

    include_pressure = False
    if include_pressure:
        force_factor = calculate_collapse_factor(pressure_array, density_array)
    else:
        force_factor = 0.0

    density = init["density"]

    # compute the new density using the modified
    # free-fall collapse as per Omukai et al. (2005)

    gravitational_constant = 4.0 * numpy.pi * 6.65259e-8 * amu_cgs.v
    freefall_time_constant = np.power(
        ((32.0 * gravitational_constant) / (3.0 * numpy.pi)), 0.5
    )

    dt = safety_factor * np.power(
        (3.0 * np.pi) / (32.0 * gravitational_constant * density), 0.5
    )

    # calculate new density from altered free-fall solution

    new_density = np.power(
        (
            np.power(density, -0.5)
            - (0.5 * freefall_time_constant * dt * np.power((1 - force_factor), 0.5))
        ),
        -2.0,
    )

    # multiply this with the elemental abundances
    density_ratio = new_density / density

    # update densities
    # only update the species array only
    for sp in primordial.required_species:
        if sp.name != "ge":
            init[sp.name] *= density_ratio

    current_temperature = calculate_temperature(init, primordial)
    Gamma = calculate_gamma(init, primordial, current_temperature)

    # update internal energy
    init["ge"] += (
        (Gamma - 1.0) * init["ge"] * freefall_time_constant * new_density**0.5 * dt
    )

    # update density
    init["density"] = new_density
    density_array = numpy.append(density_array, new_density)

    # update temperature with the updated internal energy
    # init['T'] = calculate_temperature(init, primordial)

    return init, pressure_array, density_array, dt, force_factor


def generate_init_from_results(rv_int, primordial, old_init):
    flag = rv_int["successful"]
    init = {}
    for sp in primordial.required_species:
        init[sp.name] = rv_int[sp.name][0][flag][-1]  # *sp.weight
    density = old_init["density"]
    init["density"] = density
    init["T"] = numpy.array([rv_int["T"][0][flag][-1]])
    return init


def convert_from_grackle_to_dengo(grackle_dict, init=True):
    dengo_dict = {}
    for key in grackle_dict:
        key = str(key)

        ele = key.split("I")[0]
        charge = key.count("I")
        if charge > 0:
            dengo_name = ele + "_" + str(charge)
            dengo_dict[dengo_name] = numpy.array(grackle_dict[key]) / amu_cgs.v
        elif "M" in key:
            ele = key.split("M")[0]
            dengo_name = ele + "_" + str("m0")
            dengo_dict[dengo_name] = numpy.array(grackle_dict[key]) / amu_cgs.v
        elif key == "temperature":
            dengo_name = "T"
            dengo_dict[dengo_name] = numpy.array(grackle_dict[key])
        elif key == "de":
            dengo_name = "de"
            dengo_dict[dengo_name] = numpy.array(grackle_dict[key]) / amu_cgs.v
        elif key == "density":
            dengo_dict[key] = numpy.array(grackle_dict[key]) / amu_cgs.v
        else:
            dengo_dict[key] = numpy.array(grackle_dict[key])
    if init:
        for k, v in dengo_dict.items():
            dengo_dict[k] = v[0]
    return dengo_dict


def convert_from_grackle_to_dengo_all(grackle_dict):
    dengo_dict = {}
    for key in grackle_dict:
        key = str(key)

        ele = key.split("I")[0]
        charge = key.count("I")
        if charge > 0:
            dengo_name = ele + "_" + str(charge)
            if ele == "H":
                dengo_dict[dengo_name] = (
                    numpy.array(grackle_dict[key]) / amu_cgs.v
                )  # / 1.00794
            elif ele == "He":
                dengo_dict[dengo_name] = (
                    numpy.array(grackle_dict[key]) / amu_cgs.v
                )  # / 4.002602
            elif ele == "H2":
                dengo_dict[dengo_name] = (
                    numpy.array(grackle_dict[key]) / amu_cgs.v
                )  # / 1.00794 / 2.0
        elif "M" in key:
            ele = key.split("M")[0]
            dengo_name = ele + "_" + str("m0")
            dengo_dict[dengo_name] = numpy.array(grackle_dict[key]) / amu_cgs.v
        elif key == "temperature":
            dengo_name = "T"
            dengo_dict[dengo_name] = numpy.array(grackle_dict[key])
        elif key == "de":
            dengo_name = "de"
            dengo_dict[dengo_name] = numpy.array(grackle_dict[key]) / amu_cgs.v
    return dengo_dict


def run_dengo_freefall(
        update_options,
        final_density = 1.0e12
    ):
    solver_options = {
        "output_dir": "temp_freefall",
        "solver_name": "test_freefall",
        "use_omp": True,
        "use_cvode": True,
        "use_suitesparse": True,
        "niters": 1,
        "NCELLS": 1,
        "reltol": 1e-5,
    }
    solver_options.update(update_options)
    # Initial conditions
    temperature = 20000.0  # K
    density = 1.0e-1  # cm^-3
    h2frac = 1.0e-3

    solver_name = solver_options["solver_name"]
    network = setup_primordial_network()
    write_network(network, solver_options=solver_options)

    os.chdir(pytest_dir)
    os.chdir(solver_options["output_dir"])
    pyximport.install(
        setup_args={"include_dirs": np.get_include()}, reload_support=True, inplace=True
    )

    chemistry_run = pyximport.load_module(
        "{}_solver_run".format(solver_name),
        "{}_solver_run.pyx".format(solver_name),
        build_inplace=True,
        pyxbuild_dir="_dengo_temp",
    )
    ncells = 1
    init_values = setup_initial_conditions(
        network, density, temperature, h2frac, ncells
    )

    total_t = 0.0
    final_density *= 1.00794
    density_array = numpy.array([init_values["density"]])
    pressure_array = numpy.array([])
    ttt = []
    run_time = []
    current_density = density_array[-1]

    all_data = {}
    for key in init_values.keys():
        all_data[key] = []
    all_data["force_factor"] = []

    import h5py

    print("os.getcwd()")
    print(os.getcwd())
    dir_ff_grackle = os.path.join(
        os.path.dirname(__file__), "freefall_solution/freefall.h5"
    )
    f = h5py.File(dir_ff_grackle)
    fdata = f["data"]
    grackle_init = convert_from_grackle_to_dengo(fdata)

    new_init = setup_initial_conditions(network, density, temperature, h2frac, ncells)
    # new_init, primordial = Init_values(np.array([temperature]),
    #                                   np.array([density]) ,
    #                                   n_species = 9)
    # print(grackle_init)
    for i in new_init.keys():
        if i not in ["density", "ge"]:
            # print(i, grackle_init[i])
            new_init[i] = numpy.array([grackle_init[i]])

    f.close()

    new_init["de"] = network.calculate_free_electrons(new_init)
    new_init["ge"] = calculate_energy(new_init, network)
    rv, rv_int = eval(f"chemistry_run.run_{solver_name}")(new_init, 1e-4, niter=1e0)
    count = 0
    time0 = time.time()
    while current_density < final_density:

        # keep track of time in here
        new_init = generate_init_from_results(rv_int, network, new_init)
        (
            init,
            pressure_array,
            density_array,
            dt,
            force_factor,
        ) = update_initial_condition(
            new_init, network, pressure_array, density_array, safety_factor=1.0e-2
        )
        tic = time.time()
        rv, rv_int = eval(f"chemistry_run.run_{solver_name}")(
            init, dt, niter=1, intermediate=1, verbose=False
        )
        toc = time.time()
        total_t += dt
        ttt.append(float(total_t))
        run_time.append(toc - tic)
        flag = rv_int["successful"]
        for key in init.keys():
            if key not in ["density"]:
                data = rv_int[key][0][flag][-1]
                all_data[key].append(data)

        all_data["force_factor"].append(float(force_factor))
        current_density = density_array[-1]
        if count % 500 == 0:
            current_time = time.time()
            print("Time Lapsed: {}".format(current_time - time0))
            print(
                "density = {0:.2E}, percentage: {1:0.2g}".format(
                    current_density, current_density / final_density
                )
            )
        count += 1
        # if count > 5:
        #   break
    all_data["density"] = density_array
    all_data["run_time"] = run_time
    for k, v in all_data.items():
        all_data[k] = np.array(v)
    save_obj(all_data, "freefall_dengo")
    os.chdir("../")


def compare_dengo_grackle(solver_dir):

    # load in the results from dengo
    import h5py

    dir_ff_grackle = os.path.join(
        os.path.dirname(__file__), "freefall_solution/freefall.h5"
    )
    f = h5py.File(dir_ff_grackle)

    fdata = f["data"]
    grackle_results = convert_from_grackle_to_dengo(fdata, init=False)
    f.close()

    dengo_results = load_obj("{}/freefall_dengo".format(solver_dir))

    for k, v in dengo_results.items():
        dengo_results[k] = np.array(v)

    f, ax = plt.subplots()
    ax.loglog(
        grackle_results["density"] * 1.67e-24,
        grackle_results["T"],
        label="Grackle",
        color="C0",
    )
    ax.loglog(
        dengo_results["density"][1:] * 1.67e-24,
        dengo_results["T"],
        label="Dengo",
        color="C0",
        ls="--",
    )
    ax.set_ylabel(r"Temperature ($\mathrm{K}$)")

    ax.legend()
    ax2 = ax.twinx()
    ax2.loglog(
        grackle_results["density"] * 1.67e-24,
        grackle_results["H2_1"] / (grackle_results["H_1"] + grackle_results["H2_1"]),
        color="C1",
        label="Grackle",
    )
    ax2.loglog(
        dengo_results["density"][1:] * 1.67e-24,
        dengo_results["H2_1"] / (dengo_results["H_1"] + dengo_results["H2_1"]),
        color="C1",
        ls="--",
        label="Dengo",
    )
    ax2.set_ylabel("Molecular Hydrogen Fraction")
    ax.set_xlabel(r"density ($\mathrm{g / cm^{-3}}$)")
    plt.tight_layout()
    f.savefig("freefall_grackle_dengo_comparison.png", dpi=300)

    h2frac_dengo = dengo_results["H2_1"] / (
        dengo_results["H_1"] + dengo_results["H2_1"]
    )
    h2frac_grackle = grackle_results["H2_1"] / (
        grackle_results["H_1"] + grackle_results["H2_1"]
    )

    f, ax = plt.subplots()
    tdiff = (dengo_results["T"][:-1] - grackle_results["T"]) / dengo_results["T"][:-1]
    fh2diff = (h2frac_dengo[:-1] - h2frac_grackle) / h2frac_dengo[:-1]
    ax.semilogx(
        grackle_results["density"] * 1.67e-24, tdiff, label="temperature", color="C0"
    )
    ax.semilogx(
        grackle_results["density"] * 1.67e-24,
        fh2diff,
        label=r"$\rm{f_{H_2}}$",
        color="C1",
    )
    ax.set_ylabel("Fractional Difference")
    ax.set_xlabel(r"density ($\mathrm{g / cm^{-3}}$)")
    ax.legend()
    plt.tight_layout()
    f.savefig("freefall_grackle_dengo_difference.png", dpi=300)

    ones = np.ones_like(h2frac_grackle)
    # within 10 % error
    np.testing.assert_array_almost_equal(
        ones, h2frac_dengo[1:] / h2frac_grackle, decimal=1
    )
    np.testing.assert_array_almost_equal(
        ones, dengo_results["T"][1:] / grackle_results["T"], decimal=1
    )

import argparse
if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--output_dir", default="temp_freefall", type=str)
    parser.add_argument("-n", "--name", default="test_freefall", type=str)
    parser.add_argument("-s", "--solver", default="bechem", type=str, choices=['bechem','cvode_dls', 'cvode_klu'])
    parser.add_argument('-f', '--final_density', default=1e12, type=float)

    args = parser.parse_args()
    use_cvode=False
    use_suitesparse=False
    if args.solver.startswith('cvode'):
        use_cvode == True
    if args.solver.endswith('_klu'):
        use_suitesparse=True

    solver_options = {
        "output_dir": args.output_dir,
        "solver_name": args.name,
        "use_omp": False,
        "use_cvode": use_cvode,
        "use_suitesparse": use_suitesparse,
        "niters": 1,
        "NCELLS": 1,
        "reltol": 1e-5,
    }
    run_dengo_freefall(solver_options, final_density= args.final_density)
    compare_dengo_grackle(args.output_dir)
