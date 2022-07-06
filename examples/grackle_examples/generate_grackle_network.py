import os
import pickle
import timeit
from collections import defaultdict

import h5py
import numpy
import numpy as np
import pygrackle
import yt
import yt.units as u
from matplotlib import pyplot
from pygrackle import FluidContainer, chemistry_data, evolve_constant_density
from pygrackle.utilities.physical_constants import (
    cm_per_mpc,
    mass_hydrogen_cgs,
    sec_per_year,
)


def save_obj(obj, name):
    loc_pwd = ""
    with open(loc_pwd + name + ".pkl", "wb") as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)


# my monkey-patch part
def _new_evolve_constant_density(
    fc, final_temperature=None, final_time=None, safety_factor=0.01
):
    my_chemistry = fc.chemistry_data

    if final_temperature is None and final_time is None:
        raise RuntimeError("Must specify either final_temperature " + "or final_time.")

    data = defaultdict(list)
    current_time = 0.0

    fc.calculate_cooling_time()
    dt = safety_factor * np.abs(fc["cooling_time"][0])
    fc.calculate_temperature()

    while True:
        if final_temperature is not None and fc["temperature"][0] <= final_temperature:
            break
        if final_time is not None and current_time >= final_time:
            break

        fc.calculate_temperature()
        fc.calculate_cooling_time()
        dt = safety_factor * np.abs(fc["cooling_time"][0])
        if current_time + dt > final_time:
            dt = final_time - current_time

        #         print("Evolve constant density - t: %e s, rho: %e g/cm^3, T: %e K." %
        #                (current_time * my_chemistry.time_units ,
        #                 fc["density"][0] * my_chemistry.density_units,
        #                 fc["temperature"][0]))

        fc.solve_chemistry(dt)

        for field in fc.density_fields:
            data[field].append(fc[field][0] * my_chemistry.density_units)
        data["energy"].append(fc["energy"][0])
        fc.calculate_temperature()
        data["temperature"].append(fc["temperature"][0])
        fc.calculate_pressure()
        data["pressure"].append(fc["pressure"][0])
        data["time"].append(current_time * my_chemistry.time_units)
        current_time += dt

    for field in data:
        if field in fc.density_fields:
            data[field] = yt.YTArray(data[field], "g/cm**3")
        elif field == "energy":
            data[field] = yt.YTArray(data[field], "erg/g")
        elif field == "time":
            data[field] = yt.YTArray(data[field], "s")
        elif field == "temperature":
            data[field] = yt.YTArray(data[field], "K")
        elif field == "pressure":
            data[field] = yt.YTArray(data[field], "dyne/cm**2")
        else:
            data[field] = np.array(data[field])
    return data


def run_grackle(
    density=1e10, initial_temperature=1e3, final_time=1e5, safety_factor=1e-2
):
    current_redshift = 0.0
    tiny_number = 1.0e-20

    pygrackle.evolve_constant_density = _new_evolve_constant_density

    # Set solver parameters
    my_chemistry = chemistry_data()
    my_chemistry.three_body_rate = 4
    my_chemistry.use_grackle = 1
    my_chemistry.with_radiative_cooling = 1
    my_chemistry.primordial_chemistry = 2
    my_chemistry.metal_cooling = 0
    my_chemistry.UVbackground = 0
    my_chemistry.self_shielding_method = 0
    my_chemistry.H2_self_shielding = 0
    my_chemistry.cie_cooling = 0
    grackle_dir = "/home/kwoksun2/grackle"
    my_chemistry.grackle_data_file = os.sep.join(
        [grackle_dir, "input", "CloudyData_UVB=HM2012.h5"]
    )

    # Set units
    my_chemistry.comoving_coordinates = 0  # proper units
    my_chemistry.a_units = 1.0
    my_chemistry.three_body_rate = 4
    my_chemistry.a_value = 1.0 / (1.0 + current_redshift) / my_chemistry.a_units
    my_chemistry.density_units = 1.66053904e-24  # 1u: amu # rho = 1.0 is 1.67e-24 g
    my_chemistry.length_units = 1.0  #  cm
    my_chemistry.time_units = 1.0  #  s
    my_chemistry.velocity_units = (
        my_chemistry.a_units
        * (my_chemistry.length_units / my_chemistry.a_value)
        / my_chemistry.time_units
    )

    rval = my_chemistry.initialize()

    fc = FluidContainer(my_chemistry, 1)
    fc["density"][:] = density

    tiny_number = 1.0e-20
    if my_chemistry.primordial_chemistry > 0:
        # these are in mass_density
        fc["HI"][:] = 0.76 * fc["density"] / 3.0
        fc["HII"][:] = 0.76 * fc["density"] / 3.0
        fc["HeI"][:] = (1 - 0.76) * fc["density"] / 2.0
        fc["HeII"][:] = (1 - 0.76) * fc["density"] / 2.0
        fc["HeIII"][:] = tiny_number

    if my_chemistry.primordial_chemistry > 1:

        fc["H2I"][:] = 0.76 * fc["density"] / 3.0
        fc["H2II"][:] = tiny_number
        fc["HM"][:] = tiny_number

        fc["de"][:] = fc["HII"][:] + fc["HeII"][:] / 4.0 + fc["HeIII"][:] / 4.0 * 2.0

    if my_chemistry.primordial_chemistry > 2:
        fc["DI"][:] = 2.0 * 3.4e-5 * fc["density"]
        fc["DII"][:] = tiny_number * fc["density"]
        fc["HDI"][:] = tiny_number * fc["density"]
    if my_chemistry.metal_cooling == 1:
        fc["metal"][:] = 0.1 * fc["density"] * my_chemistry.SolarMetalFractionByMass

    fc["x-velocity"][:] = 0.0
    fc["y-velocity"][:] = 0.0
    fc["z-velocity"][:] = 0.0

    fc["energy"][:] = initial_temperature / fc.chemistry_data.temperature_units
    fc.calculate_temperature()
    fc["energy"][:] *= initial_temperature / fc["temperature"]

    tic = timeit.default_timer()
    # let gas cool at constant density
    data = pygrackle.evolve_constant_density(
        fc, final_time=final_time, safety_factor=safety_factor
    )
    toc = timeit.default_timer()
    run_time = toc - tic
    print("time lapsed", run_time)

    # convert grackle language to dengo language
    name_map_dict = {
        "HI": "H_1",
        "HII": "H_2",
        "HeI": "He_1",
        "HeII": "He_2",
        "HeIII": "He_3",
        "H2I": "H2_1",
        "H2II": "H2_2",
        "HM": "H_m0",
        "de": "de",
        "temperature": "T",
        "time": "t",
        "energy": "ge",
    }

    weight_map_dict = {
        "HI": 1.00794,
        "HII": 1.00794,
        "HeI": 4.002602,
        "HeII": 4.002602,
        "HeIII": 4.002602,
        "H2II": 2.01588,
        "H2I": 2.01588,
        "HM": 1.00794,
        "de": 1.00794,
        "time": 1.0,
        "temperature": 1.0,
        "energy": 1.0,
    }

    data_grackle = {}
    for key in name_map_dict.keys():
        key_dengo = name_map_dict[key]
        if key_dengo in ["t", "T", "ge"]:
            data_grackle[key_dengo] = data[key].value
        else:
            data_grackle[key_dengo] = (
                data[key].value / weight_map_dict[key] / u.amu_cgs.value
            )
    save_obj(data_grackle, "grackle_data")
    return data_grackle, run_time


def solver_performance(Tdim=20, Ddim=20, solver_name="grackle", safety_factor=1e-3):
    time_taken = []

    total_t_arr = []

    temp_list = numpy.logspace(2, 4, Tdim)
    den_list = numpy.logspace(0, 18, Ddim)

    temp_2d, den_2d = numpy.meshgrid(temp_list, den_list)
    den_temp_pair2d = numpy.dstack((temp_2d, den_2d))

    den_temp_pair2d = numpy.reshape(den_temp_pair2d, (Tdim * Ddim, 2))

    temp_arr = []
    den_arr = []

    # initialize data dictionary
    data_grackle, run_time = run_grackle(
        density=1e10, initial_temperature=1e3, final_time=1e5, safety_factor=1e-2
    )
    data_dict = {}

    for sp in data_grackle.keys():
        data_dict[sp] = []

    for temp, den in den_temp_pair2d:
        print("temp: {}, density: {}".format(temp, den))
        # initial conditions for a given temperature and density
        total_t = calc_fftime(den).v
        # Calculate the time it takes
        # run grackle

        data_grackle, run_time = run_grackle(
            density=den,
            initial_temperature=temp,
            final_time=total_t,
            safety_factor=safety_factor,
        )
        time_taken.append(run_time)

        temp_arr.append(temp)
        den_arr.append(den)
        total_t_arr.append(total_t)

        t_interp_arr = numpy.logspace(-2, numpy.log10(total_t), 300)
        t_arr = data_grackle["t"]
        for sp in data_grackle.keys():
            data = data_grackle[sp]
            interp_data = numpy.interp(t_interp_arr, t_arr, data)
            data_dict[sp].append(interp_data)

    time_taken = numpy.array(time_taken)
    total_t_arr = numpy.array(total_t_arr)
    temp_arr = numpy.array(temp_arr)
    den_arr = numpy.array(den_arr)

    filename = "data_{}.hdf5".format(solver_name)
    try:
        f = h5py.File(filename, "w")
    except:
        f = h5py.File(filename)
        f.close()
        f = h5py.File(filename, "w")
    for key in data_dict.keys():
        data = numpy.array(data_dict[key])
        dset = f.create_dataset(key, data=data)
    dset = f.create_dataset("density", data=den_arr)
    dset = f.create_dataset("temperature", data=temp_arr)
    dset = f.create_dataset("time_taken", data=time_taken)
    dset = f.create_dataset("total_t", data=total_t_arr)

    dset = f.create_dataset("Tdim", data=(Tdim))
    dset = f.create_dataset("Ddim", data=(Ddim))

    f.close()

    return time_taken, temp_arr, den_arr, filename


def calc_fftime(den):
    mu = 1.0
    rho = mu * u.mass_hydrogen * den * (u.cm**-3)
    tff = numpy.sqrt(1.0 / u.G / rho).in_units("s")
    return tff


# initialize data dictionary
solver_name = "grackle"
time_taken, temp_arr, den_arr, filename = solver_performance(
    Tdim=20, Ddim=20, solver_name=solver_name
)
