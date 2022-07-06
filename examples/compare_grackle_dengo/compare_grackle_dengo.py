########################################################################
#
# Cooling cell example script
#
#  This will initialize a single cell at a given temperature,
#  iterate the cooling solver for a fixed time, and output the
#  temperature vs. time.
#
#
# Copyright (c) 2015-2016, Grackle Development Team.
#
# Distributed under the terms of the Enzo Public Licence.
#
# The full license is in the file LICENSE, distributed with this
# software.
########################################################################
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

    fc.calculate_temperature()

    # use the dt trial-and-error with dengo BE_CHEM_SOLVE

    iii = 1

    while True:
        if final_temperature is not None and fc["temperature"][0] <= final_temperature:
            break
        if final_time is not None and current_time >= final_time:
            break

        fc.calculate_cooling_time()
        dt = safety_factor * np.abs(fc["cooling_time"][-1])
        fc.calculate_temperature()

        # print("Evolve constant density - t: %e s, rho: %e g/cm^3, T: %e K." %
        #       (current_time * my_chemistry.time_units ,
        #        fc["density"][0] * my_chemistry.density_units,
        #        fc["temperature"][0]))

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

        iii += 1

    for field in data:
        if field in fc.density_fields:
            data[field] = np.array(data[field])
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


def run_grackle(initial_temperature, density, final_time, safety_factor=1e-2):
    current_redshift = 0.0
    tiny_number = 1.0e-20

    pygrackle.evolve_constant_density = _new_evolve_constant_density

    # Set solver parameters
    my_chemistry = chemistry_data()
    my_chemistry.three_body_rate = 4
    my_chemistry.use_grackle = 1
    my_chemistry.with_radiative_cooling = 0
    my_chemistry.primordial_chemistry = 2
    my_chemistry.metal_cooling = 0
    my_chemistry.UVbackground = 0
    my_chemistry.self_shielding_method = 0
    my_chemistry.H2_self_shielding = 0
    grackle_dir = os.path.dirname(
        os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    )
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
        fc["H2II"][:] = tiny_number  # 0.76 * fc["density"] /4.
        fc["HM"][:] = tiny_number  # 0.76 * fc["density"] /4.

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

    # let gas cool at constant density
    data = pygrackle.evolve_constant_density(
        fc, final_time=final_time, safety_factor=safety_factor
    )

    return data


def solver_performance(temperature, density, grackle_safety_factor=1e-3):
    time_taken = []
    success = []

    total_t_arr = []

    temp_list = numpy.logspace(2, 4, Tdim)
    den_list = numpy.logspace(0, 18, Ddim)

    temp_2d, den_2d = numpy.meshgrid(temp_list, den_list)
    den_temp_pair2d = numpy.dstack((temp_2d, den_2d))

    den_temp_pair2d = numpy.reshape(den_temp_pair2d, (Tdim * Ddim, 2))

    temp_arr = []
    den_arr = []

    # initialize data dictionary
    data = run_grackle(1000, 1e5, 1e5, safety_factor=1e-3)

    data_dict = {}
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
        "energy": "ge",
        "temperature": "T",
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
    }

    for sp in name_map_dict.keys():
        data_dict[name_map_dict[sp]] = []

    for temp, den in den_temp_pair2d:

        # initial conditions for a given temperature and density
        total_t = calc_fftime(den).v
        # Calculate the time it takes
        tic = timeit.default_timer()
        # run grackle
        data_grackle = run_grackle(temp, den, total_t, safety_factor=1e-2)

        toc = timeit.default_timer()
        time_taken.append(toc - tic)
        temp_arr.append(temp)
        den_arr.append(den)
        total_t_arr.append(total_t)

        t_interp_arr = numpy.logspace(-2, numpy.log10(total_t), 300)
        t_arr = data_grackle["time"]
        for sp in data_grackle.keys():
            if sp in weight_map_dict.keys():
                data = data_grackle[sp] / u.amu_cgs.value / weight_map_dict[sp]
                interp_data = numpy.interp(t_interp_arr, t_arr, data)
                data_dict[name_map_dict[sp]].append(interp_data)

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


# success, time_taken, temp_arr, den_arr, filename = solver_performance(Tdim = 50 ,Ddim = 50 , n_species=9, solver_name = "cvdls_9species_cooling", cooling=True)


# initialize data dictionary
solver_name = "grackle"
time_taken, temp_arr, den_arr, filename = solver_performance(
    Tdim=20, Ddim=20, solver_name=solver_name
)
