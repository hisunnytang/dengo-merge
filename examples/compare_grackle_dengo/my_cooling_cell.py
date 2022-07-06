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

import numpy as np
import pygrackle
import yt
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
    dt = safety_factor * 0.0005 * final_time

    fc.calculate_temperature()

    fc.calculate_cooling_time()
    dt = safety_factor * np.abs(fc["cooling_time"][0])

    while True:
        if final_temperature is not None and fc["temperature"][0] <= final_temperature:
            break
        if final_time is not None and current_time >= final_time:
            break

        fc.calculate_temperature()

        print(
            "Evolve constant density - t: %e s, rho: %e g/cm^3, T: %e K."
            % (
                current_time * my_chemistry.time_units,
                fc["density"][0] * my_chemistry.density_units,
                fc["temperature"][0],
            )
        )

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


if __name__ == "__main__":
    current_redshift = 0.0
    tiny_number = 1.0e-20

    pygrackle.evolve_constant_density = _new_evolve_constant_density
    # Set initial values
    density = 1.0e8  # 1 amu /cm^3 // mass of hydrogen atom
    initial_temperature = 1.0e3  # K
    final_time = 1.0e5  # s

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

    # timestepping safety factor
    safety_factor = 0.0001

    tic = timeit.default_timer()
    # let gas cool at constant density
    data = pygrackle.evolve_constant_density(
        fc, final_time=final_time, safety_factor=safety_factor
    )
    toc = timeit.default_timer()
    print("time lapsed", toc - tic)

    save_obj(data, "grackle_data")
    yt.save_as_dataset({}, "cooling_cell.h5", data)
