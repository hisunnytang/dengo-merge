
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

from matplotlib import pyplot
import os
import yt
import numpy as np
from collections import defaultdict
from pygrackle import \
    FluidContainer, \
    chemistry_data, \
    evolve_constant_density
import pygrackle
from pygrackle.utilities.physical_constants import \
    mass_hydrogen_cgs, \
    sec_per_year, \
    cm_per_mpc

import pickle

# my monkey-patch part
def _new_evolve_constant_density(fc, final_temperature=None,
                            final_time=None, safety_factor=0.01):
    my_chemistry = fc.chemistry_data


    if final_temperature is None and final_time is None:
        raise RuntimeError("Must specify either final_temperature " +
                           "or final_time.")

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
        dt = safety_factor * np.abs(fc["cooling_time"][0])
        fc.calculate_temperature()

        print("Evolve constant density - t: %e s, rho: %e g/cm^3, T: %e K." %
               (current_time * my_chemistry.time_units ,
                fc["density"][0] * my_chemistry.density_units,
                fc["temperature"][0]))

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


def run_grackle(initial_temperature, density, final_time, safety_factor = 1e-2 ):
    current_redshift = 0.
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
    grackle_dir = os.path.dirname(os.path.dirname(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
    my_chemistry.grackle_data_file = os.sep.join(
        [grackle_dir, "input", "CloudyData_UVB=HM2012.h5"])

    # Set units
    my_chemistry.comoving_coordinates = 0 # proper units
    my_chemistry.a_units = 1.0
    my_chemistry.three_body_rate = 4
    my_chemistry.a_value = 1. / (1. + current_redshift) / \
        my_chemistry.a_units
    my_chemistry.density_units = 1.66053904e-24 # 1u: amu # rho = 1.0 is 1.67e-24 g
    my_chemistry.length_units = 1.0         #  cm
    my_chemistry.time_units = 1.0           #  s
    my_chemistry.velocity_units = my_chemistry.a_units * \
        (my_chemistry.length_units / my_chemistry.a_value) / \
        my_chemistry.time_units

    rval = my_chemistry.initialize()

    fc = FluidContainer(my_chemistry, 1)
    fc["density"][:] = density

    tiny_number = 1.0e-20
    if my_chemistry.primordial_chemistry > 0:
        # these are in mass_density
        fc["HI"][:] = 0.76 * fc["density"] /3.
        fc["HII"][:] =  0.76 * fc["density"] /3.
        fc["HeI"][:] =  (1 - 0.76) * fc["density"] /2.
        fc["HeII"][:] =  (1 - 0.76) * fc["density"] /2.
        fc["HeIII"][:] =  tiny_number

    if my_chemistry.primordial_chemistry > 1:

        fc["H2I"][:]  = 0.76 * fc["density"] /3.
        fc["H2II"][:] = tiny_number #0.76 * fc["density"] /4.
        fc["HM"][:]   = tiny_number #0.76 * fc["density"] /4.

    if my_chemistry.primordial_chemistry > 2:
        fc["DI"][:] = 2.0 * 3.4e-5 * fc["density"]
        fc["DII"][:] = tiny_number * fc["density"]
        fc["HDI"][:] = tiny_number * fc["density"]
    if my_chemistry.metal_cooling == 1:
        fc["metal"][:] = 0.1 * fc["density"] * \
          my_chemistry.SolarMetalFractionByMass

    fc["x-velocity"][:] = 0.0
    fc["y-velocity"][:] = 0.0
    fc["z-velocity"][:] = 0.0

    fc["energy"][:] = initial_temperature / \
        fc.chemistry_data.temperature_units
    fc.calculate_temperature()
    fc["energy"][:] *= initial_temperature / fc["temperature"]

    # let gas cool at constant density
    data = pygrackle.evolve_constant_density(
        fc, final_time=final_time,
        safety_factor=safety_factor)

    return data


data = run_grackle( 1000.0, 1e10, 1e10 , safety_factor = 1e-3)
print(data)
