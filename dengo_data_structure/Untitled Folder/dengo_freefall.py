#!/usr/bin/python
import os
os.environ['YT_DEST'] = '/home/kwoksun2/anaconda2/'
import numpy as np
from dengo.chemical_network import \
    ChemicalNetwork, \
    reaction_registry, \
    cooling_registry, species_registry

import sys
import dengo.primordial_rates, dengo.primordial_cooling
from dengo.chemistry_constants import tiny, kboltz, mh
import yt
import yt.units as u
import numpy
import pickle
import pyximport
import h5py
import timeit
import time
import matplotlib.pyplot as plt
from sympy import lambdify

def calculate_pressure(init, primordial):
    P = numpy.zeros((1))
    T = init['T']
    for sp in primordial.required_species:
        if sp.name != 'ge':
            n_sp = init[sp.name]/sp.weight
            P += n_sp * u.boltzmann_constant_cgs.v * T
    return P

def calculate_collapse_factor(pressure, density):
    # Calculate the effective adiabatic index, dlog(p)/dlog(rho).
    if len(pressure) < 3:
        return 0.

    # compute dlog(p) / dlog(rho) using last two timesteps
    gamma_eff = np.log10(pressure[-1] / pressure[-2]) / \
        np.log10(density[-1] / density[-2])

    # compute a higher order derivative if more than two points available
    if len(pressure) > 2:
        gamma_eff += 0.5 * ((np.log10(pressure[-2] / pressure[-3]) /
                             np.log10(density[-2] / density[-3])) - gamma_eff)

    gamma_eff = min(gamma_eff, 4./3.)

    # Equation 9 of Omukai et al. (2005)
    if gamma_eff < 0.83:
        force_factor = 0.0
    elif gamma_eff < 1.0:
        force_factor = 0.6 + 2.5 * (gamma_eff - 1) - \
            6.0 * np.power((gamma_eff - 1.0), 2.)
    else:
        force_factor = 1.0 + 0.2 * (gamma_eff - (4./3.)) - \
            2.9 * np.power((gamma_eff - (4./3.)), 2.)

    force_factor = max(force_factor, 0.0)
    force_factor = min(force_factor, 0.95)
    return force_factor

def calculate_gamma(init, primordial):
    gamma = 5.0/3.0
    temperature = init['T']

    num_den = {}
    for sp in primordial.required_species:
        try:
            num_den[sp.name] = init[sp.name]/ sp.weight
        except:
            pass

    for sp in primordial.required_species:
        if sp.name == 'H2_1':
            sp_H2 = sp
            break
    gammaH2 = primordial.species_gamma(sp_H2, temp=True, name=False).subs({'T':temperature})

    gamma_fac = primordial.gamma_factor()
    gamma_factor = gamma_fac.subs(num_den).subs({'gamma':gamma}).subs({'gammaH2':gammaH2})

    n_density = 0.0
    for sp in primordial.required_species:
        if sp.name != 'ge':
            n_density += init[sp.name] / sp.weight

    gamma_ad = n_density/gamma_factor + 1
    gamma_ad = float(gamma_ad)
    return gamma_ad



def calculate_temperature(init, primordial):
    dT = 10.0
    temperature = init['T']

    num_den = {}
    for sp in primordial.required_species:
        try:
            num_den[sp.name] = init[sp.name]/ sp.weight
        except:
            pass
    while dT > 0.1:

        for sp in primordial.required_species:
            if sp.name == 'H2_1':
                sp_H2 = sp
                break
        gammaH2 = primordial.species_gamma(sp_H2, temp=True, name=False).subs({'T':temperature})

        gamma_factor = primordial.gamma_factor().subs(num_den).subs({'gammaH2': gammaH2 , 'gamma': 5./3.,'T':temperature })

        # with ge updated from compressional heating
        ge = init['ge']

        new_T = numpy.array([numpy.float64(init['density']*ge*mh / kboltz / gamma_factor)])
        dT = numpy.abs(new_T - temperature)
        temperature = new_T

    return new_T

def calculate_energy(init, primordial):
    """Calculate energy from the abundance and temperature
    """
    num_den = {}
    for sp in primordial.required_species:
        try:
            num_den[sp.name] = init[sp.name]/ sp.weight
        except:
            pass

    # set up initial temperatures values used to define ge
    temperature = init['T']

    # calculate gammaH2
    x = 6100.0/temperature
    gammaH2 = 2.0 / (5.0 + 2.0*x*x*numpy.exp(x) / (numpy.exp(x) - 1 )**2.0 ) + 1

    gamma_factor = primordial.gamma_factor().subs(num_den).subs({'gammaH2': gammaH2 , 'gamma': 5./3.,'T': temperature })

    ge  = ((temperature *  kboltz) *gamma_factor
                         / (init['density'] * mh  ))

    T = init['density']*ge*mh / kboltz / gamma_factor

    print(T-temperature)

    return numpy.array( [numpy.float64(ge)] )



def update_initial_condition(init, primordial, pressure_array, density_array, safety_factor=0.01):

    # should be in cgs units
    # dyne / cm^-2
    current_pressure = calculate_pressure(init, primordial)
    pressure_array = numpy.append(pressure_array, current_pressure)


    include_pressure = False
    if include_pressure:
        force_factor = calculate_collapse_factor(pressure_array, density_array)
    else:
        force_factor = 0.0
    print("force_factor: {}".format(force_factor))

    density = init['density']

    # compute the new density using the modified
    # free-fall collapse as per Omukai et al. (2005)

    gravitational_constant = 4.0*numpy.pi*6.65259e-8 *  u.amu_cgs.v
    freefall_time_constant = np.power((( 32.0*gravitational_constant)/ (3.0*numpy.pi)), 0.5)

    dt = safety_factor* np.power( (3.0*np.pi)/ (32.0* gravitational_constant *density ), 0.5 )

    # calculate new density from altered free-fall solution

    new_density = np.power((np.power(density, -0.5) -
                                (0.5 * freefall_time_constant * dt *
                                 np.power((1 - force_factor), 0.5))), -2.)

    # multiply this with the elemental abundances
    density_ratio = new_density/density

    # update densities
    # only update the species array only
    for sp in primordial.required_species:
        if sp.name != 'ge':
            init[sp.name] *= density_ratio

    Gamma = calculate_gamma(init, primordial)

    # update internal energy
    init['ge'] += (Gamma - 1.0) * init['ge'] * \
                        freefall_time_constant* \
                        new_density**0.5 * dt

    print( "gammma - 1: {}".format((Gamma - 1.0)))
    # update density
    init['density'] = new_density
    density_array = numpy.append(density_array, new_density)


    # update temperature with the updated internal energy
    init['T'] = calculate_temperature(init, primordial)

    init['de'] = primordial.calculate_free_electrons(init)

    return init, pressure_array, density_array, dt, force_factor


def generate_init_from_results(rv_int, primordial, old_init):
    flag = rv_int['successful']
    init = {}
    for sp in primordial.required_species:
        init[sp.name] = rv_int[sp.name][0][flag][-1]
    density = old_init['density']
    init['density'] = density
    init['T'] = numpy.array([rv_int['T'][0][flag][-1]])
    return init


def convert_from_grackle_to_dengo(grackle_dict):
    dengo_dict = {}
    for key in grackle_dict:
        key = str(key)

        ele = key.split('I')[0]
        charge = key.count('I')
        if charge > 0:
            dengo_name = ele+ '_' + str(charge)
            dengo_dict[dengo_name] = numpy.array(grackle_dict[key][0])/u.amu_cgs.v
        elif 'M' in key:
            ele = key.split('M')[0]
            dengo_name = ele + '_' + str("m0")
            dengo_dict[dengo_name] = numpy.array(grackle_dict[key][0])/u.amu_cgs.v
        elif key == 'temperature':
            dengo_name = 'T'
            dengo_dict[dengo_name] = numpy.array(grackle_dict[key][0])
        elif key == 'de':
            dengo_name = 'de'
            dengo_dict[dengo_name] = numpy.array(grackle_dict[key][0])/u.amu_cgs.v
    return dengo_dict


def convert_from_grackle_to_dengo_all(grackle_dict):
    dengo_dict = {}
    for key in grackle_dict:
        key = str(key)

        ele = key.split('I')[0]
        charge = key.count('I')
        if charge > 0:
            dengo_name = ele+ '_' + str(charge)
            if ele == 'H':
                dengo_dict[dengo_name] = numpy.array(grackle_dict[key])/u.amu_cgs.v
            elif ele == 'He':
                dengo_dict[dengo_name] = numpy.array(grackle_dict[key])/u.amu_cgs.v
            elif ele == 'H2':
                dengo_dict[dengo_name] = numpy.array(grackle_dict[key])/u.amu_cgs.v
        elif 'M' in key:
            ele = key.split('M')[0]
            dengo_name = ele + '_' + str("m0")
            dengo_dict[dengo_name] = numpy.array(grackle_dict[key])/u.amu_cgs.v
        elif key == 'temperature':
            dengo_name = 'T'
            dengo_dict[dengo_name] = numpy.array(grackle_dict[key])
        elif key == 'de':
            dengo_name = 'de'
            dengo_dict[dengo_name] = numpy.array(grackle_dict[key])/u.amu_cgs.v
    return dengo_dict
