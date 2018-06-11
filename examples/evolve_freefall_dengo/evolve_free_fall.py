import os
os.environ['YT_DEST'] = '/home/kwoksun2/anaconda2/'
import numpy as np
from dengo.chemical_network import \
    ChemicalNetwork, \
    reaction_registry, \
    cooling_registry, species_registry
import dengo.primordial_rates, dengo.primordial_cooling
from dengo.chemistry_constants import tiny, kboltz, mh
import yt
import yt.units as u
import numpy
import pickle
import time

import sys
sys.path.append("cvdls_examples/9species_cooling")
sys.path.append("cvspils_examples/9species_cooling")
sys.path.append("be_chem_solve_examples/9species_cooling/")

from generate_dls_cooling_network import create_cvdls_solver, Init_values
import matplotlib.pyplot as plt

def save_obj(obj, name ):
    with open(name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_obj(name ):
    with open(name + '.pkl', 'rb') as f:
        return pickle.load(f)

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
    for sp in primordial.required_species:
        if sp.name == 'H2_1':
            sp_H2 = sp
            break
    gammaH2 = primordial.species_gamma(sp, temp=True, name=False).subs({'T':temperature})

    gamma_fac = primordial.gamma_factor()
    gamma_factor = gamma_fac.subs(init).subs({'gamma':gamma}).subs({'gammaH2':gammaH2})

    n_density = 0.0
    for sp in primordial.required_species:
        if sp.name != 'ge':
            n_density += init[sp.name]

    gamma_ad = n_density/gamma_factor + 1
    gamma_ad = float(gamma_ad)
    return gamma_ad



def calculate_temperature(init, primordial):
    dT = 10.0
    temperature = init['T']

    while dT > 0.1:
        x = 6100.0/temperature
        # update the gammaH2 which is dependent on temperature
        gammaH2 = 2.0 / (5.0 + 2.0*x*x*numpy.exp(x) / (numpy.exp(x) - 1 )**2.0 ) + 1

        gamma_factor = primordial.gamma_factor().subs(init).subs({'gammaH2': gammaH2 , 'gamma': 5./3.,'T':temperature })

        # with ge updated from compressional heating
        ge = init['ge']

        new_T = numpy.array([float(init['density']*ge*mh / kboltz / gamma_factor)])
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

    # update density
    init['density'] = new_density
    density_array = numpy.append(density_array, new_density)


    # update temperature with the updated internal energy
    init['T'] = calculate_temperature(init, primordial)

    return init, pressure_array, density_array, dt, force_factor


def generate_init_from_results(rv_int, primordial, old_init):
    flag = rv_int['successful']
    init = {}
    for sp in primordial.required_species:
        init[sp.name] = rv_int[sp.name][0][flag][-1]*sp.weight
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
                dengo_dict[dengo_name] = numpy.array(grackle_dict[key])/u.amu_cgs.v / 1.00794
            elif ele == 'He':
                dengo_dict[dengo_name] = numpy.array(grackle_dict[key])/u.amu_cgs.v / 4.002602
            elif ele == 'H2':
                dengo_dict[dengo_name] = numpy.array(grackle_dict[key])/u.amu_cgs.v / 1.00794 /2.0
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

if __name__ == "__main__":
    # Initial conditions
    temperature = 20000.0# K
    density = 1.0e-1 # cm^-3

    solver_name = 'cvdls_9species'
    init, primordial = Init_values(np.array([temperature]), np.array([density]) ,
                                   n_species = 9, cooling=True)

    chemistry_run = create_cvdls_solver(init,primordial, solver_name, cooling=True)


    total_t = 0.0
    final_density = 5.0e11*1.00794
    density_array = numpy.array([ init['density'] ])
    pressure_array = numpy.array([])
    ttt = []
    run_time = []
    current_density = density_array[-1]

    all_data = {}
    for key in init.keys():
        all_data[key] = []
    all_data['force_factor'] = []


    dir_ff_grackle = "/home/kwoksun2/grackle/src/python/examples/freefall.h5"
    f = h5py.File(dir_ff_grackle)
    fdata = f['data']
    grackle_init = convert_from_grackle_to_dengo(fdata)


    new_init, primordial = Init_values(np.array([temperature]),
                                       np.array([density]) ,
                                       n_species = 9)
    for i in new_init.keys():
        if i not in ['density','ge']:
            print(i, grackle_init[i])
            new_init[i] = numpy.array([grackle_init[i]])

    new_init['de'] = primordial.calculate_free_electrons(new_init)
    new_init['ge'] =  calculate_energy(new_init, primordial)
    rv, rv_int = chemistry_run.run_cvdls_9species(new_init, 1e-5,niter=1e0)


    while current_density < final_density:

    # keep track of time in here

        new_init = generate_init_from_results(rv_int,primordial, new_init)
        init, pressure_array, \
        density_array, dt, force_factor = update_initial_condition(new_init,
                                                primordial, pressure_array,
                                                density_array, safety_factor=0.01)
        tic = time.time()
        rv, rv_int = chemistry_run.run_cvdls_9species(init, dt,niter=1e4)
        toc = time.time()
        total_t += dt
        ttt.append(float(total_t))
        run_time.append(toc-tic)

        flag = rv_int['successful']
        for key in init.keys():
            if key not in ['density']:
                data = rv_int[key][0][flag][-1]
                all_data[key].append(data)
                print(key, data)

        all_data['force_factor'].append( float(force_factor))
        current_density = density_array[-1]
        print("percentage: {}".format(current_density/final_density) )
    all_data['density']  = density_array
    all_data['run_time'] = run_time
    save_obj(all_data, 'freefall_dengo')

