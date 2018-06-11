import os
os.environ['YT_DEST'] = '/home/kwoksun2/anaconda2/pkgs/yt-3.3.5-np111py27_2/'
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
import pyximport
import h5py
import timeit



def Init_values(temperature, density, n_species = 9, cooling=True):
    """ Create a initial value dictionary,
        for a given temperature, density, number of species
    Args:
        temperature -- in Kelvin
        density     -- in amu / cm**3
        n_species   -- number of species (6/9)
        cooling

    Returns:
        init_values: initial value dictionary with
                     self-consistent energy/ electron density
        primordial : chemical_network classes

    """

    # initialize and setup the network
    dengo.primordial_rates.setup_primordial()
    primordial = ChemicalNetwork()

    if n_species == 9:
        for i in range(23):
            try:
                primordial.add_reaction("k{0:02d}".format(i+1))
            except:
                pass
    else:
        for i in range(6):
            try:
                primordial.add_reaction("k{0:02d}".format(i+1))
            except:
                pass

    # the temperature array required to interpolates the rates
    primordial.init_temperature((1e0, 1e5))

    tiny = 1.0e-20

    # init_array are is in fractional abundances
    init_values = dict()

    if n_species == 6:
        # 6-species model
        init_values["He_1"]    = density * (1.0 - 0.76) /2.
        init_values["He_2"]    = density * (1.0 - 0.76) /2.
        init_values["He_3"]    = np.array([tiny])
        init_values["H_1"]     = density *  (0.76)  /2.
        init_values['H_2']     = density *  (0.76)  /2.
    else:
        # 9-species model
        init_values["He_1"]    = density * (1.0 - 0.76)
        init_values["He_2"]    = density * np.array([tiny])
        init_values["He_3"]    = density * np.array([tiny])
        init_values["H_1"]     = density *  (0.76)
        init_values['H_2']     = density * np.array([tiny])

        init_values["H_m0"]    = density * np.array([tiny])
        init_values["H2_1"]    = density * np.array([tiny])
        init_values["H2_2"]    = density * np.array([tiny])

    # now everything in mass density
    init_values['de'] = primordial.calculate_free_electrons(init_values)
    # one signle value: again mass density
    init_values['density'] = primordial.calculate_total_density(init_values)

    num_den = {}
    for sp in primordial.required_species:
        try:
            num_den[sp.name] = init_values[sp.name]/ sp.weight
        except:
            pass

    # set up initial temperatures values used to define ge
    init_values['T'] = temperature

    # calculate gammaH2
    x = 6100.0/temperature
    gammaH2 = 2.0 / (5.0 + 2.0*x*x*numpy.exp(x) / (numpy.exp(x) - 1 )**2.0 ) + 1

    gamma_factor = primordial.gamma_factor().subs(num_den).subs({'gammaH2': gammaH2 , 'gamma': 5./3.,'T': temperature })

    ge  = ((temperature *  kboltz) *gamma_factor
                         / (init_values['density'] * mh  ))

    T = init_values['density']*ge*mh / kboltz / gamma_factor
    print("difference in temperature:", T - temperature)
    init_values['ge'] = numpy.array( [numpy.float64(ge)] )

    if cooling:
        for cooling_action in cooling_registry:
            k = cooling_registry[cooling_action]
            if (k.species).issubset( primordial.required_species ):
                print("adding: {} {}".format(k.name, k.equation) )
                primordial.add_cooling(cooling_action)
                print('---------------------------')
    return init_values, primordial



def create_cvdls_solver(init, primordial, solver_name, cooling):
    print("three body rate: {}".format(primordial.threebody))
    # name of the solver
    pyximport.install(setup_args={"include_dirs":np.get_include()},
                      reload_support=True, inplace=True)

    # write the network
    primordial.write_solver(solver_name, output_dir = ".",
        solver_template = "cvdls/sundials_CVDls",
        ode_solver_source = "cvodes_solver_CVDls.C",
        init_values=init,
        input_is_number=False)

    # import the pyx module
    sundials_cvdls_run = pyximport.load_module("{}_run".format(solver_name),
                    "{}_solver_run.pyx".format(solver_name),
                     build_inplace = True, pyxbuild_dir = "_dengo_temp")
    return sundials_cvdls_run


def solver_performance(Tdim = 50 ,Ddim = 50 , n_species=9, solver_name = "cvdls_9species", cooling=False):
    time_taken = []
    success = []

    total_t_arr = []

    temp_list = numpy.logspace(2,4,Tdim)
    den_list  = numpy.logspace(0,18,Ddim)

    temp_2d, den_2d = numpy.meshgrid(temp_list,den_list)
    den_temp_pair2d = numpy.dstack( (temp_2d,den_2d)  )

    den_temp_pair2d = ( numpy.reshape(den_temp_pair2d, (Tdim*Ddim,2) )  )

    temp_arr = []
    den_arr = []

    # initialize data dictionary
    init, primordial = Init_values(np.array([1000]), np.array([1e10]), n_species = n_species)

    chemistry_run = create_cvdls_solver(init, primordial, solver_name, cooling=cooling)
    rv, rv_int = eval("chemistry_run.run_"+solver_name+"(init, 1e1, niter=1e4)")

    data_dict = {}
    for sp in rv_int.keys():
        data_dict[sp] = []




    for temp, den in den_temp_pair2d:

        # initial conditions for a given temperature and density
        init, primordial = Init_values(np.array([temp]), np.array([den]), n_species = n_species)
        total_t = calc_fftime(den).v
        # Calculate the time it takes
        tic=timeit.default_timer()
        rv, rv_int = eval("chemistry_run.run_"+solver_name+"(init, total_t, niter=1e4)")
        toc=timeit.default_timer()
        time_taken.append(toc - tic)
        temp_arr.append(temp)
        den_arr.append(den)

        flag = rv_int['successful']
        try:
            final_t = rv_int['t'][flag][-1]
        except:
            final_t = 0.0
        success.append(final_t)
        total_t_arr.append(total_t)

        t_interp_arr = numpy.logspace( -2, numpy.log10(total_t), 300 )
        t_arr = rv_int['t'][flag]
        for sp in rv_int.keys():
            if sp not in ["successful", 'dt', 't']:
                data = rv_int[sp][0][flag]
                interp_data = numpy.interp( t_interp_arr, t_arr, data )
                data_dict[sp].append( interp_data )

    success = numpy.array(success)
    time_taken = numpy.array(time_taken)
    total_t_arr = numpy.array(total_t_arr)
    temp_arr = numpy.array(temp_arr)
    den_arr = numpy.array(den_arr)

    filename = "data_{}.hdf5".format(solver_name)
    try:
        f = h5py.File(filename, 'w')
    except:
        f = h5py.File(filename)
        f.close()
        f = h5py.File(filename, 'w')
    for key in data_dict.keys():
        data = numpy.array(data_dict[key])
        dset = f.create_dataset(key, data = data )
    dset = f.create_dataset( "density", data = den_arr )
    dset = f.create_dataset( "temperature", data = temp_arr )
    dset = f.create_dataset( "success", data = success )
    dset = f.create_dataset( "time_taken", data = time_taken)
    dset = f.create_dataset( "total_t", data = total_t_arr)

    dset = f.create_dataset( "Tdim", data = (Tdim))
    dset = f.create_dataset( "Ddim", data = (Ddim))

    f.close()

    return success, time_taken, temp_arr, den_arr, filename

def calc_fftime(den):
    mu = 1.0
    rho = mu* u.mass_hydrogen *den *(u.cm**-3)
    tff = numpy.sqrt(1.0 / u.G / rho).in_units('s')
    return tff


def run_cvdls( temperature, density , total_t, init=None, primordial=None, max_iter=1e4):

    solver_name = 'cvdls_9species'
    if init == None:
        init, primordial = Init_values(np.array([temperature]), np.array([density]), n_species = 9)
    chemistry_run = create_cvdls_solver(init, primordial, solver_name, cooling=True)
    import timeit

    tic = timeit.default_timer()
    rv, rv_int = eval("chemistry_run.run_"+solver_name+"(init, total_t, niter={})".format(max_iter))
    toc = timeit.default_timer()

    run_time = toc - tic
    return rv,rv_int, run_time

if __name__ == "__main__":
    # initialize data dictionary
    solver_name = 'cvdls_9species'
    init, primordial = Init_values(np.array([2000]), np.array([1e10]), n_species = 9)
    print(init);
    chemistry_run = create_cvdls_solver(init, primordial, solver_name, cooling=False)

    rv, rv_int = chemistry_run.run_cvdls_9species(init, 1e-1)

    #success, time_taken, temp_arr, den_arr, filename = solver_performance(Tdim = 20 ,
    #                                                                  Ddim = 20 ,
    #                                                                  n_species=9,
    #                                                                  solver_name = "cvdls_9species",
    #                                                                  cooling=False)


def run_cvdls_evolve_freefall( temperature, density, init=None, primordial=None, max_iter=1e4 ):

    solver_name = 'cvdls_9species'
    if init == None:
        init, primordial = Init_values(np.array([temperature]), np.array([density]) , n_species = 9)
    # evov
    chemistry_run = create_cvdls_solver(init, primordial, solver_name, cooling=False)

    chemistry_run_cooling = create_cvdls_solver(init,primordial, solver_name, cooling=True)


    pressure = init['pressure']
    include_pressure = True

    # compute the new density using the modified
    # free-fall collapse as per Omukai et al. (2005)
    if include_pressure:
        force_factor = calculate_collapse_factor(pressure, density)
    else:
        force_factor = 0.0

    gravitational_constant = 6.6526e-8
    freefall_time_constant = np.power((( 32.0*gravitational_constant)/ (3.0*numpy.pi)), 0.5)

    safety_factor = 0.01
    dt = safety_factor* np.power( (3.0*np.pi)/ (32.0* gravitational_constant *density ), 0.5 )

    # calculate new density from altered free-fall solution
    new_density = np.power( (np.power(density, -0.5) -
                            (0.5*freefall_time_constant * dt*
                            np.power((1.0 - force_factor), 0.5) ) ), -2.0 )

    # multiply this with the elemental abundances
    density_ratio = new_density/density

    # update densities
    # only update the species array only
    for key in init.keys():
        init[key] *= density_ratio

def calculate_pressure(init):
    pressure = 0.0
    return pressure

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
