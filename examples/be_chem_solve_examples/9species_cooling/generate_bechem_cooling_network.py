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
        init_values["He_1"]    = density * (1.0 - 0.76) /2.0
        init_values["He_2"]    = density * (1.0 - 0.76) /2.0
        init_values["He_3"]    = np.array([tiny])
        init_values["H_1"]     = density *  (0.76)  /3.
        init_values['H_2']     = density *  (0.76)  /3.

        init_values["H_m0"]    = np.array([tiny])
        init_values["H2_1"]    = density *  (0.76)  /3.
        init_values["H2_2"]    = np.array([tiny])

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
    gammaH2 = 2.0 / (5.0 + 2.0*x*x / (numpy.exp(x) - 1 )**2.0 ) + 1

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
                if k.name != "cie_cooling":
                    print("adding:", k.name, k.equation)
                    primordial.add_cooling(cooling_action)
                    print('---------------------------')
    return init_values, primordial



def create_bechem_solver(init, primordial, solver_name, cooling):
    # name of the solver
    pyximport.install(setup_args={"include_dirs":np.get_include()},
                      reload_support=True, inplace=True)

    # write the network
    primordial.write_solver(solver_name, output_dir = ".",
        solver_template = "rates_and_rate_tables",
        ode_solver_source = "BE_chem_solve.C",
        init_values=init,
        input_is_number=False)

    # import the pyx module
    bechem_run = pyximport.load_module("{}_run".format(solver_name),
                    "{}_solver_run.pyx".format(solver_name),
                     build_inplace = True, pyxbuild_dir = "_dengo_temp")
    return bechem_run


def solver_performance(Tdim = 50 ,Ddim = 50 , n_species=9, solver_name = "bechem_9species", cooling=True):
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

    chemistry_run = create_bechem_solver(init, primordial, solver_name, cooling=cooling)
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
        print(flag)
        try:
            final_t = rv_int['t'][flag][-1]
        except:
            final_t = 0.0

        success.append(final_t)
        total_t_arr.append(total_t)

        if numpy.sum(flag) > 0:
            t_interp_arr = numpy.logspace( -2, numpy.log10(total_t), 300 )
            t_arr = rv_int['t'][flag]
            for sp in rv_int.keys():
                if sp not in ["successful", 'dt', 't']:
                    data = rv_int[sp][0][flag]
                    interp_data = numpy.interp( t_interp_arr, t_arr, data )
                    data_dict[sp].append( interp_data )
        else:
            nans = numpy.empty((300))
            nans[:] = numpy.nan
            for sp in rv_int.keys():
                if sp not in ["successful", 'dt', 't']:
                    data_dict[sp].append( nans )


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


def run_bechem( temperature, density , total_t, init=None, primordial=None):

    solver_name = 'cvspils_9species'
    if init == None:
        init, primordial = Init_values(np.array([temperature]), np.array([density]), n_species = 9)
    chemistry_run = create_bechem_solver(init, primordial, solver_name, cooling=True)
    import timeit

    tic = timeit.default_timer()
    rv, rv_int = eval("chemistry_run.run_"+solver_name+"(init, total_t, niter=1e4)")
    toc = timeit.default_timer()

    run_time = toc - tic
    return rv,rv_int, run_time


if __name__ == '__main__':
    # initialize data dictionary
    solver_name = 'bechem_9species'
    init, primordial = Init_values(np.array([2000]), np.array([1e10]), n_species = 9)
    chemistry_run = create_bechem_solver(init, primordial, solver_name, cooling=True)
    success, time_taken, temp_arr, den_arr, filename = solver_performance(Tdim = 20 ,
                                                                      Ddim = 20 ,
                                                                      n_species=9,
                                                                      solver_name = "bechem_9species",
                                                                      cooling=True)

