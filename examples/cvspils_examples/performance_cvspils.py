

def time_taken_cvspils(Tdim = 50 ,Ddim = 50 ,total_t = 1e10, n_species=6):
    time_taken = []
    success = []

    HI_arr = []
    HII_arr = []
    HeI_arr = []
    HeII_arr = []
    HeIII_arr = []
    de_arr = []
    ge_arr = []
    T_arr = []
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
    rv, rv_int_sundials = sundials_cvspils_run.run_sundials_CVSpils(init, total_t, niter=1e5)
    data_dict = {}
    for sp in rv_int_sundials.keys():
        data_dict[sp] = []


    for temp, den in den_temp_pair2d:

        # initial conditions for a given temperature and density
        init, primordial = Init_values(np.array([temp]), np.array([den]), n_species = n_species)
        total_t = calc_fftime(den).v
        # Calculate the time it takes
        tic=timeit.default_timer()
        rv, rv_int_sundials = sundials_cvspils_run.run_sundials_CVSpils(init, total_t, niter=2e5)
        toc=timeit.default_timer()
        time_taken.append(toc - tic)

        temp_arr.append(temp)
        den_arr.append(den)

        flag = rv_int_sundials['successful']
        try:
            final_t = rv_int_sundials['t'][flag][-1]
        except:
            final_t = 0.0
        success.append(final_t)
        total_t_arr.append(total_t)
        for sp in rv_int_sundials.keys():
            if sp not in ["successful", 'dt', 't']:
                data_dict[sp].append( rv_int_sundials[sp][0][flag][-1] )
            else:
                data_dict[sp].append( rv_int_sundials[sp][flag][-1] )
    success = numpy.array(success)
    time_taken = numpy.array(time_taken)
    total_t_arr = numpy.array(total_t_arr)
    temp_arr = numpy.array(temp_arr)
    den_arr = numpy.array(den_arr)

    filename = 'cvspils_final_%d.hdf5' %(n_species)
    try:
        f = h5py.File(filename, 'w')
    except:
        f = h5py.File(filename)
        f.close()
        f = h5py.File(filename, 'w')
    for key in data_dict.keys():
        dset = f.create_dataset(key, data = data_dict[key])
    dset = f.create_dataset( "density", data = den_arr )
    dset = f.create_dataset( "temperature", data = temp_arr )
    dset = f.create_dataset( "success", data = success )
    dset = f.create_dataset( "time_taken", data = time_taken)
    dset = f.create_dataset( "total_t", data = total_t_arr)

    dset = f.create_dataset( "Tdim", data = (Tdim))
    dset = f.create_dataset( "Ddim", data = (Ddim))

    f.close()

    return success, time_taken, temp_arr, den_arr, filename
