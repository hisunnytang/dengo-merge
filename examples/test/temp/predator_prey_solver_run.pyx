#cython: language_level=3
cimport numpy as np
import numpy as np
import time
from libc.stdlib cimport malloc, free
from libcpp cimport bool
from cpython cimport array

# NSPECIES here is N in the .C.template file
DEF NSPECIES = 5
DEF MAX_NCELLS=1

cdef extern from "predator_prey_solver.h":
    cdef int _MAX_NCELLS  "MAX_NCELLS"
    cdef int _NSPECIES "NSPECIES"

    ctypedef struct dengo_field_data:
        unsigned long int nstrip;
        unsigned long int ncells;
        double *dead_predator_density;
        double *dead_prey_density;
        double *ge_density;
        double *predator_density;
        double *prey_density;
        double *density;
        double *CoolingTime;
        double *MolecularWeight;
        double *temperature;
        double *Gamma;
        double reltol;
        double floor_value;
        const char *dengo_data_file;
        int *grid_start;
        int *grid_end;
        int *grid_dimension;

    ctypedef struct code_units:
        int comoving_coordinates
        double density_units
        double length_units
        double time_units
        double velocity_units
        double a_units
        double a_value


    ctypedef struct predator_prey_data:
        double dbin
        double idbin
        double bounds[2]
        int nbins

        double d_zbin
        double id_zbin
        double z_bounds[2]
        int n_zbins

        double current_z
        double zdef
        double dz

        double Ts[MAX_NCELLS]
        double Tdef[MAX_NCELLS]
        double dT[MAX_NCELLS]
        double logTs[MAX_NCELLS]
        double dTs_ge[MAX_NCELLS]
        double r_exp_growth_prey[1024]
        double rs_exp_growth_prey[MAX_NCELLS]
        double drs_exp_growth_prey[MAX_NCELLS]
        double r_natural_death_predator[1024]
        double rs_natural_death_predator[MAX_NCELLS]
        double drs_natural_death_predator[MAX_NCELLS]
        double r_predation[1024]
        double rs_predation[MAX_NCELLS]
        double drs_predation[MAX_NCELLS]
        int bin_id[MAX_NCELLS]
        int ncells

    double dengo_evolve_predator_prey (double dtf, double dt, double z,
                                         double *input, double *rtol,
                                         double *atol, int dims,
                                         predator_prey_data *data, double *temp)
    int predator_prey_solve_chemistry_dt( code_units *units, dengo_field_data *field_data, double* rtol, double* atol, double dt)
    int predator_prey_solve_chemistry_enzo( code_units *units, dengo_field_data *field_data, double dt );
    int predator_prey_main(int argc, char **argv);



def main_run_predator_prey():
    t1 = time.time()
    predator_prey_main(0, NULL)
    t2 = time.time()
    print("Total elapsed time: %0.3e" % (t2-t1))

def run_predator_prey(ics, double tf, int niter = 10000,
                        int intermediate = 1, float z = -1.0, bool adaptive_step = True,
			float reltol = 1.0e-5, float floor_value = 1.0e-20):
    #assert(_MAX_NCELLS == MAX_NCELLS)
    assert(_NSPECIES == NSPECIES)
    cdef np.ndarray[np.float64_t, ndim=1] dead_predator_arr = ics["dead_predator"]
    if not dead_predator_arr.flags['C_CONTIGUOUS']:
        dead_predator_arr = np.ascontiguousarray( dead_predator_arr )
    cdef double[::1] dead_predator_memview = dead_predator_arr
    cdef np.ndarray[np.float64_t, ndim=2] dead_predator_int
    cdef np.ndarray[np.float64_t, ndim=1] dead_prey_arr = ics["dead_prey"]
    if not dead_prey_arr.flags['C_CONTIGUOUS']:
        dead_prey_arr = np.ascontiguousarray( dead_prey_arr )
    cdef double[::1] dead_prey_memview = dead_prey_arr
    cdef np.ndarray[np.float64_t, ndim=2] dead_prey_int
    cdef np.ndarray[np.float64_t, ndim=1] ge_arr = ics["ge"]
    if not ge_arr.flags['C_CONTIGUOUS']:
        ge_arr = np.ascontiguousarray( ge_arr )
    cdef double[::1] ge_memview = ge_arr
    cdef np.ndarray[np.float64_t, ndim=2] ge_int
    cdef np.ndarray[np.float64_t, ndim=1] predator_arr = ics["predator"]
    if not predator_arr.flags['C_CONTIGUOUS']:
        predator_arr = np.ascontiguousarray( predator_arr )
    cdef double[::1] predator_memview = predator_arr
    cdef np.ndarray[np.float64_t, ndim=2] predator_int
    cdef np.ndarray[np.float64_t, ndim=1] prey_arr = ics["prey"]
    if not prey_arr.flags['C_CONTIGUOUS']:
        prey_arr = np.ascontiguousarray( prey_arr )
    cdef double[::1] prey_memview = prey_arr
    cdef np.ndarray[np.float64_t, ndim=2] prey_int
    cdef np.ndarray[np.float64_t, ndim=1] CoolingTime;
    cdef np.ndarray[np.float64_t, ndim=1] Gamma;
    cdef np.ndarray[np.float64_t, ndim=1] Temperature;
    cdef np.ndarray[np.float64_t, ndim=1] MolecularWeight;

    cdef np.ndarray[np.float64_t, ndim=2] CoolingTime_int;
    cdef np.ndarray[np.float64_t, ndim=2] Gamma_int;
    cdef np.ndarray[np.float64_t, ndim=2] Temperature_int;
    cdef np.ndarray[np.float64_t, ndim=2] MolecularWeight_int;

    cdef np.ndarray[np.float64_t, ndim=1] result_int
    cdef np.ndarray[np.float64_t, ndim=2] temp_int
    cdef np.ndarray[np.float64_t, ndim=1] t_int
    cdef np.ndarray[np.float64_t, ndim=1] dt_int

    cdef int i, j, k, iter
    cdef int dims = ge_arr.shape[0]
    cdef int NTOT = NSPECIES * dims
    cdef double *input = <double *> malloc(NTOT * sizeof(double))
    cdef double *prev = <double *> malloc(NTOT * sizeof(double))

    cdef dengo_field_data *field_data = <dengo_field_data *> malloc(sizeof(dengo_field_data));
    cdef code_units       * units     = <code_units *> malloc(sizeof(code_units));
    units.density_units = 1.67e-24;
    units.length_units  = 1.0;
    units.time_units    = 1.0;
    units.velocity_units = 1.0;
    CoolingTime = np.ascontiguousarray(np.zeros((dims), dtype=np.double))
    Gamma       = np.ascontiguousarray(np.zeros((dims), dtype=np.double))
    Temperature = np.ascontiguousarray(np.zeros((dims), dtype=np.double))
    MolecularWeight = np.ascontiguousarray(np.zeros((dims), dtype=np.double))
    cdef double[::1] ctime_memview = CoolingTime
    cdef double[::1] gamma_memview = Gamma
    cdef double[::1] temp_memview  = Temperature
    cdef double[::1] mweight_memview = MolecularWeight
    field_data.dead_predator_density = &dead_predator_memview[0];
    field_data.dead_prey_density = &dead_prey_memview[0];
    field_data.ge_density = &ge_memview[0];
    field_data.predator_density = &predator_memview[0];
    field_data.prey_density = &prey_memview[0];
    field_data.CoolingTime = &ctime_memview[0]
    field_data.Gamma       = &gamma_memview[0]
    field_data.temperature = &temp_memview[0]
    field_data.MolecularWeight = &mweight_memview[0]
    field_data.ncells  = dims;

    cdef np.ndarray[np.float64_t, ndim=1] atol = np.ones( (NTOT) )*floor_value*reltol
    cdef np.ndarray[np.float64_t, ndim=1] rtol = np.array([reltol])

    if intermediate == 1:
        # allocate memory for the intermediate results
        dead_predator_int = np.zeros((dims, niter), "float64")
        dead_prey_int = np.zeros((dims, niter), "float64")
        ge_int = np.zeros((dims, niter), "float64")
        predator_int = np.zeros((dims, niter), "float64")
        prey_int = np.zeros((dims, niter), "float64")
        CoolingTime_int = np.zeros((dims, niter), "float64")
        Gamma_int       = np.zeros((dims, niter), "float64")
        Temperature_int = np.zeros((dims, niter), "float64")
        MolecularWeight_int = np.zeros((dims, niter), "float64")

        temp_int = np.zeros((dims, niter), "float64")
        result_int = np.zeros(niter, "float64")
        t_int = np.zeros(niter, "float64")
        dt_int = np.zeros(niter, "float64")

    j = 0
    for i in range(dims):
        input[j] = prev[j] = dead_predator_arr[i]
        j += 1
        input[j] = prev[j] = dead_prey_arr[i]
        j += 1
        input[j] = prev[j] = ge_arr[i]
        j += 1
        input[j] = prev[j] = predator_arr[i]
        j += 1
        input[j] = prev[j] = prey_arr[i]
        j += 1

    cdef double ttot = 0.0
    cdef int status
    cdef int Enzo_Success = 1
    cdef int Enzo_Fail    = 0
    # Allocate some temporary data
    # Now we manually evolve
    #ttot = dengo_evolve_predator_prey(tf, dt, input, rtol, atol, dims, data)
    cdef double *t_now = <double *> malloc( sizeof(double) )
    cdef double *dt_arr = <double *> malloc(sizeof(double) * dims * niter)
    cdef double *success_arr = <double *> malloc(sizeof(double) * dims * niter)
    cdef double *ttot_arr = <double *> malloc(sizeof(double) * dims * niter)
    cdef double *temp = <double *> malloc(sizeof(double) * niter)
    cdef int    *grid_dimension = <int *> malloc(sizeof(int)*3)
    cdef int    *grid_start     = <int *> malloc(sizeof(int)*3)
    cdef int    *grid_end       = <int *> malloc(sizeof(int)*3)

    
    grid_dimension[0] = dims;
    grid_start[0]     = 0
    grid_end  [0]     = dims-1;

    grid_dimension[1] = 1;
    grid_start[1]     = 0
    grid_end  [1]     = 0;

    grid_dimension[2] = 1;
    grid_start[2]     = 0
    grid_end  [2]     = 0;

    field_data.grid_dimension = grid_dimension
    field_data.grid_start     = grid_start
    field_data.grid_end       = grid_end

    cdef int niter_cvodes = niter
    cdef double dt  = tf / float(niter)
    field_data.dengo_data_file = "predator_prey_tables.h5"

    field_data.reltol = reltol;
    field_data.floor_value = floor_value;

    for iter in range(niter):
        status = predator_prey_solve_chemistry_enzo( units, field_data, dt);
        j = 0;

        for i in range(dims):
            dead_predator_int[i, iter] = field_data.dead_predator_density[i]
            j += 1
            dead_prey_int[i, iter] = field_data.dead_prey_density[i]
            j += 1
            ge_int[i, iter] = field_data.ge_density[i]
            j += 1
            predator_int[i, iter] = field_data.predator_density[i]
            j += 1
            prey_int[i, iter] = field_data.prey_density[i]
            j += 1
            temp_int[ i, iter ] = field_data.temperature[i]

        if status == Enzo_Success:
            result_int[iter] = 1
            ttot += dt
        elif status == Enzo_Fail:
            result_int[iter] = 0

        t_int[iter] = ttot
        dt_int[iter] = dt

        if status == Enzo_Success:
            if iter % 100 == 0:
                print("Successful iteration[% 5i]: (%0.3e) %0.3e / %0.3e" % (iter, dt, ttot, tf))

            if adaptive_step: 
                 dt *= 1.1

            if tf - ttot < dt:
                dt = tf - ttot
        elif status == Enzo_Fail:
            dt /= 2.0;
            # copy_array(prev, input, NTOT)
            # Reset the scaling array to match the new values
            # copy_array(input, scale, NTOT)
            if dt < 1e-20 * tf:
                print("dt too small (%0.3e / %0.3e) so breaking" % (dt, tf))
                break
            continue
        if ttot >= tf: break

    free(input)
    free(prev)
    free(t_now)
    free(temp)
    free(dt_arr)
    free(ttot_arr)
    free(success_arr)
    free(field_data)
    free(units)
    free(grid_dimension)
    free(grid_start)
    free(grid_end)

    print("End in %s iterations: %0.5e / %0.5e (%0.5e)" % (iter + 1, ttot, tf, tf - ttot))

    rv, rv_t = {}, {}
    dead_predator_arr = rv["dead_predator"] = np.zeros(dims, "float64")
    dead_prey_arr = rv["dead_prey"] = np.zeros(dims, "float64")
    ge_arr = rv["ge"] = np.zeros(dims, "float64")
    predator_arr = rv["predator"] = np.zeros(dims, "float64")
    prey_arr = rv["prey"] = np.zeros(dims, "float64")
    if intermediate:
        rv_t["dead_predator"] = dead_predator_int[:niter]
        rv_t["dead_prey"] = dead_prey_int[:niter]
        rv_t["ge"] = ge_int[:niter]
        rv_t["predator"] = predator_int[:niter]
        rv_t["prey"] = prey_int[:niter]
        rv_t["successful"] = result_int.astype("bool")
        rv_t['T'] = temp_int
        rv_t['t'] = t_int
        rv_t['dt'] = dt_int

#    j = 0
#    for i in range(dims):
#        dead_predator_arr[i] = input[j] #* 1.0
#        j += 1
#        dead_prey_arr[i] = input[j] #* 1.0
#        j += 1
#        ge_arr[i] = input[j] #* 1.0
#        j += 1
#        predator_arr[i] = input[j] #* 1.0
#        j += 1
#        prey_arr[i] = input[j] #* 1.0
#        j += 1
    return rv, rv_t

cdef copy_array(double *input, double *output, int dims):
    cdef int i
    for i in range(dims):
        output[i] = input[i]

cdef floor_values(double *input, int dims, double floor):
    cdef int i
    for i in range(dims):
        if input[i] < floor:
            input[i] = floor
"""
cimport numpy as np
import numpy as np
import time
from libc.stdlib cimport malloc, free

cdef extern from "alloca.h":
    void *alloca(int)

# NSPECIES here is N in the .C.template file
DEF NSPECIES = 5
DEF MAX_NCELLS=1024

cdef extern from "predator_prey_solver.h":
    cdef int _MAX_NCELLS  "MAX_NCELLS"
    cdef int _NSPECIES "NSPECIES"
    ctypedef struct predator_prey_data:
        double dbin
        double idbin
        double bounds[2]
        int nbins

        double d_zbin
        double id_zbin
        double z_bounds[2]
        int n_zbins

        double current_z
        double zdef
        double dz

        double Ts[MAX_NCELLS]
        double Tdef[MAX_NCELLS]
        double dT[MAX_NCELLS]
        double logTs[MAX_NCELLS]
        double dTs_ge[MAX_NCELLS]
        double r_exp_growth_prey[1024]
        double rs_exp_growth_prey[MAX_NCELLS]
        double drs_exp_growth_prey[MAX_NCELLS]
        double r_natural_death_predator[1024]
        double rs_natural_death_predator[MAX_NCELLS]
        double drs_natural_death_predator[MAX_NCELLS]
        double r_predation[1024]
        double rs_predation[MAX_NCELLS]
        double drs_predation[MAX_NCELLS]
        int bin_id[MAX_NCELLS]
        int ncells

    ctypedef int(*rhs_f)(double *, double *, int, int, void *)
    ctypedef int(*jac_f)(double *, double *, int, int, void *)

    int predator_prey_main(int argc, char **argv)
    predator_prey_data *predator_prey_setup_data(int *NumberOfFields,
            char ***FieldNames)
    void predator_prey_read_rate_tables(predator_prey_data*)
    void predator_prey_read_cooling_tables(predator_prey_data*)
    double dengo_evolve_predator_prey (double dtf, double &dt, double z,
                                         double *input, double *rtol,
                                         double *atol, int dims,
                                         predator_prey_data *data)
    int BE_chem_solve(rhs_f f, jac_f J,
		    double *u, double dt, double *rtol, 
                    double *atol, int nstrip, int nchem, 
		    double *scaling, void *sdata, double *u0, double *s,
            double *gu, double *Ju
           )
    int calculate_jacobian_predator_prey(double *input, double *Joutput,
            int nstrip, int nchem, void *sdata)
    int calculate_rhs_predator_prey(double *input, double *rhs, int nstrip,
                      int nchem, void *sdata)
    int ensure_electron_consistency(double *input, int nstrip, int nchem)
    void setting_up_extra_variables( simpleBE_data * data, double * input, int nstrip  )

def main_run_predator_prey():
    t1 = time.time()
    predator_prey_main(0, NULL)
    t2 = time.time()
    print "Total elapsed time: %0.3e" % (t2-t1)

def run_predator_prey(ics, double tf, int niter = 10000,
                        int intermediate = 1, z = -1.0, reltol=1.0e-5, floor_value=1.0e-20):
    assert(_MAX_NCELLS == MAX_NCELLS)
    assert(_NSPECIES == NSPECIES)
    cdef np.ndarray[np.float64_t, ndim=1] dead_predator_arr = ics["dead_predator"]
    # All of the intermediate variables get declared, but not necessarily assigned
    cdef np.ndarray[np.float64_t, ndim=2] dead_predator_int
    cdef np.ndarray[np.float64_t, ndim=1] dead_prey_arr = ics["dead_prey"]
    # All of the intermediate variables get declared, but not necessarily assigned
    cdef np.ndarray[np.float64_t, ndim=2] dead_prey_int
    cdef np.ndarray[np.float64_t, ndim=1] ge_arr = ics["ge"]
    # All of the intermediate variables get declared, but not necessarily assigned
    cdef np.ndarray[np.float64_t, ndim=2] ge_int
    cdef np.ndarray[np.float64_t, ndim=1] predator_arr = ics["predator"]
    # All of the intermediate variables get declared, but not necessarily assigned
    cdef np.ndarray[np.float64_t, ndim=2] predator_int
    cdef np.ndarray[np.float64_t, ndim=1] prey_arr = ics["prey"]
    # All of the intermediate variables get declared, but not necessarily assigned
    cdef np.ndarray[np.float64_t, ndim=2] prey_int
    cdef np.ndarray[np.uint8_t, ndim=1] result_int
    cdef np.ndarray[np.float64_t, ndim=2] temp_int
    cdef np.ndarray[np.float64_t, ndim=1] t_int
    cdef np.ndarray[np.float64_t, ndim=1] dt_int

    cdef int i, j, k, iter
    cdef int dims = ge_arr.shape[0]
    cdef int NTOT = NSPECIES * dims
    cdef double *input = <double *> malloc(NTOT * sizeof(double))
    cdef double *prev = <double *> malloc(NTOT * sizeof(double))
    cdef double *atol = <double *> malloc(NTOT * sizeof(double))
    cdef double *rtol = <double *> malloc(NTOT * sizeof(double))
    cdef double *scale = <double *> malloc(NTOT * sizeof(double))

    if intermediate == 1:
        dead_predator_int = np.zeros((dims, niter), "float64")
        dead_prey_int = np.zeros((dims, niter), "float64")
        ge_int = np.zeros((dims, niter), "float64")
        predator_int = np.zeros((dims, niter), "float64")
        prey_int = np.zeros((dims, niter), "float64")
        temp_int = np.zeros((dims, niter), "float64")
        result_int = np.zeros(niter, "uint8")
        t_int = np.zeros(niter, "float64")
        dt_int = np.zeros(niter, "float64")

    j = 0
    for i in range(dims):
        input[j] = prev[j] = dead_predator_arr[i] / 1.0
        atol[j] = input[j] * reltol
        rtol[j] = reltol
        scale[j] = input[j]
        j += 1
        input[j] = prev[j] = dead_prey_arr[i] / 1.0
        atol[j] = input[j] * reltol
        rtol[j] = reltol
        scale[j] = input[j]
        j += 1
        input[j] = prev[j] = ge_arr[i] / 1.0
        atol[j] = input[j] * reltol
        rtol[j] = reltol
        scale[j] = input[j]
        j += 1
        input[j] = prev[j] = predator_arr[i] / 1.0
        atol[j] = input[j] * reltol
        rtol[j] = reltol
        scale[j] = input[j]
        j += 1
        input[j] = prev[j] = prey_arr[i] / 1.0
        atol[j] = input[j] * reltol
        rtol[j] = reltol
        scale[j] = input[j]
        j += 1

    ensure_electron_consistency(input, dims, NSPECIES);

    cdef predator_prey_data *data = predator_prey_setup_data(NULL, NULL)
    cdef rhs_f f = calculate_rhs_predator_prey
    cdef jac_f jf = calculate_jacobian_predator_prey

    cdef double dt = tf / niter
    cdef double ttot = 0.0
    cdef int status
    cdef int Enzo_Success = 1;
    cdef int Enzo_Fail    = 0;
    # Allocate some temporary data
    # Now we manually evolve
    #ttot = dengo_evolve_predator_prey(tf, dt, input, rtol, atol, dims, data)
    data.current_z = z
    cdef double *u0 = <double *> malloc(sizeof(double) * dims * NSPECIES)
    cdef double *s = <double *> malloc(sizeof(double) * NSPECIES)
    cdef double *gu = <double *> malloc(sizeof(double) * dims * NSPECIES)
    cdef double *Ju = <double *> malloc(sizeof(double) * dims * NSPECIES * NSPECIES)
    setting_up_extra_variables( <predator_prey_data *> data, input, dims )
    for iter in range(niter):
        status = BE_chem_solve(f, jf, input, dt, rtol, atol, dims, NSPECIES, scale,
                               <void *> data, u0, s, gu, Ju)
        if intermediate == 1:
            j = 0
            for i in range(dims):
                dead_predator_int[i, iter] = input[j]
                j += 1
                dead_prey_int[i, iter] = input[j]
                j += 1
                ge_int[i, iter] = input[j]
                j += 1
                predator_int[i, iter] = input[j]
                j += 1
                prey_int[i, iter] = input[j]
                j += 1
                # temp_int[i, iter] = data.Ts[i]
            if status == Enzo_Success:
                result_int[iter] = 1
                ttot += dt
            elif status == Enzo_Fail:
                result_int[iter] = 0
            t_int[iter] = ttot
            dt_int[iter] = dt
        if status == Enzo_Success:
            if iter % 1000 == 0:
                print "Successful iteration[% 5i]: (%0.3e) %0.3e / %0.3e" % (
                    iter, dt, ttot, tf)
            copy_array(input, prev, NTOT)
            # Reset the scaling array to match the new values
            copy_array(input, scale, NTOT)
            dt *= 1.1
            if tf - ttot < dt:
                dt = tf - ttot
        elif status == Enzo_Fail:
            dt /= 2.0
            copy_array(prev, input, NTOT)
            # Reset the scaling array to match the new values
            copy_array(input, scale, NTOT)
            if dt < 1e-30 * tf:
                print "dt too small (%0.3e / %0.3e) so breaking" % (dt, tf)
                break
            continue
        if ttot >= tf: break
    free(u0)
    free(s)
    free(gu)
    free(Ju)
    free(input)
    free(prev)
    free(atol)
    free(rtol)
    free(scale)
    free(data)
    print "End in %s iterations: %0.5e / %0.5e (%0.5e)" % (iter + 1, ttot, tf, tf - ttot)

    rv, rv_t = {}, {}
    dead_predator_arr = rv["dead_predator"] = np.zeros(dims, "float64")
    dead_prey_arr = rv["dead_prey"] = np.zeros(dims, "float64")
    ge_arr = rv["ge"] = np.zeros(dims, "float64")
    predator_arr = rv["predator"] = np.zeros(dims, "float64")
    prey_arr = rv["prey"] = np.zeros(dims, "float64")
    if intermediate:
        rv_t["dead_predator"] = dead_predator_int[:niter]* 1.0
        rv_t["dead_prey"] = dead_prey_int[:niter]* 1.0
        rv_t["ge"] = ge_int[:niter]* 1.0
        rv_t["predator"] = predator_int[:niter]* 1.0
        rv_t["prey"] = prey_int[:niter]* 1.0
        rv_t["successful"] = result_int.astype("bool")
        rv_t['T'] = temp_int
        rv_t['t'] = t_int
        rv_t['dt'] = dt_int

    j = 0
    for i in range(dims):
        dead_predator_arr[i] = input[j] * 1.0
        j += 1
        dead_prey_arr[i] = input[j] * 1.0
        j += 1
        ge_arr[i] = input[j] * 1.0
        j += 1
        predator_arr[i] = input[j] * 1.0
        j += 1
        prey_arr[i] = input[j] * 1.0
        j += 1
    return rv, rv_t

cdef copy_array(double *input, double *output, int dims):
    cdef int i
    for i in range(dims):
        output[i] = input[i]

cdef floor_values(double *input, int dims, double floor):
    cdef int i
    for i in range(dims):
        if input[i] < floor:
            input[i] = floor

"""
