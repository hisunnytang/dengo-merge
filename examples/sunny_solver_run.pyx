
cimport numpy as np
import numpy as np
import time
from libc.stdlib cimport malloc, free

cdef extern from "alloca.h":
    void *alloca(int)

# NSPECIES here is N in the .C.template file
DEF NSPECIES = 10
DEF MAX_NCELLS=1024

cdef extern from "sunny_solver.h":
    cdef int _MAX_NCELLS  "MAX_NCELLS"
    cdef int _NSPECIES "NSPECIES"
    ctypedef struct sunny_data:
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
        double r_k01[1024]
        double rs_k01[MAX_NCELLS]
        double drs_k01[MAX_NCELLS]
        double r_k02[1024]
        double rs_k02[MAX_NCELLS]
        double drs_k02[MAX_NCELLS]
        double r_k03[1024]
        double rs_k03[MAX_NCELLS]
        double drs_k03[MAX_NCELLS]
        double r_k04[1024]
        double rs_k04[MAX_NCELLS]
        double drs_k04[MAX_NCELLS]
        double r_k05[1024]
        double rs_k05[MAX_NCELLS]
        double drs_k05[MAX_NCELLS]
        double r_k06[1024]
        double rs_k06[MAX_NCELLS]
        double drs_k06[MAX_NCELLS]
        double r_k07[1024]
        double rs_k07[MAX_NCELLS]
        double drs_k07[MAX_NCELLS]
        double r_k08[1024]
        double rs_k08[MAX_NCELLS]
        double drs_k08[MAX_NCELLS]
        double r_k09[1024]
        double rs_k09[MAX_NCELLS]
        double drs_k09[MAX_NCELLS]
        double r_k10[1024]
        double rs_k10[MAX_NCELLS]
        double drs_k10[MAX_NCELLS]
        double r_k11[1024]
        double rs_k11[MAX_NCELLS]
        double drs_k11[MAX_NCELLS]
        double r_k12[1024]
        double rs_k12[MAX_NCELLS]
        double drs_k12[MAX_NCELLS]
        double r_k13[1024]
        double rs_k13[MAX_NCELLS]
        double drs_k13[MAX_NCELLS]
        double r_k14[1024]
        double rs_k14[MAX_NCELLS]
        double drs_k14[MAX_NCELLS]
        double r_k15[1024]
        double rs_k15[MAX_NCELLS]
        double drs_k15[MAX_NCELLS]
        double r_k16[1024]
        double rs_k16[MAX_NCELLS]
        double drs_k16[MAX_NCELLS]
        double r_k17[1024]
        double rs_k17[MAX_NCELLS]
        double drs_k17[MAX_NCELLS]
        double r_k18[1024]
        double rs_k18[MAX_NCELLS]
        double drs_k18[MAX_NCELLS]
        double r_k19[1024]
        double rs_k19[MAX_NCELLS]
        double drs_k19[MAX_NCELLS]
        double r_k21[1024]
        double rs_k21[MAX_NCELLS]
        double drs_k21[MAX_NCELLS]
        double r_k22[1024]
        double rs_k22[MAX_NCELLS]
        double drs_k22[MAX_NCELLS]
        double r_k23[1024]
        double rs_k23[MAX_NCELLS]
        double drs_k23[MAX_NCELLS]
        int bin_id[MAX_NCELLS]
        int ncells

    ctypedef int(*rhs_f)(double *, double *, int, int, void *)
    ctypedef int(*jac_f)(double *, double *, int, int, void *)

    int sunny_main(int argc, char **argv)
    sunny_data *sunny_setup_data(int *NumberOfFields,
            char ***FieldNames)
    void sunny_read_rate_tables(sunny_data*)
    void sunny_read_cooling_tables(sunny_data*)
    double dengo_evolve_sunny (double dtf, double &dt, double z,
                                         double *input, double *rtol,
                                         double *atol, int dims,
                                         sunny_data *data)
    int BE_chem_solve(rhs_f f, jac_f J,
		    double *u, double dt, double *rtol, 
                    double *atol, int nstrip, int nchem, 
		    double *scaling, void *sdata, double *u0, double *s,
            double *gu, double *Ju
           )
    int calculate_jacobian_sunny(double *input, double *Joutput,
            int nstrip, int nchem, void *sdata)
    int calculate_rhs_sunny(double *input, double *rhs, int nstrip,
                      int nchem, void *sdata)
    int ensure_electron_consistency(double *input, int nstrip, int nchem)

def main_run_sunny():
    t1 = time.time()
    sunny_main(0, NULL)
    t2 = time.time()
    print "Total elapsed time: %0.3e" % (t2-t1)

def run_sunny(ics, double tf, int niter = 10000,
                        int intermediate = 1, z = -1.0):
    assert(_MAX_NCELLS == MAX_NCELLS)
    assert(_NSPECIES == NSPECIES)
    cdef np.ndarray[np.float64_t, ndim=1] H2_1_arr = ics["H2_1"]
    # All of the intermediate variables get declared, but not necessarily assigned
    cdef np.ndarray[np.float64_t, ndim=2] H2_1_int
    cdef np.ndarray[np.float64_t, ndim=1] H2_2_arr = ics["H2_2"]
    # All of the intermediate variables get declared, but not necessarily assigned
    cdef np.ndarray[np.float64_t, ndim=2] H2_2_int
    cdef np.ndarray[np.float64_t, ndim=1] H_1_arr = ics["H_1"]
    # All of the intermediate variables get declared, but not necessarily assigned
    cdef np.ndarray[np.float64_t, ndim=2] H_1_int
    cdef np.ndarray[np.float64_t, ndim=1] H_2_arr = ics["H_2"]
    # All of the intermediate variables get declared, but not necessarily assigned
    cdef np.ndarray[np.float64_t, ndim=2] H_2_int
    cdef np.ndarray[np.float64_t, ndim=1] H_m0_arr = ics["H_m0"]
    # All of the intermediate variables get declared, but not necessarily assigned
    cdef np.ndarray[np.float64_t, ndim=2] H_m0_int
    cdef np.ndarray[np.float64_t, ndim=1] He_1_arr = ics["He_1"]
    # All of the intermediate variables get declared, but not necessarily assigned
    cdef np.ndarray[np.float64_t, ndim=2] He_1_int
    cdef np.ndarray[np.float64_t, ndim=1] He_2_arr = ics["He_2"]
    # All of the intermediate variables get declared, but not necessarily assigned
    cdef np.ndarray[np.float64_t, ndim=2] He_2_int
    cdef np.ndarray[np.float64_t, ndim=1] He_3_arr = ics["He_3"]
    # All of the intermediate variables get declared, but not necessarily assigned
    cdef np.ndarray[np.float64_t, ndim=2] He_3_int
    cdef np.ndarray[np.float64_t, ndim=1] de_arr = ics["de"]
    # All of the intermediate variables get declared, but not necessarily assigned
    cdef np.ndarray[np.float64_t, ndim=2] de_int
    cdef np.ndarray[np.float64_t, ndim=1] ge_arr = ics["ge"]
    # All of the intermediate variables get declared, but not necessarily assigned
    cdef np.ndarray[np.float64_t, ndim=2] ge_int
    cdef np.ndarray[np.uint8_t, ndim=1] result_int
    cdef np.ndarray[np.float64_t, ndim=2] temp_int
    cdef np.ndarray[np.float64_t, ndim=1] t_int
    cdef np.ndarray[np.float64_t, ndim=1] dt_int
    cdef np.ndarray[np.float64_t, ndim=2] Ju_int
    
    
    cdef int i, j, k, iter
    cdef int dims = ge_arr.shape[0]
    cdef int NTOT = NSPECIES * dims
    cdef double *input = <double *> alloca(NTOT * sizeof(double))
    cdef double *prev = <double *> alloca(NTOT * sizeof(double))
    cdef double *atol = <double *> alloca(NTOT * sizeof(double))
    cdef double *rtol = <double *> alloca(NTOT * sizeof(double))
    cdef double *scale = <double *> alloca(NTOT * sizeof(double))
    cdef double v
    cdef double *total_density = <double *> alloca(dims * sizeof(double))
    
    

    if intermediate == 1:
        H2_1_int = np.zeros((dims, niter), "float64")
        H2_2_int = np.zeros((dims, niter), "float64")
        H_1_int = np.zeros((dims, niter), "float64")
        H_2_int = np.zeros((dims, niter), "float64")
        H_m0_int = np.zeros((dims, niter), "float64")
        He_1_int = np.zeros((dims, niter), "float64")
        He_2_int = np.zeros((dims, niter), "float64")
        He_3_int = np.zeros((dims, niter), "float64")
        de_int = np.zeros((dims, niter), "float64")
        ge_int = np.zeros((dims, niter), "float64")
        temp_int = np.zeros((dims, niter), "float64")
        result_int = np.zeros(niter, "uint8")
        t_int = np.zeros(niter, "float64")
        dt_int = np.zeros(niter, "float64")
        
        Ju_int = np.zeros( (dims * NSPECIES * NSPECIES, niter), "float64" )



    j = 0
    for i in range(dims):
        input[j] = prev[j] = H2_1_arr[i] / 2.01588
        atol[j] = input[j] * 1e-12
        rtol[j] = 1e-12
        scale[j] = input[j] 
        print "H2_1", scale[j], H2_1_arr[i]
        j += 1
        input[j] = prev[j] = H2_2_arr[i] / 2.01588
        atol[j] = input[j] * 1e-12
        rtol[j] = 1e-12
        scale[j] = input[j] 
        print "H2_2", scale[j], H2_2_arr[i]
        j += 1
        input[j] = prev[j] = H_1_arr[i] / 1.00794
        atol[j] = input[j] * 1e-12
        rtol[j] = 1e-12
        scale[j] = input[j] 
        print "H_1", scale[j], H_1_arr[i]
        j += 1
        input[j] = prev[j] = H_2_arr[i] / 1.00794
        atol[j] = input[j] * 1e-12
        rtol[j] = 1e-12
        scale[j] = input[j] 
        print "H_2", scale[j], H_2_arr[i]
        j += 1
        input[j] = prev[j] = H_m0_arr[i] / 1.00794
        atol[j] = input[j] * 1e-12
        rtol[j] = 1e-12
        scale[j] = input[j] 
        print "H_m0", scale[j], H_m0_arr[i]
        j += 1
        input[j] = prev[j] = He_1_arr[i] / 4.002602
        atol[j] = input[j] * 1e-12
        rtol[j] = 1e-12
        scale[j] = input[j] 
        print "He_1", scale[j], He_1_arr[i]
        j += 1
        input[j] = prev[j] = He_2_arr[i] / 4.002602
        atol[j] = input[j] * 1e-12
        rtol[j] = 1e-12
        scale[j] = input[j] 
        print "He_2", scale[j], He_2_arr[i]
        j += 1
        input[j] = prev[j] = He_3_arr[i] / 4.002602
        atol[j] = input[j] * 1e-12
        rtol[j] = 1e-12
        scale[j] = input[j] 
        print "He_3", scale[j], He_3_arr[i]
        j += 1
        input[j] = prev[j] = de_arr[i] / 1.0
        atol[j] = input[j] * 1e-12
        rtol[j] = 1e-12
        scale[j] = input[j] 
        print "de", scale[j], de_arr[i]
        j += 1
        input[j] = prev[j] = ge_arr[i] / 1.0
        atol[j] = input[j] * 1e-12
        rtol[j] = 1e-12
        scale[j] = input[j] 
        print "ge", scale[j], ge_arr[i]
        j += 1


    
    
    cdef sunny_data *data = sunny_setup_data(NULL, NULL)
    cdef rhs_f f = calculate_rhs_sunny
    cdef jac_f jf = calculate_jacobian_sunny
    

    cdef double floor_val = 1e-30
    cdef double dt = tf / 1e5
    cdef double ttot = 0.0
    cdef int status
    # Allocate some temporary data
    # Now we manually evolve
    #ttot = dengo_evolve_sunny(tf, dt, input, rtol, atol, dims, data)
    data.current_z = z
    cdef double *u0 = <double *> malloc(sizeof(double) * dims * NSPECIES)
    cdef double *s = <double *> malloc(sizeof(double) * NSPECIES)
    cdef double *gu = <double *> malloc(sizeof(double) * dims * NSPECIES)
    cdef double *Ju = <double *> malloc(sizeof(double) * dims * NSPECIES * NSPECIES)
    
    
    
    
    ensure_electron_consistency(input, dims, NSPECIES); 
    for iter in range(niter):
        # ensure_electron_consistency(input, dims, NSPECIES);
        status = BE_chem_solve(f, jf, input, dt, rtol, atol, dims, NSPECIES, scale,
                               <void *> data, u0, s, gu, Ju)
        
        
        
      

        floor_values( input , NTOT, floor_val )     
        floor_values( prev  , NTOT, floor_val )             

        if intermediate == 1:
            j = 0
            for i in range(dims):
                H2_1_int[i, iter] = input[j]
                j += 1
                H2_2_int[i, iter] = input[j]
                j += 1
                H_1_int[i, iter] = input[j]
                j += 1
                H_2_int[i, iter] = input[j]
                j += 1
                H_m0_int[i, iter] = input[j]
                j += 1
                He_1_int[i, iter] = input[j]
                j += 1
                He_2_int[i, iter] = input[j]
                j += 1
                He_3_int[i, iter] = input[j]
                j += 1
                de_int[i, iter] = input[j]
                j += 1
                ge_int[i, iter] = input[j]
                j += 1                
                temp_int[i, iter] = data.Ts[i]
    
            i = 0
            for Ju_idx in range( dims * NSPECIES * NSPECIES ):
                Ju_ele = Ju[Ju_idx]
                Ju_int[i, iter] = Ju_ele
                i += 1
    
            if status == 0:
                result_int[iter] = 1
            elif status == 1:
                result_int[iter] = 0
            t_int[iter] = ttot
            dt_int[iter] = dt
        if status == 0:
            if iter % 1000 == 0:
                print "Successful iteration[% 5i]: (%0.3e) %0.3e / %0.3e" % (
                    iter, dt, ttot, tf)
            copy_array(input, prev, NTOT)
            # Reset the scaling array to match the new values
            copy_array(input, scale, NTOT)
            ttot += dt
            dt *= 1.1
            if tf - ttot < dt:
                dt = tf - ttot
        elif status == 1:
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

    print "End in %s iterations: %0.5e / %0.5e (%0.5e)" % (iter + 1, ttot, tf, tf - ttot)

    rv, rv_t = {}, {}
    H2_1_arr = rv["H2_1"] = np.zeros(dims, "float64")
    H2_2_arr = rv["H2_2"] = np.zeros(dims, "float64")
    H_1_arr = rv["H_1"] = np.zeros(dims, "float64")
    H_2_arr = rv["H_2"] = np.zeros(dims, "float64")
    H_m0_arr = rv["H_m0"] = np.zeros(dims, "float64")
    He_1_arr = rv["He_1"] = np.zeros(dims, "float64")
    He_2_arr = rv["He_2"] = np.zeros(dims, "float64")
    He_3_arr = rv["He_3"] = np.zeros(dims, "float64")
    de_arr = rv["de"] = np.zeros(dims, "float64")
    ge_arr = rv["ge"] = np.zeros(dims, "float64")
    if intermediate:
        rv_t["H2_1"] = H2_1_int[:niter]
        rv_t["H2_2"] = H2_2_int[:niter]
        rv_t["H_1"] = H_1_int[:niter]
        rv_t["H_2"] = H_2_int[:niter]
        rv_t["H_m0"] = H_m0_int[:niter]
        rv_t["He_1"] = He_1_int[:niter]
        rv_t["He_2"] = He_2_int[:niter]
        rv_t["He_3"] = He_3_int[:niter]
        rv_t["de"] = de_int[:niter]
        rv_t["ge"] = ge_int[:niter]
        rv_t["successful"] = result_int.astype("bool")
        rv_t['T'] = temp_int
        rv_t['t'] = t_int
        rv_t['dt'] = dt_int
    
        rv_t['Ju'] = Ju_int

    j = 0
    for i in range(dims):
        H2_1_arr[i] = input[j] * 2.01588
        j += 1
        H2_2_arr[i] = input[j] * 2.01588
        j += 1
        H_1_arr[i] = input[j] * 1.00794
        j += 1
        H_2_arr[i] = input[j] * 1.00794
        j += 1
        H_m0_arr[i] = input[j] * 1.00794
        j += 1
        He_1_arr[i] = input[j] * 4.002602
        j += 1
        He_2_arr[i] = input[j] * 4.002602
        j += 1
        He_3_arr[i] = input[j] * 4.002602
        j += 1
        de_arr[i] = input[j] * 1.0
        j += 1
        ge_arr[i] = input[j] * 1.0
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