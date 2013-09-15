cimport numpy as np
import numpy as np

cdef extern from "alloca.h":
    void *alloca(int)

DEF NSPECIES = 11
DEF MAX_NCELLS = 1024

cdef extern from "oxygen_solver.h":
    ctypedef struct oxygen_data:
        double dbin
        double idbin
        double bounds[2]
        int nbins

        double Ts[MAX_NCELLS]
        double Tdef[MAX_NCELLS]
        double dT[MAX_NCELLS]
        double logTs[MAX_NCELLS]
        double dTs_ge[MAX_NCELLS]
        double r_OI_i[1024]
        double rs_OI_i[MAX_NCELLS]
        double drs_OI_i[MAX_NCELLS]
        double r_OII_i[1024]
        double rs_OII_i[MAX_NCELLS]
        double drs_OII_i[MAX_NCELLS]
        double r_OII_r[1024]
        double rs_OII_r[MAX_NCELLS]
        double drs_OII_r[MAX_NCELLS]
        double r_OIII_i[1024]
        double rs_OIII_i[MAX_NCELLS]
        double drs_OIII_i[MAX_NCELLS]
        double r_OIII_r[1024]
        double rs_OIII_r[MAX_NCELLS]
        double drs_OIII_r[MAX_NCELLS]
        double r_OIV_i[1024]
        double rs_OIV_i[MAX_NCELLS]
        double drs_OIV_i[MAX_NCELLS]
        double r_OIV_r[1024]
        double rs_OIV_r[MAX_NCELLS]
        double drs_OIV_r[MAX_NCELLS]
        double r_OIX_r[1024]
        double rs_OIX_r[MAX_NCELLS]
        double drs_OIX_r[MAX_NCELLS]
        double r_OV_i[1024]
        double rs_OV_i[MAX_NCELLS]
        double drs_OV_i[MAX_NCELLS]
        double r_OV_r[1024]
        double rs_OV_r[MAX_NCELLS]
        double drs_OV_r[MAX_NCELLS]
        double r_OVI_i[1024]
        double rs_OVI_i[MAX_NCELLS]
        double drs_OVI_i[MAX_NCELLS]
        double r_OVI_r[1024]
        double rs_OVI_r[MAX_NCELLS]
        double drs_OVI_r[MAX_NCELLS]
        double r_OVII_i[1024]
        double rs_OVII_i[MAX_NCELLS]
        double drs_OVII_i[MAX_NCELLS]
        double r_OVII_r[1024]
        double rs_OVII_r[MAX_NCELLS]
        double drs_OVII_r[MAX_NCELLS]
        double r_OVIII_i[1024]
        double rs_OVIII_i[MAX_NCELLS]
        double drs_OVIII_i[MAX_NCELLS]
        double r_OVIII_r[1024]
        double rs_OVIII_r[MAX_NCELLS]
        double drs_OVIII_r[MAX_NCELLS]
        double c_OI_c_OI_c[1024]
        double cs_OI_c_OI_c[MAX_NCELLS]
        double dcs_OI_c_OI_c[MAX_NCELLS]
        
        double c_OII_c_OII_c[1024]
        double cs_OII_c_OII_c[MAX_NCELLS]
        double dcs_OII_c_OII_c[MAX_NCELLS]
        
        double c_OIII_c_OIII_c[1024]
        double cs_OIII_c_OIII_c[MAX_NCELLS]
        double dcs_OIII_c_OIII_c[MAX_NCELLS]
        
        double c_OIV_c_OIV_c[1024]
        double cs_OIV_c_OIV_c[MAX_NCELLS]
        double dcs_OIV_c_OIV_c[MAX_NCELLS]
        
        double c_OIX_c_OIX_c[1024]
        double cs_OIX_c_OIX_c[MAX_NCELLS]
        double dcs_OIX_c_OIX_c[MAX_NCELLS]
        
        double c_OV_c_OV_c[1024]
        double cs_OV_c_OV_c[MAX_NCELLS]
        double dcs_OV_c_OV_c[MAX_NCELLS]
        
        double c_OVI_c_OVI_c[1024]
        double cs_OVI_c_OVI_c[MAX_NCELLS]
        double dcs_OVI_c_OVI_c[MAX_NCELLS]
        
        double c_OVII_c_OVII_c[1024]
        double cs_OVII_c_OVII_c[MAX_NCELLS]
        double dcs_OVII_c_OVII_c[MAX_NCELLS]
        
        double c_OVIII_c_OVIII_c[1024]
        double cs_OVIII_c_OVIII_c[MAX_NCELLS]
        double dcs_OVIII_c_OVIII_c[MAX_NCELLS]
        
        int bin_id[MAX_NCELLS]
        int ncells

    ctypedef int(*rhs_f)(double *, double *, int, int, void *)
    ctypedef int(*jac_f)(double *, double *, int, int, void *)

    int oxygen_main(int argc, char **argv)
    oxygen_data *oxygen_setup_data()
    void oxygen_read_rate_tables(oxygen_data*)
    void oxygen_read_cooling_tables(oxygen_data*)
    double dengo_evolve_oxygen (double dtf, double &dt, double *input,
                double *rtol, double *atol, int dims,
                oxygen_data *data)
    int BE_chem_solve(rhs_f f, jac_f J,
		    double *u, double dt, double *rtol, 
                    double *atol, int nstrip, int nchem, 
		    double *scaling, void *sdata)
    int calculate_jacobian_oxygen(double *input, double *Joutput,
            int nstrip, int nchem, void *sdata)
    int calculate_rhs_oxygen(double *input, double *rhs, int nstrip,
                      int nchem, void *sdata)

def main_run_oxygen():
    oxygen_main(0, NULL)

def run_oxygen(ics, double tf, int niter = 10000, int intermediate = 1):
    cdef np.ndarray[np.float64_t, ndim=1] de_arr = ics["de"]
    # All of the intermediate variables get declared, but not necessarily assigned
    cdef np.ndarray[np.float64_t, ndim=2] de_int
    cdef np.ndarray[np.float64_t, ndim=1] ge_arr = ics["ge"]
    # All of the intermediate variables get declared, but not necessarily assigned
    cdef np.ndarray[np.float64_t, ndim=2] ge_int
    cdef np.ndarray[np.float64_t, ndim=1] OI_arr = ics["OI"]
    # All of the intermediate variables get declared, but not necessarily assigned
    cdef np.ndarray[np.float64_t, ndim=2] OI_int
    cdef np.ndarray[np.float64_t, ndim=1] OII_arr = ics["OII"]
    # All of the intermediate variables get declared, but not necessarily assigned
    cdef np.ndarray[np.float64_t, ndim=2] OII_int
    cdef np.ndarray[np.float64_t, ndim=1] OIII_arr = ics["OIII"]
    # All of the intermediate variables get declared, but not necessarily assigned
    cdef np.ndarray[np.float64_t, ndim=2] OIII_int
    cdef np.ndarray[np.float64_t, ndim=1] OIV_arr = ics["OIV"]
    # All of the intermediate variables get declared, but not necessarily assigned
    cdef np.ndarray[np.float64_t, ndim=2] OIV_int
    cdef np.ndarray[np.float64_t, ndim=1] OV_arr = ics["OV"]
    # All of the intermediate variables get declared, but not necessarily assigned
    cdef np.ndarray[np.float64_t, ndim=2] OV_int
    cdef np.ndarray[np.float64_t, ndim=1] OVI_arr = ics["OVI"]
    # All of the intermediate variables get declared, but not necessarily assigned
    cdef np.ndarray[np.float64_t, ndim=2] OVI_int
    cdef np.ndarray[np.float64_t, ndim=1] OVII_arr = ics["OVII"]
    # All of the intermediate variables get declared, but not necessarily assigned
    cdef np.ndarray[np.float64_t, ndim=2] OVII_int
    cdef np.ndarray[np.float64_t, ndim=1] OVIII_arr = ics["OVIII"]
    # All of the intermediate variables get declared, but not necessarily assigned
    cdef np.ndarray[np.float64_t, ndim=2] OVIII_int
    cdef np.ndarray[np.float64_t, ndim=1] OIX_arr = ics["OIX"]
    # All of the intermediate variables get declared, but not necessarily assigned
    cdef np.ndarray[np.float64_t, ndim=2] OIX_int
    cdef np.ndarray[np.uint8_t, ndim=1] result_int
    cdef np.ndarray[np.float64_t, ndim=2] temp_int
    cdef np.ndarray[np.float64_t, ndim=1] t_int
    cdef np.ndarray[np.float64_t, ndim=1] dt_int

    cdef int i, j, k, iter
    cdef int N = ge_arr.shape[0]
    cdef int NTOT = NSPECIES * N
    cdef double *input = <double *> alloca(NTOT * sizeof(double))
    cdef double *prev = <double *> alloca(NTOT * sizeof(double))
    cdef double *atol = <double *> alloca(NTOT * sizeof(double))
    cdef double *rtol = <double *> alloca(NTOT * sizeof(double))
    cdef double *scale = <double *> alloca(NTOT * sizeof(double))
    cdef double v

    if intermediate == 1:
        de_int = np.zeros((N, niter), "float64")
        ge_int = np.zeros((N, niter), "float64")
        OI_int = np.zeros((N, niter), "float64")
        OII_int = np.zeros((N, niter), "float64")
        OIII_int = np.zeros((N, niter), "float64")
        OIV_int = np.zeros((N, niter), "float64")
        OV_int = np.zeros((N, niter), "float64")
        OVI_int = np.zeros((N, niter), "float64")
        OVII_int = np.zeros((N, niter), "float64")
        OVIII_int = np.zeros((N, niter), "float64")
        OIX_int = np.zeros((N, niter), "float64")
        temp_int = np.zeros((N, niter), "float64")
        result_int = np.zeros(niter, "uint8")
        t_int = np.zeros(niter, "float64")
        dt_int = np.zeros(niter, "float64")

    j = 0
    for i in range(N):
        input[j] = prev[j] = de_arr[i] / 1.0
        atol[j] = input[j] * 1e-11
        rtol[j] = 1e-11
        scale[j] = 1.0
        j += 1
        input[j] = prev[j] = ge_arr[i] / 1.0
        atol[j] = input[j] * 1e-11
        rtol[j] = 1e-11
        scale[j] = 1.0
        j += 1
        input[j] = prev[j] = OI_arr[i] / 16
        atol[j] = input[j] * 1e-11
        rtol[j] = 1e-11
        scale[j] = 1.0
        j += 1
        input[j] = prev[j] = OII_arr[i] / 16
        atol[j] = input[j] * 1e-11
        rtol[j] = 1e-11
        scale[j] = 1.0
        j += 1
        input[j] = prev[j] = OIII_arr[i] / 16
        atol[j] = input[j] * 1e-11
        rtol[j] = 1e-11
        scale[j] = 1.0
        j += 1
        input[j] = prev[j] = OIV_arr[i] / 16
        atol[j] = input[j] * 1e-11
        rtol[j] = 1e-11
        scale[j] = 1.0
        j += 1
        input[j] = prev[j] = OV_arr[i] / 16
        atol[j] = input[j] * 1e-11
        rtol[j] = 1e-11
        scale[j] = 1.0
        j += 1
        input[j] = prev[j] = OVI_arr[i] / 16
        atol[j] = input[j] * 1e-11
        rtol[j] = 1e-11
        scale[j] = 1.0
        j += 1
        input[j] = prev[j] = OVII_arr[i] / 16
        atol[j] = input[j] * 1e-11
        rtol[j] = 1e-11
        scale[j] = 1.0
        j += 1
        input[j] = prev[j] = OVIII_arr[i] / 16
        atol[j] = input[j] * 1e-11
        rtol[j] = 1e-11
        scale[j] = 1.0
        j += 1
        input[j] = prev[j] = OIX_arr[i] / 16
        atol[j] = input[j] * 1e-11
        rtol[j] = 1e-11
        scale[j] = 1.0
        j += 1

    cdef oxygen_data *data = oxygen_setup_data()
    cdef rhs_f f = calculate_rhs_oxygen
    cdef jac_f jf = calculate_jacobian_oxygen

    cdef double dt = tf / 1e5
    cdef double ttot = 0.0
    cdef int status
    # Now we manually evolve
    #ttot = dengo_evolve_oxygen(tf, dt, input, rtol, atol, N, data)
    for iter in range(niter):
        status = BE_chem_solve(f, jf, input, dt, rtol, atol, N, NSPECIES, scale,
                               <void *> data)
        if intermediate == 1:
            j = 0
            for i in range(N):
                de_int[i, iter] = input[j]
                j += 1
                ge_int[i, iter] = input[j]
                j += 1
                OI_int[i, iter] = input[j]
                j += 1
                OII_int[i, iter] = input[j]
                j += 1
                OIII_int[i, iter] = input[j]
                j += 1
                OIV_int[i, iter] = input[j]
                j += 1
                OV_int[i, iter] = input[j]
                j += 1
                OVI_int[i, iter] = input[j]
                j += 1
                OVII_int[i, iter] = input[j]
                j += 1
                OVIII_int[i, iter] = input[j]
                j += 1
                OIX_int[i, iter] = input[j]
                j += 1
                temp_int[i, iter] = data.Ts[i]
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
            ttot += dt
            dt *= 1.1
            if tf - ttot < dt:
                dt = tf- ttot
        elif status == 1:
            dt /= 2.0
            copy_array(prev, input, NTOT)
            if dt < 1e-30 * tf:
                print "dt too small (%0.3e / %0.3e) so breaking" % (dt, tf)
                break
            continue
        if ttot >= tf: break

    print "End in %s iterations: %0.5e / %0.5e (%0.5e)" % (iter + 1, ttot, tf, tf - ttot)

    rv, rv_t = {}, {}
    de_arr = rv["de"] = np.zeros(N, "float64")
    ge_arr = rv["ge"] = np.zeros(N, "float64")
    OI_arr = rv["OI"] = np.zeros(N, "float64")
    OII_arr = rv["OII"] = np.zeros(N, "float64")
    OIII_arr = rv["OIII"] = np.zeros(N, "float64")
    OIV_arr = rv["OIV"] = np.zeros(N, "float64")
    OV_arr = rv["OV"] = np.zeros(N, "float64")
    OVI_arr = rv["OVI"] = np.zeros(N, "float64")
    OVII_arr = rv["OVII"] = np.zeros(N, "float64")
    OVIII_arr = rv["OVIII"] = np.zeros(N, "float64")
    OIX_arr = rv["OIX"] = np.zeros(N, "float64")
    if intermediate:
        rv_t["de"] = de_int[:niter]
        rv_t["ge"] = ge_int[:niter]
        rv_t["OI"] = OI_int[:niter]
        rv_t["OII"] = OII_int[:niter]
        rv_t["OIII"] = OIII_int[:niter]
        rv_t["OIV"] = OIV_int[:niter]
        rv_t["OV"] = OV_int[:niter]
        rv_t["OVI"] = OVI_int[:niter]
        rv_t["OVII"] = OVII_int[:niter]
        rv_t["OVIII"] = OVIII_int[:niter]
        rv_t["OIX"] = OIX_int[:niter]
        rv_t["successful"] = result_int.astype("bool")
        rv_t['T'] = temp_int
        rv_t['t'] = t_int
        rv_t['dt'] = dt_int

    j = 0
    for i in range(N):
        de_arr[i] = input[j] * 1.0
        j += 1
        ge_arr[i] = input[j] * 1.0
        j += 1
        OI_arr[i] = input[j] * 16
        j += 1
        OII_arr[i] = input[j] * 16
        j += 1
        OIII_arr[i] = input[j] * 16
        j += 1
        OIV_arr[i] = input[j] * 16
        j += 1
        OV_arr[i] = input[j] * 16
        j += 1
        OVI_arr[i] = input[j] * 16
        j += 1
        OVII_arr[i] = input[j] * 16
        j += 1
        OVIII_arr[i] = input[j] * 16
        j += 1
        OIX_arr[i] = input[j] * 16
        j += 1
    return rv, rv_t

cdef copy_array(double *input, double *output, int N):
    cdef int i
    for i in range(N):
        output[i] = input[i]