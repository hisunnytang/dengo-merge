cimport numpy as np
import numpy as np

cdef extern from "alloca.h":
    void *alloca(int)

DEF NSPECIES = 17
DEF MAX_NCELLS = 1024

cdef extern from "umist_solver.h":
    ctypedef struct umist_data:
        double dbin
        double idbin
        double bounds[2]
        int nbins

        double Ts[MAX_NCELLS]
        double Tdef[MAX_NCELLS]
        double dT[MAX_NCELLS]
        double logTs[MAX_NCELLS]
        double dTs_ge[MAX_NCELLS]
        double r_us_C_plus_us_Om[1024]
        double rs_us_C_plus_us_Om[MAX_NCELLS]
        double drs_us_C_plus_us_Om[MAX_NCELLS]
        double r_us_Cm_plus_us_Cp[1024]
        double rs_us_Cm_plus_us_Cp[MAX_NCELLS]
        double drs_us_Cm_plus_us_Cp[MAX_NCELLS]
        double r_us_Cm_plus_us_Hp[1024]
        double rs_us_Cm_plus_us_Hp[MAX_NCELLS]
        double drs_us_Cm_plus_us_Hp[MAX_NCELLS]
        double r_us_Cm_plus_us_O[1024]
        double rs_us_Cm_plus_us_O[MAX_NCELLS]
        double drs_us_Cm_plus_us_O[MAX_NCELLS]
        double r_us_Cm_plus_us_Op[1024]
        double rs_us_Cm_plus_us_Op[MAX_NCELLS]
        double drs_us_Cm_plus_us_Op[MAX_NCELLS]
        double r_us_H2_plus_us_Op[1024]
        double rs_us_H2_plus_us_Op[MAX_NCELLS]
        double drs_us_H2_plus_us_Op[MAX_NCELLS]
        double r_us_H2p_plus_us_O[1024]
        double rs_us_H2p_plus_us_O[MAX_NCELLS]
        double drs_us_H2p_plus_us_O[MAX_NCELLS]
        double r_us_H2p_plus_us_OH[1024]
        double rs_us_H2p_plus_us_OH[MAX_NCELLS]
        double drs_us_H2p_plus_us_OH[MAX_NCELLS]
        double r_us_H_plus_us_H2p[1024]
        double rs_us_H_plus_us_H2p[MAX_NCELLS]
        double drs_us_H_plus_us_H2p[MAX_NCELLS]
        double r_us_H_plus_us_Om[1024]
        double rs_us_H_plus_us_Om[MAX_NCELLS]
        double drs_us_H_plus_us_Om[MAX_NCELLS]
        double r_us_H_plus_us_Op[1024]
        double rs_us_H_plus_us_Op[MAX_NCELLS]
        double drs_us_H_plus_us_Op[MAX_NCELLS]
        double r_us_Hm_plus_us_Cp[1024]
        double rs_us_Hm_plus_us_Cp[MAX_NCELLS]
        double drs_us_Hm_plus_us_Cp[MAX_NCELLS]
        double r_us_Hm_plus_us_Hp[1024]
        double rs_us_Hm_plus_us_Hp[MAX_NCELLS]
        double drs_us_Hm_plus_us_Hp[MAX_NCELLS]
        double r_us_Hm_plus_us_O[1024]
        double rs_us_Hm_plus_us_O[MAX_NCELLS]
        double drs_us_Hm_plus_us_O[MAX_NCELLS]
        double r_us_Hm_plus_us_Op[1024]
        double rs_us_Hm_plus_us_Op[MAX_NCELLS]
        double drs_us_Hm_plus_us_Op[MAX_NCELLS]
        double r_us_Hp_plus_us_O[1024]
        double rs_us_Hp_plus_us_O[MAX_NCELLS]
        double drs_us_Hp_plus_us_O[MAX_NCELLS]
        double r_us_Hp_plus_us_OH[1024]
        double rs_us_Hp_plus_us_OH[MAX_NCELLS]
        double drs_us_Hp_plus_us_OH[MAX_NCELLS]
        double r_us_OHm_plus_us_Cp[1024]
        double rs_us_OHm_plus_us_Cp[MAX_NCELLS]
        double drs_us_OHm_plus_us_Cp[MAX_NCELLS]
        double r_us_OHm_plus_us_Hp[1024]
        double rs_us_OHm_plus_us_Hp[MAX_NCELLS]
        double drs_us_OHm_plus_us_Hp[MAX_NCELLS]
        double r_us_OHm_plus_us_Op[1024]
        double rs_us_OHm_plus_us_Op[MAX_NCELLS]
        double drs_us_OHm_plus_us_Op[MAX_NCELLS]
        double r_us_Om_plus_us_Cp[1024]
        double rs_us_Om_plus_us_Cp[MAX_NCELLS]
        double drs_us_Om_plus_us_Cp[MAX_NCELLS]
        double r_us_Om_plus_us_Hp[1024]
        double rs_us_Om_plus_us_Hp[MAX_NCELLS]
        double drs_us_Om_plus_us_Hp[MAX_NCELLS]
        double r_us_Om_plus_us_Op[1024]
        double rs_us_Om_plus_us_Op[MAX_NCELLS]
        double drs_us_Om_plus_us_Op[MAX_NCELLS]
        int bin_id[MAX_NCELLS]
        int ncells

    ctypedef int(*rhs_f)(double *, double *, int, int, void *)
    ctypedef int(*jac_f)(double *, double *, int, int, void *)

    int umist_main(int argc, char **argv)
    umist_data *umist_setup_data()
    void umist_read_rate_tables(umist_data*)
    void umist_read_cooling_tables(umist_data*)
    double dengo_evolve_umist (double dtf, double &dt, double *input,
                double *rtol, double *atol, int dims,
                umist_data *data)
    int BE_chem_solve(rhs_f f, jac_f J,
		    double *u, double dt, double *rtol, 
                    double *atol, int nstrip, int nchem, 
		    double *scaling, void *sdata)
    int calculate_jacobian_umist(double *input, double *Joutput,
            int nstrip, int nchem, void *sdata)
    int calculate_rhs_umist(double *input, double *rhs, int nstrip,
                      int nchem, void *sdata)

def main_run_umist():
    umist_main(0, NULL)

def run_umist(ics, double tf, int niter = 10000, int intermediate = 1):
    cdef np.ndarray[np.float64_t, ndim=1] ge_arr = ics["ge"]
    # All of the intermediate variables get declared, but not necessarily assigned
    cdef np.ndarray[np.float64_t, ndim=2] ge_int
    cdef np.ndarray[np.float64_t, ndim=1] us_CO_arr = ics["us_CO"]
    # All of the intermediate variables get declared, but not necessarily assigned
    cdef np.ndarray[np.float64_t, ndim=2] us_CO_int
    cdef np.ndarray[np.float64_t, ndim=1] us_Cm_arr = ics["us_Cm"]
    # All of the intermediate variables get declared, but not necessarily assigned
    cdef np.ndarray[np.float64_t, ndim=2] us_Cm_int
    cdef np.ndarray[np.float64_t, ndim=1] us_em_arr = ics["us_em"]
    # All of the intermediate variables get declared, but not necessarily assigned
    cdef np.ndarray[np.float64_t, ndim=2] us_em_int
    cdef np.ndarray[np.float64_t, ndim=1] us_O_arr = ics["us_O"]
    # All of the intermediate variables get declared, but not necessarily assigned
    cdef np.ndarray[np.float64_t, ndim=2] us_O_int
    cdef np.ndarray[np.float64_t, ndim=1] us_C_arr = ics["us_C"]
    # All of the intermediate variables get declared, but not necessarily assigned
    cdef np.ndarray[np.float64_t, ndim=2] us_C_int
    cdef np.ndarray[np.float64_t, ndim=1] us_Om_arr = ics["us_Om"]
    # All of the intermediate variables get declared, but not necessarily assigned
    cdef np.ndarray[np.float64_t, ndim=2] us_Om_int
    cdef np.ndarray[np.float64_t, ndim=1] us_OHm_arr = ics["us_OHm"]
    # All of the intermediate variables get declared, but not necessarily assigned
    cdef np.ndarray[np.float64_t, ndim=2] us_OHm_int
    cdef np.ndarray[np.float64_t, ndim=1] us_Hm_arr = ics["us_Hm"]
    # All of the intermediate variables get declared, but not necessarily assigned
    cdef np.ndarray[np.float64_t, ndim=2] us_Hm_int
    cdef np.ndarray[np.float64_t, ndim=1] us_Cp_arr = ics["us_Cp"]
    # All of the intermediate variables get declared, but not necessarily assigned
    cdef np.ndarray[np.float64_t, ndim=2] us_Cp_int
    cdef np.ndarray[np.float64_t, ndim=1] us_H2_arr = ics["us_H2"]
    # All of the intermediate variables get declared, but not necessarily assigned
    cdef np.ndarray[np.float64_t, ndim=2] us_H2_int
    cdef np.ndarray[np.float64_t, ndim=1] us_H2p_arr = ics["us_H2p"]
    # All of the intermediate variables get declared, but not necessarily assigned
    cdef np.ndarray[np.float64_t, ndim=2] us_H2p_int
    cdef np.ndarray[np.float64_t, ndim=1] us_H_arr = ics["us_H"]
    # All of the intermediate variables get declared, but not necessarily assigned
    cdef np.ndarray[np.float64_t, ndim=2] us_H_int
    cdef np.ndarray[np.float64_t, ndim=1] us_Hp_arr = ics["us_Hp"]
    # All of the intermediate variables get declared, but not necessarily assigned
    cdef np.ndarray[np.float64_t, ndim=2] us_Hp_int
    cdef np.ndarray[np.float64_t, ndim=1] us_Op_arr = ics["us_Op"]
    # All of the intermediate variables get declared, but not necessarily assigned
    cdef np.ndarray[np.float64_t, ndim=2] us_Op_int
    cdef np.ndarray[np.float64_t, ndim=1] us_OHp_arr = ics["us_OHp"]
    # All of the intermediate variables get declared, but not necessarily assigned
    cdef np.ndarray[np.float64_t, ndim=2] us_OHp_int
    cdef np.ndarray[np.float64_t, ndim=1] us_OH_arr = ics["us_OH"]
    # All of the intermediate variables get declared, but not necessarily assigned
    cdef np.ndarray[np.float64_t, ndim=2] us_OH_int
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
        ge_int = np.zeros((N, niter), "float64")
        us_CO_int = np.zeros((N, niter), "float64")
        us_Cm_int = np.zeros((N, niter), "float64")
        us_em_int = np.zeros((N, niter), "float64")
        us_O_int = np.zeros((N, niter), "float64")
        us_C_int = np.zeros((N, niter), "float64")
        us_Om_int = np.zeros((N, niter), "float64")
        us_OHm_int = np.zeros((N, niter), "float64")
        us_Hm_int = np.zeros((N, niter), "float64")
        us_Cp_int = np.zeros((N, niter), "float64")
        us_H2_int = np.zeros((N, niter), "float64")
        us_H2p_int = np.zeros((N, niter), "float64")
        us_H_int = np.zeros((N, niter), "float64")
        us_Hp_int = np.zeros((N, niter), "float64")
        us_Op_int = np.zeros((N, niter), "float64")
        us_OHp_int = np.zeros((N, niter), "float64")
        us_OH_int = np.zeros((N, niter), "float64")
        temp_int = np.zeros((N, niter), "float64")
        result_int = np.zeros(niter, "uint8")
        t_int = np.zeros(niter, "float64")
        dt_int = np.zeros(niter, "float64")

    j = 0
    for i in range(N):
        input[j] = prev[j] = ge_arr[i] / 1.0
        atol[j] = input[j] * 1e-11
        rtol[j] = 1e-11
        scale[j] = 1.0
        j += 1
        input[j] = prev[j] = us_CO_arr[i] / 28
        atol[j] = input[j] * 1e-11
        rtol[j] = 1e-11
        scale[j] = 1.0
        j += 1
        input[j] = prev[j] = us_Cm_arr[i] / 1.0
        atol[j] = input[j] * 1e-11
        rtol[j] = 1e-11
        scale[j] = 1.0
        j += 1
        input[j] = prev[j] = us_em_arr[i] / 1.0
        atol[j] = input[j] * 1e-11
        rtol[j] = 1e-11
        scale[j] = 1.0
        j += 1
        input[j] = prev[j] = us_O_arr[i] / 1.0
        atol[j] = input[j] * 1e-11
        rtol[j] = 1e-11
        scale[j] = 1.0
        j += 1
        input[j] = prev[j] = us_C_arr[i] / 1.0
        atol[j] = input[j] * 1e-11
        rtol[j] = 1e-11
        scale[j] = 1.0
        j += 1
        input[j] = prev[j] = us_Om_arr[i] / 1.0
        atol[j] = input[j] * 1e-11
        rtol[j] = 1e-11
        scale[j] = 1.0
        j += 1
        input[j] = prev[j] = us_OHm_arr[i] / 1.0
        atol[j] = input[j] * 1e-11
        rtol[j] = 1e-11
        scale[j] = 1.0
        j += 1
        input[j] = prev[j] = us_Hm_arr[i] / 1.0
        atol[j] = input[j] * 1e-11
        rtol[j] = 1e-11
        scale[j] = 1.0
        j += 1
        input[j] = prev[j] = us_Cp_arr[i] / 1.0
        atol[j] = input[j] * 1e-11
        rtol[j] = 1e-11
        scale[j] = 1.0
        j += 1
        input[j] = prev[j] = us_H2_arr[i] / 1.0
        atol[j] = input[j] * 1e-11
        rtol[j] = 1e-11
        scale[j] = 1.0
        j += 1
        input[j] = prev[j] = us_H2p_arr[i] / 1.0
        atol[j] = input[j] * 1e-11
        rtol[j] = 1e-11
        scale[j] = 1.0
        j += 1
        input[j] = prev[j] = us_H_arr[i] / 1.0
        atol[j] = input[j] * 1e-11
        rtol[j] = 1e-11
        scale[j] = 1.0
        j += 1
        input[j] = prev[j] = us_Hp_arr[i] / 1.0
        atol[j] = input[j] * 1e-11
        rtol[j] = 1e-11
        scale[j] = 1.0
        j += 1
        input[j] = prev[j] = us_Op_arr[i] / 1.0
        atol[j] = input[j] * 1e-11
        rtol[j] = 1e-11
        scale[j] = 1.0
        j += 1
        input[j] = prev[j] = us_OHp_arr[i] / 1.0
        atol[j] = input[j] * 1e-11
        rtol[j] = 1e-11
        scale[j] = 1.0
        j += 1
        input[j] = prev[j] = us_OH_arr[i] / 1.0
        atol[j] = input[j] * 1e-11
        rtol[j] = 1e-11
        scale[j] = 1.0
        j += 1

    cdef umist_data *data = umist_setup_data()
    cdef rhs_f f = calculate_rhs_umist
    cdef jac_f jf = calculate_jacobian_umist

    cdef double dt = tf / 1e5
    cdef double ttot = 0.0
    cdef int status
    # Now we manually evolve
    #ttot = dengo_evolve_umist(tf, dt, input, rtol, atol, N, data)
    for iter in range(niter):
        status = BE_chem_solve(f, jf, input, dt, rtol, atol, N, NSPECIES, scale,
                               <void *> data)
        if intermediate == 1:
            j = 0
            for i in range(N):
                ge_int[i, iter] = input[j]
                j += 1
                us_CO_int[i, iter] = input[j]
                j += 1
                us_Cm_int[i, iter] = input[j]
                j += 1
                us_em_int[i, iter] = input[j]
                j += 1
                us_O_int[i, iter] = input[j]
                j += 1
                us_C_int[i, iter] = input[j]
                j += 1
                us_Om_int[i, iter] = input[j]
                j += 1
                us_OHm_int[i, iter] = input[j]
                j += 1
                us_Hm_int[i, iter] = input[j]
                j += 1
                us_Cp_int[i, iter] = input[j]
                j += 1
                us_H2_int[i, iter] = input[j]
                j += 1
                us_H2p_int[i, iter] = input[j]
                j += 1
                us_H_int[i, iter] = input[j]
                j += 1
                us_Hp_int[i, iter] = input[j]
                j += 1
                us_Op_int[i, iter] = input[j]
                j += 1
                us_OHp_int[i, iter] = input[j]
                j += 1
                us_OH_int[i, iter] = input[j]
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
    ge_arr = rv["ge"] = np.zeros(N, "float64")
    us_CO_arr = rv["us_CO"] = np.zeros(N, "float64")
    us_Cm_arr = rv["us_Cm"] = np.zeros(N, "float64")
    us_em_arr = rv["us_em"] = np.zeros(N, "float64")
    us_O_arr = rv["us_O"] = np.zeros(N, "float64")
    us_C_arr = rv["us_C"] = np.zeros(N, "float64")
    us_Om_arr = rv["us_Om"] = np.zeros(N, "float64")
    us_OHm_arr = rv["us_OHm"] = np.zeros(N, "float64")
    us_Hm_arr = rv["us_Hm"] = np.zeros(N, "float64")
    us_Cp_arr = rv["us_Cp"] = np.zeros(N, "float64")
    us_H2_arr = rv["us_H2"] = np.zeros(N, "float64")
    us_H2p_arr = rv["us_H2p"] = np.zeros(N, "float64")
    us_H_arr = rv["us_H"] = np.zeros(N, "float64")
    us_Hp_arr = rv["us_Hp"] = np.zeros(N, "float64")
    us_Op_arr = rv["us_Op"] = np.zeros(N, "float64")
    us_OHp_arr = rv["us_OHp"] = np.zeros(N, "float64")
    us_OH_arr = rv["us_OH"] = np.zeros(N, "float64")
    if intermediate:
        rv_t["ge"] = ge_int[:niter]
        rv_t["us_CO"] = us_CO_int[:niter]
        rv_t["us_Cm"] = us_Cm_int[:niter]
        rv_t["us_em"] = us_em_int[:niter]
        rv_t["us_O"] = us_O_int[:niter]
        rv_t["us_C"] = us_C_int[:niter]
        rv_t["us_Om"] = us_Om_int[:niter]
        rv_t["us_OHm"] = us_OHm_int[:niter]
        rv_t["us_Hm"] = us_Hm_int[:niter]
        rv_t["us_Cp"] = us_Cp_int[:niter]
        rv_t["us_H2"] = us_H2_int[:niter]
        rv_t["us_H2p"] = us_H2p_int[:niter]
        rv_t["us_H"] = us_H_int[:niter]
        rv_t["us_Hp"] = us_Hp_int[:niter]
        rv_t["us_Op"] = us_Op_int[:niter]
        rv_t["us_OHp"] = us_OHp_int[:niter]
        rv_t["us_OH"] = us_OH_int[:niter]
        rv_t["successful"] = result_int.astype("bool")
        rv_t['T'] = temp_int
        rv_t['t'] = t_int
        rv_t['dt'] = dt_int

    j = 0
    for i in range(N):
        ge_arr[i] = input[j] * 1.0
        j += 1
        us_CO_arr[i] = input[j] * 28
        j += 1
        us_Cm_arr[i] = input[j] * 1.0
        j += 1
        us_em_arr[i] = input[j] * 1.0
        j += 1
        us_O_arr[i] = input[j] * 1.0
        j += 1
        us_C_arr[i] = input[j] * 1.0
        j += 1
        us_Om_arr[i] = input[j] * 1.0
        j += 1
        us_OHm_arr[i] = input[j] * 1.0
        j += 1
        us_Hm_arr[i] = input[j] * 1.0
        j += 1
        us_Cp_arr[i] = input[j] * 1.0
        j += 1
        us_H2_arr[i] = input[j] * 1.0
        j += 1
        us_H2p_arr[i] = input[j] * 1.0
        j += 1
        us_H_arr[i] = input[j] * 1.0
        j += 1
        us_Hp_arr[i] = input[j] * 1.0
        j += 1
        us_Op_arr[i] = input[j] * 1.0
        j += 1
        us_OHp_arr[i] = input[j] * 1.0
        j += 1
        us_OH_arr[i] = input[j] * 1.0
        j += 1
    return rv, rv_t

cdef copy_array(double *input, double *output, int N):
    cdef int i
    for i in range(N):
        output[i] = input[i]