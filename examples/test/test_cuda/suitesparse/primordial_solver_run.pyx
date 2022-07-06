cimport numpy as np

import time

import numpy as np

from cpython cimport array
from libc.stdlib cimport free, malloc
from libcpp cimport bool

# NSPECIES here is N in the .C.template file
DEF NSPECIES = 10
DEF MAX_NCELLS=1

cdef extern from "primordial_solver.h":
    cdef int _MAX_NCELLS  "MAX_NCELLS"
    cdef int _NSPECIES "NSPECIES"

    ctypedef struct dengo_field_data:
        unsigned long int nstrip;
        unsigned long int ncells;
        double *H2_1_density;
        double *H2_2_density;
        double *H_1_density;
        double *H_2_density;
        double *H_m0_density;
        double *He_1_density;
        double *He_2_density;
        double *He_3_density;
        double *de_density;
        double *ge_density;
        double *density;
        double *CoolingTime;
        double *MolecularWeight;
        double *temperature;
        double *Gamma;
        double reltol;
        double floor_value;
        const char *dengo_data_file;

    ctypedef struct code_units:
        int comoving_coordinates
        double density_units
        double length_units
        double time_units
        double velocity_units
        double a_units
        double a_value


    ctypedef struct primordial_data:
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
        double c_brem_brem[1024]
        double cs_brem_brem[MAX_NCELLS]
        double dcs_brem_brem[MAX_NCELLS]

        double c_ceHeI_ceHeI[1024]
        double cs_ceHeI_ceHeI[MAX_NCELLS]
        double dcs_ceHeI_ceHeI[MAX_NCELLS]

        double c_ceHeII_ceHeII[1024]
        double cs_ceHeII_ceHeII[MAX_NCELLS]
        double dcs_ceHeII_ceHeII[MAX_NCELLS]

        double c_ceHI_ceHI[1024]
        double cs_ceHI_ceHI[MAX_NCELLS]
        double dcs_ceHI_ceHI[MAX_NCELLS]

        double c_cie_cooling_cieco[1024]
        double cs_cie_cooling_cieco[MAX_NCELLS]
        double dcs_cie_cooling_cieco[MAX_NCELLS]

        double c_ciHeI_ciHeI[1024]
        double cs_ciHeI_ciHeI[MAX_NCELLS]
        double dcs_ciHeI_ciHeI[MAX_NCELLS]

        double c_ciHeII_ciHeII[1024]
        double cs_ciHeII_ciHeII[MAX_NCELLS]
        double dcs_ciHeII_ciHeII[MAX_NCELLS]

        double c_ciHeIS_ciHeIS[1024]
        double cs_ciHeIS_ciHeIS[MAX_NCELLS]
        double dcs_ciHeIS_ciHeIS[MAX_NCELLS]

        double c_ciHI_ciHI[1024]
        double cs_ciHI_ciHI[MAX_NCELLS]
        double dcs_ciHI_ciHI[MAX_NCELLS]

        double c_compton_comp_[1024]
        double cs_compton_comp_[MAX_NCELLS]
        double dcs_compton_comp_[MAX_NCELLS]

        double c_gammah_gammah[1024]
        double cs_gammah_gammah[MAX_NCELLS]
        double dcs_gammah_gammah[MAX_NCELLS]

        double c_gloverabel08_gael[1024]
        double cs_gloverabel08_gael[MAX_NCELLS]
        double dcs_gloverabel08_gael[MAX_NCELLS]
        double c_gloverabel08_gaH2[1024]
        double cs_gloverabel08_gaH2[MAX_NCELLS]
        double dcs_gloverabel08_gaH2[MAX_NCELLS]
        double c_gloverabel08_gaHe[1024]
        double cs_gloverabel08_gaHe[MAX_NCELLS]
        double dcs_gloverabel08_gaHe[MAX_NCELLS]
        double c_gloverabel08_gaHI[1024]
        double cs_gloverabel08_gaHI[MAX_NCELLS]
        double dcs_gloverabel08_gaHI[MAX_NCELLS]
        double c_gloverabel08_gaHp[1024]
        double cs_gloverabel08_gaHp[MAX_NCELLS]
        double dcs_gloverabel08_gaHp[MAX_NCELLS]
        double c_gloverabel08_gphdl[1024]
        double cs_gloverabel08_gphdl[MAX_NCELLS]
        double dcs_gloverabel08_gphdl[MAX_NCELLS]
        double c_gloverabel08_gpldl[1024]
        double cs_gloverabel08_gpldl[MAX_NCELLS]
        double dcs_gloverabel08_gpldl[MAX_NCELLS]
        double c_gloverabel08_h2lte[1024]
        double cs_gloverabel08_h2lte[MAX_NCELLS]
        double dcs_gloverabel08_h2lte[MAX_NCELLS]

        double c_h2formation_h2mcool[1024]
        double cs_h2formation_h2mcool[MAX_NCELLS]
        double dcs_h2formation_h2mcool[MAX_NCELLS]
        double c_h2formation_h2mheat[1024]
        double cs_h2formation_h2mheat[MAX_NCELLS]
        double dcs_h2formation_h2mheat[MAX_NCELLS]
        double c_h2formation_ncrd1[1024]
        double cs_h2formation_ncrd1[MAX_NCELLS]
        double dcs_h2formation_ncrd1[MAX_NCELLS]
        double c_h2formation_ncrd2[1024]
        double cs_h2formation_ncrd2[MAX_NCELLS]
        double dcs_h2formation_ncrd2[MAX_NCELLS]
        double c_h2formation_ncrn[1024]
        double cs_h2formation_ncrn[MAX_NCELLS]
        double dcs_h2formation_ncrn[MAX_NCELLS]

        double c_reHeII1_reHeII1[1024]
        double cs_reHeII1_reHeII1[MAX_NCELLS]
        double dcs_reHeII1_reHeII1[MAX_NCELLS]

        double c_reHeII2_reHeII2[1024]
        double cs_reHeII2_reHeII2[MAX_NCELLS]
        double dcs_reHeII2_reHeII2[MAX_NCELLS]

        double c_reHeIII_reHeIII[1024]
        double cs_reHeIII_reHeIII[MAX_NCELLS]
        double dcs_reHeIII_reHeIII[MAX_NCELLS]

        double c_reHII_reHII[1024]
        double cs_reHII_reHII[MAX_NCELLS]
        double dcs_reHII_reHII[MAX_NCELLS]

        int bin_id[MAX_NCELLS]
        int ncells

    double dengo_evolve_primordial (double dtf, double &dt, double z,
                                         double *input, double *rtol,
                                         double *atol, int dims,
                                         primordial_data *data, double *temp)
    int primordial_solve_chemistry_dt( code_units *units, dengo_field_data *field_data, double* rtol, double* atol, double dt)
    int primordial_solve_chemistry( code_units *units, dengo_field_data *field_data, double dt );
    int primordial_main(int argc, char **argv);



def main_run_primordial():
    t1 = time.time()
    primordial_main(0, NULL)
    t2 = time.time()
    print("Total elapsed time: %0.3e" % (t2-t1))

def run_primordial(ics, double tf, int niter = 10000,
                        int intermediate = 1, float z = -1.0, bool adaptive_step = True,
			float reltol = 1.0e-5, float floor_value = 1.0e-20):
    #assert(_MAX_NCELLS == MAX_NCELLS)
    assert(_NSPECIES == NSPECIES)
    cdef np.ndarray[np.float64_t, ndim=1] H2_1_arr = ics["H2_1"]
    if not H2_1_arr.flags['C_CONTIGUOUS']:
        H2_1_arr = np.ascontiguousarray( H2_1_arr )
    cdef double[::1] H2_1_memview = H2_1_arr
    cdef np.ndarray[np.float64_t, ndim=2] H2_1_int
    cdef np.ndarray[np.float64_t, ndim=1] H2_2_arr = ics["H2_2"]
    if not H2_2_arr.flags['C_CONTIGUOUS']:
        H2_2_arr = np.ascontiguousarray( H2_2_arr )
    cdef double[::1] H2_2_memview = H2_2_arr
    cdef np.ndarray[np.float64_t, ndim=2] H2_2_int
    cdef np.ndarray[np.float64_t, ndim=1] H_1_arr = ics["H_1"]
    if not H_1_arr.flags['C_CONTIGUOUS']:
        H_1_arr = np.ascontiguousarray( H_1_arr )
    cdef double[::1] H_1_memview = H_1_arr
    cdef np.ndarray[np.float64_t, ndim=2] H_1_int
    cdef np.ndarray[np.float64_t, ndim=1] H_2_arr = ics["H_2"]
    if not H_2_arr.flags['C_CONTIGUOUS']:
        H_2_arr = np.ascontiguousarray( H_2_arr )
    cdef double[::1] H_2_memview = H_2_arr
    cdef np.ndarray[np.float64_t, ndim=2] H_2_int
    cdef np.ndarray[np.float64_t, ndim=1] H_m0_arr = ics["H_m0"]
    if not H_m0_arr.flags['C_CONTIGUOUS']:
        H_m0_arr = np.ascontiguousarray( H_m0_arr )
    cdef double[::1] H_m0_memview = H_m0_arr
    cdef np.ndarray[np.float64_t, ndim=2] H_m0_int
    cdef np.ndarray[np.float64_t, ndim=1] He_1_arr = ics["He_1"]
    if not He_1_arr.flags['C_CONTIGUOUS']:
        He_1_arr = np.ascontiguousarray( He_1_arr )
    cdef double[::1] He_1_memview = He_1_arr
    cdef np.ndarray[np.float64_t, ndim=2] He_1_int
    cdef np.ndarray[np.float64_t, ndim=1] He_2_arr = ics["He_2"]
    if not He_2_arr.flags['C_CONTIGUOUS']:
        He_2_arr = np.ascontiguousarray( He_2_arr )
    cdef double[::1] He_2_memview = He_2_arr
    cdef np.ndarray[np.float64_t, ndim=2] He_2_int
    cdef np.ndarray[np.float64_t, ndim=1] He_3_arr = ics["He_3"]
    if not He_3_arr.flags['C_CONTIGUOUS']:
        He_3_arr = np.ascontiguousarray( He_3_arr )
    cdef double[::1] He_3_memview = He_3_arr
    cdef np.ndarray[np.float64_t, ndim=2] He_3_int
    cdef np.ndarray[np.float64_t, ndim=1] de_arr = ics["de"]
    if not de_arr.flags['C_CONTIGUOUS']:
        de_arr = np.ascontiguousarray( de_arr )
    cdef double[::1] de_memview = de_arr
    cdef np.ndarray[np.float64_t, ndim=2] de_int
    cdef np.ndarray[np.float64_t, ndim=1] ge_arr = ics["ge"]
    if not ge_arr.flags['C_CONTIGUOUS']:
        ge_arr = np.ascontiguousarray( ge_arr )
    cdef double[::1] ge_memview = ge_arr
    cdef np.ndarray[np.float64_t, ndim=2] ge_int
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
    field_data.H2_1_density = &H2_1_memview[0];
    field_data.H2_2_density = &H2_2_memview[0];
    field_data.H_1_density = &H_1_memview[0];
    field_data.H_2_density = &H_2_memview[0];
    field_data.H_m0_density = &H_m0_memview[0];
    field_data.He_1_density = &He_1_memview[0];
    field_data.He_2_density = &He_2_memview[0];
    field_data.He_3_density = &He_3_memview[0];
    field_data.de_density = &de_memview[0];
    field_data.ge_density = &ge_memview[0];
    field_data.CoolingTime = &ctime_memview[0]
    field_data.Gamma       = &gamma_memview[0]
    field_data.temperature = &temp_memview[0]
    field_data.MolecularWeight = &mweight_memview[0]
    field_data.ncells  = dims;

    cdef np.ndarray[np.float64_t, ndim=1] atol = np.ones( (NTOT) )*floor_value*reltol
    cdef np.ndarray[np.float64_t, ndim=1] rtol = np.array([reltol])

    if intermediate == 1:
        # allocate memory for the intermediate results
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
        input[j] = prev[j] = H2_1_arr[i]
        j += 1
        input[j] = prev[j] = H2_2_arr[i]
        j += 1
        input[j] = prev[j] = H_1_arr[i]
        j += 1
        input[j] = prev[j] = H_2_arr[i]
        j += 1
        input[j] = prev[j] = H_m0_arr[i]
        j += 1
        input[j] = prev[j] = He_1_arr[i]
        j += 1
        input[j] = prev[j] = He_2_arr[i]
        j += 1
        input[j] = prev[j] = He_3_arr[i]
        j += 1
        input[j] = prev[j] = de_arr[i]
        j += 1
        input[j] = prev[j] = ge_arr[i]
        j += 1

    cdef double ttot = 0.0
    cdef int status
    # Allocate some temporary data
    # Now we manually evolve
    #ttot = dengo_evolve_primordial(tf, dt, input, rtol, atol, dims, data)
    cdef double *t_now = <double *> malloc( sizeof(double) )
    cdef double *dt_arr = <double *> malloc(sizeof(double) * dims * niter)
    cdef double *success_arr = <double *> malloc(sizeof(double) * dims * niter)
    cdef double *ttot_arr = <double *> malloc(sizeof(double) * dims * niter)
    cdef double *temp = <double *> malloc(sizeof(double) * niter)

    cdef int niter_cvodes = niter
    cdef double dt  = tf / float(niter)
    field_data.dengo_data_file = "primordial_tables.h5"

    field_data.reltol = reltol;
    field_data.floor_value = floor_value;

    for iter in range(niter):
        status = primordial_solve_chemistry( units, field_data, dt);
        j = 0;

        for i in range(dims):
            H2_1_int[i, iter] = field_data.H2_1_density[i]
            j += 1
            H2_2_int[i, iter] = field_data.H2_2_density[i]
            j += 1
            H_1_int[i, iter] = field_data.H_1_density[i]
            j += 1
            H_2_int[i, iter] = field_data.H_2_density[i]
            j += 1
            H_m0_int[i, iter] = field_data.H_m0_density[i]
            j += 1
            He_1_int[i, iter] = field_data.He_1_density[i]
            j += 1
            He_2_int[i, iter] = field_data.He_2_density[i]
            j += 1
            He_3_int[i, iter] = field_data.He_3_density[i]
            j += 1
            de_int[i, iter] = field_data.de_density[i]
            j += 1
            ge_int[i, iter] = field_data.ge_density[i]
            j += 1
            temp_int[ i, iter ] = field_data.temperature[i]

        if status == 0:
            result_int[iter] = 1
            ttot += dt
        elif status == 1:
            result_int[iter] = 0

        t_int[iter] = ttot
        dt_int[iter] = dt

        if status == 0:
            if iter % 100 == 0:
                print("Successful iteration[% 5i]: (%0.3e) %0.3e / %0.3e" % (iter, dt, ttot, tf))

            if adaptive_step:
                 dt *= 1.1

            if tf - ttot < dt:
                dt = tf - ttot
        elif status == 1:
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

    print("End in %s iterations: %0.5e / %0.5e (%0.5e)" % (iter + 1, ttot, tf, tf - ttot))

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

    j = 0
    for i in range(dims):
        H2_1_arr[i] = input[j] #* 2.0
        j += 1
        H2_2_arr[i] = input[j] #* 2.0
        j += 1
        H_1_arr[i] = input[j] #* 1.00794
        j += 1
        H_2_arr[i] = input[j] #* 1.00794
        j += 1
        H_m0_arr[i] = input[j] #* 1.00794
        j += 1
        He_1_arr[i] = input[j] #* 4.002602
        j += 1
        He_2_arr[i] = input[j] #* 4.002602
        j += 1
        He_3_arr[i] = input[j] #* 4.002602
        j += 1
        de_arr[i] = input[j] #* 1.0
        j += 1
        ge_arr[i] = input[j] #* 1.0
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
