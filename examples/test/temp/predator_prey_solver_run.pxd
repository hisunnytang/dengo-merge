DEF MAX_NCELLS = 1
DEF NSPECIES = 5
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
    double dengo_evolve_predator_prey (double dtf, double dt, double z,
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
