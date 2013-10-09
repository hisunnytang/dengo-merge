/*

The generalized rate data type holders.

*/


/* stdlib, hdf5, local includes */

#include "time.h"
#include "sys/time.h"
#include "stdlib.h"
#include "math.h"
#include "hdf5.h"
#include "hdf5_hl.h"
#include "stdio.h"

#define MAX_NCELLS 1024
#define DMAX(A,B) ((A) > (B) ? (A) : (B))
#define DMIN(A,B) ((A) < (B) ? (A) : (B))

 

int umist_main(int argc, char **argv);



typedef struct umist_data_ {
    /* All of the network bins will be the same width */
    double dbin;
    double idbin;
    double bounds[2];
    int nbins;

    double Ts[MAX_NCELLS];
    double Tdef[MAX_NCELLS]; /* t1, t2, tdef */
    double dT[MAX_NCELLS]; /* t1, t2, tdef */
    double logTs[MAX_NCELLS];
    double dTs_ge[MAX_NCELLS];
    /* Now we do all of our cooling and chemical tables */
    double r_us_C_plus_us_Om[1024];
    double rs_us_C_plus_us_Om[MAX_NCELLS];
    double drs_us_C_plus_us_Om[MAX_NCELLS];
    
    double r_us_Cm_plus_us_Cp[1024];
    double rs_us_Cm_plus_us_Cp[MAX_NCELLS];
    double drs_us_Cm_plus_us_Cp[MAX_NCELLS];
    
    double r_us_Cm_plus_us_Hp[1024];
    double rs_us_Cm_plus_us_Hp[MAX_NCELLS];
    double drs_us_Cm_plus_us_Hp[MAX_NCELLS];
    
    double r_us_Cm_plus_us_O[1024];
    double rs_us_Cm_plus_us_O[MAX_NCELLS];
    double drs_us_Cm_plus_us_O[MAX_NCELLS];
    
    double r_us_Cm_plus_us_Op[1024];
    double rs_us_Cm_plus_us_Op[MAX_NCELLS];
    double drs_us_Cm_plus_us_Op[MAX_NCELLS];
    
    double r_us_H2_plus_us_Op[1024];
    double rs_us_H2_plus_us_Op[MAX_NCELLS];
    double drs_us_H2_plus_us_Op[MAX_NCELLS];
    
    double r_us_H2p_plus_us_O[1024];
    double rs_us_H2p_plus_us_O[MAX_NCELLS];
    double drs_us_H2p_plus_us_O[MAX_NCELLS];
    
    double r_us_H2p_plus_us_OH[1024];
    double rs_us_H2p_plus_us_OH[MAX_NCELLS];
    double drs_us_H2p_plus_us_OH[MAX_NCELLS];
    
    double r_us_H_plus_us_H2p[1024];
    double rs_us_H_plus_us_H2p[MAX_NCELLS];
    double drs_us_H_plus_us_H2p[MAX_NCELLS];
    
    double r_us_H_plus_us_Om[1024];
    double rs_us_H_plus_us_Om[MAX_NCELLS];
    double drs_us_H_plus_us_Om[MAX_NCELLS];
    
    double r_us_H_plus_us_Op[1024];
    double rs_us_H_plus_us_Op[MAX_NCELLS];
    double drs_us_H_plus_us_Op[MAX_NCELLS];
    
    double r_us_Hm_plus_us_Cp[1024];
    double rs_us_Hm_plus_us_Cp[MAX_NCELLS];
    double drs_us_Hm_plus_us_Cp[MAX_NCELLS];
    
    double r_us_Hm_plus_us_Hp[1024];
    double rs_us_Hm_plus_us_Hp[MAX_NCELLS];
    double drs_us_Hm_plus_us_Hp[MAX_NCELLS];
    
    double r_us_Hm_plus_us_O[1024];
    double rs_us_Hm_plus_us_O[MAX_NCELLS];
    double drs_us_Hm_plus_us_O[MAX_NCELLS];
    
    double r_us_Hm_plus_us_Op[1024];
    double rs_us_Hm_plus_us_Op[MAX_NCELLS];
    double drs_us_Hm_plus_us_Op[MAX_NCELLS];
    
    double r_us_Hp_plus_us_O[1024];
    double rs_us_Hp_plus_us_O[MAX_NCELLS];
    double drs_us_Hp_plus_us_O[MAX_NCELLS];
    
    double r_us_Hp_plus_us_OH[1024];
    double rs_us_Hp_plus_us_OH[MAX_NCELLS];
    double drs_us_Hp_plus_us_OH[MAX_NCELLS];
    
    double r_us_OHm_plus_us_Cp[1024];
    double rs_us_OHm_plus_us_Cp[MAX_NCELLS];
    double drs_us_OHm_plus_us_Cp[MAX_NCELLS];
    
    double r_us_OHm_plus_us_Hp[1024];
    double rs_us_OHm_plus_us_Hp[MAX_NCELLS];
    double drs_us_OHm_plus_us_Hp[MAX_NCELLS];
    
    double r_us_OHm_plus_us_Op[1024];
    double rs_us_OHm_plus_us_Op[MAX_NCELLS];
    double drs_us_OHm_plus_us_Op[MAX_NCELLS];
    
    double r_us_Om_plus_us_Cp[1024];
    double rs_us_Om_plus_us_Cp[MAX_NCELLS];
    double drs_us_Om_plus_us_Cp[MAX_NCELLS];
    
    double r_us_Om_plus_us_Hp[1024];
    double rs_us_Om_plus_us_Hp[MAX_NCELLS];
    double drs_us_Om_plus_us_Hp[MAX_NCELLS];
    
    double r_us_Om_plus_us_Op[1024];
    double rs_us_Om_plus_us_Op[MAX_NCELLS];
    double drs_us_Om_plus_us_Op[MAX_NCELLS];
    
    int bin_id[MAX_NCELLS];
    int ncells;
} umist_data;

umist_data *umist_setup_data(void);
void umist_read_rate_tables(umist_data*);
void umist_read_cooling_tables(umist_data*);
double dengo_evolve_umist (double dtf, double &dt, double *input,
            double *rtol, double *atol, int dims,
            umist_data *data);
 

typedef int(*rhs_f)(double *, double *, int, int, void *);
typedef int(*jac_f)(double *, double *, int, int, void *);
int BE_chem_solve(rhs_f f, jac_f J, double *u, double dt, double *rtol, 
                  double *atol, int nstrip, int nchem, double *scaling, void *sdata);



int calculate_jacobian_umist(double *input, double *Joutput,
        int nstrip, int nchem, void *sdata);
int calculate_rhs_umist(double *input, double *rhs, int nstrip,
                  int nchem, void *sdata);

 
