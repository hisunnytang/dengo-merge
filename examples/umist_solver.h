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
#include "string.h"

#define MAX_NCELLS 1
#define DMAX(A,B) ((A) > (B) ? (A) : (B))
#define DMIN(A,B) ((A) < (B) ? (A) : (B))

 

int umist_main(int argc, char **argv);



typedef struct umist_data_ {
    /* All of the network bins will be the same width */
    double dbin;
    double idbin;
    double bounds[2];
    int nbins;

    /* These will be for bins in redshift space */
    double d_zbin;
    double id_zbin;
    double z_bounds[2];
    int n_zbins;

    /* For storing and passing around
       redshift information */
    double current_z;
    double zdef;
    double dz;

    double Ts[MAX_NCELLS];
    double Tdef[MAX_NCELLS]; /* t1, t2, tdef */
    double dT[MAX_NCELLS]; /* t1, t2, tdef */
    double logTs[MAX_NCELLS];
    double invTs[MAX_NCELLS];
    double dTs_ge[MAX_NCELLS];

    /* Now we do all of our cooling and chemical tables */
    double r_us_H2_1_plus_us_H2_1[1024];
    double rs_us_H2_1_plus_us_H2_1[MAX_NCELLS];
    double drs_us_H2_1_plus_us_H2_1[MAX_NCELLS];
    
    double r_us_H_1_plus_us_H2_1[1024];
    double rs_us_H_1_plus_us_H2_1[MAX_NCELLS];
    double drs_us_H_1_plus_us_H2_1[MAX_NCELLS];
    
    int bin_id[MAX_NCELLS];
    int ncells;
} umist_data;

umist_data *umist_setup_data(int *, char***);
void umist_read_rate_tables(umist_data*);
void umist_read_cooling_tables(umist_data*);
double dengo_evolve_umist (double dtf, double &dt, double z,
                                     double *input, double *rtol,
                                     double *atol, long long dims,
                                     umist_data *data);
 

typedef int(*rhs_f)(double *, double *, int, int, void *);
typedef int(*jac_f)(double *, double *, int, int, void *);
int BE_chem_solve(rhs_f f, jac_f J, double *u, double dt, double *rtol, 
                  double *atol, int nstrip, int nchem, double *scaling, void *sdata,
                  double *, double *, double *, double *);



int calculate_jacobian_umist(double *input, double *Joutput,
        int nstrip, int nchem, void *sdata);
int calculate_rhs_umist(double *input, double *rhs, int nstrip,
                  int nchem, void *sdata);
void ensure_electron_consistency(double *input, int nstrip, int nchem);
void temperature_from_mass_density(double *input, int nstrip, int nchem, 
                                   double *strip_temperature);

 
