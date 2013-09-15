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

 

int oxygen_main(int argc, char **argv);



typedef struct oxygen_data_ {
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
    double r_OI_i[1024];
    double rs_OI_i[MAX_NCELLS];
    double drs_OI_i[MAX_NCELLS];
    
    double r_OII_i[1024];
    double rs_OII_i[MAX_NCELLS];
    double drs_OII_i[MAX_NCELLS];
    
    double r_OII_r[1024];
    double rs_OII_r[MAX_NCELLS];
    double drs_OII_r[MAX_NCELLS];
    
    double r_OIII_i[1024];
    double rs_OIII_i[MAX_NCELLS];
    double drs_OIII_i[MAX_NCELLS];
    
    double r_OIII_r[1024];
    double rs_OIII_r[MAX_NCELLS];
    double drs_OIII_r[MAX_NCELLS];
    
    double r_OIV_i[1024];
    double rs_OIV_i[MAX_NCELLS];
    double drs_OIV_i[MAX_NCELLS];
    
    double r_OIV_r[1024];
    double rs_OIV_r[MAX_NCELLS];
    double drs_OIV_r[MAX_NCELLS];
    
    double r_OIX_r[1024];
    double rs_OIX_r[MAX_NCELLS];
    double drs_OIX_r[MAX_NCELLS];
    
    double r_OV_i[1024];
    double rs_OV_i[MAX_NCELLS];
    double drs_OV_i[MAX_NCELLS];
    
    double r_OV_r[1024];
    double rs_OV_r[MAX_NCELLS];
    double drs_OV_r[MAX_NCELLS];
    
    double r_OVI_i[1024];
    double rs_OVI_i[MAX_NCELLS];
    double drs_OVI_i[MAX_NCELLS];
    
    double r_OVI_r[1024];
    double rs_OVI_r[MAX_NCELLS];
    double drs_OVI_r[MAX_NCELLS];
    
    double r_OVII_i[1024];
    double rs_OVII_i[MAX_NCELLS];
    double drs_OVII_i[MAX_NCELLS];
    
    double r_OVII_r[1024];
    double rs_OVII_r[MAX_NCELLS];
    double drs_OVII_r[MAX_NCELLS];
    
    double r_OVIII_i[1024];
    double rs_OVIII_i[MAX_NCELLS];
    double drs_OVIII_i[MAX_NCELLS];
    
    double r_OVIII_r[1024];
    double rs_OVIII_r[MAX_NCELLS];
    double drs_OVIII_r[MAX_NCELLS];
    
    double c_OI_c_OI_c[1024];
    double cs_OI_c_OI_c[MAX_NCELLS];
    double dcs_OI_c_OI_c[MAX_NCELLS];
    
    
    double c_OII_c_OII_c[1024];
    double cs_OII_c_OII_c[MAX_NCELLS];
    double dcs_OII_c_OII_c[MAX_NCELLS];
    
    
    double c_OIII_c_OIII_c[1024];
    double cs_OIII_c_OIII_c[MAX_NCELLS];
    double dcs_OIII_c_OIII_c[MAX_NCELLS];
    
    
    double c_OIV_c_OIV_c[1024];
    double cs_OIV_c_OIV_c[MAX_NCELLS];
    double dcs_OIV_c_OIV_c[MAX_NCELLS];
    
    
    double c_OIX_c_OIX_c[1024];
    double cs_OIX_c_OIX_c[MAX_NCELLS];
    double dcs_OIX_c_OIX_c[MAX_NCELLS];
    
    
    double c_OV_c_OV_c[1024];
    double cs_OV_c_OV_c[MAX_NCELLS];
    double dcs_OV_c_OV_c[MAX_NCELLS];
    
    
    double c_OVI_c_OVI_c[1024];
    double cs_OVI_c_OVI_c[MAX_NCELLS];
    double dcs_OVI_c_OVI_c[MAX_NCELLS];
    
    
    double c_OVII_c_OVII_c[1024];
    double cs_OVII_c_OVII_c[MAX_NCELLS];
    double dcs_OVII_c_OVII_c[MAX_NCELLS];
    
    
    double c_OVIII_c_OVIII_c[1024];
    double cs_OVIII_c_OVIII_c[MAX_NCELLS];
    double dcs_OVIII_c_OVIII_c[MAX_NCELLS];
    
    
    int bin_id[MAX_NCELLS];
    int ncells;
} oxygen_data;

oxygen_data *oxygen_setup_data(void);
void oxygen_read_rate_tables(oxygen_data*);
void oxygen_read_cooling_tables(oxygen_data*);
double dengo_evolve_oxygen (double dtf, double &dt, double *input,
            double *rtol, double *atol, int dims,
            oxygen_data *data);
 

typedef int(*rhs_f)(double *, double *, int, int, void *);
typedef int(*jac_f)(double *, double *, int, int, void *);
int BE_chem_solve(rhs_f f, jac_f J, double *u, double dt, double *rtol, 
                  double *atol, int nstrip, int nchem, double *scaling, void *sdata);



int calculate_jacobian_oxygen(double *input, double *Joutput,
        int nstrip, int nchem, void *sdata);
int calculate_rhs_oxygen(double *input, double *rhs, int nstrip,
                  int nchem, void *sdata);

 
