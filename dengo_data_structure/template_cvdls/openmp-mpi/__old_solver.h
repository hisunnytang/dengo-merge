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

#define MAX_NCELLS 1024
#define NSPECIES 10
#define DMAX(A,B) ((A) > (B) ? (A) : (B))
#define DMIN(A,B) ((A) < (B) ? (A) : (B))

 

int __old_main(int argc, char **argv);



typedef struct __old_data_ {
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
    double r_k01[1024];
    double rs_k01[MAX_NCELLS];
    double drs_k01[MAX_NCELLS];
    
    double r_k02[1024];
    double rs_k02[MAX_NCELLS];
    double drs_k02[MAX_NCELLS];
    
    double r_k03[1024];
    double rs_k03[MAX_NCELLS];
    double drs_k03[MAX_NCELLS];
    
    double r_k04[1024];
    double rs_k04[MAX_NCELLS];
    double drs_k04[MAX_NCELLS];
    
    double r_k05[1024];
    double rs_k05[MAX_NCELLS];
    double drs_k05[MAX_NCELLS];
    
    double r_k06[1024];
    double rs_k06[MAX_NCELLS];
    double drs_k06[MAX_NCELLS];
    
    double r_k07[1024];
    double rs_k07[MAX_NCELLS];
    double drs_k07[MAX_NCELLS];
    
    double r_k08[1024];
    double rs_k08[MAX_NCELLS];
    double drs_k08[MAX_NCELLS];
    
    double r_k09[1024];
    double rs_k09[MAX_NCELLS];
    double drs_k09[MAX_NCELLS];
    
    double r_k10[1024];
    double rs_k10[MAX_NCELLS];
    double drs_k10[MAX_NCELLS];
    
    double r_k11[1024];
    double rs_k11[MAX_NCELLS];
    double drs_k11[MAX_NCELLS];
    
    double r_k12[1024];
    double rs_k12[MAX_NCELLS];
    double drs_k12[MAX_NCELLS];
    
    double r_k13[1024];
    double rs_k13[MAX_NCELLS];
    double drs_k13[MAX_NCELLS];
    
    double r_k14[1024];
    double rs_k14[MAX_NCELLS];
    double drs_k14[MAX_NCELLS];
    
    double r_k15[1024];
    double rs_k15[MAX_NCELLS];
    double drs_k15[MAX_NCELLS];
    
    double r_k16[1024];
    double rs_k16[MAX_NCELLS];
    double drs_k16[MAX_NCELLS];
    
    double r_k17[1024];
    double rs_k17[MAX_NCELLS];
    double drs_k17[MAX_NCELLS];
    
    double r_k18[1024];
    double rs_k18[MAX_NCELLS];
    double drs_k18[MAX_NCELLS];
    
    double r_k19[1024];
    double rs_k19[MAX_NCELLS];
    double drs_k19[MAX_NCELLS];
    
    double r_k21[1024];
    double rs_k21[MAX_NCELLS];
    double drs_k21[MAX_NCELLS];
    
    double r_k22[1024];
    double rs_k22[MAX_NCELLS];
    double drs_k22[MAX_NCELLS];
    
    double c_brem_brem[1024];
    double cs_brem_brem[MAX_NCELLS];
    double dcs_brem_brem[MAX_NCELLS];
    
    double c_ceHeI_ceHeI[1024];
    double cs_ceHeI_ceHeI[MAX_NCELLS];
    double dcs_ceHeI_ceHeI[MAX_NCELLS];
    
    double c_ceHeII_ceHeII[1024];
    double cs_ceHeII_ceHeII[MAX_NCELLS];
    double dcs_ceHeII_ceHeII[MAX_NCELLS];
    
    double c_ceHI_ceHI[1024];
    double cs_ceHI_ceHI[MAX_NCELLS];
    double dcs_ceHI_ceHI[MAX_NCELLS];
    
    double c_ciHeI_ciHeI[1024];
    double cs_ciHeI_ciHeI[MAX_NCELLS];
    double dcs_ciHeI_ciHeI[MAX_NCELLS];
    
    double c_ciHeII_ciHeII[1024];
    double cs_ciHeII_ciHeII[MAX_NCELLS];
    double dcs_ciHeII_ciHeII[MAX_NCELLS];
    
    double c_ciHeIS_ciHeIS[1024];
    double cs_ciHeIS_ciHeIS[MAX_NCELLS];
    double dcs_ciHeIS_ciHeIS[MAX_NCELLS];
    
    double c_ciHI_ciHI[1024];
    double cs_ciHI_ciHI[MAX_NCELLS];
    double dcs_ciHI_ciHI[MAX_NCELLS];
    
    double c_compton_comp_[1024];
    double cs_compton_comp_[MAX_NCELLS];
    double dcs_compton_comp_[MAX_NCELLS];
    
    double c_gammah_gammah[1024];
    double cs_gammah_gammah[MAX_NCELLS];
    double dcs_gammah_gammah[MAX_NCELLS];
    
    double c_gloverabel08_gael[1024];
    double cs_gloverabel08_gael[MAX_NCELLS];
    double dcs_gloverabel08_gael[MAX_NCELLS];
    double c_gloverabel08_gaH2[1024];
    double cs_gloverabel08_gaH2[MAX_NCELLS];
    double dcs_gloverabel08_gaH2[MAX_NCELLS];
    double c_gloverabel08_gaHe[1024];
    double cs_gloverabel08_gaHe[MAX_NCELLS];
    double dcs_gloverabel08_gaHe[MAX_NCELLS];
    double c_gloverabel08_gaHI[1024];
    double cs_gloverabel08_gaHI[MAX_NCELLS];
    double dcs_gloverabel08_gaHI[MAX_NCELLS];
    double c_gloverabel08_gaHp[1024];
    double cs_gloverabel08_gaHp[MAX_NCELLS];
    double dcs_gloverabel08_gaHp[MAX_NCELLS];
    double c_gloverabel08_gphdl[1024];
    double cs_gloverabel08_gphdl[MAX_NCELLS];
    double dcs_gloverabel08_gphdl[MAX_NCELLS];
    double c_gloverabel08_gpldl[1024];
    double cs_gloverabel08_gpldl[MAX_NCELLS];
    double dcs_gloverabel08_gpldl[MAX_NCELLS];
    double c_gloverabel08_h2lte[1024];
    double cs_gloverabel08_h2lte[MAX_NCELLS];
    double dcs_gloverabel08_h2lte[MAX_NCELLS];
    
    double c_h2formation_h2mcool[1024];
    double cs_h2formation_h2mcool[MAX_NCELLS];
    double dcs_h2formation_h2mcool[MAX_NCELLS];
    double c_h2formation_h2mheat[1024];
    double cs_h2formation_h2mheat[MAX_NCELLS];
    double dcs_h2formation_h2mheat[MAX_NCELLS];
    double c_h2formation_ncrd1[1024];
    double cs_h2formation_ncrd1[MAX_NCELLS];
    double dcs_h2formation_ncrd1[MAX_NCELLS];
    double c_h2formation_ncrd2[1024];
    double cs_h2formation_ncrd2[MAX_NCELLS];
    double dcs_h2formation_ncrd2[MAX_NCELLS];
    double c_h2formation_ncrn[1024];
    double cs_h2formation_ncrn[MAX_NCELLS];
    double dcs_h2formation_ncrn[MAX_NCELLS];
    
    double c_reHeII1_reHeII1[1024];
    double cs_reHeII1_reHeII1[MAX_NCELLS];
    double dcs_reHeII1_reHeII1[MAX_NCELLS];
    
    double c_reHeII2_reHeII2[1024];
    double cs_reHeII2_reHeII2[MAX_NCELLS];
    double dcs_reHeII2_reHeII2[MAX_NCELLS];
    
    double c_reHeIII_reHeIII[1024];
    double cs_reHeIII_reHeIII[MAX_NCELLS];
    double dcs_reHeIII_reHeIII[MAX_NCELLS];
    
    double c_reHII_reHII[1024];
    double cs_reHII_reHII[MAX_NCELLS];
    double dcs_reHII_reHII[MAX_NCELLS];
    
    int bin_id[MAX_NCELLS];
    int ncells;
} __old_data;

__old_data *__old_setup_data(int *, char***);
void __old_read_rate_tables(__old_data*);
void __old_read_cooling_tables(__old_data*);
double dengo_evolve___old (double dtf, double &dt, double z,
                                     double *input, double *rtol,
                                     double *atol, long long dims,
                                     __old_data *data);
 

typedef int(*rhs_f)(double *, double *, int, int, void *);
typedef int(*jac_f)(double *, double *, int, int, void *);
int BE_chem_solve(rhs_f f, jac_f J, double *u, double dt, double *rtol, 
                  double *atol, int nstrip, int nchem, double *scaling, void *sdata,
                  double *, double *, double *, double *);



int calculate_jacobian___old(double *input, double *Joutput,
        int nstrip, int nchem, void *sdata);
int calculate_rhs___old(double *input, double *rhs, int nstrip,
                  int nchem, void *sdata);
void ensure_electron_consistency(double *input, int nstrip, int nchem);
void temperature_from_mass_density(double *input, int nstrip, int nchem, 
                                   double *strip_temperature);

 
