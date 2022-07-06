/*

The generalized rate data type holders.

*/


/* stdlib, hdf5, local includes */

#include "omp.h"

#include "time.h"
#include "sys/time.h"
#include "stdlib.h"
#include "math.h"
#include "hdf5.h"
#include "hdf5_hl.h"
#include "stdio.h"
#include "string.h"

/* header files for CVODES/SUNDIALS */
#include <cvode/cvode.h>               /* prototypes for CVODE fcts., consts.  */
#include <nvector/nvector_serial.h>    /* access to serial N_Vector            */
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix            */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver      */
#include <cvode/cvode_direct.h>        /* access to CVDls interface            */
#include <sundials/sundials_types.h>   /* defs. of realtype, sunindextype      */

#ifdef CVKLU
#include <sunmatrix/sunmatrix_sparse.h>
#include <sunlinsol/sunlinsol_klu.h>
#endif


/* User-defined vector and matrix accessor macros: Ith, IJth */

/* These macros are defined in order to write code which exactly matches
   the mathematical problem description given above.

   Ith(v,i) references the ith component of the vector v, where i is in
   the range [1..NEQ] and NEQ is defined below. The Ith macro is defined
   using the N_VIth macro in nvector.h. N_VIth numbers the components of
   a vector starting from 0.

   IJth(A,i,j) references the (i,j)th element of the dense matrix A, where
   i and j are in the range [1..NEQ]. The IJth macro is defined using the
   DENSE_ELEM macro in dense.h. DENSE_ELEM numbers rows and columns of a
   dense matrix starting from 0. */

#define Ith(v,i)    NV_Ith_S(v,i-1)       /* Ith numbers components 1..NEQ */
#define IJth(A,i,j) SM_ELEMENT_D(A,i-1,j-1) /* IJth numbers rows,cols 1..NEQ */


#ifndef MAX_NCELLS
#define MAX_NCELLS 1024
#endif

#ifndef NTHREADS
#define NTHREADS 8
#endif

#define NSPECIES 10
#define DMAX(A,B) ((A) > (B) ? (A) : (B))
#define DMIN(A,B) ((A) < (B) ? (A) : (B))



int primordial_main(int argc, char **argv);



typedef struct primordial_data {
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

    double Ts[NTHREADS][MAX_NCELLS];
    double Tdef[NTHREADS][MAX_NCELLS]; /* t1, t2, tdef */
    double dT[NTHREADS][MAX_NCELLS]; /* t1, t2, tdef */
    double logTs[NTHREADS][MAX_NCELLS];
    double invTs[NTHREADS][MAX_NCELLS];
    double dTs_ge[NTHREADS][MAX_NCELLS];

    /* Now we do all of our cooling and chemical tables */
    double r_k01[1024];
    double rs_k01[NTHREADS][MAX_NCELLS];
    double drs_k01[NTHREADS][MAX_NCELLS];

    double r_k02[1024];
    double rs_k02[NTHREADS][MAX_NCELLS];
    double drs_k02[NTHREADS][MAX_NCELLS];

    double r_k03[1024];
    double rs_k03[NTHREADS][MAX_NCELLS];
    double drs_k03[NTHREADS][MAX_NCELLS];

    double r_k04[1024];
    double rs_k04[NTHREADS][MAX_NCELLS];
    double drs_k04[NTHREADS][MAX_NCELLS];

    double r_k05[1024];
    double rs_k05[NTHREADS][MAX_NCELLS];
    double drs_k05[NTHREADS][MAX_NCELLS];

    double r_k06[1024];
    double rs_k06[NTHREADS][MAX_NCELLS];
    double drs_k06[NTHREADS][MAX_NCELLS];

    double r_k07[1024];
    double rs_k07[NTHREADS][MAX_NCELLS];
    double drs_k07[NTHREADS][MAX_NCELLS];

    double r_k08[1024];
    double rs_k08[NTHREADS][MAX_NCELLS];
    double drs_k08[NTHREADS][MAX_NCELLS];

    double r_k09[1024];
    double rs_k09[NTHREADS][MAX_NCELLS];
    double drs_k09[NTHREADS][MAX_NCELLS];

    double r_k10[1024];
    double rs_k10[NTHREADS][MAX_NCELLS];
    double drs_k10[NTHREADS][MAX_NCELLS];

    double r_k11[1024];
    double rs_k11[NTHREADS][MAX_NCELLS];
    double drs_k11[NTHREADS][MAX_NCELLS];

    double r_k12[1024];
    double rs_k12[NTHREADS][MAX_NCELLS];
    double drs_k12[NTHREADS][MAX_NCELLS];

    double r_k13[1024];
    double rs_k13[NTHREADS][MAX_NCELLS];
    double drs_k13[NTHREADS][MAX_NCELLS];

    double r_k14[1024];
    double rs_k14[NTHREADS][MAX_NCELLS];
    double drs_k14[NTHREADS][MAX_NCELLS];

    double r_k15[1024];
    double rs_k15[NTHREADS][MAX_NCELLS];
    double drs_k15[NTHREADS][MAX_NCELLS];

    double r_k16[1024];
    double rs_k16[NTHREADS][MAX_NCELLS];
    double drs_k16[NTHREADS][MAX_NCELLS];

    double r_k17[1024];
    double rs_k17[NTHREADS][MAX_NCELLS];
    double drs_k17[NTHREADS][MAX_NCELLS];

    double r_k18[1024];
    double rs_k18[NTHREADS][MAX_NCELLS];
    double drs_k18[NTHREADS][MAX_NCELLS];

    double r_k19[1024];
    double rs_k19[NTHREADS][MAX_NCELLS];
    double drs_k19[NTHREADS][MAX_NCELLS];

    double r_k21[1024];
    double rs_k21[NTHREADS][MAX_NCELLS];
    double drs_k21[NTHREADS][MAX_NCELLS];

    double r_k22[1024];
    double rs_k22[NTHREADS][MAX_NCELLS];
    double drs_k22[NTHREADS][MAX_NCELLS];

    double r_k23[1024];
    double rs_k23[NTHREADS][MAX_NCELLS];
    double drs_k23[NTHREADS][MAX_NCELLS];

    double c_brem_brem[1024];
    double cs_brem_brem[NTHREADS][MAX_NCELLS];
    double dcs_brem_brem[NTHREADS][MAX_NCELLS];

    double c_ceHeI_ceHeI[1024];
    double cs_ceHeI_ceHeI[NTHREADS][MAX_NCELLS];
    double dcs_ceHeI_ceHeI[NTHREADS][MAX_NCELLS];

    double c_ceHeII_ceHeII[1024];
    double cs_ceHeII_ceHeII[NTHREADS][MAX_NCELLS];
    double dcs_ceHeII_ceHeII[NTHREADS][MAX_NCELLS];

    double c_ceHI_ceHI[1024];
    double cs_ceHI_ceHI[NTHREADS][MAX_NCELLS];
    double dcs_ceHI_ceHI[NTHREADS][MAX_NCELLS];

    double c_cie_cooling_cieco[1024];
    double cs_cie_cooling_cieco[NTHREADS][MAX_NCELLS];
    double dcs_cie_cooling_cieco[NTHREADS][MAX_NCELLS];

    double c_ciHeI_ciHeI[1024];
    double cs_ciHeI_ciHeI[NTHREADS][MAX_NCELLS];
    double dcs_ciHeI_ciHeI[NTHREADS][MAX_NCELLS];

    double c_ciHeII_ciHeII[1024];
    double cs_ciHeII_ciHeII[NTHREADS][MAX_NCELLS];
    double dcs_ciHeII_ciHeII[NTHREADS][MAX_NCELLS];

    double c_ciHeIS_ciHeIS[1024];
    double cs_ciHeIS_ciHeIS[NTHREADS][MAX_NCELLS];
    double dcs_ciHeIS_ciHeIS[NTHREADS][MAX_NCELLS];

    double c_ciHI_ciHI[1024];
    double cs_ciHI_ciHI[NTHREADS][MAX_NCELLS];
    double dcs_ciHI_ciHI[NTHREADS][MAX_NCELLS];

    double c_compton_comp_[1024];
    double cs_compton_comp_[NTHREADS][MAX_NCELLS];
    double dcs_compton_comp_[NTHREADS][MAX_NCELLS];

    double c_gammah_gammah[1024];
    double cs_gammah_gammah[NTHREADS][MAX_NCELLS];
    double dcs_gammah_gammah[NTHREADS][MAX_NCELLS];

    double c_gloverabel08_gael[1024];
    double cs_gloverabel08_gael[NTHREADS][MAX_NCELLS];
    double dcs_gloverabel08_gael[NTHREADS][MAX_NCELLS];
    double c_gloverabel08_gaH2[1024];
    double cs_gloverabel08_gaH2[NTHREADS][MAX_NCELLS];
    double dcs_gloverabel08_gaH2[NTHREADS][MAX_NCELLS];
    double c_gloverabel08_gaHe[1024];
    double cs_gloverabel08_gaHe[NTHREADS][MAX_NCELLS];
    double dcs_gloverabel08_gaHe[NTHREADS][MAX_NCELLS];
    double c_gloverabel08_gaHI[1024];
    double cs_gloverabel08_gaHI[NTHREADS][MAX_NCELLS];
    double dcs_gloverabel08_gaHI[NTHREADS][MAX_NCELLS];
    double c_gloverabel08_gaHp[1024];
    double cs_gloverabel08_gaHp[NTHREADS][MAX_NCELLS];
    double dcs_gloverabel08_gaHp[NTHREADS][MAX_NCELLS];
    double c_gloverabel08_gphdl[1024];
    double cs_gloverabel08_gphdl[NTHREADS][MAX_NCELLS];
    double dcs_gloverabel08_gphdl[NTHREADS][MAX_NCELLS];
    double c_gloverabel08_gpldl[1024];
    double cs_gloverabel08_gpldl[NTHREADS][MAX_NCELLS];
    double dcs_gloverabel08_gpldl[NTHREADS][MAX_NCELLS];
    double c_gloverabel08_h2lte[1024];
    double cs_gloverabel08_h2lte[NTHREADS][MAX_NCELLS];
    double dcs_gloverabel08_h2lte[NTHREADS][MAX_NCELLS];

    double c_h2formation_h2mcool[1024];
    double cs_h2formation_h2mcool[NTHREADS][MAX_NCELLS];
    double dcs_h2formation_h2mcool[NTHREADS][MAX_NCELLS];
    double c_h2formation_h2mheat[1024];
    double cs_h2formation_h2mheat[NTHREADS][MAX_NCELLS];
    double dcs_h2formation_h2mheat[NTHREADS][MAX_NCELLS];
    double c_h2formation_ncrd1[1024];
    double cs_h2formation_ncrd1[NTHREADS][MAX_NCELLS];
    double dcs_h2formation_ncrd1[NTHREADS][MAX_NCELLS];
    double c_h2formation_ncrd2[1024];
    double cs_h2formation_ncrd2[NTHREADS][MAX_NCELLS];
    double dcs_h2formation_ncrd2[NTHREADS][MAX_NCELLS];
    double c_h2formation_ncrn[1024];
    double cs_h2formation_ncrn[NTHREADS][MAX_NCELLS];
    double dcs_h2formation_ncrn[NTHREADS][MAX_NCELLS];

    double c_reHeII1_reHeII1[1024];
    double cs_reHeII1_reHeII1[NTHREADS][MAX_NCELLS];
    double dcs_reHeII1_reHeII1[NTHREADS][MAX_NCELLS];

    double c_reHeII2_reHeII2[1024];
    double cs_reHeII2_reHeII2[NTHREADS][MAX_NCELLS];
    double dcs_reHeII2_reHeII2[NTHREADS][MAX_NCELLS];

    double c_reHeIII_reHeIII[1024];
    double cs_reHeIII_reHeIII[NTHREADS][MAX_NCELLS];
    double dcs_reHeIII_reHeIII[NTHREADS][MAX_NCELLS];

    double c_reHII_reHII[1024];
    double cs_reHII_reHII[NTHREADS][MAX_NCELLS];
    double dcs_reHII_reHII[NTHREADS][MAX_NCELLS];

    int bin_id[NTHREADS][MAX_NCELLS];
    int ncells;



    // gamma as a function of temperature
    double g_gammaH2_1[1024];
    double g_dgammaH2_1_dT[1024];

    // store the gamma for that particular step
    double gammaH2_1[NTHREADS][MAX_NCELLS];
    double dgammaH2_1_dT[NTHREADS][MAX_NCELLS];

    // store 1 / (gamma - 1)
    double _gammaH2_1_m1[NTHREADS][MAX_NCELLS];

    double g_gammaH2_2[1024];
    double g_dgammaH2_2_dT[1024];

    // store the gamma for that particular step
    double gammaH2_2[NTHREADS][MAX_NCELLS];
    double dgammaH2_2_dT[NTHREADS][MAX_NCELLS];

    // store 1 / (gamma - 1)
    double _gammaH2_2_m1[NTHREADS][MAX_NCELLS];


    double scale[NTHREADS][10 * MAX_NCELLS ];
    double inv_scale[NTHREADS][10 * MAX_NCELLS];


    int nstrip;
    double mdensity[NTHREADS][MAX_NCELLS];
    double inv_mdensity[NTHREADS][MAX_NCELLS];


    double cie_optical_depth_approx[NTHREADS][MAX_NCELLS];


    double h2_optical_depth_approx[NTHREADS][MAX_NCELLS];


    const char *dengo_data_file;

    double reltol;
    double floor_value;
} primordial_data;


/* Declare ctype RHS and Jacobian */
typedef int(*rhs_f)( realtype, N_Vector , N_Vector , void * );
#ifndef CVSPILS
typedef int(*jac_f)( realtype, N_Vector  , N_Vector , SUNMatrix , void *, N_Vector, N_Vector, N_Vector);
#endif
#ifdef CVSPILS
typedef int(*jac_f)(N_Vector , N_Vector , realtype,
             N_Vector, N_Vector,
             void *user_data, N_Vector);
#endif


void *setup_cvode_solver( rhs_f f, jac_f Jac,  int NEQ,
        primordial_data *data, SUNLinearSolver LS, SUNMatrix A, N_Vector y, double reltol, N_Vector abstol);

int cvode_solver( void *cvode_mem, double *output, int NEQ, double *dt, primordial_data * data, N_Vector y, double reltol, N_Vector abstol );


primordial_data *primordial_setup_data( const char *FileLocation, int *NumberOfFields, char ***FieldNames);
void primordial_read_rate_tables(primordial_data*);
void primordial_read_cooling_tables(primordial_data*);
void primordial_read_gamma(primordial_data*);
void primordial_interpolate_gamma(primordial_data*, int );

void setting_up_extra_variables( primordial_data * data, double * input, int nstrip );

int dengo_evolve_primordial (double dtf, double &dt, double z,
                                     double *input, double *rtol,
                                     double *atol, unsigned long dims,
                                     primordial_data *data, double *temp);

double evolve_in_batches( void * cvode_mem, N_Vector y_vec, N_Vector abstol,
                          double reltol,double *input, int v_size, int d, int start_idx,
                          int MAX_ITERATION, double dtf, primordial_data *data );







#ifndef CVSPILS
#ifdef  CVKLU
int calculate_sparse_jacobian_primordial( realtype t,
                                        N_Vector y, N_Vector fy,
                                        SUNMatrix J, void *user_data,
                                        N_Vector tmp1, N_Vector tmp2,
                                        N_Vector tmp3);
#else
int calculate_jacobian_primordial( realtype t,
               N_Vector y, N_Vector fy, SUNMatrix J, void *user_data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
#endif
#endif

#ifdef CVSPILS
int calculate_JacTimesVec_primordial
            (N_Vector v, N_Vector Jv, realtype t,
             N_Vector y, N_Vector fy,
             void *user_data, N_Vector tmp);
#endif

int calculate_rhs_primordial(realtype t, N_Vector y, N_Vector ydot, void *user_data);
void ensure_electron_consistency(double *input, double *equil_array, unsigned long nstrip, int nchem);
void temperature_from_mass_density(double *input, int nstrip, int nchem,
                                   double *strip_temperature);

int primordial_calculate_temperature(primordial_data *data, double *input, int nstrip, int nchem);




typedef struct code_units
{

  int comoving_coordinates;
  double density_units;
  double length_units;
  double time_units;
  double velocity_units;
  double a_units;
  double a_value;

} code_units;


typedef struct dengo_field_data
{

  unsigned long int nstrip;
  unsigned long int ncells;
  // let's just pass them passively through field_data
  double reltol;
  double floor_value;
  // This should be updated dynamically
  // with dengo
  double *density;
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

  double *CoolingTime;
  double *MolecularWeight;
  double *temperature;
  double *Gamma;

  int *grid_start;
  int *grid_end;
  int *grid_dimension;

  const char *dengo_data_file;
  code_units *units;
} dengo_field_data;


// we can call this and pass reltol and floor_value to solve_chemistry
int primordial_solve_chemistry( code_units *units, dengo_field_data *field_data, double dt );
int primordial_solve_chemistry_dt( code_units *units, dengo_field_data *field_data, double* rtol, double *atol, double dt );

int dengo_estimate_cooling_time( code_units* units, dengo_field_data * field_data );

int primordial_calculate_cooling_timescale( double *cooling_time, double *input, int nstrip, primordial_data *data);

int dengo_calculate_temperature( code_units*, dengo_field_data* );
int dengo_calculate_gamma( double* gamma_eff, primordial_data*, double*, int );
int dengo_calculate_mean_molecular_weight( code_units*, dengo_field_data * );

int calculate_equilibrium_abundance( primordial_data *data, double *input,
	int nstrip, unsigned long d, unsigned long dims, double *equil_array);
// Enzo interface
int dengo_estimate_cooling_time_enzo( code_units* units, dengo_field_data *field_data );

int primordial_solve_chemistry_enzo( code_units *units, dengo_field_data *field_data, double dt );


int reshape_to_dengo_field_data( code_units* units, dengo_field_data *field_data, double* input );
int flatten_dengo_field_data( code_units* units, dengo_field_data *field_data, double *input );
int reshape_to_dengo_field_data_enzo( code_units* units, dengo_field_data *field_data, double* input, double *temp );
int flatten_dengo_field_data_enzo( code_units* units, dengo_field_data *field_data, double *input );


void ensure_species_conservation(double *input, double *mdensity, double *equil_array, primordial_data *data, int nstrip, unsigned long d, unsigned long dims, int nchem);
