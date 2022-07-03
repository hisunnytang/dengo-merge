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
#include <omp.h>

#define NTHREADS 40
#define MAX_NCELLS 256
#define NSPECIES 5
#define DMAX(A,B) ((A) > (B) ? (A) : (B))
#define DMIN(A,B) ((A) < (B) ? (A) : (B))

 

int predator_prey_main(int argc, char **argv);



typedef struct predator_prey_data {
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
    double r_exp_growth_prey[1024];
    double rs_exp_growth_prey[NTHREADS][MAX_NCELLS];
    double drs_exp_growth_prey[NTHREADS][MAX_NCELLS];
    
    double r_natural_death_predator[1024];
    double rs_natural_death_predator[NTHREADS][MAX_NCELLS];
    double drs_natural_death_predator[NTHREADS][MAX_NCELLS];
    
    double r_predation[1024];
    double rs_predation[NTHREADS][MAX_NCELLS];
    double drs_predation[NTHREADS][MAX_NCELLS];
    
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
    

    double scale[NTHREADS][5 * MAX_NCELLS ];
    double inv_scale[NTHREADS][5 * MAX_NCELLS];
    

    int nstrip;
    double mdensity[NTHREADS][MAX_NCELLS];
    double inv_mdensity[NTHREADS][MAX_NCELLS];

    
    
    
    const char *dengo_data_file;
    
    double reltol;
    double floor_value;
} predator_prey_data;
predator_prey_data *predator_prey_setup_data(const char *, int *, char***);
void predator_prey_read_rate_tables(predator_prey_data*);
void predator_prey_read_cooling_tables(predator_prey_data*);
void predator_prey_read_gamma(predator_prey_data*);
int dengo_evolve_predator_prey (double dtf, double dt, double z,
                                     double *input, double *rtol,
                                     double *atol, long long dims,
                                     predator_prey_data *data);
 

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
  double *dead_predator_density;
  double *dead_prey_density;
  double *ge_density;
  double *predator_density;
  double *prey_density;
    
  double *CoolingTime;
  double *MolecularWeight;
  double *temperature;
  double *Gamma;
  double *Pressure;

  int *grid_start;
  int *grid_end;
  int *grid_dimension;

  const char *dengo_data_file;
  code_units *units;
} dengo_field_data;

// Enzo interface
int dengo_estimate_cooling_time_enzo( code_units* units, dengo_field_data *field_data );

int predator_prey_solve_chemistry_enzo( code_units *units, dengo_field_data *field_data, double dt );
int reshape_to_dengo_field_data_enzo( code_units* units, dengo_field_data *field_data, double* input, double *temp );
int flatten_dengo_field_data_enzo( code_units* units, dengo_field_data *field_data, double *input );

typedef int(*rhs_f)(double *, double *, int, int, void *);
typedef int(*jac_f)(double *, double *, int, int, void *);
int BE_chem_solve(rhs_f f, jac_f J, double *u, double dt, double *rtol, 
                  double *atol, int nstrip, int nchem, double *scaling, void *sdata,
                  double *, double *, double *, double *);



int calculate_jacobian_predator_prey(double *input, double *Joutput,
        int nstrip, int nchem, void *sdata);
int calculate_rhs_predator_prey(double *input, double *rhs, int nstrip,
                  int nchem, void *sdata);
void ensure_electron_consistency(double *input, long long nstrip, int nchem);
void temperature_from_mass_density(double *input, int nstrip, int nchem, 
                                   double *strip_temperature);
void setting_up_extra_variables( predator_prey_data * data, double * input, int nstrip );
    
 


int predator_prey_calculate_cooling_timescale( double *cooling_time, double *input, int nstrip, predator_prey_data *data);

int dengo_estimate_cooling_time_enzo( code_units* units, dengo_field_data *field_data );
