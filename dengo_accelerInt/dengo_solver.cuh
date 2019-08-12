#ifndef DENGO_SOLVER_H
#define DENGO_SOLVER_H

#include "omp.h"
#include "time.h"
#include "sys/time.h"
#include "stdlib.h"
#include "math.h"
#include "hdf5.h"
#include "hdf5_hl.h"
#include "stdio.h"
#include "string.h" 

// Defining constants
#define RDATA_NBINS  (1024)    // number of bins beween the bounds
#define LOG_RDATA_LB (0.0)         // min Temp = 1 K
#define LOG_RDATA_UB (11.5129255)  // max Temp = 1.0e5 K

#define RDATA_DBIN  ( (LOG_RDATA_UB - LOG_RDATA_LB)/ (RDATA_NBINS-1) )
#define RDATA_IDBIN ( 1.0 / RDATA_DBIN )

#define NRATE (23)
#define NCOOL (26)

struct cuArray_ratedata {
  /* Now we do all of our cooling and chemical tables */
  cudaArray *r_k01;
  cudaArray *r_k02;
  cudaArray *r_k03;
  cudaArray *r_k04;
  cudaArray *r_k05;
  cudaArray *r_k06;
  cudaArray *r_k07;
  cudaArray *r_k08;
  cudaArray *r_k09;
  cudaArray *r_k10;
  cudaArray *r_k11;
  cudaArray *r_k12;
  cudaArray *r_k13;
  cudaArray *r_k14;
  cudaArray *r_k15;
  cudaArray *r_k16;
  cudaArray *r_k17;
  cudaArray *r_k18;
  cudaArray *r_k19;
  cudaArray *r_k20;
  cudaArray *r_k21;
  cudaArray *r_k22;
  cudaArray *r_k23;
  cudaArray *c_brem_brem;
  cudaArray *c_ceHI_ceHI;
  cudaArray *c_ceHeII_ceHeII;
  cudaArray *c_ceHeI_ceHeI;
  cudaArray *c_ciHI_ciHI;
  cudaArray *c_ciHeII_ciHeII;
  cudaArray *c_ciHeIS_ciHeIS;
  cudaArray *c_ciHeI_ciHeI;
  cudaArray *c_cie_cooling_cieco;
  cudaArray *c_cie_optical_depth_approx;
  cudaArray *c_compton_comp_;
  cudaArray *c_gammah_gammah;
  cudaArray *c_gloverabel08_gaH2;
  cudaArray *c_gloverabel08_gaHI;
  cudaArray *c_gloverabel08_gaHe;
  cudaArray *c_gloverabel08_gaHp;
  cudaArray *c_gloverabel08_gael;
  cudaArray *c_gloverabel08_h2lte;
  cudaArray *c_h2formation_h2mcool;
  cudaArray *c_h2formation_h2mheat;
  cudaArray *c_h2formation_ncrd1;
  cudaArray *c_h2formation_ncrd2;
  cudaArray *c_h2formation_ncrn;
  cudaArray *c_reHII_reHII;
  cudaArray *c_reHeII1_reHeII1;
  cudaArray *c_reHeII2_reHeII2;
  cudaArray *c_reHeIII_reHeIII;  
  
  // gamma as a function of temperature
  cudaArray *g_gammaH2_1;
  cudaArray *g_dgammaH2_1dT;
    
  cudaArray *g_gammaH2_2;
  float *g_dgammaH2_2_dT;


  const char *dengo_data_file;
};

struct cvklu_data {
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

    /* Now we do all of our cooling and chemical tables */
    float r_k01[1024];
    float r_k02[1024];
    float r_k03[1024];
    float r_k04[1024];
    float r_k05[1024];
    float r_k06[1024];
    float r_k07[1024];
    float r_k08[1024];
    float r_k09[1024];
    float r_k10[1024];
    float r_k11[1024];
    float r_k12[1024];
    float r_k13[1024];
    float r_k14[1024];
    float r_k15[1024];
    float r_k16[1024];
    float r_k17[1024];
    float r_k18[1024];
    float r_k19[1024];
    float r_k21[1024];
    float r_k22[1024];
    float c_brem_brem[1024];
    float c_ceHeI_ceHeI[1024];
    float c_ceHeII_ceHeII[1024];
    float c_ceHI_ceHI[1024];
    float c_cie_cooling_cieco[1024];
    float c_ciHeI_ciHeI[1024];
    float c_ciHeII_ciHeII[1024];
    float c_ciHeIS_ciHeIS[1024];
    float c_ciHI_ciHI[1024];
    float c_compton_comp_[1024];
    float c_gammah_gammah[1024];
    float c_gloverabel08_gael[1024];
    float c_gloverabel08_gaH2[1024];
    float c_gloverabel08_gaHe[1024];
    float c_gloverabel08_gaHI[1024];
    float c_gloverabel08_gaHp[1024];
    float c_gloverabel08_gphdl[1024];
    float c_gloverabel08_gpldl[1024];
    float c_gloverabel08_h2lte[1024];
    
    float c_h2formation_h2mcool[1024];
    float c_h2formation_h2mheat[1024];
    float c_h2formation_ncrd1[1024];
    float c_h2formation_ncrd2[1024];
    float c_h2formation_ncrn[1024];
    
    float c_reHeII1_reHeII1[1024];
    float c_reHeII2_reHeII2[1024];
    float c_reHeIII_reHeIII[1024];
    float c_reHII_reHII[1024];
    
    int ncells;

    // gamma as a function of temperature
    float g_gammaH2_1[1024];
    float g_dgammaH2_1_dT[1024];

    float g_gammaH2_2[1024];
    float g_dgammaH2_2_dT[1024];

    int nstrip;
    const char *dengo_data_file;
};



__host__ void cvklu_read_gamma(cvklu_data *data);
__host__ void cvklu_read_cooling_tables(cvklu_data *data);
__host__ void cvklu_read_rate_tables(cvklu_data *data);
__host__ cvklu_data *cvklu_setup_data( const char *FileLocation, int *NumberOfFields, char ***FieldNames);
__host__ void dengo_set_initial_conditions( double density, double T0, double fH2, int NUM, double **y_host, double** var_host );
__host__ void dengo_set_additional_constant( double denisty, double temperature, int NUM, double **y_host, double **temperature_array, double **density_array, double **h2_optical_depth_approx );

__host__ void bind_rate_to_texture( cvklu_data *data );
__device__ void get_normalized_T( float *Tdef, float temp_out );
__global__ void interpolate_rates_texture( double *reaction_rate_out, float *temp_out );
__global__ void interpolate_rates_double( double *reaction_rates_out, float * temp_out, cvklu_data *rate_data);


#endif


