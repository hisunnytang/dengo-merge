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
    double r_k01[1024];
    double r_k02[1024];
    double r_k03[1024];
    double r_k04[1024];
    double r_k05[1024];
    double r_k06[1024];
    double r_k07[1024];
    double r_k08[1024];
    double r_k09[1024];
    double r_k10[1024];
    double r_k11[1024];
    double r_k12[1024];
    double r_k13[1024];
    double r_k14[1024];
    double r_k15[1024];
    double r_k16[1024];
    double r_k17[1024];
    double r_k18[1024];
    double r_k19[1024];
    double r_k21[1024];
    double r_k22[1024];
    double c_brem_brem[1024];
    double c_ceHeI_ceHeI[1024];
    double c_ceHeII_ceHeII[1024];
    double c_ceHI_ceHI[1024];
    double c_cie_cooling_cieco[1024];
    double c_ciHeI_ciHeI[1024];
    double c_ciHeII_ciHeII[1024];
    double c_ciHeIS_ciHeIS[1024];
    double c_ciHI_ciHI[1024];
    double c_compton_comp_[1024];
    double c_gammah_gammah[1024];
    double c_gloverabel08_gael[1024];
    double c_gloverabel08_gaH2[1024];
    double c_gloverabel08_gaHe[1024];
    double c_gloverabel08_gaHI[1024];
    double c_gloverabel08_gaHp[1024];
    double c_gloverabel08_gphdl[1024];
    double c_gloverabel08_gpldl[1024];
    double c_gloverabel08_h2lte[1024];
    
    double c_h2formation_h2mcool[1024];
    double c_h2formation_h2mheat[1024];
    double c_h2formation_ncrd1[1024];
    double c_h2formation_ncrd2[1024];
    double c_h2formation_ncrn[1024];
    
    double c_reHeII1_reHeII1[1024];
    double c_reHeII2_reHeII2[1024];
    double c_reHeIII_reHeIII[1024];
    double c_reHII_reHII[1024];
    
    int ncells;

    // gamma as a function of temperature
    double g_gammaH2_1[1024];
    double g_dgammaH2_1_dT[1024];

    double g_gammaH2_2[1024];
    double g_dgammaH2_2_dT[1024];

    int nstrip;
    const char *dengo_data_file;
};



void cvklu_read_gamma(cvklu_data *data);
void cvklu_read_cooling_tables(cvklu_data *data);
void cvklu_read_rate_tables(cvklu_data *data);
cvklu_data *cvklu_setup_data( const char *FileLocation, int *NumberOfFields, char ***FieldNames);
void dengo_set_initial_conditions( double density, double T0, double fH2, int NUM, double **y_host, double** var_host );
void dengo_set_additional_constant( double denisty, double temperature, int NUM, double **y_host, double **temperature_array, double **density_array, double **h2_optical_depth_approx );

#endif


