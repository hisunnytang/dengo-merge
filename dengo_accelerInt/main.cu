#include <cuda.h>

/* Include common code. */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>
#include <complex.h>
#include <sys/types.h>
#include <sys/stat.h>
#ifdef DEBUG
#include <fenv.h>
#endif

/* Include CUDA libraries. */
#include <cuda.h>
#include <cuda_runtime.h>
#include <helper_cuda.h>
#include <cuComplex.h>
#include "dengo_solver.cuh"
#include "mechanism.cuh"

texture<float, 1, cudaReadModeElementType> rk01_tex;
texture<float, 1, cudaReadModeElementType> rk02_tex;
texture<float, 1, cudaReadModeElementType> rk03_tex;
texture<float, 1, cudaReadModeElementType> rk04_tex;
texture<float, 1, cudaReadModeElementType> rk05_tex;
texture<float, 1, cudaReadModeElementType> rk06_tex;
texture<float, 1, cudaReadModeElementType> rk07_tex;
texture<float, 1, cudaReadModeElementType> rk08_tex;
texture<float, 1, cudaReadModeElementType> rk09_tex;
texture<float, 1, cudaReadModeElementType> rk10_tex;
texture<float, 1, cudaReadModeElementType> rk11_tex;
texture<float, 1, cudaReadModeElementType> rk12_tex;
texture<float, 1, cudaReadModeElementType> rk13_tex;
texture<float, 1, cudaReadModeElementType> rk14_tex;
texture<float, 1, cudaReadModeElementType> rk15_tex;
texture<float, 1, cudaReadModeElementType> rk16_tex;
texture<float, 1, cudaReadModeElementType> rk17_tex;
texture<float, 1, cudaReadModeElementType> rk18_tex;
texture<float, 1, cudaReadModeElementType> rk19_tex;
texture<float, 1, cudaReadModeElementType> rk20_tex;
texture<float, 1, cudaReadModeElementType> rk21_tex;
texture<float, 1, cudaReadModeElementType> rk22_tex;
texture<float, 1, cudaReadModeElementType> rk23_tex;

// https://stackoverflow.com/questions/37080835/how-to-bind-a-float-array-to-a-1d-texture-in-cuda




__global__ void interpolate_rates_texture( double *reaction_rate_out, float *temp_out )
{
  int tid = threadIdx.x + blockDim.x * blockIdx.x;

  // norm_T runs between
  float norm_T = (log(temp_out[tid]) - LOG_RDATA_LB) / RDATA_DBIN / 1024.0;
  reaction_rate_out[INDEX( 0)] = tex1D<float>(rk01_tex, norm_T);
  reaction_rate_out[INDEX( 1)] = tex1D<float>(rk02_tex, norm_T);
  reaction_rate_out[INDEX( 2)] = tex1D<float>(rk03_tex, norm_T);
  reaction_rate_out[INDEX( 3)] = tex1D<float>(rk04_tex, norm_T);
  reaction_rate_out[INDEX( 4)] = tex1D<float>(rk05_tex, norm_T);
  reaction_rate_out[INDEX( 5)] = tex1D<float>(rk06_tex, norm_T);
  reaction_rate_out[INDEX( 6)] = tex1D<float>(rk07_tex, norm_T);
  reaction_rate_out[INDEX( 7)] = tex1D<float>(rk08_tex, norm_T);
  reaction_rate_out[INDEX( 8)] = tex1D<float>(rk09_tex, norm_T);
  reaction_rate_out[INDEX( 9)] = tex1D<float>(rk10_tex, norm_T);
  reaction_rate_out[INDEX(10)] = tex1D<float>(rk11_tex, norm_T);
  reaction_rate_out[INDEX(11)] = tex1D<float>(rk12_tex, norm_T);
  reaction_rate_out[INDEX(12)] = tex1D<float>(rk13_tex, norm_T);
  reaction_rate_out[INDEX(13)] = tex1D<float>(rk14_tex, norm_T);
  reaction_rate_out[INDEX(14)] = tex1D<float>(rk15_tex, norm_T);
  reaction_rate_out[INDEX(15)] = tex1D<float>(rk16_tex, norm_T);
  reaction_rate_out[INDEX(16)] = tex1D<float>(rk17_tex, norm_T);
  reaction_rate_out[INDEX(17)] = tex1D<float>(rk18_tex, norm_T);
  reaction_rate_out[INDEX(18)] = tex1D<float>(rk19_tex, norm_T);
//  reaction_rate_out[INDEX(19)] = tex1D<float>(rk20_tex, norm_T);
  reaction_rate_out[INDEX(20)] = tex1D<float>(rk21_tex, norm_T);
  reaction_rate_out[INDEX(21)] = tex1D<float>(rk22_tex, norm_T);
//  reaction_rate_out[INDEX(22)] = tex1D<float>(rk23_tex, norm_T);
  if (tid == 0){
    for (int i = 0; i< 21; i++){
      printf("norm_T = %0.5g, rate[INDEX(%d)]= %0.5g\n", norm_T, i, reaction_rate_out[INDEX(i)]);
    }
}

}

__global__ void interpolate_rates_double( double *reaction_rates_out, float *temp_out, cvklu_data *rate_data)
{
    
    int tid, bin_id, zbin_id;
    double t1, t2;
    double Tdef, dT, invTs, log_temp_out;
    int no_photo = 0;
    double lb = log(rate_data->bounds[0]);

    tid = threadIdx.x + blockDim.x * blockIdx.x;
 
    log_temp_out = log(temp_out[tid]);
    bin_id = (int) ( rate_data->idbin * ( log_temp_out -  lb ) );
    if ( bin_id <= 0) {
        bin_id = 0;
    } else if ( bin_id >= rate_data->nbins) {
        bin_id = rate_data->nbins - 1;
    }
    
    //printf( "bin_id = %d; temp_out = %0.5g \n", bin_id, temp_out[tid]);
    t1 = (lb + (bin_id    ) * rate_data->dbin);
    t2 = (lb + (bin_id + 1) * rate_data->dbin);
    Tdef = (log_temp_out - t1)/(t2 - t1);
    dT    = t2 - t1;
    invTs = 1.0 / temp_out[tid];

 
    // rate_out is a long 1D array
    // NRATE is the number of rate required by the solver network

    reaction_rates_out[INDEX( 0)] = (double) (rate_data->r_k01[bin_id] + Tdef * (rate_data->r_k01[bin_id+1] - rate_data->r_k01[bin_id]) );
    reaction_rates_out[INDEX( 1)] = (double) (rate_data->r_k02[bin_id] + Tdef * (rate_data->r_k02[bin_id+1] - rate_data->r_k02[bin_id]) );
    reaction_rates_out[INDEX( 2)] = (double) (rate_data->r_k03[bin_id] + Tdef * (rate_data->r_k03[bin_id+1] - rate_data->r_k03[bin_id]) );
    reaction_rates_out[INDEX( 3)] = (double) (rate_data->r_k04[bin_id] + Tdef * (rate_data->r_k04[bin_id+1] - rate_data->r_k04[bin_id]) );
    reaction_rates_out[INDEX( 4)] = (double) (rate_data->r_k05[bin_id] + Tdef * (rate_data->r_k05[bin_id+1] - rate_data->r_k05[bin_id]) );
    reaction_rates_out[INDEX( 5)] = (double) (rate_data->r_k06[bin_id] + Tdef * (rate_data->r_k06[bin_id+1] - rate_data->r_k06[bin_id]) );
    reaction_rates_out[INDEX( 6)] = (double) (rate_data->r_k07[bin_id] + Tdef * (rate_data->r_k07[bin_id+1] - rate_data->r_k07[bin_id]) );
    reaction_rates_out[INDEX( 7)] = (double) (rate_data->r_k08[bin_id] + Tdef * (rate_data->r_k08[bin_id+1] - rate_data->r_k08[bin_id]) );
    reaction_rates_out[INDEX( 8)] = (double) (rate_data->r_k09[bin_id] + Tdef * (rate_data->r_k09[bin_id+1] - rate_data->r_k09[bin_id]) );
    reaction_rates_out[INDEX( 9)] = (double) (rate_data->r_k10[bin_id] + Tdef * (rate_data->r_k10[bin_id+1] - rate_data->r_k10[bin_id]) );
    reaction_rates_out[INDEX(10)] = (double) (rate_data->r_k11[bin_id] + Tdef * (rate_data->r_k11[bin_id+1] - rate_data->r_k11[bin_id]) );
    reaction_rates_out[INDEX(11)] = (double) (rate_data->r_k12[bin_id] + Tdef * (rate_data->r_k12[bin_id+1] - rate_data->r_k12[bin_id]) );
    reaction_rates_out[INDEX(12)] = (double) (rate_data->r_k13[bin_id] + Tdef * (rate_data->r_k13[bin_id+1] - rate_data->r_k13[bin_id]) );
    reaction_rates_out[INDEX(13)] = (double) (rate_data->r_k14[bin_id] + Tdef * (rate_data->r_k14[bin_id+1] - rate_data->r_k14[bin_id]) );
    reaction_rates_out[INDEX(14)] = (double) (rate_data->r_k15[bin_id] + Tdef * (rate_data->r_k15[bin_id+1] - rate_data->r_k15[bin_id]) );
    reaction_rates_out[INDEX(15)] = (double) (rate_data->r_k16[bin_id] + Tdef * (rate_data->r_k16[bin_id+1] - rate_data->r_k16[bin_id]) );
    reaction_rates_out[INDEX(16)] = (double) (rate_data->r_k17[bin_id] + Tdef * (rate_data->r_k17[bin_id+1] - rate_data->r_k17[bin_id]) );
    reaction_rates_out[INDEX(17)] = (double) (rate_data->r_k18[bin_id] + Tdef * (rate_data->r_k18[bin_id+1] - rate_data->r_k18[bin_id]) );
    reaction_rates_out[INDEX(18)] = (double) (rate_data->r_k19[bin_id] + Tdef * (rate_data->r_k19[bin_id+1] - rate_data->r_k19[bin_id]) );
//    reaction_rates_out[INDEX(19)] = (double) (rate_data->r_k20[bin_id] + Tdef * (rate_data->r_k20[bin_id+1] - rate_data->r_k20[bin_id]) );
    reaction_rates_out[INDEX(20)] = (double) (rate_data->r_k21[bin_id] + Tdef * (rate_data->r_k21[bin_id+1] - rate_data->r_k21[bin_id]) );
    reaction_rates_out[INDEX(21)] = (double) (rate_data->r_k22[bin_id] + Tdef * (rate_data->r_k22[bin_id+1] - rate_data->r_k22[bin_id]) );
//    reaction_rates_out[INDEX(22)] = (double) (rate_data->r_k23[bin_id] + Tdef * (rate_data->r_k23[bin_id+1] - rate_data->r_k23[bin_id]) );
  if (tid == 0){
    for (int i = 0; i< 21; i++){
      printf("Tdef = %0.5g; rate[INDEX(%d)]= %0.5g\n",Tdef,  i, reaction_rates_out[INDEX(i)]);
    }
  }
}




void cvklu_read_rate_tables(cvklu_data *data)
{
    const char * filedir;
    if (data->dengo_data_file != NULL){
        filedir =  data->dengo_data_file; 
    } else{
        filedir = "cvklu_tables.h5";   
    }

    hid_t file_id = H5Fopen( filedir , H5F_ACC_RDONLY, H5P_DEFAULT);
    /* Allocate the correct number of rate tables */
    H5LTread_dataset_float(file_id, "/k01", data->r_k01);
    H5LTread_dataset_float(file_id, "/k02", data->r_k02);
    H5LTread_dataset_float(file_id, "/k03", data->r_k03);
    H5LTread_dataset_float(file_id, "/k04", data->r_k04);
    H5LTread_dataset_float(file_id, "/k05", data->r_k05);
    H5LTread_dataset_float(file_id, "/k06", data->r_k06);
    H5LTread_dataset_float(file_id, "/k07", data->r_k07);
    H5LTread_dataset_float(file_id, "/k08", data->r_k08);
    H5LTread_dataset_float(file_id, "/k09", data->r_k09);
    H5LTread_dataset_float(file_id, "/k10", data->r_k10);
    H5LTread_dataset_float(file_id, "/k11", data->r_k11);
    H5LTread_dataset_float(file_id, "/k12", data->r_k12);
    H5LTread_dataset_float(file_id, "/k13", data->r_k13);
    H5LTread_dataset_float(file_id, "/k14", data->r_k14);
    H5LTread_dataset_float(file_id, "/k15", data->r_k15);
    H5LTread_dataset_float(file_id, "/k16", data->r_k16);
    H5LTread_dataset_float(file_id, "/k17", data->r_k17);
    H5LTread_dataset_float(file_id, "/k18", data->r_k18);
    H5LTread_dataset_float(file_id, "/k19", data->r_k19);
    H5LTread_dataset_float(file_id, "/k21", data->r_k21);
    H5LTread_dataset_float(file_id, "/k22", data->r_k22);
    
    H5Fclose(file_id);
}


void cvklu_read_cooling_tables(cvklu_data *data)
{

    const char * filedir;
    if (data->dengo_data_file != NULL){
        filedir =  data->dengo_data_file; 
    } else{
        filedir = "cvklu_tables.h5";   
    }
    hid_t file_id = H5Fopen( filedir , H5F_ACC_RDONLY, H5P_DEFAULT);
    /* Allocate the correct number of rate tables */
    H5LTread_dataset_float(file_id, "/brem_brem",
                            data->c_brem_brem);
    H5LTread_dataset_float(file_id, "/ceHeI_ceHeI",
                            data->c_ceHeI_ceHeI);
    H5LTread_dataset_float(file_id, "/ceHeII_ceHeII",
                            data->c_ceHeII_ceHeII);
    H5LTread_dataset_float(file_id, "/ceHI_ceHI",
                            data->c_ceHI_ceHI);
    H5LTread_dataset_float(file_id, "/cie_cooling_cieco",
                            data->c_cie_cooling_cieco);
    H5LTread_dataset_float(file_id, "/ciHeI_ciHeI",
                            data->c_ciHeI_ciHeI);
    H5LTread_dataset_float(file_id, "/ciHeII_ciHeII",
                            data->c_ciHeII_ciHeII);
    H5LTread_dataset_float(file_id, "/ciHeIS_ciHeIS",
                            data->c_ciHeIS_ciHeIS);
    H5LTread_dataset_float(file_id, "/ciHI_ciHI",
                            data->c_ciHI_ciHI);
    H5LTread_dataset_float(file_id, "/compton_comp_",
                            data->c_compton_comp_);
    H5LTread_dataset_float(file_id, "/gammah_gammah",
                            data->c_gammah_gammah);
    H5LTread_dataset_float(file_id, "/gloverabel08_gael",
                            data->c_gloverabel08_gael);
    H5LTread_dataset_float(file_id, "/gloverabel08_gaH2",
                            data->c_gloverabel08_gaH2);
    H5LTread_dataset_float(file_id, "/gloverabel08_gaHe",
                            data->c_gloverabel08_gaHe);
    H5LTread_dataset_float(file_id, "/gloverabel08_gaHI",
                            data->c_gloverabel08_gaHI);
    H5LTread_dataset_float(file_id, "/gloverabel08_gaHp",
                            data->c_gloverabel08_gaHp);
    H5LTread_dataset_float(file_id, "/gloverabel08_gphdl",
                            data->c_gloverabel08_gphdl);
    H5LTread_dataset_float(file_id, "/gloverabel08_gpldl",
                            data->c_gloverabel08_gpldl);
    H5LTread_dataset_float(file_id, "/gloverabel08_h2lte",
                            data->c_gloverabel08_h2lte);
    H5LTread_dataset_float(file_id, "/h2formation_h2mcool",
                            data->c_h2formation_h2mcool);
    H5LTread_dataset_float(file_id, "/h2formation_h2mheat",
                            data->c_h2formation_h2mheat);
    H5LTread_dataset_float(file_id, "/h2formation_ncrd1",
                            data->c_h2formation_ncrd1);
    H5LTread_dataset_float(file_id, "/h2formation_ncrd2",
                            data->c_h2formation_ncrd2);
    H5LTread_dataset_float(file_id, "/h2formation_ncrn",
                            data->c_h2formation_ncrn);
    H5LTread_dataset_float(file_id, "/reHeII1_reHeII1",
                            data->c_reHeII1_reHeII1);
    H5LTread_dataset_float(file_id, "/reHeII2_reHeII2",
                            data->c_reHeII2_reHeII2);
    H5LTread_dataset_float(file_id, "/reHeIII_reHeIII",
                            data->c_reHeIII_reHeIII);
    H5LTread_dataset_float(file_id, "/reHII_reHII",
                            data->c_reHII_reHII);

    H5Fclose(file_id);
}

void cvklu_read_gamma(cvklu_data *data)
{

    const char * filedir;
    if (data->dengo_data_file != NULL){
        filedir =  data->dengo_data_file; 
    } else{
        filedir = "cvklu_tables.h5";   
    }
    
    hid_t file_id = H5Fopen( filedir , H5F_ACC_RDONLY, H5P_DEFAULT);
    /* Allocate the correct number of rate tables */
    H5LTread_dataset_float(file_id, "/gammaH2_1",
                            data->g_gammaH2_1 );
    H5LTread_dataset_float(file_id, "/dgammaH2_1_dT",
                            data->g_dgammaH2_1_dT );   
    
    H5LTread_dataset_float(file_id, "/gammaH2_2",
                            data->g_gammaH2_2 );
    H5LTread_dataset_float(file_id, "/dgammaH2_2_dT",
                            data->g_dgammaH2_2_dT );   
    

    H5Fclose(file_id);

}


cvklu_data *cvklu_setup_data( const char *FileLocation, int *NumberOfFields, char ***FieldNames)
{

    //-----------------------------------------------------
    // Function : cvklu_setup_data
    // Description: Initialize a data object that stores the reaction/ cooling rate data 
    //-----------------------------------------------------

    int i, n;
    
    cvklu_data *data = (cvklu_data *) malloc(sizeof(cvklu_data));
    
    // point the module to look for cvklu_tables.h5
    data->dengo_data_file = FileLocation;

    /* allocate space for the scale related pieces */

    /* Temperature-related pieces */
    data->bounds[0] = 1.0;
    data->bounds[1] = 100000.0;
    data->nbins = 1024 - 1;
    data->dbin = (log(data->bounds[1]) - log(data->bounds[0])) / data->nbins;
    data->idbin = 1.0L / data->dbin;

    /* Redshift-related pieces */
    data->z_bounds[0] = 0.0;
    data->z_bounds[1] = 0.0;
    data->n_zbins = 0 - 1;
    data->d_zbin = (log(data->z_bounds[1] + 1.0) - log(data->z_bounds[0] + 1.0)) / data->n_zbins;
    data->id_zbin = 1.0L / data->d_zbin;
    
    cvklu_read_rate_tables(data);
    //fprintf(stderr, "Successfully read in rate tables.\n");

    cvklu_read_cooling_tables(data);
    //fprintf(stderr, "Successfully read in cooling rate tables.\n");
    
    cvklu_read_gamma(data);
    //fprintf(stderr, "Successfully read in gamma tables. \n");

    if (FieldNames != NULL && NumberOfFields != NULL) {
        NumberOfFields[0] = 10;
        FieldNames[0] = new char*[10];
        i = 0;
        
        FieldNames[0][i++] = strdup("H2_1");
        
        FieldNames[0][i++] = strdup("H2_2");
        
        FieldNames[0][i++] = strdup("H_1");
        
        FieldNames[0][i++] = strdup("H_2");
        
        FieldNames[0][i++] = strdup("H_m0");
        
        FieldNames[0][i++] = strdup("He_1");
        
        FieldNames[0][i++] = strdup("He_2");
        
        FieldNames[0][i++] = strdup("He_3");
        
        FieldNames[0][i++] = strdup("de");
        
        FieldNames[0][i++] = strdup("ge");
        
    }

    data->dengo_data_file = NULL;

    return data;

}

void dengo_set_initial_conditions( double density, double T0, double fH2, int NUM, double **y_host, double** var_host ){

/*
{'H2_1': array([10000.]),
 'H2_2': array([10000.]),
 'H_1': array([7.6e+13]),
 'H_2': array([10000.]),
 'H_m0': array([10000.]),
 'He_1': array([2.4e+13]),
 'He_2': array([10000.]),
 'He_3': array([10000.]),
 'T': array([3000.]),
 'de': array([12495.12442156]),
 'density': array([1.e+14]),
 'ge': array([3.02681398e+11])}
*/

    double mH = 1.67e-24;
    double k  = 1.3806488e-16;
    double tiny = 1.0e-20;

    (*y_host) = (double*)malloc(NUM * NSP * sizeof(double));
    (*var_host) = (double*)malloc(NUM * sizeof(double));
    //load temperature and mass fractions for all threads (cells)
    printf("NUM = %d; NSP = %d \n", NUM, NSP ); 
    
    
    double m_amu = 1.66053904e-24;
    density *= mH/ m_amu;

    int j = 1;

    for (int i = 0; i < NUM; ++i) {
        //loop through species
	j = 0;
    	// H2I
	(*y_host)[i + NUM * j] = 0.76 * fH2 * density / 2.0;
        j += 1;
        // H2II
	(*y_host)[i + NUM * j] = density * tiny / 2.0;
	j += 1;
        // HI
	(*y_host)[i + NUM * j] = 0.76 * (1.0 - fH2)* density / 1.00794 + density * tiny / 1.00794;
	j += 1;
        // HII
	(*y_host)[i + NUM * j] = density * tiny / 1.00794;
	j += 1;
        // H-
	(*y_host)[i + NUM * j] = density * tiny / 1.00794;
	j += 1;
        // HeI
	(*y_host)[i + NUM * j] = 0.24 * density / 4.002602;
	j += 1;
        // HeII
	(*y_host)[i + NUM * j] = tiny * density / 4.002602;
	j += 1;
        // HeIII
	(*y_host)[i + NUM * j] = tiny * density / 4.002602;
	j += 1;
        // electron (de)
	(*y_host)[i + NUM * j] = tiny * density;
	j += 1;
        // internal energy (ge)
	(*y_host)[i + NUM * j] = 3.0 / 2.0 * k * T0 / mH;
	j += 1;
    }

}



void dengo_set_additional_constant( double density, double temperature, int NUM, double **y_host, double **temperature_array, double **density_array, double **h2_optical_depth_approx ){

    double mH = 1.67e-24;
    double k  = 1.3806488e-16;
    double m_amu = 1.66053904e-24;
    density *= mH/ m_amu;


    (*temperature_array) = (double*)malloc(NUM * sizeof(double));
    (*density_array)     = (double*)malloc(NUM * sizeof(double));
    (*h2_optical_depth_approx) = (double*)malloc(NUM *sizeof(double));

    for (int i = 0; i < NUM; ++i) {
        (*temperature_array)[i] = temperature;
        (*density_array)    [i] = 1.0 * density * mH;
        (*h2_optical_depth_approx)[i] = fmin( 1.0, pow(( mH *density / (1.34e-14)), -0.45) ); 
    }


}


void initialize_cuArray_ratedata( cuArray_ratedata *H_data , cvklu_data * ratedata )
{
  // this describes the format of a texture element
  // for float texels we could create a channel with:
  cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();

  // allocate the device arrays on the host pointer
  cudaMallocArray( &( H_data->r_k01 ), &channelDesc, RDATA_NBINS);
  cudaMallocArray( &( H_data->r_k02 ), &channelDesc, RDATA_NBINS);
  cudaMallocArray( &( H_data->r_k03 ), &channelDesc, RDATA_NBINS);
  cudaMallocArray( &( H_data->r_k04 ), &channelDesc, RDATA_NBINS);
  cudaMallocArray( &( H_data->r_k05 ), &channelDesc, RDATA_NBINS);
  cudaMallocArray( &( H_data->r_k06 ), &channelDesc, RDATA_NBINS);
  cudaMallocArray( &( H_data->r_k07 ), &channelDesc, RDATA_NBINS);
  cudaMallocArray( &( H_data->r_k08 ), &channelDesc, RDATA_NBINS);
  cudaMallocArray( &( H_data->r_k09 ), &channelDesc, RDATA_NBINS);
  cudaMallocArray( &( H_data->r_k10 ), &channelDesc, RDATA_NBINS);
  cudaMallocArray( &( H_data->r_k11 ), &channelDesc, RDATA_NBINS);
  cudaMallocArray( &( H_data->r_k12 ), &channelDesc, RDATA_NBINS);
  cudaMallocArray( &( H_data->r_k13 ), &channelDesc, RDATA_NBINS);
  cudaMallocArray( &( H_data->r_k14 ), &channelDesc, RDATA_NBINS);
  cudaMallocArray( &( H_data->r_k15 ), &channelDesc, RDATA_NBINS);
  cudaMallocArray( &( H_data->r_k16 ), &channelDesc, RDATA_NBINS);
  cudaMallocArray( &( H_data->r_k17 ), &channelDesc, RDATA_NBINS);
  cudaMallocArray( &( H_data->r_k18 ), &channelDesc, RDATA_NBINS);
  cudaMallocArray( &( H_data->r_k19 ), &channelDesc, RDATA_NBINS);
  cudaMallocArray( &( H_data->r_k20 ), &channelDesc, RDATA_NBINS);
  cudaMallocArray( &( H_data->r_k21 ), &channelDesc, RDATA_NBINS);
  cudaMallocArray( &( H_data->r_k22 ), &channelDesc, RDATA_NBINS);
  cudaMallocArray( &( H_data->r_k23 ), &channelDesc, RDATA_NBINS);

  cudaMemcpyToArray( H_data->r_k01, 0, 0, ratedata->r_k01, sizeof(float)*RDATA_NBINS, cudaMemcpyHostToDevice);
  cudaMemcpyToArray( H_data->r_k02, 0, 0, ratedata->r_k02, sizeof(float)*RDATA_NBINS, cudaMemcpyHostToDevice);
  cudaMemcpyToArray( H_data->r_k03, 0, 0, ratedata->r_k03, sizeof(float)*RDATA_NBINS, cudaMemcpyHostToDevice);
  cudaMemcpyToArray( H_data->r_k04, 0, 0, ratedata->r_k04, sizeof(float)*RDATA_NBINS, cudaMemcpyHostToDevice);
  cudaMemcpyToArray( H_data->r_k05, 0, 0, ratedata->r_k05, sizeof(float)*RDATA_NBINS, cudaMemcpyHostToDevice);
  cudaMemcpyToArray( H_data->r_k06, 0, 0, ratedata->r_k06, sizeof(float)*RDATA_NBINS, cudaMemcpyHostToDevice);
  cudaMemcpyToArray( H_data->r_k07, 0, 0, ratedata->r_k07, sizeof(float)*RDATA_NBINS, cudaMemcpyHostToDevice);
  cudaMemcpyToArray( H_data->r_k08, 0, 0, ratedata->r_k08, sizeof(float)*RDATA_NBINS, cudaMemcpyHostToDevice);
  cudaMemcpyToArray( H_data->r_k09, 0, 0, ratedata->r_k09, sizeof(float)*RDATA_NBINS, cudaMemcpyHostToDevice);
  cudaMemcpyToArray( H_data->r_k10, 0, 0, ratedata->r_k10, sizeof(float)*RDATA_NBINS, cudaMemcpyHostToDevice);
  cudaMemcpyToArray( H_data->r_k11, 0, 0, ratedata->r_k11, sizeof(float)*RDATA_NBINS, cudaMemcpyHostToDevice);
  cudaMemcpyToArray( H_data->r_k12, 0, 0, ratedata->r_k12, sizeof(float)*RDATA_NBINS, cudaMemcpyHostToDevice);
  cudaMemcpyToArray( H_data->r_k13, 0, 0, ratedata->r_k13, sizeof(float)*RDATA_NBINS, cudaMemcpyHostToDevice);
  cudaMemcpyToArray( H_data->r_k14, 0, 0, ratedata->r_k14, sizeof(float)*RDATA_NBINS, cudaMemcpyHostToDevice);
  cudaMemcpyToArray( H_data->r_k15, 0, 0, ratedata->r_k15, sizeof(float)*RDATA_NBINS, cudaMemcpyHostToDevice);
  cudaMemcpyToArray( H_data->r_k16, 0, 0, ratedata->r_k16, sizeof(float)*RDATA_NBINS, cudaMemcpyHostToDevice);
  cudaMemcpyToArray( H_data->r_k17, 0, 0, ratedata->r_k17, sizeof(float)*RDATA_NBINS, cudaMemcpyHostToDevice);
  cudaMemcpyToArray( H_data->r_k18, 0, 0, ratedata->r_k18, sizeof(float)*RDATA_NBINS, cudaMemcpyHostToDevice);
  cudaMemcpyToArray( H_data->r_k19, 0, 0, ratedata->r_k19, sizeof(float)*RDATA_NBINS, cudaMemcpyHostToDevice);
//  cudaMemcpyToArray( H_data->r_k20, 0, 0, ratedata->r_k20, sizeof(float)*RDATA_NBINS, cudaMemcpyHostToDevice);
  cudaMemcpyToArray( H_data->r_k21, 0, 0, ratedata->r_k21, sizeof(float)*RDATA_NBINS, cudaMemcpyHostToDevice);
  cudaMemcpyToArray( H_data->r_k22, 0, 0, ratedata->r_k22, sizeof(float)*RDATA_NBINS, cudaMemcpyHostToDevice);
//  cudaMemcpyToArray( H_data->r_k23, 0, 0, ratedata->r_k23, sizeof(float)*RDATA_NBINS, cudaMemcpyHostToDevice);


}

void bind_ratedata_to_texture( cuArray_ratedata *H_data )
{
  cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
  // Set texture parameters
  rk01_tex.addressMode[0] = cudaAddressModeClamp;
  rk01_tex.addressMode[1] = cudaAddressModeClamp;
  rk01_tex.filterMode = cudaFilterModeLinear;
  rk01_tex.normalized = true;
  cudaBindTextureToArray(rk01_tex, H_data->r_k01, channelDesc);
 
  rk02_tex.addressMode[0] = cudaAddressModeClamp;
  rk02_tex.addressMode[1] = cudaAddressModeClamp;
  rk02_tex.filterMode = cudaFilterModeLinear;
  rk02_tex.normalized = true;
  cudaBindTextureToArray(rk02_tex, H_data->r_k02, channelDesc);
 
  rk03_tex.addressMode[0] = cudaAddressModeClamp;
  rk03_tex.addressMode[1] = cudaAddressModeClamp;
  rk03_tex.filterMode = cudaFilterModeLinear;
  rk03_tex.normalized = true;
  cudaBindTextureToArray(rk03_tex, H_data->r_k03, channelDesc);
 
  rk04_tex.addressMode[0] = cudaAddressModeClamp;
  rk04_tex.addressMode[1] = cudaAddressModeClamp;
  rk04_tex.filterMode = cudaFilterModeLinear;
  rk04_tex.normalized = true;
  cudaBindTextureToArray(rk04_tex, H_data->r_k04, channelDesc);
 
  rk05_tex.addressMode[0] = cudaAddressModeClamp;
  rk05_tex.addressMode[1] = cudaAddressModeClamp;
  rk05_tex.filterMode = cudaFilterModeLinear;
  rk05_tex.normalized = true;
  cudaBindTextureToArray(rk05_tex, H_data->r_k05, channelDesc);
 
  rk06_tex.addressMode[0] = cudaAddressModeClamp;
  rk06_tex.addressMode[1] = cudaAddressModeClamp;
  rk06_tex.filterMode = cudaFilterModeLinear;
  rk06_tex.normalized = true;
  cudaBindTextureToArray(rk06_tex, H_data->r_k06, channelDesc);
 
  rk07_tex.addressMode[0] = cudaAddressModeClamp;
  rk07_tex.addressMode[1] = cudaAddressModeClamp;
  rk07_tex.filterMode = cudaFilterModeLinear;
  rk07_tex.normalized = true;
  cudaBindTextureToArray(rk07_tex, H_data->r_k07, channelDesc);
 
  rk08_tex.addressMode[0] = cudaAddressModeClamp;
  rk08_tex.addressMode[1] = cudaAddressModeClamp;
  rk08_tex.filterMode = cudaFilterModeLinear;
  rk08_tex.normalized = true;
  cudaBindTextureToArray(rk08_tex, H_data->r_k08, channelDesc);
 
  rk09_tex.addressMode[0] = cudaAddressModeClamp;
  rk09_tex.addressMode[1] = cudaAddressModeClamp;
  rk09_tex.filterMode = cudaFilterModeLinear;
  rk09_tex.normalized = true;
  cudaBindTextureToArray(rk09_tex, H_data->r_k09, channelDesc);
 
  rk10_tex.addressMode[0] = cudaAddressModeClamp;
  rk10_tex.addressMode[1] = cudaAddressModeClamp;
  rk10_tex.filterMode = cudaFilterModeLinear;
  rk10_tex.normalized = true;
  cudaBindTextureToArray(rk10_tex, H_data->r_k10, channelDesc);
 
  rk11_tex.addressMode[0] = cudaAddressModeClamp;
  rk11_tex.addressMode[1] = cudaAddressModeClamp;
  rk11_tex.filterMode = cudaFilterModeLinear;
  rk11_tex.normalized = true;
  cudaBindTextureToArray(rk11_tex, H_data->r_k11, channelDesc);
 
  rk12_tex.addressMode[0] = cudaAddressModeClamp;
  rk12_tex.addressMode[1] = cudaAddressModeClamp;
  rk12_tex.filterMode = cudaFilterModeLinear;
  rk12_tex.normalized = true;
  cudaBindTextureToArray(rk12_tex, H_data->r_k12, channelDesc);
 
  rk13_tex.addressMode[0] = cudaAddressModeClamp;
  rk13_tex.addressMode[1] = cudaAddressModeClamp;
  rk13_tex.filterMode = cudaFilterModeLinear;
  rk13_tex.normalized = true;
  cudaBindTextureToArray(rk13_tex, H_data->r_k13, channelDesc);
 
  rk14_tex.addressMode[0] = cudaAddressModeClamp;
  rk14_tex.addressMode[1] = cudaAddressModeClamp;
  rk14_tex.filterMode = cudaFilterModeLinear;
  rk14_tex.normalized = true;
  cudaBindTextureToArray(rk14_tex, H_data->r_k14, channelDesc);
 
  rk15_tex.addressMode[0] = cudaAddressModeClamp;
  rk15_tex.addressMode[1] = cudaAddressModeClamp;
  rk15_tex.filterMode = cudaFilterModeLinear;
  rk15_tex.normalized = true;
  cudaBindTextureToArray(rk15_tex, H_data->r_k15, channelDesc);
 
  rk16_tex.addressMode[0] = cudaAddressModeClamp;
  rk16_tex.addressMode[1] = cudaAddressModeClamp;
  rk16_tex.filterMode = cudaFilterModeLinear;
  rk16_tex.normalized = true;
  cudaBindTextureToArray(rk16_tex, H_data->r_k16, channelDesc);
 
  rk17_tex.addressMode[0] = cudaAddressModeClamp;
  rk17_tex.addressMode[1] = cudaAddressModeClamp;
  rk17_tex.filterMode = cudaFilterModeLinear;
  rk17_tex.normalized = true;
  cudaBindTextureToArray(rk17_tex, H_data->r_k17, channelDesc);
 
  rk18_tex.addressMode[0] = cudaAddressModeClamp;
  rk18_tex.addressMode[1] = cudaAddressModeClamp;
  rk18_tex.filterMode = cudaFilterModeLinear;
  rk18_tex.normalized = true;
  cudaBindTextureToArray(rk18_tex, H_data->r_k18, channelDesc);
 
  rk19_tex.addressMode[0] = cudaAddressModeClamp;
  rk19_tex.addressMode[1] = cudaAddressModeClamp;
  rk19_tex.filterMode = cudaFilterModeLinear;
  rk19_tex.normalized = true;
  cudaBindTextureToArray(rk19_tex, H_data->r_k19, channelDesc);
 
/*
  rk20_tex.addressMode[0] = cudaAddressModeClamp;
  rk20_tex.addressMode[1] = cudaAddressModeClamp;
  rk20_tex.filterMode = cudaFilterModeLinear;
  rk20_tex.normalized = true;
  cudaBindTextureToArray(rk20_tex, H_data->r_k20, channelDesc);
*/
  rk21_tex.addressMode[0] = cudaAddressModeClamp;
  rk21_tex.addressMode[1] = cudaAddressModeClamp;
  rk21_tex.filterMode = cudaFilterModeLinear;
  rk21_tex.normalized = true;
  cudaBindTextureToArray(rk21_tex, H_data->r_k21, channelDesc);
 
  rk22_tex.addressMode[0] = cudaAddressModeClamp;
  rk22_tex.addressMode[1] = cudaAddressModeClamp;
  rk22_tex.filterMode = cudaFilterModeLinear;
  rk22_tex.normalized = true;
  cudaBindTextureToArray(rk22_tex, H_data->r_k22, channelDesc);
 
/*
  rk23_tex.addressMode[0] = cudaAddressModeClamp;
  rk23_tex.addressMode[1] = cudaAddressModeClamp;
  rk23_tex.filterMode = cudaFilterModeLinear;
  rk23_tex.normalized = true;
  cudaBindTextureToArray(rk23_tex, H_data->r_k23, channelDesc);
 */

}

int main(){

  // read the ratedata from the hdf5 file
  // as a default cvklu_data
  const char * FileLocation = "cvklu_tables.h5";
  cvklu_data *host_rateData = cvklu_setup_data( FileLocation, NULL, NULL);

  // now initialize a cuda rate object
  // which each pointer points to a cudaArray
  // that is then binded to the globally defined texture object
  cuArray_ratedata *cuda_ratedata = NULL;
  cuda_ratedata = (cuArray_ratedata *) malloc(sizeof(cuArray_ratedata));
  initialize_cuArray_ratedata( cuda_ratedata, host_rateData ); 
  bind_ratedata_to_texture( cuda_ratedata );

  // copy rate data to device memory
  cvklu_data* device_rateData= NULL;
  cudaMalloc( (void**) &device_rateData, sizeof(cvklu_data) );
  cudaMemcpy( device_rateData, host_rateData, sizeof(cvklu_data), cudaMemcpyHostToDevice );  

  // Allocate CUDA array in device memoryi

  /*
  // this describes the format of a texture element
  // for float texels we could create a channel with:
  cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc(32, 0,0,0, cudaChannelFormatKindFloat);

  // copy to device memory that are located at r_k01
  cudaArray* cu_rk01;
  cudaMallocArray(&cu_rk01, &channelDesc, RDATA_NBINS);
  cudaMemcpyToArray( cu_rk01, 0, 0, host_rateData->r_k01, sizeof(float)*RDATA_NBINS, cudaMemcpyHostToDevice );

  // Set texture parameters
  rk01_tex.addressMode[0] = cudaAddressModeClamp;
  rk01_tex.addressMode[1] = cudaAddressModeWrap;
  rk01_tex.filterMode = cudaFilterModeLinear;
  rk01_tex.normalized = true;
  cudaBindTextureToArray( rk01_tex, cu_rk01, channelDesc );
  */  



/*
  float *cu_rk01 = NULL;
  cudaMalloc( (void**) &cu_rk01,  sizeof(float) * 1024  );
  cudaMemcpy( cu_rk01, host_rateData->r_k01, sizeof(float)*RDATA_NBINS, cudaMemcpyHostToDevice );
*/

 
  printf("finished reading memory \n ");
  // bind the rate data to individual texture object
//  bind_rate_to_texture( cu_rateData );
//  cudaBindTexture(NULL, rk01_tex, cu_rk01, sizeof(float) * 1024 );
  printf("finshed bindinf texture\n" ); 
 
  int NBLOCK  = 512;
  int NTHREAD = 256;
  int N = NBLOCK * NTHREAD;


  float  *temp_out = (float *) malloc(NBLOCK * NTHREAD * sizeof(float));
  double *rrate_out_tex = (double *) malloc( NRATE * NBLOCK * NTHREAD* sizeof(double) );
  double *rrate_out     = (double *) malloc( NRATE * NBLOCK * NTHREAD * sizeof(double) );

  float  *d_temp_out;
  double *d_rrate_tex;
  double *d_rrate;
 
  cudaMalloc( (void**)&d_temp_out,  sizeof(float)*N );
  cudaMalloc( (void**)&d_rrate_tex, sizeof(double)*N*NRATE);
  cudaMalloc( (void**)&d_rrate,     sizeof(double)*N*NRATE);

  float PI = 3.14159;

  float Tmax = 5.0e4;
  float Tmin = 10.0;
  for (int i = 0; i < N; i++){
    temp_out[i] = Tmax *( (float) i / (float) N * PI/ 2.0) + Tmin;
  }      

  cudaMemcpy(d_temp_out, temp_out, sizeof(float)*N , cudaMemcpyHostToDevice);

cudaEvent_t start, stop;
cudaEventCreate(&start);
cudaEventCreate(&stop); 
  float milliseconds = 0;

  cudaEventRecord(start);
for (int i = 0; i< 100; i++)  interpolate_rates_double <<< NBLOCK, NTHREAD >>>( d_rrate, d_temp_out, device_rateData );
  cudaEventRecord(stop);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&milliseconds, start, stop);
  printf("interpolate time = %0.5g ms\n", milliseconds );


  cudaEventRecord(start);
 for (int i = 0; i< 100; i++)  interpolate_rates_texture<<< NBLOCK, NTHREAD >>>( d_rrate_tex, d_temp_out);
  cudaEventRecord(stop);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&milliseconds, start, stop);
  printf("texture time = %0.5g ms\n", milliseconds );



  
  cudaMemcpy( rrate_out_tex, d_rrate_tex, sizeof(double)*N*NRATE, cudaMemcpyDeviceToHost );
  cudaMemcpy( rrate_out,     d_rrate,     sizeof(double)*N*NRATE, cudaMemcpyDeviceToHost );

  printf("rtex[132] = %0.5g ; r[132] = %0.5g\n", rrate_out_tex[132], rrate_out[132]);

  for (int i = 0; i < N; i++ ){
    if ( fabs( rrate_out_tex[i] - rrate_out[i] ) > 1.0e-6 ){
      printf("at %d, T = %0.3g K, rate_tex = %0.5g, rate = %0.5g\n") ; 
    }
  } 

  cudaDeviceReset() ;
  cudaPeekAtLastError() ;

  return 0; 
}

/*
__global__ void interpolate_reaction_rates( double *reaction_rates_out, double temp_out ){

  int bin_id;
  double t1, t2;
  double Tdef, log_temp_out;
  double RDATA_LB = 0.0; // this should be set as global constant
 
  log_temp_out = log(temp_out);
  bin_id = ( int ) RDATA_IDBIN * ( log_temp_out - LOG_RDATA_LB );
    
  if ( bin_id <= 0) {
    bin_id = 0;
  } else if ( bin_id >= RDATA_NBINS) {
    bin_id = RDATA_NBINS - 1;
  }
  
  t1   = LOG_RDATA_LB + (bin_id) * RDATA_DBIN;
  Tdef = ( log_temp_out - t1) * RDATA_IDBIN; 

  interp1D( reaction_rates_out, bin_id, Tdef );


}


template <typename T>
__host__ __device__
inline T lerp( T v0, T v1, T t){
  //https://devblogs.nvidia.com/lerp-faster-cuda/ 
  return fma(t, v1, fma(-t, v0, v0));
}

__device__ void interp1D_double( double *du, double bin_id, double Tdef  )
{

  // fetch the data masquereding as int2
  int2 v0_int2 = tex1Dfetch(tex_u, bin_id);
  int2 v1_int2 = tex1Dfetch(tex_u, bin_id + 1);

  // convert the int2 data into double
  double v0 = __hiloint2double( v0_int2.y, v0_int2.x );
  double v1 = __hiloint2double( v1_int2.y, v1_int2.x );

  du[0] = lerp( v0, v1, Tdef ); 

}

*/
