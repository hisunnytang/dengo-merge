#include "dengo_solver.h"
#include "mechanism.cuh"

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
    H5LTread_dataset_double(file_id, "/k01", data->r_k01);
    H5LTread_dataset_double(file_id, "/k02", data->r_k02);
    H5LTread_dataset_double(file_id, "/k03", data->r_k03);
    H5LTread_dataset_double(file_id, "/k04", data->r_k04);
    H5LTread_dataset_double(file_id, "/k05", data->r_k05);
    H5LTread_dataset_double(file_id, "/k06", data->r_k06);
    H5LTread_dataset_double(file_id, "/k07", data->r_k07);
    H5LTread_dataset_double(file_id, "/k08", data->r_k08);
    H5LTread_dataset_double(file_id, "/k09", data->r_k09);
    H5LTread_dataset_double(file_id, "/k10", data->r_k10);
    H5LTread_dataset_double(file_id, "/k11", data->r_k11);
    H5LTread_dataset_double(file_id, "/k12", data->r_k12);
    H5LTread_dataset_double(file_id, "/k13", data->r_k13);
    H5LTread_dataset_double(file_id, "/k14", data->r_k14);
    H5LTread_dataset_double(file_id, "/k15", data->r_k15);
    H5LTread_dataset_double(file_id, "/k16", data->r_k16);
    H5LTread_dataset_double(file_id, "/k17", data->r_k17);
    H5LTread_dataset_double(file_id, "/k18", data->r_k18);
    H5LTread_dataset_double(file_id, "/k19", data->r_k19);
    H5LTread_dataset_double(file_id, "/k21", data->r_k21);
    H5LTread_dataset_double(file_id, "/k22", data->r_k22);
    
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
    H5LTread_dataset_double(file_id, "/brem_brem",
                            data->c_brem_brem);
    H5LTread_dataset_double(file_id, "/ceHeI_ceHeI",
                            data->c_ceHeI_ceHeI);
    H5LTread_dataset_double(file_id, "/ceHeII_ceHeII",
                            data->c_ceHeII_ceHeII);
    H5LTread_dataset_double(file_id, "/ceHI_ceHI",
                            data->c_ceHI_ceHI);
    H5LTread_dataset_double(file_id, "/cie_cooling_cieco",
                            data->c_cie_cooling_cieco);
    H5LTread_dataset_double(file_id, "/ciHeI_ciHeI",
                            data->c_ciHeI_ciHeI);
    H5LTread_dataset_double(file_id, "/ciHeII_ciHeII",
                            data->c_ciHeII_ciHeII);
    H5LTread_dataset_double(file_id, "/ciHeIS_ciHeIS",
                            data->c_ciHeIS_ciHeIS);
    H5LTread_dataset_double(file_id, "/ciHI_ciHI",
                            data->c_ciHI_ciHI);
    H5LTread_dataset_double(file_id, "/compton_comp_",
                            data->c_compton_comp_);
    H5LTread_dataset_double(file_id, "/gammah_gammah",
                            data->c_gammah_gammah);
    H5LTread_dataset_double(file_id, "/gloverabel08_gael",
                            data->c_gloverabel08_gael);
    H5LTread_dataset_double(file_id, "/gloverabel08_gaH2",
                            data->c_gloverabel08_gaH2);
    H5LTread_dataset_double(file_id, "/gloverabel08_gaHe",
                            data->c_gloverabel08_gaHe);
    H5LTread_dataset_double(file_id, "/gloverabel08_gaHI",
                            data->c_gloverabel08_gaHI);
    H5LTread_dataset_double(file_id, "/gloverabel08_gaHp",
                            data->c_gloverabel08_gaHp);
    H5LTread_dataset_double(file_id, "/gloverabel08_gphdl",
                            data->c_gloverabel08_gphdl);
    H5LTread_dataset_double(file_id, "/gloverabel08_gpldl",
                            data->c_gloverabel08_gpldl);
    H5LTread_dataset_double(file_id, "/gloverabel08_h2lte",
                            data->c_gloverabel08_h2lte);
    H5LTread_dataset_double(file_id, "/h2formation_h2mcool",
                            data->c_h2formation_h2mcool);
    H5LTread_dataset_double(file_id, "/h2formation_h2mheat",
                            data->c_h2formation_h2mheat);
    H5LTread_dataset_double(file_id, "/h2formation_ncrd1",
                            data->c_h2formation_ncrd1);
    H5LTread_dataset_double(file_id, "/h2formation_ncrd2",
                            data->c_h2formation_ncrd2);
    H5LTread_dataset_double(file_id, "/h2formation_ncrn",
                            data->c_h2formation_ncrn);
    H5LTread_dataset_double(file_id, "/reHeII1_reHeII1",
                            data->c_reHeII1_reHeII1);
    H5LTread_dataset_double(file_id, "/reHeII2_reHeII2",
                            data->c_reHeII2_reHeII2);
    H5LTread_dataset_double(file_id, "/reHeIII_reHeIII",
                            data->c_reHeIII_reHeIII);
    H5LTread_dataset_double(file_id, "/reHII_reHII",
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
    H5LTread_dataset_double(file_id, "/gammaH2_1",
                            data->g_gammaH2_1 );
    H5LTread_dataset_double(file_id, "/dgammaH2_1_dT",
                            data->g_dgammaH2_1_dT );   
    
    H5LTread_dataset_double(file_id, "/gammaH2_2",
                            data->g_gammaH2_2 );
    H5LTread_dataset_double(file_id, "/dgammaH2_2_dT",
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


