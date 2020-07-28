
/* THIS FILE HAS BEEN AUTO-GENERATED.  DO NOT EDIT. */

/* This is C++ code to read HDF5 files for
   reaction rates, cooling rates, and initial
   conditions for the chemical network defined
   by the user.  In addition, this contains
   code for calculating temperature from the
   gas energy and computing the RHS and the
   Jacobian of the system of equations which
   will be fed into the solver.
*/


#include "primordial_solver.h"

///////////////////////////////////////////////////////////////////////////////
/////////// Setup the reaction, cooling rate data table ///////////////////////
///////////////////////////////////////////////////////////////////////////////
primordial_data *primordial_setup_data( const char *FileLocation, int *NumberOfFields, char ***FieldNames)
{

    //-----------------------------------------------------
    // Function : primordial_setup_data
    // Description: Initialize a data object that stores the reaction/ cooling rate data 
    //-----------------------------------------------------

    int i, n;
    
    primordial_data *data = (primordial_data *) malloc(sizeof(primordial_data));
    
    // point the module to look for primordial_tables.h5
    data->dengo_data_file = FileLocation;

    /* allocate space for the scale related pieces */

    // Number of cells to be solved in a batch 
    data->nstrip = MAX_NCELLS;
    /*initialize temperature so it wont crash*/
    for ( i = 0; i < MAX_NCELLS; i++ ){
        for( n = 0; n < NTHREADS; n++ ){
            data->Ts[n][i]    = 1000.0;
            data->logTs[n][i] = log(1000.0);
        }
    }

    /* Temperature-related pieces */
    data->bounds[0] = 1.0;
    data->bounds[1] = 100000000.0;
    data->nbins = 1024 - 1;
    data->dbin = (log(data->bounds[1]) - log(data->bounds[0])) / data->nbins;
    data->idbin = 1.0L / data->dbin;

    /* Redshift-related pieces */
    data->z_bounds[0] = 0.0;
    data->z_bounds[1] = 10.0;
    data->n_zbins = 0 - 1;
    data->d_zbin = (log(data->z_bounds[1] + 1.0) - log(data->z_bounds[0] + 1.0)) / data->n_zbins;
    data->id_zbin = 1.0L / data->d_zbin;
    
    primordial_read_rate_tables(data);
    //fprintf(stderr, "Successfully read in rate tables.\n");

    primordial_read_cooling_tables(data);
    //fprintf(stderr, "Successfully read in cooling rate tables.\n");
    
    primordial_read_gamma(data);
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


void primordial_read_rate_tables(primordial_data *data)
{
    const char * filedir;
    if (data->dengo_data_file != NULL){
        filedir =  data->dengo_data_file; 
    } else{
        filedir = "primordial_tables.h5";   
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
    H5LTread_dataset_double(file_id, "/k23", data->r_k23);
    
    H5Fclose(file_id);
}


void primordial_read_cooling_tables(primordial_data *data)
{

    const char * filedir;
    if (data->dengo_data_file != NULL){
        filedir =  data->dengo_data_file; 
    } else{
        filedir = "primordial_tables.h5";   
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

void primordial_read_gamma(primordial_data *data)
{

    const char * filedir;
    if (data->dengo_data_file != NULL){
        filedir =  data->dengo_data_file; 
    } else{
        filedir = "primordial_tables.h5";   
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
 


/*
   This setup may be different than the user may anticipate, as a result
   of the lockstep timestep we use for a pencil beam through the grid.
   As such, it accepts the number of things to interpolate and makes
   assumptions about the sizes of the rates.
*/

/* This also requires no templating other than for the solver name...*/
void primordial_interpolate_rates(primordial_data *data,
                    int nstrip)
{
    int i, bin_id, zbin_id;
    double lb, t1, t2;
    double lbz, z1, z2, Tdef, zdef;
    int no_photo = 0;
    lb = log(data->bounds[0]);
    lbz = log(data->z_bounds[0] + 1.0);


    i = 0;
    #ifdef _OPENMP
    int threadID = omp_get_thread_num();
    #else
    int threadID = 0;
    #endif

    for ( i = 0; i < nstrip; i++ ){
        data->bin_id[threadID][i] = bin_id = (int) (data->idbin * (data->logTs[threadID][i] - lb));
        if (data->bin_id[threadID][i] <= 0) {
            data->bin_id[threadID][i] = 0;
        } else if (data->bin_id[threadID][i] >= data->nbins) {
            data->bin_id[threadID][i] = data->nbins - 1;
        }
        t1 = (lb + (bin_id    ) * data->dbin);
        t2 = (lb + (bin_id + 1) * data->dbin);
        data->Tdef[threadID][i] = (data->logTs[threadID][i] - t1)/(t2 - t1);
        data->dT[threadID][i] = (t2 - t1);
        /*fprintf(stderr, "INTERP: %d, bin_id = %d, dT = % 0.16g, T = % 0.16g, logT = % 0.16g\n",
                i, data->bin_id[i], data->dT[i], data->Ts[i],
                data->logTs[i]);*/
    
    if ((data->current_z >= data->z_bounds[0]) && (data->current_z < data->z_bounds[1])) {
        zbin_id = (int) (data->id_zbin * (log(data->current_z + 1.0) - lbz));
        if (zbin_id <= 0) {
            zbin_id = 0;
        } else if (zbin_id >= data->n_zbins) {
            zbin_id = data->n_zbins - 1;
        }
        z1 = (lbz + (zbin_id    ) * data->d_zbin);
        z2 = (lbz + (zbin_id + 1) * data->d_zbin);
        data->zdef = (log(data->current_z + 1.0) - z1)/(z2 - z1);
        data->dz = (exp(z2) - exp(z1)); //note: given this, we don't have to divide rate of change by z
    } else {
        no_photo = 1;
    }
    }

    zdef   = data->zdef;
    
    
    
    for ( i = 0; i < nstrip; i++ ){
        bin_id = data->bin_id[threadID][i];
        Tdef   = data->Tdef[threadID][i];
        data->rs_k01[threadID][i] = data->r_k01[bin_id] +
            Tdef * (data->r_k01[bin_id+1] - data->r_k01[bin_id]);
        data->drs_k01[threadID][i] = (data->r_k01[bin_id+1] - data->r_k01[bin_id]);
        data->drs_k01[threadID][i] /= data->dT[threadID][i];
        data->drs_k01[threadID][i] *= data->invTs[threadID][i];
    }
    
    
    for ( i = 0; i < nstrip; i++ ){
        bin_id = data->bin_id[threadID][i];
        Tdef   = data->Tdef[threadID][i];
        data->rs_k02[threadID][i] = data->r_k02[bin_id] +
            Tdef * (data->r_k02[bin_id+1] - data->r_k02[bin_id]);
        data->drs_k02[threadID][i] = (data->r_k02[bin_id+1] - data->r_k02[bin_id]);
        data->drs_k02[threadID][i] /= data->dT[threadID][i];
        data->drs_k02[threadID][i] *= data->invTs[threadID][i];
    }
    
    
    for ( i = 0; i < nstrip; i++ ){
        bin_id = data->bin_id[threadID][i];
        Tdef   = data->Tdef[threadID][i];
        data->rs_k03[threadID][i] = data->r_k03[bin_id] +
            Tdef * (data->r_k03[bin_id+1] - data->r_k03[bin_id]);
        data->drs_k03[threadID][i] = (data->r_k03[bin_id+1] - data->r_k03[bin_id]);
        data->drs_k03[threadID][i] /= data->dT[threadID][i];
        data->drs_k03[threadID][i] *= data->invTs[threadID][i];
    }
    
    
    for ( i = 0; i < nstrip; i++ ){
        bin_id = data->bin_id[threadID][i];
        Tdef   = data->Tdef[threadID][i];
        data->rs_k04[threadID][i] = data->r_k04[bin_id] +
            Tdef * (data->r_k04[bin_id+1] - data->r_k04[bin_id]);
        data->drs_k04[threadID][i] = (data->r_k04[bin_id+1] - data->r_k04[bin_id]);
        data->drs_k04[threadID][i] /= data->dT[threadID][i];
        data->drs_k04[threadID][i] *= data->invTs[threadID][i];
    }
    
    
    for ( i = 0; i < nstrip; i++ ){
        bin_id = data->bin_id[threadID][i];
        Tdef   = data->Tdef[threadID][i];
        data->rs_k05[threadID][i] = data->r_k05[bin_id] +
            Tdef * (data->r_k05[bin_id+1] - data->r_k05[bin_id]);
        data->drs_k05[threadID][i] = (data->r_k05[bin_id+1] - data->r_k05[bin_id]);
        data->drs_k05[threadID][i] /= data->dT[threadID][i];
        data->drs_k05[threadID][i] *= data->invTs[threadID][i];
    }
    
    
    for ( i = 0; i < nstrip; i++ ){
        bin_id = data->bin_id[threadID][i];
        Tdef   = data->Tdef[threadID][i];
        data->rs_k06[threadID][i] = data->r_k06[bin_id] +
            Tdef * (data->r_k06[bin_id+1] - data->r_k06[bin_id]);
        data->drs_k06[threadID][i] = (data->r_k06[bin_id+1] - data->r_k06[bin_id]);
        data->drs_k06[threadID][i] /= data->dT[threadID][i];
        data->drs_k06[threadID][i] *= data->invTs[threadID][i];
    }
    
    
    for ( i = 0; i < nstrip; i++ ){
        bin_id = data->bin_id[threadID][i];
        Tdef   = data->Tdef[threadID][i];
        data->rs_k07[threadID][i] = data->r_k07[bin_id] +
            Tdef * (data->r_k07[bin_id+1] - data->r_k07[bin_id]);
        data->drs_k07[threadID][i] = (data->r_k07[bin_id+1] - data->r_k07[bin_id]);
        data->drs_k07[threadID][i] /= data->dT[threadID][i];
        data->drs_k07[threadID][i] *= data->invTs[threadID][i];
    }
    
    
    for ( i = 0; i < nstrip; i++ ){
        bin_id = data->bin_id[threadID][i];
        Tdef   = data->Tdef[threadID][i];
        data->rs_k08[threadID][i] = data->r_k08[bin_id] +
            Tdef * (data->r_k08[bin_id+1] - data->r_k08[bin_id]);
        data->drs_k08[threadID][i] = (data->r_k08[bin_id+1] - data->r_k08[bin_id]);
        data->drs_k08[threadID][i] /= data->dT[threadID][i];
        data->drs_k08[threadID][i] *= data->invTs[threadID][i];
    }
    
    
    for ( i = 0; i < nstrip; i++ ){
        bin_id = data->bin_id[threadID][i];
        Tdef   = data->Tdef[threadID][i];
        data->rs_k09[threadID][i] = data->r_k09[bin_id] +
            Tdef * (data->r_k09[bin_id+1] - data->r_k09[bin_id]);
        data->drs_k09[threadID][i] = (data->r_k09[bin_id+1] - data->r_k09[bin_id]);
        data->drs_k09[threadID][i] /= data->dT[threadID][i];
        data->drs_k09[threadID][i] *= data->invTs[threadID][i];
    }
    
    
    for ( i = 0; i < nstrip; i++ ){
        bin_id = data->bin_id[threadID][i];
        Tdef   = data->Tdef[threadID][i];
        data->rs_k10[threadID][i] = data->r_k10[bin_id] +
            Tdef * (data->r_k10[bin_id+1] - data->r_k10[bin_id]);
        data->drs_k10[threadID][i] = (data->r_k10[bin_id+1] - data->r_k10[bin_id]);
        data->drs_k10[threadID][i] /= data->dT[threadID][i];
        data->drs_k10[threadID][i] *= data->invTs[threadID][i];
    }
    
    
    for ( i = 0; i < nstrip; i++ ){
        bin_id = data->bin_id[threadID][i];
        Tdef   = data->Tdef[threadID][i];
        data->rs_k11[threadID][i] = data->r_k11[bin_id] +
            Tdef * (data->r_k11[bin_id+1] - data->r_k11[bin_id]);
        data->drs_k11[threadID][i] = (data->r_k11[bin_id+1] - data->r_k11[bin_id]);
        data->drs_k11[threadID][i] /= data->dT[threadID][i];
        data->drs_k11[threadID][i] *= data->invTs[threadID][i];
    }
    
    
    for ( i = 0; i < nstrip; i++ ){
        bin_id = data->bin_id[threadID][i];
        Tdef   = data->Tdef[threadID][i];
        data->rs_k12[threadID][i] = data->r_k12[bin_id] +
            Tdef * (data->r_k12[bin_id+1] - data->r_k12[bin_id]);
        data->drs_k12[threadID][i] = (data->r_k12[bin_id+1] - data->r_k12[bin_id]);
        data->drs_k12[threadID][i] /= data->dT[threadID][i];
        data->drs_k12[threadID][i] *= data->invTs[threadID][i];
    }
    
    
    for ( i = 0; i < nstrip; i++ ){
        bin_id = data->bin_id[threadID][i];
        Tdef   = data->Tdef[threadID][i];
        data->rs_k13[threadID][i] = data->r_k13[bin_id] +
            Tdef * (data->r_k13[bin_id+1] - data->r_k13[bin_id]);
        data->drs_k13[threadID][i] = (data->r_k13[bin_id+1] - data->r_k13[bin_id]);
        data->drs_k13[threadID][i] /= data->dT[threadID][i];
        data->drs_k13[threadID][i] *= data->invTs[threadID][i];
    }
    
    
    for ( i = 0; i < nstrip; i++ ){
        bin_id = data->bin_id[threadID][i];
        Tdef   = data->Tdef[threadID][i];
        data->rs_k14[threadID][i] = data->r_k14[bin_id] +
            Tdef * (data->r_k14[bin_id+1] - data->r_k14[bin_id]);
        data->drs_k14[threadID][i] = (data->r_k14[bin_id+1] - data->r_k14[bin_id]);
        data->drs_k14[threadID][i] /= data->dT[threadID][i];
        data->drs_k14[threadID][i] *= data->invTs[threadID][i];
    }
    
    
    for ( i = 0; i < nstrip; i++ ){
        bin_id = data->bin_id[threadID][i];
        Tdef   = data->Tdef[threadID][i];
        data->rs_k15[threadID][i] = data->r_k15[bin_id] +
            Tdef * (data->r_k15[bin_id+1] - data->r_k15[bin_id]);
        data->drs_k15[threadID][i] = (data->r_k15[bin_id+1] - data->r_k15[bin_id]);
        data->drs_k15[threadID][i] /= data->dT[threadID][i];
        data->drs_k15[threadID][i] *= data->invTs[threadID][i];
    }
    
    
    for ( i = 0; i < nstrip; i++ ){
        bin_id = data->bin_id[threadID][i];
        Tdef   = data->Tdef[threadID][i];
        data->rs_k16[threadID][i] = data->r_k16[bin_id] +
            Tdef * (data->r_k16[bin_id+1] - data->r_k16[bin_id]);
        data->drs_k16[threadID][i] = (data->r_k16[bin_id+1] - data->r_k16[bin_id]);
        data->drs_k16[threadID][i] /= data->dT[threadID][i];
        data->drs_k16[threadID][i] *= data->invTs[threadID][i];
    }
    
    
    for ( i = 0; i < nstrip; i++ ){
        bin_id = data->bin_id[threadID][i];
        Tdef   = data->Tdef[threadID][i];
        data->rs_k17[threadID][i] = data->r_k17[bin_id] +
            Tdef * (data->r_k17[bin_id+1] - data->r_k17[bin_id]);
        data->drs_k17[threadID][i] = (data->r_k17[bin_id+1] - data->r_k17[bin_id]);
        data->drs_k17[threadID][i] /= data->dT[threadID][i];
        data->drs_k17[threadID][i] *= data->invTs[threadID][i];
    }
    
    
    for ( i = 0; i < nstrip; i++ ){
        bin_id = data->bin_id[threadID][i];
        Tdef   = data->Tdef[threadID][i];
        data->rs_k18[threadID][i] = data->r_k18[bin_id] +
            Tdef * (data->r_k18[bin_id+1] - data->r_k18[bin_id]);
        data->drs_k18[threadID][i] = (data->r_k18[bin_id+1] - data->r_k18[bin_id]);
        data->drs_k18[threadID][i] /= data->dT[threadID][i];
        data->drs_k18[threadID][i] *= data->invTs[threadID][i];
    }
    
    
    for ( i = 0; i < nstrip; i++ ){
        bin_id = data->bin_id[threadID][i];
        Tdef   = data->Tdef[threadID][i];
        data->rs_k19[threadID][i] = data->r_k19[bin_id] +
            Tdef * (data->r_k19[bin_id+1] - data->r_k19[bin_id]);
        data->drs_k19[threadID][i] = (data->r_k19[bin_id+1] - data->r_k19[bin_id]);
        data->drs_k19[threadID][i] /= data->dT[threadID][i];
        data->drs_k19[threadID][i] *= data->invTs[threadID][i];
    }
    
    
    for ( i = 0; i < nstrip; i++ ){
        bin_id = data->bin_id[threadID][i];
        Tdef   = data->Tdef[threadID][i];
        data->rs_k21[threadID][i] = data->r_k21[bin_id] +
            Tdef * (data->r_k21[bin_id+1] - data->r_k21[bin_id]);
        data->drs_k21[threadID][i] = (data->r_k21[bin_id+1] - data->r_k21[bin_id]);
        data->drs_k21[threadID][i] /= data->dT[threadID][i];
        data->drs_k21[threadID][i] *= data->invTs[threadID][i];
    }
    
    
    for ( i = 0; i < nstrip; i++ ){
        bin_id = data->bin_id[threadID][i];
        Tdef   = data->Tdef[threadID][i];
        data->rs_k22[threadID][i] = data->r_k22[bin_id] +
            Tdef * (data->r_k22[bin_id+1] - data->r_k22[bin_id]);
        data->drs_k22[threadID][i] = (data->r_k22[bin_id+1] - data->r_k22[bin_id]);
        data->drs_k22[threadID][i] /= data->dT[threadID][i];
        data->drs_k22[threadID][i] *= data->invTs[threadID][i];
    }
    
    
    for ( i = 0; i < nstrip; i++ ){
        bin_id = data->bin_id[threadID][i];
        Tdef   = data->Tdef[threadID][i];
        data->rs_k23[threadID][i] = data->r_k23[bin_id] +
            Tdef * (data->r_k23[bin_id+1] - data->r_k23[bin_id]);
        data->drs_k23[threadID][i] = (data->r_k23[bin_id+1] - data->r_k23[bin_id]);
        data->drs_k23[threadID][i] /= data->dT[threadID][i];
        data->drs_k23[threadID][i] *= data->invTs[threadID][i];
    }
    
    for ( i = 0; i < nstrip; i++ ){
        bin_id = data->bin_id[threadID][i];
        Tdef   = data->Tdef[threadID][i];
        data->cs_brem_brem[threadID][i] = data->c_brem_brem[bin_id] +
            Tdef * (data->c_brem_brem[bin_id+1] - data->c_brem_brem[bin_id]);
        data->dcs_brem_brem[threadID][i] = (data->c_brem_brem[bin_id+1] - data->c_brem_brem[bin_id]);
        data->dcs_brem_brem[threadID][i] /= data->dT[threadID][i];
        data->dcs_brem_brem[threadID][i] *= data->invTs[threadID][i];
    }
    
    for ( i = 0; i < nstrip; i++ ){
        bin_id = data->bin_id[threadID][i];
        Tdef   = data->Tdef[threadID][i];
        data->cs_ceHeI_ceHeI[threadID][i] = data->c_ceHeI_ceHeI[bin_id] +
            Tdef * (data->c_ceHeI_ceHeI[bin_id+1] - data->c_ceHeI_ceHeI[bin_id]);
        data->dcs_ceHeI_ceHeI[threadID][i] = (data->c_ceHeI_ceHeI[bin_id+1] - data->c_ceHeI_ceHeI[bin_id]);
        data->dcs_ceHeI_ceHeI[threadID][i] /= data->dT[threadID][i];
        data->dcs_ceHeI_ceHeI[threadID][i] *= data->invTs[threadID][i];
    }
    
    for ( i = 0; i < nstrip; i++ ){
        bin_id = data->bin_id[threadID][i];
        Tdef   = data->Tdef[threadID][i];
        data->cs_ceHeII_ceHeII[threadID][i] = data->c_ceHeII_ceHeII[bin_id] +
            Tdef * (data->c_ceHeII_ceHeII[bin_id+1] - data->c_ceHeII_ceHeII[bin_id]);
        data->dcs_ceHeII_ceHeII[threadID][i] = (data->c_ceHeII_ceHeII[bin_id+1] - data->c_ceHeII_ceHeII[bin_id]);
        data->dcs_ceHeII_ceHeII[threadID][i] /= data->dT[threadID][i];
        data->dcs_ceHeII_ceHeII[threadID][i] *= data->invTs[threadID][i];
    }
    
    for ( i = 0; i < nstrip; i++ ){
        bin_id = data->bin_id[threadID][i];
        Tdef   = data->Tdef[threadID][i];
        data->cs_ceHI_ceHI[threadID][i] = data->c_ceHI_ceHI[bin_id] +
            Tdef * (data->c_ceHI_ceHI[bin_id+1] - data->c_ceHI_ceHI[bin_id]);
        data->dcs_ceHI_ceHI[threadID][i] = (data->c_ceHI_ceHI[bin_id+1] - data->c_ceHI_ceHI[bin_id]);
        data->dcs_ceHI_ceHI[threadID][i] /= data->dT[threadID][i];
        data->dcs_ceHI_ceHI[threadID][i] *= data->invTs[threadID][i];
    }
    
    for ( i = 0; i < nstrip; i++ ){
        bin_id = data->bin_id[threadID][i];
        Tdef   = data->Tdef[threadID][i];
        data->cs_cie_cooling_cieco[threadID][i] = data->c_cie_cooling_cieco[bin_id] +
            Tdef * (data->c_cie_cooling_cieco[bin_id+1] - data->c_cie_cooling_cieco[bin_id]);
        data->dcs_cie_cooling_cieco[threadID][i] = (data->c_cie_cooling_cieco[bin_id+1] - data->c_cie_cooling_cieco[bin_id]);
        data->dcs_cie_cooling_cieco[threadID][i] /= data->dT[threadID][i];
        data->dcs_cie_cooling_cieco[threadID][i] *= data->invTs[threadID][i];
    }
    
    for ( i = 0; i < nstrip; i++ ){
        bin_id = data->bin_id[threadID][i];
        Tdef   = data->Tdef[threadID][i];
        data->cs_ciHeI_ciHeI[threadID][i] = data->c_ciHeI_ciHeI[bin_id] +
            Tdef * (data->c_ciHeI_ciHeI[bin_id+1] - data->c_ciHeI_ciHeI[bin_id]);
        data->dcs_ciHeI_ciHeI[threadID][i] = (data->c_ciHeI_ciHeI[bin_id+1] - data->c_ciHeI_ciHeI[bin_id]);
        data->dcs_ciHeI_ciHeI[threadID][i] /= data->dT[threadID][i];
        data->dcs_ciHeI_ciHeI[threadID][i] *= data->invTs[threadID][i];
    }
    
    for ( i = 0; i < nstrip; i++ ){
        bin_id = data->bin_id[threadID][i];
        Tdef   = data->Tdef[threadID][i];
        data->cs_ciHeII_ciHeII[threadID][i] = data->c_ciHeII_ciHeII[bin_id] +
            Tdef * (data->c_ciHeII_ciHeII[bin_id+1] - data->c_ciHeII_ciHeII[bin_id]);
        data->dcs_ciHeII_ciHeII[threadID][i] = (data->c_ciHeII_ciHeII[bin_id+1] - data->c_ciHeII_ciHeII[bin_id]);
        data->dcs_ciHeII_ciHeII[threadID][i] /= data->dT[threadID][i];
        data->dcs_ciHeII_ciHeII[threadID][i] *= data->invTs[threadID][i];
    }
    
    for ( i = 0; i < nstrip; i++ ){
        bin_id = data->bin_id[threadID][i];
        Tdef   = data->Tdef[threadID][i];
        data->cs_ciHeIS_ciHeIS[threadID][i] = data->c_ciHeIS_ciHeIS[bin_id] +
            Tdef * (data->c_ciHeIS_ciHeIS[bin_id+1] - data->c_ciHeIS_ciHeIS[bin_id]);
        data->dcs_ciHeIS_ciHeIS[threadID][i] = (data->c_ciHeIS_ciHeIS[bin_id+1] - data->c_ciHeIS_ciHeIS[bin_id]);
        data->dcs_ciHeIS_ciHeIS[threadID][i] /= data->dT[threadID][i];
        data->dcs_ciHeIS_ciHeIS[threadID][i] *= data->invTs[threadID][i];
    }
    
    for ( i = 0; i < nstrip; i++ ){
        bin_id = data->bin_id[threadID][i];
        Tdef   = data->Tdef[threadID][i];
        data->cs_ciHI_ciHI[threadID][i] = data->c_ciHI_ciHI[bin_id] +
            Tdef * (data->c_ciHI_ciHI[bin_id+1] - data->c_ciHI_ciHI[bin_id]);
        data->dcs_ciHI_ciHI[threadID][i] = (data->c_ciHI_ciHI[bin_id+1] - data->c_ciHI_ciHI[bin_id]);
        data->dcs_ciHI_ciHI[threadID][i] /= data->dT[threadID][i];
        data->dcs_ciHI_ciHI[threadID][i] *= data->invTs[threadID][i];
    }
    
    for ( i = 0; i < nstrip; i++ ){
        bin_id = data->bin_id[threadID][i];
        Tdef   = data->Tdef[threadID][i];
        data->cs_compton_comp_[threadID][i] = data->c_compton_comp_[bin_id] +
            Tdef * (data->c_compton_comp_[bin_id+1] - data->c_compton_comp_[bin_id]);
        data->dcs_compton_comp_[threadID][i] = (data->c_compton_comp_[bin_id+1] - data->c_compton_comp_[bin_id]);
        data->dcs_compton_comp_[threadID][i] /= data->dT[threadID][i];
        data->dcs_compton_comp_[threadID][i] *= data->invTs[threadID][i];
    }
    
    for ( i = 0; i < nstrip; i++ ){
        bin_id = data->bin_id[threadID][i];
        Tdef   = data->Tdef[threadID][i];
        data->cs_gammah_gammah[threadID][i] = data->c_gammah_gammah[bin_id] +
            Tdef * (data->c_gammah_gammah[bin_id+1] - data->c_gammah_gammah[bin_id]);
        data->dcs_gammah_gammah[threadID][i] = (data->c_gammah_gammah[bin_id+1] - data->c_gammah_gammah[bin_id]);
        data->dcs_gammah_gammah[threadID][i] /= data->dT[threadID][i];
        data->dcs_gammah_gammah[threadID][i] *= data->invTs[threadID][i];
    }
    
    for ( i = 0; i < nstrip; i++ ){
        bin_id = data->bin_id[threadID][i];
        Tdef   = data->Tdef[threadID][i];
        data->cs_gloverabel08_gael[threadID][i] = data->c_gloverabel08_gael[bin_id] +
            Tdef * (data->c_gloverabel08_gael[bin_id+1] - data->c_gloverabel08_gael[bin_id]);
        data->dcs_gloverabel08_gael[threadID][i] = (data->c_gloverabel08_gael[bin_id+1] - data->c_gloverabel08_gael[bin_id]);
        data->dcs_gloverabel08_gael[threadID][i] /= data->dT[threadID][i];
        data->dcs_gloverabel08_gael[threadID][i] *= data->invTs[threadID][i];
    }
    for ( i = 0; i < nstrip; i++ ){
        bin_id = data->bin_id[threadID][i];
        Tdef   = data->Tdef[threadID][i];
        data->cs_gloverabel08_gaH2[threadID][i] = data->c_gloverabel08_gaH2[bin_id] +
            Tdef * (data->c_gloverabel08_gaH2[bin_id+1] - data->c_gloverabel08_gaH2[bin_id]);
        data->dcs_gloverabel08_gaH2[threadID][i] = (data->c_gloverabel08_gaH2[bin_id+1] - data->c_gloverabel08_gaH2[bin_id]);
        data->dcs_gloverabel08_gaH2[threadID][i] /= data->dT[threadID][i];
        data->dcs_gloverabel08_gaH2[threadID][i] *= data->invTs[threadID][i];
    }
    for ( i = 0; i < nstrip; i++ ){
        bin_id = data->bin_id[threadID][i];
        Tdef   = data->Tdef[threadID][i];
        data->cs_gloverabel08_gaHe[threadID][i] = data->c_gloverabel08_gaHe[bin_id] +
            Tdef * (data->c_gloverabel08_gaHe[bin_id+1] - data->c_gloverabel08_gaHe[bin_id]);
        data->dcs_gloverabel08_gaHe[threadID][i] = (data->c_gloverabel08_gaHe[bin_id+1] - data->c_gloverabel08_gaHe[bin_id]);
        data->dcs_gloverabel08_gaHe[threadID][i] /= data->dT[threadID][i];
        data->dcs_gloverabel08_gaHe[threadID][i] *= data->invTs[threadID][i];
    }
    for ( i = 0; i < nstrip; i++ ){
        bin_id = data->bin_id[threadID][i];
        Tdef   = data->Tdef[threadID][i];
        data->cs_gloverabel08_gaHI[threadID][i] = data->c_gloverabel08_gaHI[bin_id] +
            Tdef * (data->c_gloverabel08_gaHI[bin_id+1] - data->c_gloverabel08_gaHI[bin_id]);
        data->dcs_gloverabel08_gaHI[threadID][i] = (data->c_gloverabel08_gaHI[bin_id+1] - data->c_gloverabel08_gaHI[bin_id]);
        data->dcs_gloverabel08_gaHI[threadID][i] /= data->dT[threadID][i];
        data->dcs_gloverabel08_gaHI[threadID][i] *= data->invTs[threadID][i];
    }
    for ( i = 0; i < nstrip; i++ ){
        bin_id = data->bin_id[threadID][i];
        Tdef   = data->Tdef[threadID][i];
        data->cs_gloverabel08_gaHp[threadID][i] = data->c_gloverabel08_gaHp[bin_id] +
            Tdef * (data->c_gloverabel08_gaHp[bin_id+1] - data->c_gloverabel08_gaHp[bin_id]);
        data->dcs_gloverabel08_gaHp[threadID][i] = (data->c_gloverabel08_gaHp[bin_id+1] - data->c_gloverabel08_gaHp[bin_id]);
        data->dcs_gloverabel08_gaHp[threadID][i] /= data->dT[threadID][i];
        data->dcs_gloverabel08_gaHp[threadID][i] *= data->invTs[threadID][i];
    }
    for ( i = 0; i < nstrip; i++ ){
        bin_id = data->bin_id[threadID][i];
        Tdef   = data->Tdef[threadID][i];
        data->cs_gloverabel08_gphdl[threadID][i] = data->c_gloverabel08_gphdl[bin_id] +
            Tdef * (data->c_gloverabel08_gphdl[bin_id+1] - data->c_gloverabel08_gphdl[bin_id]);
        data->dcs_gloverabel08_gphdl[threadID][i] = (data->c_gloverabel08_gphdl[bin_id+1] - data->c_gloverabel08_gphdl[bin_id]);
        data->dcs_gloverabel08_gphdl[threadID][i] /= data->dT[threadID][i];
        data->dcs_gloverabel08_gphdl[threadID][i] *= data->invTs[threadID][i];
    }
    for ( i = 0; i < nstrip; i++ ){
        bin_id = data->bin_id[threadID][i];
        Tdef   = data->Tdef[threadID][i];
        data->cs_gloverabel08_gpldl[threadID][i] = data->c_gloverabel08_gpldl[bin_id] +
            Tdef * (data->c_gloverabel08_gpldl[bin_id+1] - data->c_gloverabel08_gpldl[bin_id]);
        data->dcs_gloverabel08_gpldl[threadID][i] = (data->c_gloverabel08_gpldl[bin_id+1] - data->c_gloverabel08_gpldl[bin_id]);
        data->dcs_gloverabel08_gpldl[threadID][i] /= data->dT[threadID][i];
        data->dcs_gloverabel08_gpldl[threadID][i] *= data->invTs[threadID][i];
    }
    for ( i = 0; i < nstrip; i++ ){
        bin_id = data->bin_id[threadID][i];
        Tdef   = data->Tdef[threadID][i];
        data->cs_gloverabel08_h2lte[threadID][i] = data->c_gloverabel08_h2lte[bin_id] +
            Tdef * (data->c_gloverabel08_h2lte[bin_id+1] - data->c_gloverabel08_h2lte[bin_id]);
        data->dcs_gloverabel08_h2lte[threadID][i] = (data->c_gloverabel08_h2lte[bin_id+1] - data->c_gloverabel08_h2lte[bin_id]);
        data->dcs_gloverabel08_h2lte[threadID][i] /= data->dT[threadID][i];
        data->dcs_gloverabel08_h2lte[threadID][i] *= data->invTs[threadID][i];
    }
    
    for ( i = 0; i < nstrip; i++ ){
        bin_id = data->bin_id[threadID][i];
        Tdef   = data->Tdef[threadID][i];
        data->cs_h2formation_h2mcool[threadID][i] = data->c_h2formation_h2mcool[bin_id] +
            Tdef * (data->c_h2formation_h2mcool[bin_id+1] - data->c_h2formation_h2mcool[bin_id]);
        data->dcs_h2formation_h2mcool[threadID][i] = (data->c_h2formation_h2mcool[bin_id+1] - data->c_h2formation_h2mcool[bin_id]);
        data->dcs_h2formation_h2mcool[threadID][i] /= data->dT[threadID][i];
        data->dcs_h2formation_h2mcool[threadID][i] *= data->invTs[threadID][i];
    }
    for ( i = 0; i < nstrip; i++ ){
        bin_id = data->bin_id[threadID][i];
        Tdef   = data->Tdef[threadID][i];
        data->cs_h2formation_h2mheat[threadID][i] = data->c_h2formation_h2mheat[bin_id] +
            Tdef * (data->c_h2formation_h2mheat[bin_id+1] - data->c_h2formation_h2mheat[bin_id]);
        data->dcs_h2formation_h2mheat[threadID][i] = (data->c_h2formation_h2mheat[bin_id+1] - data->c_h2formation_h2mheat[bin_id]);
        data->dcs_h2formation_h2mheat[threadID][i] /= data->dT[threadID][i];
        data->dcs_h2formation_h2mheat[threadID][i] *= data->invTs[threadID][i];
    }
    for ( i = 0; i < nstrip; i++ ){
        bin_id = data->bin_id[threadID][i];
        Tdef   = data->Tdef[threadID][i];
        data->cs_h2formation_ncrd1[threadID][i] = data->c_h2formation_ncrd1[bin_id] +
            Tdef * (data->c_h2formation_ncrd1[bin_id+1] - data->c_h2formation_ncrd1[bin_id]);
        data->dcs_h2formation_ncrd1[threadID][i] = (data->c_h2formation_ncrd1[bin_id+1] - data->c_h2formation_ncrd1[bin_id]);
        data->dcs_h2formation_ncrd1[threadID][i] /= data->dT[threadID][i];
        data->dcs_h2formation_ncrd1[threadID][i] *= data->invTs[threadID][i];
    }
    for ( i = 0; i < nstrip; i++ ){
        bin_id = data->bin_id[threadID][i];
        Tdef   = data->Tdef[threadID][i];
        data->cs_h2formation_ncrd2[threadID][i] = data->c_h2formation_ncrd2[bin_id] +
            Tdef * (data->c_h2formation_ncrd2[bin_id+1] - data->c_h2formation_ncrd2[bin_id]);
        data->dcs_h2formation_ncrd2[threadID][i] = (data->c_h2formation_ncrd2[bin_id+1] - data->c_h2formation_ncrd2[bin_id]);
        data->dcs_h2formation_ncrd2[threadID][i] /= data->dT[threadID][i];
        data->dcs_h2formation_ncrd2[threadID][i] *= data->invTs[threadID][i];
    }
    for ( i = 0; i < nstrip; i++ ){
        bin_id = data->bin_id[threadID][i];
        Tdef   = data->Tdef[threadID][i];
        data->cs_h2formation_ncrn[threadID][i] = data->c_h2formation_ncrn[bin_id] +
            Tdef * (data->c_h2formation_ncrn[bin_id+1] - data->c_h2formation_ncrn[bin_id]);
        data->dcs_h2formation_ncrn[threadID][i] = (data->c_h2formation_ncrn[bin_id+1] - data->c_h2formation_ncrn[bin_id]);
        data->dcs_h2formation_ncrn[threadID][i] /= data->dT[threadID][i];
        data->dcs_h2formation_ncrn[threadID][i] *= data->invTs[threadID][i];
    }
    
    for ( i = 0; i < nstrip; i++ ){
        bin_id = data->bin_id[threadID][i];
        Tdef   = data->Tdef[threadID][i];
        data->cs_reHeII1_reHeII1[threadID][i] = data->c_reHeII1_reHeII1[bin_id] +
            Tdef * (data->c_reHeII1_reHeII1[bin_id+1] - data->c_reHeII1_reHeII1[bin_id]);
        data->dcs_reHeII1_reHeII1[threadID][i] = (data->c_reHeII1_reHeII1[bin_id+1] - data->c_reHeII1_reHeII1[bin_id]);
        data->dcs_reHeII1_reHeII1[threadID][i] /= data->dT[threadID][i];
        data->dcs_reHeII1_reHeII1[threadID][i] *= data->invTs[threadID][i];
    }
    
    for ( i = 0; i < nstrip; i++ ){
        bin_id = data->bin_id[threadID][i];
        Tdef   = data->Tdef[threadID][i];
        data->cs_reHeII2_reHeII2[threadID][i] = data->c_reHeII2_reHeII2[bin_id] +
            Tdef * (data->c_reHeII2_reHeII2[bin_id+1] - data->c_reHeII2_reHeII2[bin_id]);
        data->dcs_reHeII2_reHeII2[threadID][i] = (data->c_reHeII2_reHeII2[bin_id+1] - data->c_reHeII2_reHeII2[bin_id]);
        data->dcs_reHeII2_reHeII2[threadID][i] /= data->dT[threadID][i];
        data->dcs_reHeII2_reHeII2[threadID][i] *= data->invTs[threadID][i];
    }
    
    for ( i = 0; i < nstrip; i++ ){
        bin_id = data->bin_id[threadID][i];
        Tdef   = data->Tdef[threadID][i];
        data->cs_reHeIII_reHeIII[threadID][i] = data->c_reHeIII_reHeIII[bin_id] +
            Tdef * (data->c_reHeIII_reHeIII[bin_id+1] - data->c_reHeIII_reHeIII[bin_id]);
        data->dcs_reHeIII_reHeIII[threadID][i] = (data->c_reHeIII_reHeIII[bin_id+1] - data->c_reHeIII_reHeIII[bin_id]);
        data->dcs_reHeIII_reHeIII[threadID][i] /= data->dT[threadID][i];
        data->dcs_reHeIII_reHeIII[threadID][i] *= data->invTs[threadID][i];
    }
    
    for ( i = 0; i < nstrip; i++ ){
        bin_id = data->bin_id[threadID][i];
        Tdef   = data->Tdef[threadID][i];
        data->cs_reHII_reHII[threadID][i] = data->c_reHII_reHII[bin_id] +
            Tdef * (data->c_reHII_reHII[bin_id+1] - data->c_reHII_reHII[bin_id]);
        data->dcs_reHII_reHII[threadID][i] = (data->c_reHII_reHII[bin_id+1] - data->c_reHII_reHII[bin_id]);
        data->dcs_reHII_reHII[threadID][i] /= data->dT[threadID][i];
        data->dcs_reHII_reHII[threadID][i] *= data->invTs[threadID][i];
    }
    
    
}
 


void primordial_interpolate_gamma(primordial_data *data,
                    int i)
{   

    /*
     * find the bin_id for the given temperature 
     * update dT for i_th strip
     */

    int bin_id, zbin_id;
    double lb, t1, t2;
    double lbz, z1, z2;
    int no_photo = 0;
    lb = log(data->bounds[0]);
    lbz = log(data->z_bounds[0] + 1.0);
    
    #ifdef _OPENMP
    int threadID = omp_get_thread_num();
    #else
    int threadID = 0;
    #endif

    data->bin_id[threadID][i] = bin_id = (int) (data->idbin * (data->logTs[threadID][i] - lb));
    if (data->bin_id[threadID][i] <= 0) {
        data->bin_id[threadID][i] = 0;
    } else if (data->bin_id[threadID][i] >= data->nbins) {
        data->bin_id[threadID][i] = data->nbins - 1;
    }
    t1 = (lb + (bin_id    ) * data->dbin);
    t2 = (lb + (bin_id + 1) * data->dbin);
    data->Tdef[threadID][i] = (data->logTs[threadID][i] - t1)/(t2 - t1);
    data->dT[threadID][i] = (t2 - t1);

    
    
    bin_id = data->bin_id[threadID][i];
    data->gammaH2_2[threadID][i] = data->g_gammaH2_2[bin_id] +
        data->Tdef[threadID][i] * (data->g_gammaH2_2[bin_id+1] - data->g_gammaH2_2[bin_id]);

    data->dgammaH2_2_dT[threadID][i] = data->g_dgammaH2_2_dT[bin_id] +
        data->Tdef[threadID][i] * (data->g_dgammaH2_2_dT[bin_id+1] 
        - data->g_dgammaH2_2_dT[bin_id]);
    
    
    bin_id = data->bin_id[threadID][i];
    data->gammaH2_1[threadID][i] = data->g_gammaH2_1[bin_id] +
        data->Tdef[threadID][i] * (data->g_gammaH2_1[bin_id+1] - data->g_gammaH2_1[bin_id]);

    data->dgammaH2_1_dT[threadID][i] = data->g_dgammaH2_1_dT[bin_id] +
        data->Tdef[threadID][i] * (data->g_dgammaH2_1_dT[bin_id+1] 
        - data->g_dgammaH2_1_dT[bin_id]);
    
       
    }




///////////////////////////////////////////////////////////////////////////////
/////////////////// Main Evolution Routines            ////////////////////////
///////////////////////////////////////////////////////////////////////////////

int primordial_main(int argc, char** argv )
{
    //-----------------------------------------------------
    // Function : primordial_main
    // Description: this will look for initial condition files from the CLI, 
    //              evolve the ODE system to dtf specified from the CLI, (if not it's set to freefall time),
    //              and write the result to primordial_solution.h5 if output file name is not specified
    // Parameter:   argv[1]   : initial condition file name (in hdf5 format)
    //              argv[2]   : output file name (in hdf5 format)
    //              argv[3]   : desired final time reached by solver (in seconds)
    //
    // Note:        Units of initial conditions/ output 
    //              Baryons: mass density in a.m.u * cm^-3 ( 1mH = 1.00794 amu )
    //              de     : number density of electrons in cm^-3 
    //              ge     : internal energy per mass density of the cell (erg / g )
    //-----------------------------------------------------
    primordial_data *data = primordial_setup_data(NULL, NULL, NULL);

    /* Initial conditions */
    hid_t file_id;
    if (argc < 2){
    file_id = H5Fopen("primordial_initial_conditions.h5", H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file_id < 0) {fprintf(stderr, "Failed to open "
        "primordial_initial_conditions.h5 so dying.\n");
        return(1);}
    } else {
        file_id = H5Fopen( argv[1], H5F_ACC_RDONLY, H5P_DEFAULT);
        if (file_id < 0) {fprintf(stderr, "Failed to open  your initial_conditions file so dying.\n");
        return(1);}
       
            
    }

    /* Allocate the correct number of cells */
    hsize_t dims; /* We have flat versus number of species */

    /* Check gas energy to get the number of cells */
    fprintf(stderr, "Getting dimensionality from ge:\n");
    herr_t status = H5LTget_dataset_info(file_id, "/ge", &dims, NULL, NULL);
    if(status == -1) {
        fprintf(stderr, "Error opening initial conditions file.\n");
        return 1;
    }
    fprintf(stderr, "  ncells = % 3i\n", (int) dims);
    data->ncells = dims;

    int N = 10;

    double *atol, *rtol;
    atol = (double *) malloc(N * dims * sizeof(double));
    rtol = (double *) malloc(N * dims * sizeof(double));

    double *tics = (double *) malloc(dims * sizeof(double));
    double *ics = (double *) malloc(dims * N * sizeof(double));
    double *input = (double *) malloc(dims * N * sizeof(double));
    double *temp  = (double *) malloc(dims * sizeof(double) );

    unsigned int i = 0, j;
    
    fprintf(stderr, "Reading I.C. for /H2_1\n");
    H5LTread_dataset_double(file_id, "/H2_1", tics);
    for (j = 0; j < dims; j++) {
        ics[j * N + i] = tics[j]; 
        atol[j * N + i] = tics[j] * 1e-09;
        rtol[j * N + i] = 1e-09;
        if(j==0) {
            fprintf(stderr, "H2_1[0] = %0.3g, atol => % 0.16g\n",
                    tics[j], atol[j]);
        }
    }
    i++;
    
    fprintf(stderr, "Reading I.C. for /H2_2\n");
    H5LTread_dataset_double(file_id, "/H2_2", tics);
    for (j = 0; j < dims; j++) {
        ics[j * N + i] = tics[j]; 
        atol[j * N + i] = tics[j] * 1e-09;
        rtol[j * N + i] = 1e-09;
        if(j==0) {
            fprintf(stderr, "H2_2[0] = %0.3g, atol => % 0.16g\n",
                    tics[j], atol[j]);
        }
    }
    i++;
    
    fprintf(stderr, "Reading I.C. for /H_1\n");
    H5LTread_dataset_double(file_id, "/H_1", tics);
    for (j = 0; j < dims; j++) {
        ics[j * N + i] = tics[j]; 
        atol[j * N + i] = tics[j] * 1e-09;
        rtol[j * N + i] = 1e-09;
        if(j==0) {
            fprintf(stderr, "H_1[0] = %0.3g, atol => % 0.16g\n",
                    tics[j], atol[j]);
        }
    }
    i++;
    
    fprintf(stderr, "Reading I.C. for /H_2\n");
    H5LTread_dataset_double(file_id, "/H_2", tics);
    for (j = 0; j < dims; j++) {
        ics[j * N + i] = tics[j]; 
        atol[j * N + i] = tics[j] * 1e-09;
        rtol[j * N + i] = 1e-09;
        if(j==0) {
            fprintf(stderr, "H_2[0] = %0.3g, atol => % 0.16g\n",
                    tics[j], atol[j]);
        }
    }
    i++;
    
    fprintf(stderr, "Reading I.C. for /H_m0\n");
    H5LTread_dataset_double(file_id, "/H_m0", tics);
    for (j = 0; j < dims; j++) {
        ics[j * N + i] = tics[j]; 
        atol[j * N + i] = tics[j] * 1e-09;
        rtol[j * N + i] = 1e-09;
        if(j==0) {
            fprintf(stderr, "H_m0[0] = %0.3g, atol => % 0.16g\n",
                    tics[j], atol[j]);
        }
    }
    i++;
    
    fprintf(stderr, "Reading I.C. for /He_1\n");
    H5LTread_dataset_double(file_id, "/He_1", tics);
    for (j = 0; j < dims; j++) {
        ics[j * N + i] = tics[j]; 
        atol[j * N + i] = tics[j] * 1e-09;
        rtol[j * N + i] = 1e-09;
        if(j==0) {
            fprintf(stderr, "He_1[0] = %0.3g, atol => % 0.16g\n",
                    tics[j], atol[j]);
        }
    }
    i++;
    
    fprintf(stderr, "Reading I.C. for /He_2\n");
    H5LTread_dataset_double(file_id, "/He_2", tics);
    for (j = 0; j < dims; j++) {
        ics[j * N + i] = tics[j]; 
        atol[j * N + i] = tics[j] * 1e-09;
        rtol[j * N + i] = 1e-09;
        if(j==0) {
            fprintf(stderr, "He_2[0] = %0.3g, atol => % 0.16g\n",
                    tics[j], atol[j]);
        }
    }
    i++;
    
    fprintf(stderr, "Reading I.C. for /He_3\n");
    H5LTread_dataset_double(file_id, "/He_3", tics);
    for (j = 0; j < dims; j++) {
        ics[j * N + i] = tics[j]; 
        atol[j * N + i] = tics[j] * 1e-09;
        rtol[j * N + i] = 1e-09;
        if(j==0) {
            fprintf(stderr, "He_3[0] = %0.3g, atol => % 0.16g\n",
                    tics[j], atol[j]);
        }
    }
    i++;
    
    fprintf(stderr, "Reading I.C. for /de\n");
    H5LTread_dataset_double(file_id, "/de", tics);
    for (j = 0; j < dims; j++) {
        ics[j * N + i] = tics[j]; 
        atol[j * N + i] = tics[j] * 1e-09;
        rtol[j * N + i] = 1e-09;
        if(j==0) {
            fprintf(stderr, "de[0] = %0.3g, atol => % 0.16g\n",
                    tics[j], atol[j]);
        }
    }
    i++;
    
    fprintf(stderr, "Reading I.C. for /ge\n");
    H5LTread_dataset_double(file_id, "/ge", tics);
    for (j = 0; j < dims; j++) {
        ics[j * N + i] = tics[j]; 
        atol[j * N + i] = tics[j] * 1e-09;
        rtol[j * N + i] = 1e-09;
        if(j==0) {
            fprintf(stderr, "ge[0] = %0.3g, atol => % 0.16g\n",
                    tics[j], atol[j]);
        }
    }
    i++;
    
    
    double *density = (double *) malloc(dims *sizeof(double) );
    H5LTread_dataset_double(file_id, "/density", density);

    H5Fclose(file_id);
    
    // double dtf = 31557000000000.0;
    double dtf, t0;
    t0 = 2.992e15;
    dtf = t0 / sqrt(density[0]);
    
    // if the output time is specified,
    // it overrides the freefall time
    if (argc > 3){
        dtf = atof( argv[3] ); 
    }

    double dt = -1.0;
    double z = -1.0;
    for (i = 0; i < dims * N; i++) input[i] = ics[i];
    double ttot;
    int flag = dengo_evolve_primordial(dtf, dt, z, input, rtol, atol, dims, data, temp);
    if (flag > 0){
        fprintf(stderr, "solver failed, Time reached by the solver: %0.5g \n", dt);    
    }

    /* Write results to HDF5 file */

    if (argc < 3){
        file_id = H5Fcreate("primordial_solution.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    } else{
        file_id = H5Fcreate( argv[2], H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    }

    hsize_t dimsarr[1];
    dimsarr[0] = dims;
    i = 0;
    
    double H2_1[dims];
    for (j = 0; j < dims; j++) {
        H2_1[j] = input[j * N + i]; 
    }
    fprintf(stderr, "Writing solution for /H2_1\n");
    H5LTmake_dataset_double(file_id, "/H2_1", 1, dimsarr, H2_1);
    i++;
    
    double H2_2[dims];
    for (j = 0; j < dims; j++) {
        H2_2[j] = input[j * N + i]; 
    }
    fprintf(stderr, "Writing solution for /H2_2\n");
    H5LTmake_dataset_double(file_id, "/H2_2", 1, dimsarr, H2_2);
    i++;
    
    double H_1[dims];
    for (j = 0; j < dims; j++) {
        H_1[j] = input[j * N + i]; 
    }
    fprintf(stderr, "Writing solution for /H_1\n");
    H5LTmake_dataset_double(file_id, "/H_1", 1, dimsarr, H_1);
    i++;
    
    double H_2[dims];
    for (j = 0; j < dims; j++) {
        H_2[j] = input[j * N + i]; 
    }
    fprintf(stderr, "Writing solution for /H_2\n");
    H5LTmake_dataset_double(file_id, "/H_2", 1, dimsarr, H_2);
    i++;
    
    double H_m0[dims];
    for (j = 0; j < dims; j++) {
        H_m0[j] = input[j * N + i]; 
    }
    fprintf(stderr, "Writing solution for /H_m0\n");
    H5LTmake_dataset_double(file_id, "/H_m0", 1, dimsarr, H_m0);
    i++;
    
    double He_1[dims];
    for (j = 0; j < dims; j++) {
        He_1[j] = input[j * N + i]; 
    }
    fprintf(stderr, "Writing solution for /He_1\n");
    H5LTmake_dataset_double(file_id, "/He_1", 1, dimsarr, He_1);
    i++;
    
    double He_2[dims];
    for (j = 0; j < dims; j++) {
        He_2[j] = input[j * N + i]; 
    }
    fprintf(stderr, "Writing solution for /He_2\n");
    H5LTmake_dataset_double(file_id, "/He_2", 1, dimsarr, He_2);
    i++;
    
    double He_3[dims];
    for (j = 0; j < dims; j++) {
        He_3[j] = input[j * N + i]; 
    }
    fprintf(stderr, "Writing solution for /He_3\n");
    H5LTmake_dataset_double(file_id, "/He_3", 1, dimsarr, He_3);
    i++;
    
    double de[dims];
    for (j = 0; j < dims; j++) {
        de[j] = input[j * N + i]; 
    }
    fprintf(stderr, "Writing solution for /de\n");
    H5LTmake_dataset_double(file_id, "/de", 1, dimsarr, de);
    i++;
    
    double ge[dims];
    for (j = 0; j < dims; j++) {
        ge[j] = input[j * N + i]; 
    }
    fprintf(stderr, "Writing solution for /ge\n");
    H5LTmake_dataset_double(file_id, "/ge", 1, dimsarr, ge);
    i++;
    


    H5LTmake_dataset_double(file_id, "/T", 1, dimsarr, temp);

    double time[1];
    time[0] = ttot;
    double timestep[1];
    timestep[0] = dt;
    H5LTset_attribute_double(file_id, "/", "time", time, 1); 
    H5LTset_attribute_double(file_id, "/", "timestep", timestep, 1);
    H5Fclose(file_id);
    
    free(temp);
    free(tics);
    free(ics);
    free(data);
    free(rtol);
    free(atol);
    free(input);
    free(density);

    return 0;
}
 



int dengo_evolve_primordial (double dtf, double &dt, double z, double *input,
            double *rtol, double *atol, unsigned long dims, primordial_data *data, double *temp_array ){

    //-----------------------------------------------------
    // Function     : dengo_evolve_primordial
    // Description  : Main ODE solver function in dengo
    
    // Parameter    :   dtf     : Desired time to be reached by the solver
    //                  dt      : Pointer to the actual time reached by the solver
    //                  z       : Current redshift
    //                  input   : Array to store the initial value of each species, 
    //                            it will be updated with the value at time dt
    //                  rtol    : relative tolerance required for convergenece in each internal CVODE timesteps
    //                  atol    : absolute tolerance required for convergence in each interanl CVODE timesteps
    //                  dims    : dimension of the input array, i.e. no. of species * no. of cells
    //                  data    : primordial_data object that relay the reaction/cooling rates, and normalizations 
    //                  temp_array: temperature of each cell by the end of the evolution
    //                           
    //-----------------------------------------------------
    unsigned long i, j;
    hid_t file_id;
    /* fprintf(stderr, "  ncells = % 3i\n", (int) dims); */
    int N = 10;

    data->reltol = rtol[0];
    double H2_1;
    double H2_2;
    double H_1;
    double H_2;
    double H_m0;
    double He_1;
    double He_2;
    double He_3;
    double de;
    double ge;
    for (i = 0; i < dims; i++) {
        j = i * N;
        input[j] /= 2.0; // H2_1
        j++;
        input[j] /= 2.0; // H2_2
        j++;
        input[j] /= 1.00794; // H_1
        j++;
        input[j] /= 1.00794; // H_2
        j++;
        input[j] /= 1.00794; // H_m0
        j++;
        input[j] /= 4.002602; // He_1
        j++;
        input[j] /= 4.002602; // He_2
        j++;
        input[j] /= 4.002602; // He_3
        j++;
        input[j] /= 1.0; // de
        j++;
        j++;
    }
    
    // when partial equilibrium solver is used
    // the equilibrium abundance is passed through the temp_array
    double *equil_array = &temp_array[dims];

    // TODO:
    // Need to consider the cases where 
    // the incoming array assumes some equilibrium species,
    // - should first calculate temperature, then rates
    // - then equilibrium species
    // - then electron conservation
    // calculate_equilibrium_abundance(data, input, dims, 0, dims, equil_array);
    // ensure_electron_consistency(input, equil_array, dims, N);

    rhs_f f = calculate_rhs_primordial;

    #ifndef CVSPILS
    #ifdef  CVKLU
    jac_f jf = calculate_sparse_jacobian_primordial;
    #else
    jac_f jf = calculate_jacobian_primordial;
    #endif
    #endif
    
    #ifdef CVSPILS
    jac_f jf = calculate_JacTimesVec_primordial;
    #endif

    if (dt < 0) dt = dtf / 1e0;
    data->current_z = z;
    int niter = 0;
    int siter = 0;

    double floor_value = 1e-20;

    // Initialize a CVODE object, memory spaces
    // and attach rhs, jac to them
    int flag;
    double reltol = rtol[0];
    void *cvode_mem;
    int MAX_ITERATION = 1;
    double mh = 1.66054e-24;
    int nstrip = data->nstrip;

    // inputs are now grouped into a batch of nstrip
    // and send into the CVode solver
    //
    // ntimes: the number of times the CVode solver will be called
    // v_size: the size of the input vector for the solver
    //       : v_size = N * nstrip
    // v_size_res: leftover strip that doesn't fit into a v_size batch
    // N     : Number of species
    
    int v_size      = N * nstrip;
    int nstrip_res  = dims % nstrip;
    int v_size_res  = N * nstrip_res;
    unsigned long ntimes      = dims / nstrip;
    
    SUNLinearSolver LS;
    SUNMatrix A;
    N_Vector y_vec, abstol;

    y_vec = NULL;   
    LS = NULL;
    A  = NULL;
    
    double *yvec_ptr;
    double *atol_ptr;

    // these objects should be initialize once !!
    // in each separate thread!!
    double *ttot = (double *) malloc( (ntimes + 1)* sizeof(double));
    int d;
    int sum;
    int NSPARSE = 64; // no. of sparse jacobian compoents
    
    int threadID;

    // this routine is also called from the cythonized routine too
    // where MAX_NCELLS > input strip length
    // to avoid error message, we will catch it here
    if (ntimes > 0){

    #pragma omp parallel private (A, LS, cvode_mem, threadID, y_vec, abstol) num_threads(NTHREADS)
    {

    y_vec  = N_VNew_Serial(v_size);
    abstol = N_VNew_Serial(v_size); 
    yvec_ptr = N_VGetArrayPointer(y_vec);
    atol_ptr = N_VGetArrayPointer(abstol);
    
    // Need to be initialized before feeding into A, LS
    for (i = 0; i < v_size; i++) {
            yvec_ptr[i]   = 1.0;
            atol_ptr[i]   = reltol;
    }
    
    #ifdef CVKLU
    A         = SUNSparseMatrix      ( v_size, v_size, nstrip * NSPARSE, CSR_MAT );
    LS        = SUNKLU( y_vec, A);
    #else
    A         = SUNDenseMatrix      ( v_size, v_size );
    LS        = SUNDenseLinearSolver( y_vec, A);
    #endif
    
    cvode_mem = setup_cvode_solver  ( f, jf, v_size , data, LS, A, y_vec, reltol, abstol);
    
    // d: d th path going into the solver
    #pragma omp for private (sum, i, d, siter, threadID) schedule(static, 1)
    for (int d = 0; d < ntimes; d++){
        #ifdef _OPENMP
        threadID = omp_get_thread_num();
        #else
        threadID = 0;
        #endif
        ttot[d] = evolve_in_batches( cvode_mem, y_vec, abstol, reltol, input, v_size, d, d*v_size, MAX_ITERATION, dtf, data );

        // re-calculate temperature at the final output
        // updated in the data->Ts[threadID]
        primordial_calculate_temperature(data, &input[d*v_size], nstrip, N);
        
        // fprintf(stderr, "%d th strip from thread %d = %0.5g\n", d, threadID, ttot[d]);
        if (ttot[d] < dtf){
            fprintf(stderr, "FAILED FROM thread %d at ttot[%d] = %0.5g \n", threadID, d, ttot[d]);    
        } 

        for ( i = 0; i < nstrip; i ++){
            temp_array[ d * nstrip + i ] = data->Ts[threadID][i];
        }
    calculate_equilibrium_abundance(data, &input[d*v_size], nstrip, d, dims, equil_array);

    } // for d dims loop
    
    CVodeFree(&cvode_mem);
    SUNLinSolFree(LS);
    SUNMatDestroy(A);
    #ifdef _OPENMP
    N_VDestroy(y_vec);
    N_VDestroy(abstol);
    #endif
    }
    } // if (ntimes > 0)

    if ( v_size_res > 0 ){
        #ifdef _OPENMP
        threadID = omp_get_thread_num();
        #else
        threadID = 0;
        #endif

        data->nstrip = nstrip_res;
        d = ntimes;
    
        y_vec  = N_VNew_Serial( v_size_res );
        abstol = N_VNew_Serial( v_size_res ); 

        yvec_ptr = N_VGetArrayPointer(y_vec);
        atol_ptr = N_VGetArrayPointer(abstol);

        for (i = 0; i < v_size_res; i++) {
            yvec_ptr[i]   = 1.0;
            atol_ptr[i]   = reltol;
        }
       
        #ifdef CVKLU
        A  = SUNSparseMatrix(v_size_res, v_size_res, nstrip_res * NSPARSE, CSR_MAT);
        LS = SUNKLU(y_vec,A); 
        #else
        A  = SUNDenseMatrix(v_size_res, v_size_res);
        LS = SUNDenseLinearSolver(y_vec, A);
        #endif
        
        cvode_mem = setup_cvode_solver( f, jf, v_size_res, data, LS, A, y_vec, reltol, abstol );
        ttot[d] = evolve_in_batches(cvode_mem, y_vec, abstol, reltol, input, v_size_res, d, d * v_size,  MAX_ITERATION, dtf, data);
        primordial_calculate_temperature(data, &input[d*v_size], nstrip_res, N);
        for ( i = 0; i < nstrip_res; i ++){
            temp_array[ d * nstrip + i ] = data->Ts[threadID][i];
        }
        
        
        calculate_equilibrium_abundance(data, &input[d*v_size], nstrip_res, d, dims, equil_array);
        //fprintf(stderr, "%d th strip = %0.5g\n", d, ttot[d]);
        
        // free all unused memory;
        CVodeFree(&cvode_mem);
        SUNLinSolFree(LS);
        SUNMatDestroy(A);
        N_VDestroy(y_vec);
        N_VDestroy(abstol);    
    } else{
      ttot[ntimes] = dtf;      
    }

    // inputs are in `number density`
    for (i = 0; i < dims; i++) {
      j = i * N;
      H2_1 = input[j];
      input[j] *= 2.0; // H2_1
      j++;
      H2_2 = input[j];
      input[j] *= 2.0; // H2_2
      j++;
      H_1 = input[j];
      input[j] *= 1.00794; // H_1
      j++;
      H_2 = input[j];
      input[j] *= 1.00794; // H_2
      j++;
      H_m0 = input[j];
      input[j] *= 1.00794; // H_m0
      j++;
      He_1 = input[j];
      input[j] *= 4.002602; // He_1
      j++;
      He_2 = input[j];
      input[j] *= 4.002602; // He_2
      j++;
      He_3 = input[j];
      input[j] *= 4.002602; // He_3
      j++;
      de = input[j];
      input[j] *= 1.0; // de
      j++;
      j++;
    }

    double dt_final = dtf;
    
    for (int d = 0; d < (ntimes + 1); d++){
        if (ttot[d] < dt_final) dt_final = ttot[d];    
    }

    // fprintf(stderr, "Fraction of completion (dt (%0.3g) / dtf (%0.3g)): %0.3g\n", dt, dt_final, dt_final/dtf);
    free(ttot);

    dt = dt_final;
    if (dt_final < dtf) return 1;
    return 0;

}
 



double evolve_in_batches( void * cvode_mem, N_Vector y_vec, N_Vector abstol,  
                          double reltol, double *input, int v_size, int d, int start_idx, 
                          int MAX_ITERATION, double dtf, primordial_data *data ){ 
    // Function   :     evolve_in_batches
    // Description:     this would evolve the ODE system in bataches of size v_size
    //                  and return the final time reached by the solver
    //
    // Parameter  :     cvode_mem   : CVODE memory object
    //                  y_vec           : Array to store/relay the data to the CVODE solver                  
    //                  abstol          : Array of absolute tolerance for each internal timestep
    //                  reltol          : Relative tolerance requested for each internal timestep
    //                  input           : Array to store the input number density/ energy of the species
    //                  v_size          : Size of the ODE system, i.e. No. of Species * No. of Strips/cells
    //                  d               : Batch count
    //                  start_idx       : Index to the first element from "input" for the current batch
    //                  MAX_ITERATION   : Maximum Call/retry to the CVODE solver
    //                  dtf             : Desired final output time
    //                  data            : primordial_data object that stores rate/ temperature/ scale data
    // Return    :      ttot            : final time reached by the solver

    int i, siter, flag;
    double dt, ttot;
    double y[v_size];
    int nstrip = data->nstrip;

    #ifdef _OPENMP
    int threadID = omp_get_thread_num();
    #else
    int threadID = 0;
    #endif

    // access the array pointer of sundials vector object
    double *yvec_ptr = N_VGetArrayPointer(y_vec);
    double *atol_ptr = N_VGetArrayPointer(abstol);

    // by default: scale the input array
    #ifdef SCALE_INPUT
    double *scale = data->scale[threadID];
    double *inv_scale = data->inv_scale[threadID];
    for (i = 0; i < v_size; i++){ 
        scale[i]      = input[ start_idx + i];
        inv_scale[i]  = 1.0 / scale[i];
        yvec_ptr[i]   = 1.0;
        atol_ptr[i]   = reltol*reltol; //atol_ptr[i]*inv_scale[i];
    }
    setting_up_extra_variables(data, data->scale[threadID], nstrip );
    #else
    for (i = 0; i < v_size; i++){ 
        yvec_ptr[i]   = input[ start_idx + i];
        atol_ptr[i]   = reltol*reltol*yvec_ptr[i]; //atol_ptr[i]*inv_scale[i];
    }
    setting_up_extra_variables(data, yvec_ptr, nstrip );
    #endif

    // initialize a dt for the solver  
    dt = dtf;
    ttot = 0.0;
    siter = 0;
            
    while (ttot < dtf) { 
        // fprintf(stderr, "%d th strip: %d iterations, time: %0.5g\n", d, siter, ttot );    
        flag = cvode_solver( cvode_mem, y, v_size , &dt, data, y_vec, reltol, abstol);

        for (i = 0; i < v_size; i++) {
            if (y[i] < 0) {
            // this catches negatives values smaller than the abstol
            // and replaces them
            // need to check the convergence of this approach
            if (y[i] + atol_ptr[i] > 0 ){
                y[i] = atol_ptr[i]; 
            } else{
                fprintf(stderr, "negative \n");
            flag = 1;
                    break;
	    }
            }
        }
            
        if (flag < 1){
            // flag = 0 => success
            // we reset the scale of each component 
            // with the solution returned by the solver
            #ifdef SCALE_INPUT
            for (i = 0; i < v_size ; i++){
                yvec_ptr[i]   = 1.0;
        	//atol_ptr[i]   = fmin(reltol, atol_ptr[i]*inv_scale[i]);
                scale[i]     = y[i] * scale[i];
                inv_scale[i] = 1.0 /  scale[i];
            }
            #else
            for (i = 0; i < v_size ; i++){
                yvec_ptr[i]   = y[i];
            }
            #endif
            ttot += dt;
            dt = DMIN(dt * 2.0, dtf - ttot);
        } else{
            dt /= 2.0;
            dt = DMIN(dt , dtf - ttot);
        }
            if (siter == MAX_ITERATION) break;
            siter++;
    } // while loop for each strip
   
    // copy the results back to the input array
    // regardless the SCALE_INPUT
    // the returned array is still in number density
    // i.e. not scaled
    #ifdef SCALE_INPUT
    for (i = 0; i < v_size; i++){ 
        input[ start_idx + i] = scale[i];
    }
    #else
    for (i = 0; i < v_size; i++){ 
        input[ start_idx + i] = y[i];
    }
    #endif


    return ttot;
}


///////////////////////////////////////////////////////////////////////////////
//////////// Evaluate Temperature /////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////


int primordial_calculate_temperature(primordial_data *data,
                        double *input, int nstrip, int nchem)
{
    int i, j;
    double density, T, Tnew;
    double kb = 1.3806504e-16; // Boltzmann constant [erg/K]
    double mh = 1.66054e-24;
    double gamma = 5.e0/3.e0;
    double _gamma_m1 = 1.0 / (gamma - 1);

    double dge_dT;
    double gammaH2 = 7.e0/5.e0; // Should be a function of temperature
    double dge;
    double gammaH2_1;
    double dgammaH2_1_dT;
    double _gammaH2_1_m1;
    double gammaH2_2;
    double dgammaH2_2_dT;
    double _gammaH2_2_m1;
        
    double Tdiff = 1.0;
    double reltol = data->reltol;
    int MAX_T_ITERATION = 100;
    int count = 0;

    /* Calculate total density */
    double H2_1;
    double H2_2;
    double H_1;
    double H_2;
    double H_m0;
    double He_1;
    double He_2;
    double He_3;
    double de;
    double ge;
    
    i = 0;
    #ifdef _OPENMP
    int threadID = omp_get_thread_num();
    #else 
    int threadID = 0;
    #endif

    for ( i = 0; i < nstrip; i++ ){
        j = i * nchem;
        H2_1 = input[j];
        j++;
        H2_2 = input[j];
        j++;
        H_1 = input[j];
        j++;
        H_2 = input[j];
        j++;
        H_m0 = input[j];
        j++;
        He_1 = input[j];
        j++;
        He_2 = input[j];
        j++;
        He_3 = input[j];
        j++;
        de = input[j];
        j++;
        ge = input[j];
        j++;
    
        /*
        */
	
        // TODO: pull the rates from primordial_data
        // these species usually contribute negligbly to the number density (?)
	// it is a little tricky here,
	// since these species depends on the temperature 
	// and the abundance of the rest of the species
	// BUT, without their abundance, temperature CANNOT be evaluated....
	// FOR NOW, a not entirely correct physically, 
	// BUT a not-too-bad surrogate is:
	// assume these species has negligible abundance....

        density = 2.0*H2_1 + 2.0*H2_2 + 1.0079400000000001*H_1 + 1.0079400000000001*H_2 + 1.0079400000000001*H_m0 + 4.0026020000000004*He_1 + 4.0026020000000004*He_2 + 4.0026020000000004*He_3;
        
        // Requires iteration on the convergence of temperature
        // since the gammaH2 is not fixed
        // Initiate the "guess" temperature
        T    = data->Ts[threadID][i];
        Tnew = T*1.1;
        Tdiff = Tnew - T;
        count = 0;

        while ( Tdiff/ Tnew > 1.0e-8){
            // We do Newton's Iteration to calculate the temperature
            // Since gammaH2 is dependent on the temperature too!

            T = data->Ts[threadID][i];
        
            primordial_interpolate_gamma(data, i);
            
            gammaH2_1 = data->gammaH2_1[threadID][i];
            dgammaH2_1_dT = data->dgammaH2_1_dT[threadID][i];
            _gammaH2_1_m1 = 1.0 / (gammaH2_1 - 1.0);
            // fprintf(stderr, ":gammaSpecies: H2_1 %0.5g , dgammaSpecies: H2_1_dT: %.5g \n", gammaSpecies: H2_1, dgammaSpecies: H2_1_dT  );
            
            gammaH2_2 = data->gammaH2_2[threadID][i];
            dgammaH2_2_dT = data->dgammaH2_2_dT[threadID][i];
            _gammaH2_2_m1 = 1.0 / (gammaH2_2 - 1.0);
            // fprintf(stderr, ":gammaSpecies: H2_2 %0.5g , dgammaSpecies: H2_2_dT: %.5g \n", gammaSpecies: H2_2, dgammaSpecies: H2_2_dT  );
            
       
        
            // update gammaH2
            // The derivatives of  sum (nkT/(gamma - 1)/mh/density) - ge
            // This is the function we want to minimize
            // which should only be dependent on the first part
            dge_dT = T*kb*(-H2_1*_gammaH2_1_m1*_gammaH2_1_m1*dgammaH2_1_dT - H2_2*_gammaH2_2_m1*_gammaH2_2_m1*dgammaH2_2_dT)/(density*mh) + kb*(H2_1*_gammaH2_1_m1 + H2_2*_gammaH2_2_m1 + H_1*_gamma_m1 + H_2*_gamma_m1 + H_m0*_gamma_m1 + He_1*_gamma_m1 + He_2*_gamma_m1 + He_3*_gamma_m1 + _gamma_m1*de)/(density*mh);
        
            //This is the change in ge for each iteration
            dge = T*kb*(H2_1*_gammaH2_1_m1 + H2_2*_gammaH2_2_m1 + H_1*_gamma_m1 + H_2*_gamma_m1 + H_m0*_gamma_m1 + He_1*_gamma_m1 + He_2*_gamma_m1 + He_3*_gamma_m1 + _gamma_m1*de)/(density*mh) - ge;

            Tnew = T - dge/dge_dT;
            data->Ts[threadID][i] = Tnew;
        
            Tdiff = fabs(T - Tnew);
            // fprintf(stderr, "T: %0.5g ; Tnew: %0.5g; dge_dT: %.5g, dge: %.5g, ge: %.5g \n", T,Tnew, dge_dT, dge, ge);
            count += 1;
            if (count > MAX_T_ITERATION){
                fprintf(stderr, "T failed to converge \n");
                return 1;
            }
        } // while loop
    
        data->Ts[threadID][i] = Tnew;


        //fprintf(stderr,"T : %0.5g, density : %0.5g, d_gammaH2: %0.5g \n", Tnew, density, gammaH2 - 7./5.);


        

        if (data->Ts[threadID][i] < data->bounds[0]) {
            data->Ts[threadID][i] = data->bounds[0];
        } else if (data->Ts[threadID][i] > data->bounds[1]) {
            data->Ts[threadID][i] = data->bounds[1];
        }
        data->logTs[threadID][i] = log(data->Ts[threadID][i]);
        data->invTs[threadID][i] = 1.0 / data->Ts[threadID][i];

        dge_dT = T*kb*(-H2_1*_gammaH2_1_m1*_gammaH2_1_m1*dgammaH2_1_dT - H2_2*_gammaH2_2_m1*_gammaH2_2_m1*dgammaH2_2_dT)/(density*mh) + kb*(H2_1*_gammaH2_1_m1 + H2_2*_gammaH2_2_m1 + H_1*_gamma_m1 + H_2*_gamma_m1 + H_m0*_gamma_m1 + He_1*_gamma_m1 + He_2*_gamma_m1 + He_3*_gamma_m1 + _gamma_m1*de)/(density*mh);
        data->dTs_ge[threadID][i] = 1.0 / dge_dT;
    } // for i in nstrip loop
    return 0;
         
}
 



///////////////////////////////////////////////////////////////////////////////
////////// Ensure Conservation of mass, charge, species ///////////////////////
///////////////////////////////////////////////////////////////////////////////


void ensure_electron_consistency(double *input, double *equil_array, unsigned long nstrip, int nchem)
{
    // inputs are assumed to be in number density
    unsigned long i, j;

    /* Now we set up some temporaries */
    double H2_1;
    double H2_2;
    double H_1;
    double H_2;
    double H_m0;
    double He_1;
    double He_2;
    double He_3;
    double de;
    double ge;

    double total_e = 0.0;
    int e_indx;

    for (i = 0; i<nstrip; i++) {
	total_e = 0.0;
        j = i * nchem;
        H2_1 = input[j];
        
        total_e += H2_1 * 0.0;
        j++;
        H2_2 = input[j];
        
        total_e += H2_2 * 1.0;
        j++;
        H_1 = input[j];
        
        total_e += H_1 * 0.0;
        j++;
        H_2 = input[j];
        
        total_e += H_2 * 1.0;
        j++;
        H_m0 = input[j];
        
        total_e += H_m0 * -1.0;
        j++;
        He_1 = input[j];
        
        total_e += He_1 * 0.0;
        j++;
        He_2 = input[j];
        
        total_e += He_2 * 1.0;
        j++;
        He_3 = input[j];
        
        total_e += He_3 * 2.0;
        j++;
        de = input[j];
        e_indx = j;
        j++;
        ge = input[j];
        
        j++;
        input[e_indx] = total_e;
    }
}



///////////////////////////////////////////////////////////////////////////////
//////////////////// Auxillary Functions that estimate state variables ////////
///////////////////////////////////////////////////////////////////////////////
// Again these are exposed to the user
// and can be called basically with units, field_data objects
///////////////////////////////////////////////////////////////////////////////

int dengo_estimate_cooling_time_enzo( code_units* units, dengo_field_data *field_data ){
    

    unsigned long int i, j, k, d, dims;
    int nchem = 10;
    
    // calculate total number of cells 
    int is, ie, js, je, ks, ke;
    // avoid ghost zones
    is = field_data->grid_start[0];
    ie = field_data->grid_end[0];

    js = field_data->grid_start[1];
    je = field_data->grid_end[1];

    ks = field_data->grid_start[2];
    ke = field_data->grid_end[2];

    // number of cells that actually required calculations
    dims = ie - is + 1;
    dims*= je - js + 1;
    dims*= ke - ks + 1;

    const char * FileLocation = field_data->dengo_data_file;
    primordial_data *data = primordial_setup_data( FileLocation, NULL, NULL);

    // fills the `input` array with field data
    // note that input are still some kind of mass density, rho/ m_amu  
    double *input = (double *) malloc(dims*nchem*sizeof(double));
    flatten_dengo_field_data_enzo(units, field_data, input);

    int nstrip      = data->nstrip;
    unsigned long ntimes      = dims / nstrip;
    int nstrip_res            = dims % nstrip;
    
    double *input_batch;
    double *cooling_time_batch;

    // update the redshift
    double a = units->a_value * units->a_units;
    double z = 1./a - 1.;
    data->current_z = z;

    #pragma omp parallel for private (i, j ,d, input_batch, cooling_time_batch) num_threads(NTHREADS) schedule(static, 1) 
    for ( d = 0; d < ntimes; d++ ){
        input_batch        = &input[d* nstrip* nchem];
        cooling_time_batch = &field_data->CoolingTime[d * nstrip];
        primordial_calculate_cooling_timescale( cooling_time_batch, input_batch, nstrip, data);
    }

    if (nstrip_res > 0){
        input_batch        = &input[ntimes* nstrip * nchem];
        cooling_time_batch = &field_data->CoolingTime[ntimes*nstrip]; 
        primordial_calculate_cooling_timescale( cooling_time_batch, input_batch, nstrip_res, data );    
    }

    for (i = 0; i < dims; i++ ){
        field_data->CoolingTime[i] /= units->time_units; 
    }
    
    free(input);
    free(data);

    // in the world of Enzo
    //FAIL = 0
    return 1;
}



int dengo_estimate_cooling_time( code_units* units, dengo_field_data *field_data ){
    
    unsigned long int i, j, d, dims;
    int nchem     = 10;
    dims = field_data->ncells;
    
    const char * FileLocation = field_data->dengo_data_file;
    primordial_data *data = primordial_setup_data( FileLocation, NULL, NULL);
    
    int nstrip      = data->nstrip;

    unsigned long ntimes      = dims / nstrip;
    int nstrip_res            = dims % nstrip;
    
    // flatten field_data entry to a 1d input array
    double *input = (double *) malloc( sizeof(double) * nchem * dims );
    flatten_dengo_field_data( units, field_data, input );
    
    double *input_batch;
    double *cooling_time_batch;
    
    #pragma omp parallel for private (i, j ,d, input_batch, cooling_time_batch) num_threads(NTHREADS) schedule(static, 1) 
    for ( d = 0; d < ntimes; d++ ){
        input_batch        = &input[d* nstrip *nchem];
        cooling_time_batch = &field_data->CoolingTime[d * nstrip];
        primordial_calculate_cooling_timescale( cooling_time_batch, input_batch, nstrip, data);
    }

    if (nstrip_res > 0){
        input_batch        = &input[ntimes* nstrip* nchem];
        cooling_time_batch = &field_data->CoolingTime[ntimes*nstrip]; 
        primordial_calculate_cooling_timescale( cooling_time_batch, input_batch, nstrip_res, data );    
    }

    for (i = 0; i < dims; i++ ){
        field_data->CoolingTime[i] /= units->time_units; 
    }
    
    free(input);
    free(data);
}




int primordial_calculate_cooling_timescale( double *cooling_time, double *input, int nstrip, primordial_data *data){
    
    #ifdef _OPENMP
    int threadID = omp_get_thread_num(); 
    #else
    int threadID = 0;
    #endif

    unsigned long i, j, dims;
    int flag;
    int nchem = 10;
    /* Now we set up some temporaries */
    // make sure the input are in number density
    for (i = 0; i < nstrip; i++) {
        j = i * nchem;
        input[j] /= 2.0; // H2_1
        j++;
        input[j] /= 2.0; // H2_2
        j++;
        input[j] /= 1.00794; // H_1
        j++;
        input[j] /= 1.00794; // H_2
        j++;
        input[j] /= 1.00794; // H_m0
        j++;
        input[j] /= 4.002602; // He_1
        j++;
        input[j] /= 4.002602; // He_2
        j++;
        input[j] /= 4.002602; // He_3
        j++;
        input[j] /= 1.0; // de
        j++;
        j++;
    }
    
    // calculate temperature and cooling rate
    setting_up_extra_variables(data, input, nstrip );
    flag = primordial_calculate_temperature(data,  input , nstrip, nchem );
    if (flag > 0){
        // check if the temperature failed to converged
        return -1;    
    }
    primordial_interpolate_rates(data, nstrip);
    double H_1;
    double He_2;
    double He_1;
    double ge;
    double H_m0;
    double H2_2;
    double de;
    double He_3;
    double H_2;
    double H2_1;

    // this might be redundant
    // should only select relavant rates
    double *k01 = data->rs_k01[threadID];
    double *k02 = data->rs_k02[threadID];
    double *k03 = data->rs_k03[threadID];
    double *k04 = data->rs_k04[threadID];
    double *k05 = data->rs_k05[threadID];
    double *k06 = data->rs_k06[threadID];
    double *k07 = data->rs_k07[threadID];
    double *k08 = data->rs_k08[threadID];
    double *k09 = data->rs_k09[threadID];
    double *k10 = data->rs_k10[threadID];
    double *k11 = data->rs_k11[threadID];
    double *k12 = data->rs_k12[threadID];
    double *k13 = data->rs_k13[threadID];
    double *k14 = data->rs_k14[threadID];
    double *k15 = data->rs_k15[threadID];
    double *k16 = data->rs_k16[threadID];
    double *k17 = data->rs_k17[threadID];
    double *k18 = data->rs_k18[threadID];
    double *k19 = data->rs_k19[threadID];
    double *k21 = data->rs_k21[threadID];
    double *k22 = data->rs_k22[threadID];
    double *k23 = data->rs_k23[threadID];
    double *brem_brem = data->cs_brem_brem[threadID];
    double *ceHeI_ceHeI = data->cs_ceHeI_ceHeI[threadID];
    double *ceHeII_ceHeII = data->cs_ceHeII_ceHeII[threadID];
    double *ceHI_ceHI = data->cs_ceHI_ceHI[threadID];
    double *cie_cooling_cieco = data->cs_cie_cooling_cieco[threadID];
    double *ciHeI_ciHeI = data->cs_ciHeI_ciHeI[threadID];
    double *ciHeII_ciHeII = data->cs_ciHeII_ciHeII[threadID];
    double *ciHeIS_ciHeIS = data->cs_ciHeIS_ciHeIS[threadID];
    double *ciHI_ciHI = data->cs_ciHI_ciHI[threadID];
    double *compton_comp_ = data->cs_compton_comp_[threadID];
    double *gammah_gammah = data->cs_gammah_gammah[threadID];
    double *gloverabel08_gael = data->cs_gloverabel08_gael[threadID];
    double *gloverabel08_gaH2 = data->cs_gloverabel08_gaH2[threadID];
    double *gloverabel08_gaHe = data->cs_gloverabel08_gaHe[threadID];
    double *gloverabel08_gaHI = data->cs_gloverabel08_gaHI[threadID];
    double *gloverabel08_gaHp = data->cs_gloverabel08_gaHp[threadID];
    double *gloverabel08_gphdl = data->cs_gloverabel08_gphdl[threadID];
    double *gloverabel08_gpldl = data->cs_gloverabel08_gpldl[threadID];
    double *gloverabel08_h2lte = data->cs_gloverabel08_h2lte[threadID];
    double *h2formation_h2mcool = data->cs_h2formation_h2mcool[threadID];
    double *h2formation_h2mheat = data->cs_h2formation_h2mheat[threadID];
    double *h2formation_ncrd1 = data->cs_h2formation_ncrd1[threadID];
    double *h2formation_ncrd2 = data->cs_h2formation_ncrd2[threadID];
    double *h2formation_ncrn = data->cs_h2formation_ncrn[threadID];
    double *reHeII1_reHeII1 = data->cs_reHeII1_reHeII1[threadID];
    double *reHeII2_reHeII2 = data->cs_reHeII2_reHeII2[threadID];
    double *reHeIII_reHeIII = data->cs_reHeIII_reHeIII[threadID];
    double *reHII_reHII = data->cs_reHII_reHII[threadID];
    
    
    double h2_optical_depth_approx;    
    
    
    double cie_optical_depth_approx;
    
    
    double z;
    double T;
    double mdensity, inv_mdensity, dge_dt;

    for ( i = 0; i < nstrip; i++ ){
        
        T            = data->Ts[threadID][i];
        z            = data->current_z;
        mdensity     = data->mdensity[threadID][i];
        inv_mdensity = data->inv_mdensity[threadID][i];
        
        h2_optical_depth_approx = data->h2_optical_depth_approx[threadID][i];
        
        
        cie_optical_depth_approx = data->cie_optical_depth_approx[threadID][i];
        

        j = i * nchem;
        H2_1 = input[j];
        j++;
        H2_2 = input[j];
        j++;
        H_1 = input[j];
        j++;
        H_2 = input[j];
        j++;
        H_m0 = input[j];
        j++;
        He_1 = input[j];
        j++;
        He_2 = input[j];
        j++;
        He_3 = input[j];
        j++;
        de = input[j];
        j++;
        ge = input[j];
        j++;
   	
        // obtain a quasis-equilibrium estimate

        //
        // Species: ge
        //
        dge_dt = -brem_brem[i]*cie_optical_depth_approx*de*(H_2 + He_2 + 4.0*He_3) - ceHI_ceHI[i]*H_1*cie_optical_depth_approx*de - ceHeII_ceHeII[i]*He_2*cie_optical_depth_approx*de - ceHeI_ceHeI[i]*He_2*cie_optical_depth_approx*pow(de, 2) - ciHI_ciHI[i]*H_1*cie_optical_depth_approx*de - ciHeII_ciHeII[i]*He_2*cie_optical_depth_approx*de - ciHeIS_ciHeIS[i]*He_2*cie_optical_depth_approx*pow(de, 2) - ciHeI_ciHeI[i]*He_1*cie_optical_depth_approx*de - 2.0158800000000001*cie_cooling_cieco[i]*H2_1*cie_optical_depth_approx*mdensity - compton_comp_[i]*cie_optical_depth_approx*de*pow(z + 1.0, 4)*(T - 2.73*z - 2.73) - gloverabel08_h2lte[i]*H2_1*cie_optical_depth_approx*h2_optical_depth_approx/(gloverabel08_h2lte[i]/(gloverabel08_gaH2[i]*H2_1 + gloverabel08_gaHI[i]*H_1 + gloverabel08_gaHe[i]*He_1 + gloverabel08_gaHp[i]*H_2 + gloverabel08_gael[i]*de) + 1.0) - reHII_reHII[i]*H_2*cie_optical_depth_approx*de - reHeII1_reHeII1[i]*He_2*cie_optical_depth_approx*de - reHeII2_reHeII2[i]*He_2*cie_optical_depth_approx*de - reHeIII_reHeIII[i]*He_3*cie_optical_depth_approx*de + 0.5*1.0/(h2formation_ncrn[i]/(h2formation_ncrd1[i]*H_1 + h2formation_ncrd2[i]*H2_1) + 1.0)*(-h2formation_h2mcool[i]*H2_1*H_1 + h2formation_h2mheat[i]*pow(H_1, 3));
        dge_dt *= inv_mdensity;
        cooling_time[i] = fabs( ge / dge_dt);
    
    //fprintf(stderr, "----------------\n");
    }
}







int dengo_calculate_temperature( code_units *units, dengo_field_data *field_data){
    
    unsigned long int i, j, d, dims;
    int nchem     = 10;
    dims = field_data->ncells;
    
    const char * FileLocation = field_data->dengo_data_file;
    primordial_data *data = primordial_setup_data( FileLocation, NULL, NULL);
    
    int nstrip      = data->nstrip;
    unsigned long ntimes      = dims / nstrip;
    int nstrip_res            = dims % nstrip;
    
    // flatten field_data entry to a 1d input array
    double *input = (double *) malloc( sizeof(double) * nchem * dims );
    flatten_dengo_field_data( units, field_data, input );

    double *input_batch;
    double *cooling_time_batch;
    double *gamma_eff_batch;

    int threadID; 
    /* Now we set up some temporaries */

    #pragma omp parallel for private (i, j ,d, input_batch, gamma_eff_batch, threadID) num_threads(NTHREADS) schedule(static, 1)
    for ( d = 0; d < ntimes; d++ ){
        #ifdef _OPENMP
        threadID = omp_get_thread_num();
        #else
        threadID = 0;
        #endif
        input_batch        = &input[d* nchem];
        gamma_eff_batch    = &field_data->Gamma[d*nstrip];

        setting_up_extra_variables( data, input_batch, nstrip );
        primordial_calculate_temperature(data,  input_batch , nstrip, nchem );
        dengo_calculate_gamma( gamma_eff_batch, data, input_batch, nstrip);
        for ( i = 0; i < nstrip; i++ ){
            field_data->temperature[d*nstrip + i] = data->Ts[threadID][i]; 
        }
    }

    if (nstrip_res > 0){
        input_batch        = &input[ntimes * nchem];
        gamma_eff_batch    = &field_data->Gamma[ntimes*nstrip];
        
        setting_up_extra_variables( data, input_batch, nstrip_res );
        primordial_calculate_temperature(data,  input_batch , nstrip_res , nchem );
        dengo_calculate_gamma( gamma_eff_batch, data, input_batch, nstrip_res);
       
        for ( i = 0; i < nstrip_res; i++ ){
            field_data->temperature[ntimes*nstrip + i] = data->Ts[threadID][i]; 
        }
    }
    
    free(input);
    free(data);
}





int dengo_calculate_gamma( double* gamma_eff, primordial_data *data, double *input, int nstrip  ){
    unsigned long int i, j, d, dims;
    int nchem     = 10;
    
    
    double H_1;
    
    double He_2;
    
    double He_1;
    
    double ge;
    
    double H_m0;
    
    double H2_2;
    
    double de;
    
    double He_3;
    
    double H_2;
    
    double H2_1;
    
    
    double gamma = 5.0/3.0;
    
    double gammaH2_1;
    double dgammaH2_1_dT;
    double _gammaH2_1_m1;
    
    double gammaH2_2;
    double dgammaH2_2_dT;
    double _gammaH2_2_m1;
    
    double n_gamma_m1, ndensity;
    
    #ifdef _OPENMP
    int threadID = omp_get_thread_num();
    #else
    int threadID = 0;
    #endif

    // both temperature and gamma are updated before calling this
    for (i = 0; i < nstrip; i++ ){
        ndensity = 0;
        j = i * nchem;
        H_1 = input[j] / 1.00794;
        ndensity += H_1;
        He_2 = input[j] / 4.002602;
        ndensity += He_2;
        He_1 = input[j] / 4.002602;
        ndensity += He_1;
        ge = input[j] / 1.0;
        H_m0 = input[j] / 1.00794;
        ndensity += H_m0;
        H2_2 = input[j] / 2.0;
        ndensity += H2_2;
        de = input[j] / 1.0;
        ndensity += de;
        He_3 = input[j] / 4.002602;
        ndensity += He_3;
        H_2 = input[j] / 1.00794;
        ndensity += H_2;
        H2_1 = input[j] / 2.0;
        ndensity += H2_1;
        gammaH2_1 = data->gammaH2_1[threadID][i];
        dgammaH2_1_dT = data->dgammaH2_1_dT[threadID][i];
        _gammaH2_1_m1 = 1.0 / (gammaH2_1 - 1.0);
        gammaH2_2 = data->gammaH2_2[threadID][i];
        dgammaH2_2_dT = data->dgammaH2_2_dT[threadID][i];
        _gammaH2_2_m1 = 1.0 / (gammaH2_2 - 1.0);
      
        n_gamma_m1 = H2_1/(gammaH2_1 - 1.0) + H2_2/(gammaH2_2 - 1.0) + H_1/(gamma - 1.0) + H_2/(gamma - 1.0) + H_m0/(gamma - 1.0) + He_1/(gamma - 1.0) + He_2/(gamma - 1.0) + He_3/(gamma - 1.0) + de/(gamma - 1.0);
        gamma_eff[i] = ndensity / n_gamma_m1 + 1.0; 
    }

}



int dengo_calculate_mean_molecular_weight( code_units *units, dengo_field_data *field_data ){
    unsigned long int i, j, d, dims;
    int nchem     = 10;
    dims = field_data->ncells;
    
    
    // flatten field_data entry to a 1d input array
    double *input = (double *) malloc( sizeof(double) * nchem * dims );
    flatten_dengo_field_data( units, field_data, input );
    double H_1;
    double He_2;
    double He_1;
    double ge;
    double H_m0;
    double H2_2;
    double de;
    double He_3;
    double H_2;
    double H2_1;
    double ndensity, mdensity;

    for (i = 0; i < dims; i++ ){
        ndensity = 0.0;
        mdensity = 0.0;
        j = i * nchem;
        // Species: H2_1 
        mdensity += input[j];
        ndensity += input[j] /2.0;
        j++;
        // Species: H2_2 
        mdensity += input[j];
        ndensity += input[j] /2.0;
        j++;
        // Species: H_1 
        mdensity += input[j];
        ndensity += input[j] /1.00794;
        j++;
        // Species: H_2 
        mdensity += input[j];
        ndensity += input[j] /1.00794;
        j++;
        // Species: H_m0 
        mdensity += input[j];
        ndensity += input[j] /1.00794;
        j++;
        // Species: He_1 
        mdensity += input[j];
        ndensity += input[j] /4.002602;
        j++;
        // Species: He_2 
        mdensity += input[j];
        ndensity += input[j] /4.002602;
        j++;
        // Species: He_3 
        mdensity += input[j];
        ndensity += input[j] /4.002602;
        j++;
        // Species: de
        ndensity += input[j] /1.0;
        j++;
        // Species: ge
        j++;
        field_data->MolecularWeight[i] = mdensity / ndensity;   
    }
    
    free(input);
}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////Solve Chemistry   ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// basically all the routines that are exposed to the user
// that can be called from dengo


int primordial_solve_chemistry( code_units *units, dengo_field_data *field_data, double dt ){
    // to use this, reltol and floor_value must be specified through 
    // dengo_field_data 
    unsigned long int i, j, d, dims;
    int N = 10;
    dims = field_data->ncells; // total number of strips to be evaluated
    
    const char * FileLocation = field_data->dengo_data_file;
    primordial_data *data = primordial_setup_data( FileLocation, NULL, NULL);

    double *input = (double *) malloc(dims * N * sizeof(double));
    double *temp  = (double *) malloc(dims * (0+1) * sizeof(double) );
    double *atol  = (double *) malloc(dims * N * sizeof(double));
    double *rtol  = (double *) malloc(sizeof(double));
    
    // code unit in terms of erg/g
    double UNIT_E_per_M = units->velocity_units * units->velocity_units;
    double m_amu = 1.66053904e-24;

    rtol[0]            = field_data->reltol;
    data->reltol       = field_data->reltol;
    double floor_value = field_data->floor_value;

    flatten_dengo_field_data(units, field_data, input);
    /*
    for (int i = 0; i < N; i++){
        fprintf(stderr, "input[%d] = %0.5g\n", i, input[i] );
    }
    */
    
    int flag;
    double z;
    double dtf;
    
    // convert code time to seconds 
    dt *= units->time_units;
    dtf = dt;
    
    // update the rate table location
    data->dengo_data_file = field_data->dengo_data_file; 

    flag = dengo_evolve_primordial(dtf, dt, z, input, rtol, atol, dims, data, temp);
    
    reshape_to_dengo_field_data(units, field_data, input);  
    
    for ( d = 0; d< dims; d++  ){
        field_data->temperature[d] = temp[d];
    }

    for ( d = 0; d< dims; d++  ){
    }

    //dengo_estimate_cooling_time( units, field_data );
    //dengo_calculate_temperature(units, field_data);
    //dengo_calculate_mean_molecular_weight( units, field_data );

    free(input);
    free(temp);
    free(data);
    free(atol);
    free(rtol);

    if (flag > 0) return 1;

    return 0;
}




int primordial_solve_chemistry_enzo( code_units *units, dengo_field_data *field_data, double dt ){
    //-------------------------------------------------------------------------
    // Function: primordial_solve_chemistry_enzo
    // Description: takes the same input as primordial_solve_chemistry
    //              BUT field_data needs to be specified with things like
    //              grid_start, grid_end, grid_dimension
    //              such that dengo can avoid evaluating ghost zones
    //-------------------------------------------------------------------------

    // to use this, reltol and floor_value must be specified through 
    // dengo_field_data 
    unsigned long int i, j, k, d, dims;
    int N = 10;
    
    // calculate total number of cells 
    int is, ie, js, je, ks, ke;
    // avoid ghost zones
    is = field_data->grid_start[0];
    ie = field_data->grid_end[0];

    js = field_data->grid_start[1];
    je = field_data->grid_end[1];

    ks = field_data->grid_start[2];
    ke = field_data->grid_end[2];

    // number of cells that actually required calculations
    dims = ie - is + 1;
    dims*= je - js + 1;
    dims*= ke - ks + 1;

    const char * FileLocation = field_data->dengo_data_file;
    primordial_data *data = primordial_setup_data( FileLocation, NULL, NULL);

    double *input = (double *) malloc(dims * N * sizeof(double));
    double *temp  = (double *) malloc(dims * (0+1) * sizeof(double) );
    double *atol  = (double *) malloc(dims * N * sizeof(double));
    double *rtol  = (double *) malloc(sizeof(double));
    
// fills the `input` array with field data
    // note that input are still some kind of mass density, rho/ m_amu  
    flatten_dengo_field_data_enzo(units, field_data, input);

    // specify the relative tolerance and floor value
    rtol[0]            = field_data->reltol;
    data->reltol       = field_data->reltol;
    double floor_value = field_data->floor_value;
    
    // specifiy the redshift
    int flag;
    double a = units->a_value * units->a_units;
    double z = 1./a - 1.;
    double dtf;
    
    // convert code time to seconds 
    dt *= units->time_units;
    dtf = dt;
    
    // update the rate table location
    data->dengo_data_file = field_data->dengo_data_file; 
    data->floor_value = floor_value;
    
    // evolve chemistry in dengo
    flag = dengo_evolve_primordial(dtf, dt, z, input, rtol, atol, dims, data, temp);

    // fill results in `input` back to field data
    // in appropriate units
    reshape_to_dengo_field_data_enzo(units, field_data, input, temp);

    // free all pointers
    free(input);
    free(temp);
    free(data);
    free(atol);
    free(rtol);

    // in Enzo; 0 == FAIL
    if (flag > 0) return 0;
    // 1 == SUCCESS
    return 1;
}



int primordial_solve_chemistry_dt( code_units *units, dengo_field_data *field_data, double* reltol, double* abstol, double dt ){
    // TODO:
    // this is not only called from the python modules 
    // but this is dumb...
    // this is almost the same as primordial_solve_chemistry!!!
    // should be replaced later....

    unsigned long int i, j, d, dims;
    int N = 10;
    dims = field_data->ncells; // total number of strips to be evaluated

    // turned the field data into a long chain of 
    // 1-D array
    //
    // N: number of species
    // d: d th number of strip to be evalulated
    // i: index for each species
    //    should be in correct order and handled 
    //    by dengo templates
    // i.e.
    // input[d*N + i] = field_data->HI_density[]
    //
    
    const char * FileLocation = field_data->dengo_data_file;
    primordial_data *data = primordial_setup_data( FileLocation, NULL, NULL);

    double *input = (double *) malloc(dims * N * sizeof(double));
    double *temp  = (double *) malloc(dims * sizeof(double) );
    double *atol  = (double *) malloc(dims * N * sizeof(double));
    double *rtol  = (double *) malloc(sizeof(double));
    
    // code unit in terms of erg/g
    double UNIT_E_per_M = units->velocity_units * units->velocity_units;
    double m_amu = 1.66053904e-24;


    #pragma omp parallel for private (i, j ,d) num_threads(NTHREADS) schedule(static,1)
    for ( d = 0; d< dims; d++  ){
        j = d*N;
        // this should be the normalized 
        // by the input units later
        // atol = input * rtol;
        // which again should be set by dengo
        // input in *mass density per amu* 
        // and energy in the units of (erg / g)
        input[j]  = field_data->H2_1_density[d] ;
        input[j] *= units->density_units / m_amu;
        j++;
        input[j]  = field_data->H2_2_density[d] ;
        input[j] *= units->density_units / m_amu;
        j++;
        input[j]  = field_data->H_1_density[d] ;
        input[j] *= units->density_units / m_amu;
        j++;
        input[j]  = field_data->H_2_density[d] ;
        input[j] *= units->density_units / m_amu;
        j++;
        input[j]  = field_data->H_m0_density[d] ;
        input[j] *= units->density_units / m_amu;
        j++;
        input[j]  = field_data->He_1_density[d] ;
        input[j] *= units->density_units / m_amu;
        j++;
        input[j]  = field_data->He_2_density[d] ;
        input[j] *= units->density_units / m_amu;
        j++;
        input[j]  = field_data->He_3_density[d] ;
        input[j] *= units->density_units / m_amu;
        j++;
        input[j]  = field_data->de_density[d] ;
        input[j] *= units->density_units / m_amu;
        j++;
        input[j]  = field_data->ge_density[d] ;
        input[j] *= UNIT_E_per_M;
        j++;
    }

    /*
    for (int i = 0; i < N; i++){
        fprintf(stderr, "input[%d] = %0.5g\n", i, input[i] );
    }
    */
    
    int flag;
    double z;
    double dtf;
    
    // convert code time to seconds 
    dt *= units->time_units;
    dtf = dt;
    
    // update the rate table location
    data->dengo_data_file = field_data->dengo_data_file; 

    flag = dengo_evolve_primordial(dtf, dt, z, input, reltol, abstol, dims, data, temp);
    
    #pragma omp parallel for private (i, j ,d) num_threads(NTHREADS) schedule(static, 1)
    for ( d = 0; d< dims; d++  ){
        j = d*N;
        // this should be the normalized 
        // by the input units later
        // atol = input * rtol;
        // which again should be set by dengo
        field_data->H2_1_density[d] = input[j];
        field_data->H2_1_density[d] /= units->density_units / m_amu;
        j++;
        field_data->H2_2_density[d] = input[j];
        field_data->H2_2_density[d] /= units->density_units / m_amu;
        j++;
        field_data->H_1_density[d] = input[j];
        field_data->H_1_density[d] /= units->density_units / m_amu;
        j++;
        field_data->H_2_density[d] = input[j];
        field_data->H_2_density[d] /= units->density_units / m_amu;
        j++;
        field_data->H_m0_density[d] = input[j];
        field_data->H_m0_density[d] /= units->density_units / m_amu;
        j++;
        field_data->He_1_density[d] = input[j];
        field_data->He_1_density[d] /= units->density_units / m_amu;
        j++;
        field_data->He_2_density[d] = input[j];
        field_data->He_2_density[d] /= units->density_units / m_amu;
        j++;
        field_data->He_3_density[d] = input[j];
        field_data->He_3_density[d] /= units->density_units / m_amu;
        j++;
        field_data->de_density[d] = input[j];
        field_data->de_density[d] /= units->density_units / m_amu;
        j++;
        field_data->ge_density[d] = input[j];
        field_data->ge_density[d] /= UNIT_E_per_M;
        j++;
    field_data->temperature[d] = temp[d];
    }
    
    dengo_estimate_cooling_time( units, field_data );
    //dengo_calculate_temperature(units, field_data);
    dengo_calculate_mean_molecular_weight( units, field_data );


    free(input);
    free(temp);
    free(data);
    free(atol);
    free(rtol);

    if (flag > 0) return 1;

    return 0;
}



int calculate_equilibrium_abundance( primordial_data *data, double *input, 
    int nstrip, unsigned long d, unsigned long dims, double *equil_array){

    //-----------------------------------------------------
    // Function     : calculate_equilibrium_abundance
    // Parameter    :   
    //                  input   : Array to store the initial value of each species, 
    //                            it will be updated with the value at time dt
    //                  data    : primordial_data object that relay the reaction/cooling rates, and normalizations 
    //                  equil_array: temperature of each cell by the end of the evolution
    //                  d          : batch count  
    //                  dims       : total number of cells to be solved
    //-----------------------------------------------------
    
    // update the rates given the updated temperature from primordial_data
    primordial_interpolate_rates(data, nstrip);
    
    #ifdef _OPENMP
    int threadID = omp_get_thread_num();
    #else
    int threadID = 0;
    #endif

    int nchem = 10;
    unsigned long i, j;
    double H2_1;
    double H2_2;
    double H_1;
    double H_2;
    double H_m0;
    double He_1;
    double He_2;
    double He_3;
    double de;
    double ge;
    double *k01 = data->rs_k01[threadID];
    double *k02 = data->rs_k02[threadID];
    double *k03 = data->rs_k03[threadID];
    double *k04 = data->rs_k04[threadID];
    double *k05 = data->rs_k05[threadID];
    double *k06 = data->rs_k06[threadID];
    double *k07 = data->rs_k07[threadID];
    double *k08 = data->rs_k08[threadID];
    double *k09 = data->rs_k09[threadID];
    double *k10 = data->rs_k10[threadID];
    double *k11 = data->rs_k11[threadID];
    double *k12 = data->rs_k12[threadID];
    double *k13 = data->rs_k13[threadID];
    double *k14 = data->rs_k14[threadID];
    double *k15 = data->rs_k15[threadID];
    double *k16 = data->rs_k16[threadID];
    double *k17 = data->rs_k17[threadID];
    double *k18 = data->rs_k18[threadID];
    double *k19 = data->rs_k19[threadID];
    double *k21 = data->rs_k21[threadID];
    double *k22 = data->rs_k22[threadID];
    double *k23 = data->rs_k23[threadID];

    for ( i = 0; i < nstrip; i++ ){
        j = i * nchem;
        H2_1 = input[j];
        j++;
        H2_2 = input[j];
        j++;
        H_1 = input[j];
        j++;
        H_2 = input[j];
        j++;
        H_m0 = input[j];
        j++;
        He_1 = input[j];
        j++;
        He_2 = input[j];
        j++;
        He_3 = input[j];
        j++;
        de = input[j];
        j++;
        ge = input[j];
        j++;
    }
    return 0;
}




/*

int calculate_JacTimesVec_primordial
            (N_Vector v, N_Vector Jv, realtype t,
             N_Vector y, N_Vector fy,
             void *user_data, N_Vector tmp)
{
    // TODO:
    // as of now this is utterly useless, 
    // cos it runs even slower than the actual dense linear solver ...
    // BUT! kept in mind that autodiff could easily replace this 
    // but some linear call to rhs/ f evauluations O(n) time
    // but not O(n^2) i think...
    
    // We iterate over all of the rates
    // Calcuate temperature first
    int nstrip = 1;
    int nchem = 10;
    primordial_data *data = (primordial_data*)user_data; 
    

    int i, j;
    j = 0;

    // change N_Vector back to an array 
    double y_arr[ 10 ];
    y_arr[0] = Ith(y , 1);
    y_arr[1] = Ith(y , 2);
    y_arr[2] = Ith(y , 3);
    y_arr[3] = Ith(y , 4);
    y_arr[4] = Ith(y , 5);
    y_arr[5] = Ith(y , 6);
    y_arr[6] = Ith(y , 7);
    y_arr[7] = Ith(y , 8);
    y_arr[8] = Ith(y , 9);
    y_arr[9] = Ith(y , 10);
    // Abundances are scaled in the calculate temperature module
    int flag;
    flag = primordial_calculate_temperature(data, y_arr, nstrip, nchem);
    if (flag > 0){
        return 1;    
    }
    
    #ifdef _OPENMP
    int threadID = omp_get_thread_num();
    #else
    int threadID = 0;
    #endif

    primordial_interpolate_rates(data, nstrip);
    
    // Now We set up some temporaries
    double *Tge = data->dTs_ge[threadID];
    
    // Define the reaction rates
    double *k01 = data->rs_k01[threadID];
    double *rk01 = data->drs_k01[threadID];
    double *k02 = data->rs_k02[threadID];
    double *rk02 = data->drs_k02[threadID];
    double *k03 = data->rs_k03[threadID];
    double *rk03 = data->drs_k03[threadID];
    double *k04 = data->rs_k04[threadID];
    double *rk04 = data->drs_k04[threadID];
    double *k05 = data->rs_k05[threadID];
    double *rk05 = data->drs_k05[threadID];
    double *k06 = data->rs_k06[threadID];
    double *rk06 = data->drs_k06[threadID];
    double *k07 = data->rs_k07[threadID];
    double *rk07 = data->drs_k07[threadID];
    double *k08 = data->rs_k08[threadID];
    double *rk08 = data->drs_k08[threadID];
    double *k09 = data->rs_k09[threadID];
    double *rk09 = data->drs_k09[threadID];
    double *k10 = data->rs_k10[threadID];
    double *rk10 = data->drs_k10[threadID];
    double *k11 = data->rs_k11[threadID];
    double *rk11 = data->drs_k11[threadID];
    double *k12 = data->rs_k12[threadID];
    double *rk12 = data->drs_k12[threadID];
    double *k13 = data->rs_k13[threadID];
    double *rk13 = data->drs_k13[threadID];
    double *k14 = data->rs_k14[threadID];
    double *rk14 = data->drs_k14[threadID];
    double *k15 = data->rs_k15[threadID];
    double *rk15 = data->drs_k15[threadID];
    double *k16 = data->rs_k16[threadID];
    double *rk16 = data->drs_k16[threadID];
    double *k17 = data->rs_k17[threadID];
    double *rk17 = data->drs_k17[threadID];
    double *k18 = data->rs_k18[threadID];
    double *rk18 = data->drs_k18[threadID];
    double *k19 = data->rs_k19[threadID];
    double *rk19 = data->drs_k19[threadID];
    double *k21 = data->rs_k21[threadID];
    double *rk21 = data->drs_k21[threadID];
    double *k22 = data->rs_k22[threadID];
    double *rk22 = data->drs_k22[threadID];
    double *k23 = data->rs_k23[threadID];
    double *rk23 = data->drs_k23[threadID];
    double *brem_brem = data->cs_brem_brem[threadID];
    double *rbrem_brem = data->dcs_brem_brem[threadID];
    double *ceHeI_ceHeI = data->cs_ceHeI_ceHeI[threadID];
    double *rceHeI_ceHeI = data->dcs_ceHeI_ceHeI[threadID];
    double *ceHeII_ceHeII = data->cs_ceHeII_ceHeII[threadID];
    double *rceHeII_ceHeII = data->dcs_ceHeII_ceHeII[threadID];
    double *ceHI_ceHI = data->cs_ceHI_ceHI[threadID];
    double *rceHI_ceHI = data->dcs_ceHI_ceHI[threadID];
    double *cie_cooling_cieco = data->cs_cie_cooling_cieco[threadID];
    double *rcie_cooling_cieco = data->dcs_cie_cooling_cieco[threadID];
    double *ciHeI_ciHeI = data->cs_ciHeI_ciHeI[threadID];
    double *rciHeI_ciHeI = data->dcs_ciHeI_ciHeI[threadID];
    double *ciHeII_ciHeII = data->cs_ciHeII_ciHeII[threadID];
    double *rciHeII_ciHeII = data->dcs_ciHeII_ciHeII[threadID];
    double *ciHeIS_ciHeIS = data->cs_ciHeIS_ciHeIS[threadID];
    double *rciHeIS_ciHeIS = data->dcs_ciHeIS_ciHeIS[threadID];
    double *ciHI_ciHI = data->cs_ciHI_ciHI[threadID];
    double *rciHI_ciHI = data->dcs_ciHI_ciHI[threadID];
    double *compton_comp_ = data->cs_compton_comp_[threadID];
    double *rcompton_comp_ = data->dcs_compton_comp_[threadID];
    double *gammah_gammah = data->cs_gammah_gammah[threadID];
    double *rgammah_gammah = data->dcs_gammah_gammah[threadID];
    double *gloverabel08_gael = data->cs_gloverabel08_gael[threadID];
    double *rgloverabel08_gael = data->dcs_gloverabel08_gael[threadID];
    double *gloverabel08_gaH2 = data->cs_gloverabel08_gaH2[threadID];
    double *rgloverabel08_gaH2 = data->dcs_gloverabel08_gaH2[threadID];
    double *gloverabel08_gaHe = data->cs_gloverabel08_gaHe[threadID];
    double *rgloverabel08_gaHe = data->dcs_gloverabel08_gaHe[threadID];
    double *gloverabel08_gaHI = data->cs_gloverabel08_gaHI[threadID];
    double *rgloverabel08_gaHI = data->dcs_gloverabel08_gaHI[threadID];
    double *gloverabel08_gaHp = data->cs_gloverabel08_gaHp[threadID];
    double *rgloverabel08_gaHp = data->dcs_gloverabel08_gaHp[threadID];
    double *gloverabel08_gphdl = data->cs_gloverabel08_gphdl[threadID];
    double *rgloverabel08_gphdl = data->dcs_gloverabel08_gphdl[threadID];
    double *gloverabel08_gpldl = data->cs_gloverabel08_gpldl[threadID];
    double *rgloverabel08_gpldl = data->dcs_gloverabel08_gpldl[threadID];
    double *gloverabel08_h2lte = data->cs_gloverabel08_h2lte[threadID];
    double *rgloverabel08_h2lte = data->dcs_gloverabel08_h2lte[threadID];
    double *h2formation_h2mcool = data->cs_h2formation_h2mcool[threadID];
    double *rh2formation_h2mcool = data->dcs_h2formation_h2mcool[threadID];
    double *h2formation_h2mheat = data->cs_h2formation_h2mheat[threadID];
    double *rh2formation_h2mheat = data->dcs_h2formation_h2mheat[threadID];
    double *h2formation_ncrd1 = data->cs_h2formation_ncrd1[threadID];
    double *rh2formation_ncrd1 = data->dcs_h2formation_ncrd1[threadID];
    double *h2formation_ncrd2 = data->cs_h2formation_ncrd2[threadID];
    double *rh2formation_ncrd2 = data->dcs_h2formation_ncrd2[threadID];
    double *h2formation_ncrn = data->cs_h2formation_ncrn[threadID];
    double *rh2formation_ncrn = data->dcs_h2formation_ncrn[threadID];
    double *reHeII1_reHeII1 = data->cs_reHeII1_reHeII1[threadID];
    double *rreHeII1_reHeII1 = data->dcs_reHeII1_reHeII1[threadID];
    double *reHeII2_reHeII2 = data->cs_reHeII2_reHeII2[threadID];
    double *rreHeII2_reHeII2 = data->dcs_reHeII2_reHeII2[threadID];
    double *reHeIII_reHeIII = data->cs_reHeIII_reHeIII[threadID];
    double *rreHeIII_reHeIII = data->dcs_reHeIII_reHeIII[threadID];
    double *reHII_reHII = data->cs_reHII_reHII[threadID];
    double *rreHII_reHII = data->dcs_reHII_reHII[threadID];
    
    
    double h2_optical_depth_approx;  
    

    double scale;
    // Define the species
    double H2_1, v0;
    double H2_2, v1;
    double H_1, v2;
    double H_2, v3;
    double H_m0, v4;
    double He_1, v5;
    double He_2, v6;
    double He_3, v7;
    double de, v8;
    double ge, v9;
    scale = data->scale[threadID][0];
    v0 = Ith( v, 1 );
    v0 *= scale;
    scale = data->scale[threadID][1];
    v1 = Ith( v, 2 );
    v1 *= scale;
    scale = data->scale[threadID][2];
    v2 = Ith( v, 3 );
    v2 *= scale;
    scale = data->scale[threadID][3];
    v3 = Ith( v, 4 );
    v3 *= scale;
    scale = data->scale[threadID][4];
    v4 = Ith( v, 5 );
    v4 *= scale;
    scale = data->scale[threadID][5];
    v5 = Ith( v, 6 );
    v5 *= scale;
    scale = data->scale[threadID][6];
    v6 = Ith( v, 7 );
    v6 *= scale;
    scale = data->scale[threadID][7];
    v7 = Ith( v, 8 );
    v7 *= scale;
    scale = data->scale[threadID][8];
    v8 = Ith( v, 9 );
    v8 *= scale;
    scale = data->scale[threadID][9];
    v9 = Ith( v, 10 );
    v9 *= scale;

    double z;
    double T;

    double mh = 1.67e-24;
    double mdensity;
    
    int jj;
    jj = 0;

    j = i*nchem;
    mdensity = 0.0;
    z = data->current_z;
        // Rescale the Species abundance
    scale = data->scale[threadID][j];
    H2_1 = Ith( y, 1  )*scale;
    
    mdensity += H2_1;
    
    j++;
    
        // Rescale the Species abundance
    scale = data->scale[threadID][j];
    H2_2 = Ith( y, 2  )*scale;
    
    mdensity += H2_2;
    
    j++;
    
        // Rescale the Species abundance
    scale = data->scale[threadID][j];
    H_1 = Ith( y, 3  )*scale;
    
    mdensity += H_1;
    
    j++;
    
        // Rescale the Species abundance
    scale = data->scale[threadID][j];
    H_2 = Ith( y, 4  )*scale;
    
    mdensity += H_2;
    
    j++;
    
        // Rescale the Species abundance
    scale = data->scale[threadID][j];
    H_m0 = Ith( y, 5  )*scale;
    
    mdensity += H_m0;
    
    j++;
    
        // Rescale the Species abundance
    scale = data->scale[threadID][j];
    He_1 = Ith( y, 6  )*scale;
    
    mdensity += He_1;
    
    j++;
    
        // Rescale the Species abundance
    scale = data->scale[threadID][j];
    He_2 = Ith( y, 7  )*scale;
    
    mdensity += He_2;
    
    j++;
    
        // Rescale the Species abundance
    scale = data->scale[threadID][j];
    He_3 = Ith( y, 8  )*scale;
    
    mdensity += He_3;
    
    j++;
    
        // Rescale the Species abundance
    scale = data->scale[threadID][j];
    de = Ith( y, 9  )*scale;
    
    j++;
    
        // Rescale the Species abundance
    scale = data->scale[threadID][j];
    ge = Ith( y, 10  )*scale;
    
    j++;
    


    mdensity *= mh;
        
    j = 0;
    //
    // Species: H2_1
    //
    Ith(Jv, 1 ) = -k11[i]*H2_1*v3 - k12[i]*H2_1*v8 + Tge[i]*v9*(-pow(H2_1, 2)*rk23[i] + H2_1*pow(H_1, 2)*rk21[i] - H2_1*H_1*rk13[i] - H2_1*H_2*rk11[i] - H2_1*de*rk12[i] + H2_2*H_1*rk10[i] + H2_2*H_m0*rk19[i] + pow(H_1, 3)*rk22[i] + H_1*H_m0*rk08[i]) + v0*(-k11[i]*H_2 - k12[i]*de - k13[i]*H_1 + k21[i]*pow(H_1, 2) - 2*k23[i]*H2_1) + v1*(k10[i]*H_1 + k19[i]*H_m0) + v2*(k08[i]*H_m0 + k10[i]*H2_2 - k13[i]*H2_1 + 2*k21[i]*H2_1*H_1 + 3*k22[i]*pow(H_1, 2)) + v4*(k08[i]*H_1 + k19[i]*H2_2);

    scale = data->scale[threadID][0];
    Ith(Jv, 1) /= scale;

    
    //
    // Species: H2_2
    //
    Ith(Jv, 2 ) = k11[i]*H_2*v0 - k18[i]*H2_2*v8 + Tge[i]*v9*(H2_1*H_2*rk11[i] - H2_2*H_1*rk10[i] - H2_2*H_m0*rk19[i] - H2_2*de*rk18[i] + H_1*H_2*rk09[i] + H_2*H_m0*rk17[i]) + v1*(-k10[i]*H_1 - k18[i]*de - k19[i]*H_m0) + v2*(k09[i]*H_2 - k10[i]*H2_2) + v3*(k09[i]*H_1 + k11[i]*H2_1 + k17[i]*H_m0) + v4*(k17[i]*H_2 - k19[i]*H2_2);

    scale = data->scale[threadID][1];
    Ith(Jv, 2) /= scale;

    
    //
    // Species: H_1
    //
    Ith(Jv, 3 ) = Tge[i]*v9*(2*pow(H2_1, 2)*rk23[i] - 2*H2_1*pow(H_1, 2)*rk21[i] + 2*H2_1*H_1*rk13[i] + H2_1*H_2*rk11[i] + 2*H2_1*de*rk12[i] - H2_2*H_1*rk10[i] + H2_2*H_m0*rk19[i] + 2*H2_2*de*rk18[i] - 2*pow(H_1, 3)*rk22[i] - H_1*H_2*rk09[i] - H_1*H_m0*rk08[i] + H_1*H_m0*rk15[i] - H_1*de*rk01[i] - H_1*de*rk07[i] + 2*H_2*H_m0*rk16[i] + H_2*de*rk02[i] + H_m0*de*rk14[i]) + v0*(k11[i]*H_2 + 2*k12[i]*de + 2*k13[i]*H_1 - 2*k21[i]*pow(H_1, 2) + 4*k23[i]*H2_1) + v1*(-k10[i]*H_1 + 2*k18[i]*de + k19[i]*H_m0) + v2*(-k01[i]*de - k07[i]*de - k08[i]*H_m0 - k09[i]*H_2 - k10[i]*H2_2 + 2*k13[i]*H2_1 + k15[i]*H_m0 - 4*k21[i]*H2_1*H_1 - 6*k22[i]*pow(H_1, 2)) + v3*(k02[i]*de - k09[i]*H_1 + k11[i]*H2_1 + 2*k16[i]*H_m0) + v4*(-k08[i]*H_1 + k14[i]*de + k15[i]*H_1 + 2*k16[i]*H_2 + k19[i]*H2_2) + v8*(-k01[i]*H_1 + k02[i]*H_2 - k07[i]*H_1 + 2*k12[i]*H2_1 + k14[i]*H_m0 + 2*k18[i]*H2_2);

    scale = data->scale[threadID][2];
    Ith(Jv, 3) /= scale;

    
    //
    // Species: H_2
    //
    Ith(Jv, 4 ) = k10[i]*H_1*v1 - k11[i]*H_2*v0 + Tge[i]*v9*(-H2_1*H_2*rk11[i] + H2_2*H_1*rk10[i] - H_1*H_2*rk09[i] + H_1*de*rk01[i] - H_2*H_m0*rk16[i] - H_2*H_m0*rk17[i] - H_2*de*rk02[i]) + v2*(k01[i]*de - k09[i]*H_2 + k10[i]*H2_2) + v3*(-k02[i]*de - k09[i]*H_1 - k11[i]*H2_1 - k16[i]*H_m0 - k17[i]*H_m0) + v4*(-k16[i]*H_2 - k17[i]*H_2) + v8*(k01[i]*H_1 - k02[i]*H_2);

    scale = data->scale[threadID][3];
    Ith(Jv, 4) /= scale;

    
    //
    // Species: H_m0
    //
    Ith(Jv, 5 ) = -k19[i]*H_m0*v1 + Tge[i]*v9*(-H2_2*H_m0*rk19[i] - H_1*H_m0*rk08[i] - H_1*H_m0*rk15[i] + H_1*de*rk07[i] - H_2*H_m0*rk16[i] - H_2*H_m0*rk17[i] - H_m0*de*rk14[i]) + v2*(k07[i]*de - k08[i]*H_m0 - k15[i]*H_m0) + v3*(-k16[i]*H_m0 - k17[i]*H_m0) + v4*(-k08[i]*H_1 - k14[i]*de - k15[i]*H_1 - k16[i]*H_2 - k17[i]*H_2 - k19[i]*H2_2) + v8*(k07[i]*H_1 - k14[i]*H_m0);

    scale = data->scale[threadID][4];
    Ith(Jv, 5) /= scale;

    
    //
    // Species: He_1
    //
    Ith(Jv, 6 ) = -k03[i]*de*v5 + k04[i]*de*v6 + Tge[i]*v9*(-He_1*de*rk03[i] + He_2*de*rk04[i]) + v8*(-k03[i]*He_1 + k04[i]*He_2);

    scale = data->scale[threadID][5];
    Ith(Jv, 6) /= scale;

    
    //
    // Species: He_2
    //
    Ith(Jv, 7 ) = k03[i]*de*v5 + k06[i]*de*v7 + Tge[i]*v9*(He_1*de*rk03[i] - He_2*de*rk04[i] - He_2*de*rk05[i] + He_3*de*rk06[i]) + v6*(-k04[i]*de - k05[i]*de) + v8*(k03[i]*He_1 - k04[i]*He_2 - k05[i]*He_2 + k06[i]*He_3);

    scale = data->scale[threadID][6];
    Ith(Jv, 7) /= scale;

    
    //
    // Species: He_3
    //
    Ith(Jv, 8 ) = k05[i]*de*v6 - k06[i]*de*v7 + Tge[i]*v9*(He_2*de*rk05[i] - He_3*de*rk06[i]) + v8*(k05[i]*He_2 - k06[i]*He_3);

    scale = data->scale[threadID][7];
    Ith(Jv, 8) /= scale;

    
    //
    // Species: de
    //
    Ith(Jv, 9 ) = k03[i]*de*v5 - k06[i]*de*v7 - k18[i]*de*v1 + Tge[i]*v9*(-H2_2*de*rk18[i] + H_1*H_m0*rk08[i] + H_1*H_m0*rk15[i] + H_1*de*rk01[i] - H_1*de*rk07[i] + H_2*H_m0*rk17[i] - H_2*de*rk02[i] + H_m0*de*rk14[i] + He_1*de*rk03[i] - He_2*de*rk04[i] + He_2*de*rk05[i] - He_3*de*rk06[i]) + v2*(k01[i]*de - k07[i]*de + k08[i]*H_m0 + k15[i]*H_m0) + v3*(-k02[i]*de + k17[i]*H_m0) + v4*(k08[i]*H_1 + k14[i]*de + k15[i]*H_1 + k17[i]*H_2) + v6*(-k04[i]*de + k05[i]*de) + v8*(k01[i]*H_1 - k02[i]*H_2 + k03[i]*He_1 - k04[i]*He_2 + k05[i]*He_2 - k06[i]*He_3 - k07[i]*H_1 + k14[i]*H_m0 - k18[i]*H2_2);

    scale = data->scale[threadID][8];
    Ith(Jv, 9) /= scale;

    
    //
    // Species: ge
    //
    Ith(Jv, 10 ) = Tge[i]*v9*(-gloverabel08_h2lte[i]*H2_1*h2_optical_depth_approx*(-gloverabel08_h2lte[i]*(-H2_1*rgloverabel08_gaH2[i] - H_1*rgloverabel08_gaHI[i] - H_2*rgloverabel08_gaHp[i] - He_1*rgloverabel08_gaHe[i] - de*rgloverabel08_gael[i])/pow(gloverabel08_gaH2[i]*H2_1 + gloverabel08_gaHI[i]*H_1 + gloverabel08_gaHe[i]*He_1 + gloverabel08_gaHp[i]*H_2 + gloverabel08_gael[i]*de, 2) - rgloverabel08_h2lte[i]/(gloverabel08_gaH2[i]*H2_1 + gloverabel08_gaHI[i]*H_1 + gloverabel08_gaHe[i]*He_1 + gloverabel08_gaHp[i]*H_2 + gloverabel08_gael[i]*de))/pow(gloverabel08_h2lte[i]/(gloverabel08_gaH2[i]*H2_1 + gloverabel08_gaHI[i]*H_1 + gloverabel08_gaHe[i]*He_1 + gloverabel08_gaHp[i]*H_2 + gloverabel08_gael[i]*de) + 1.0, 2) - H2_1*h2_optical_depth_approx*rgloverabel08_h2lte[i]/(gloverabel08_h2lte[i]/(gloverabel08_gaH2[i]*H2_1 + gloverabel08_gaHI[i]*H_1 + gloverabel08_gaHe[i]*He_1 + gloverabel08_gaHp[i]*H_2 + gloverabel08_gael[i]*de) + 1.0) - 2.0158800000000001*H2_1*mdensity*rcie_cooling_cieco[i] - H_1*de*rceHI_ceHI[i] - H_1*de*rciHI_ciHI[i] - H_2*de*rreHII_reHII[i] - He_1*de*rciHeI_ciHeI[i] - He_2*pow(de, 2)*rceHeI_ceHeI[i] - He_2*pow(de, 2)*rciHeIS_ciHeIS[i] - He_2*de*rceHeII_ceHeII[i] - He_2*de*rciHeII_ciHeII[i] - He_2*de*rreHeII1_reHeII1[i] - He_2*de*rreHeII2_reHeII2[i] - He_3*de*rreHeIII_reHeIII[i] - de*rbrem_brem[i]*(H_2 + He_2 + 4.0*He_3) - de*rcompton_comp_[i]*pow(z + 1.0, 4)*(T - 2.73*z - 2.73) + 0.5*pow(h2formation_ncrn[i]/(h2formation_ncrd1[i]*H_1 + h2formation_ncrd2[i]*H2_1) + 1.0, -2.0)*(-h2formation_h2mcool[i]*H2_1*H_1 + h2formation_h2mheat[i]*pow(H_1, 3))*(-1.0*h2formation_ncrn[i]*(-H2_1*rh2formation_ncrd2[i] - H_1*rh2formation_ncrd1[i])/pow(h2formation_ncrd1[i]*H_1 + h2formation_ncrd2[i]*H2_1, 2) - 1.0*rh2formation_ncrn[i]/(h2formation_ncrd1[i]*H_1 + h2formation_ncrd2[i]*H2_1)) + 0.5*1.0/(h2formation_ncrn[i]/(h2formation_ncrd1[i]*H_1 + h2formation_ncrd2[i]*H2_1) + 1.0)*(-H2_1*H_1*rh2formation_h2mcool[i] + pow(H_1, 3)*rh2formation_h2mheat[i]))/mdensity + v0*(-2.0158800000000001*cie_cooling_cieco[i]*mdensity - gloverabel08_gaH2[i]*pow(gloverabel08_h2lte[i], 2)*H2_1*h2_optical_depth_approx/(pow(gloverabel08_h2lte[i]/(gloverabel08_gaH2[i]*H2_1 + gloverabel08_gaHI[i]*H_1 + gloverabel08_gaHe[i]*He_1 + gloverabel08_gaHp[i]*H_2 + gloverabel08_gael[i]*de) + 1.0, 2)*pow(gloverabel08_gaH2[i]*H2_1 + gloverabel08_gaHI[i]*H_1 + gloverabel08_gaHe[i]*He_1 + gloverabel08_gaHp[i]*H_2 + gloverabel08_gael[i]*de, 2)) - gloverabel08_h2lte[i]*h2_optical_depth_approx/(gloverabel08_h2lte[i]/(gloverabel08_gaH2[i]*H2_1 + gloverabel08_gaHI[i]*H_1 + gloverabel08_gaHe[i]*He_1 + gloverabel08_gaHp[i]*H_2 + gloverabel08_gael[i]*de) + 1.0) - 0.5*h2formation_h2mcool[i]*H_1*1.0/(h2formation_ncrn[i]/(h2formation_ncrd1[i]*H_1 + h2formation_ncrd2[i]*H2_1) + 1.0) + 0.5*h2formation_ncrd2[i]*h2formation_ncrn[i]*pow(h2formation_ncrn[i]/(h2formation_ncrd1[i]*H_1 + h2formation_ncrd2[i]*H2_1) + 1.0, -2.0)*(-h2formation_h2mcool[i]*H2_1*H_1 + h2formation_h2mheat[i]*pow(H_1, 3))/pow(h2formation_ncrd1[i]*H_1 + h2formation_ncrd2[i]*H2_1, 2))/mdensity + v2*(-ceHI_ceHI[i]*de - ciHI_ciHI[i]*de - gloverabel08_gaHI[i]*pow(gloverabel08_h2lte[i], 2)*H2_1*h2_optical_depth_approx/(pow(gloverabel08_h2lte[i]/(gloverabel08_gaH2[i]*H2_1 + gloverabel08_gaHI[i]*H_1 + gloverabel08_gaHe[i]*He_1 + gloverabel08_gaHp[i]*H_2 + gloverabel08_gael[i]*de) + 1.0, 2)*pow(gloverabel08_gaH2[i]*H2_1 + gloverabel08_gaHI[i]*H_1 + gloverabel08_gaHe[i]*He_1 + gloverabel08_gaHp[i]*H_2 + gloverabel08_gael[i]*de, 2)) + 0.5*h2formation_ncrd1[i]*h2formation_ncrn[i]*pow(h2formation_ncrn[i]/(h2formation_ncrd1[i]*H_1 + h2formation_ncrd2[i]*H2_1) + 1.0, -2.0)*(-h2formation_h2mcool[i]*H2_1*H_1 + h2formation_h2mheat[i]*pow(H_1, 3))/pow(h2formation_ncrd1[i]*H_1 + h2formation_ncrd2[i]*H2_1, 2) + 0.5*(-h2formation_h2mcool[i]*H2_1 + 3*h2formation_h2mheat[i]*pow(H_1, 2))*1.0/(h2formation_ncrn[i]/(h2formation_ncrd1[i]*H_1 + h2formation_ncrd2[i]*H2_1) + 1.0))/mdensity + v3*(-brem_brem[i]*de - gloverabel08_gaHp[i]*pow(gloverabel08_h2lte[i], 2)*H2_1*h2_optical_depth_approx/(pow(gloverabel08_h2lte[i]/(gloverabel08_gaH2[i]*H2_1 + gloverabel08_gaHI[i]*H_1 + gloverabel08_gaHe[i]*He_1 + gloverabel08_gaHp[i]*H_2 + gloverabel08_gael[i]*de) + 1.0, 2)*pow(gloverabel08_gaH2[i]*H2_1 + gloverabel08_gaHI[i]*H_1 + gloverabel08_gaHe[i]*He_1 + gloverabel08_gaHp[i]*H_2 + gloverabel08_gael[i]*de, 2)) - reHII_reHII[i]*de)/mdensity + v5*(-ciHeI_ciHeI[i]*de - gloverabel08_gaHe[i]*pow(gloverabel08_h2lte[i], 2)*H2_1*h2_optical_depth_approx/(pow(gloverabel08_h2lte[i]/(gloverabel08_gaH2[i]*H2_1 + gloverabel08_gaHI[i]*H_1 + gloverabel08_gaHe[i]*He_1 + gloverabel08_gaHp[i]*H_2 + gloverabel08_gael[i]*de) + 1.0, 2)*pow(gloverabel08_gaH2[i]*H2_1 + gloverabel08_gaHI[i]*H_1 + gloverabel08_gaHe[i]*He_1 + gloverabel08_gaHp[i]*H_2 + gloverabel08_gael[i]*de, 2)))/mdensity + v6*(-brem_brem[i]*de - ceHeII_ceHeII[i]*de - ceHeI_ceHeI[i]*pow(de, 2) - ciHeII_ciHeII[i]*de - ciHeIS_ciHeIS[i]*pow(de, 2) - reHeII1_reHeII1[i]*de - reHeII2_reHeII2[i]*de)/mdensity + v7*(-4.0*brem_brem[i]*de - reHeIII_reHeIII[i]*de)/mdensity + v8*(brem_brem[i]*(-H_2 - He_2 - 4.0*He_3) - ceHI_ceHI[i]*H_1 - ceHeII_ceHeII[i]*He_2 - 2*ceHeI_ceHeI[i]*He_2*de - ciHI_ciHI[i]*H_1 - ciHeII_ciHeII[i]*He_2 - 2*ciHeIS_ciHeIS[i]*He_2*de - ciHeI_ciHeI[i]*He_1 - compton_comp_[i]*pow(z + 1.0, 4)*(T - 2.73*z - 2.73) - gloverabel08_gael[i]*pow(gloverabel08_h2lte[i], 2)*H2_1*h2_optical_depth_approx/(pow(gloverabel08_h2lte[i]/(gloverabel08_gaH2[i]*H2_1 + gloverabel08_gaHI[i]*H_1 + gloverabel08_gaHe[i]*He_1 + gloverabel08_gaHp[i]*H_2 + gloverabel08_gael[i]*de) + 1.0, 2)*pow(gloverabel08_gaH2[i]*H2_1 + gloverabel08_gaHI[i]*H_1 + gloverabel08_gaHe[i]*He_1 + gloverabel08_gaHp[i]*H_2 + gloverabel08_gael[i]*de, 2)) - reHII_reHII[i]*H_2 - reHeII1_reHeII1[i]*He_2 - reHeII2_reHeII2[i]*He_2 - reHeIII_reHeIII[i]*He_3)/mdensity;

    scale = data->scale[threadID][9];
    Ith(Jv, 10) /= scale;

    
    Ith(Jv, 10) /= mdensity;
    
    return 0;
}

*/


///////////////////////////////////////////////////////////////////////////////
////////////////// RHS functions and Jacobian Functions ///////////////////////
///////////////////////////////////////////////////////////////////////////////

int calculate_rhs_primordial(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
    primordial_data *data = (primordial_data* ) user_data;
    int i, j;

    int nchem = 10;
    int nstrip = data->nstrip;
    
    /* change N_Vector back to an array */
    double y_arr[ nchem * nstrip ];
    double H2_1;
    double H2_2;
    double H_1;
    double H_2;
    double H_m0;
    double He_1;
    double He_2;
    double He_3;
    double de;
    double ge;
    
    #ifdef _OPENMP
    int threadID = omp_get_thread_num();
    #else
    int threadID = 0;
    #endif
   
    double *yvec_ptr = N_VGetArrayPointer(y);
    double *ydot_ptr = N_VGetArrayPointer(ydot);

    #ifdef SCALE_INPUT
    double *scale     = data->scale[threadID];
    double *inv_scale = data->inv_scale[threadID];   
    for ( i = 0; i < nstrip*nchem; i++ ){
        y_arr[i] = yvec_ptr[i]*scale[i];
    }
    #else
    for ( i = 0; i < nstrip*nchem; i++ ){
        y_arr[i] = yvec_ptr[i];
    }
    #endif

    int flag;
    flag = primordial_calculate_temperature(data, y_arr , nstrip, nchem );
    if (flag > 0){
        // check if the temperature failed to converged
        return -1;    
    }
    primordial_interpolate_rates(data, nstrip);

    /* Now we set up some temporaries */
    double *k01 = data->rs_k01[threadID];
    double *k02 = data->rs_k02[threadID];
    double *k03 = data->rs_k03[threadID];
    double *k04 = data->rs_k04[threadID];
    double *k05 = data->rs_k05[threadID];
    double *k06 = data->rs_k06[threadID];
    double *k07 = data->rs_k07[threadID];
    double *k08 = data->rs_k08[threadID];
    double *k09 = data->rs_k09[threadID];
    double *k10 = data->rs_k10[threadID];
    double *k11 = data->rs_k11[threadID];
    double *k12 = data->rs_k12[threadID];
    double *k13 = data->rs_k13[threadID];
    double *k14 = data->rs_k14[threadID];
    double *k15 = data->rs_k15[threadID];
    double *k16 = data->rs_k16[threadID];
    double *k17 = data->rs_k17[threadID];
    double *k18 = data->rs_k18[threadID];
    double *k19 = data->rs_k19[threadID];
    double *k21 = data->rs_k21[threadID];
    double *k22 = data->rs_k22[threadID];
    double *k23 = data->rs_k23[threadID];
    double *brem_brem = data->cs_brem_brem[threadID];
    double *ceHeI_ceHeI = data->cs_ceHeI_ceHeI[threadID];
    double *ceHeII_ceHeII = data->cs_ceHeII_ceHeII[threadID];
    double *ceHI_ceHI = data->cs_ceHI_ceHI[threadID];
    double *cie_cooling_cieco = data->cs_cie_cooling_cieco[threadID];
    double *ciHeI_ciHeI = data->cs_ciHeI_ciHeI[threadID];
    double *ciHeII_ciHeII = data->cs_ciHeII_ciHeII[threadID];
    double *ciHeIS_ciHeIS = data->cs_ciHeIS_ciHeIS[threadID];
    double *ciHI_ciHI = data->cs_ciHI_ciHI[threadID];
    double *compton_comp_ = data->cs_compton_comp_[threadID];
    double *gammah_gammah = data->cs_gammah_gammah[threadID];
    double *gloverabel08_gael = data->cs_gloverabel08_gael[threadID];
    double *gloverabel08_gaH2 = data->cs_gloverabel08_gaH2[threadID];
    double *gloverabel08_gaHe = data->cs_gloverabel08_gaHe[threadID];
    double *gloverabel08_gaHI = data->cs_gloverabel08_gaHI[threadID];
    double *gloverabel08_gaHp = data->cs_gloverabel08_gaHp[threadID];
    double *gloverabel08_gphdl = data->cs_gloverabel08_gphdl[threadID];
    double *gloverabel08_gpldl = data->cs_gloverabel08_gpldl[threadID];
    double *gloverabel08_h2lte = data->cs_gloverabel08_h2lte[threadID];
    double *h2formation_h2mcool = data->cs_h2formation_h2mcool[threadID];
    double *h2formation_h2mheat = data->cs_h2formation_h2mheat[threadID];
    double *h2formation_ncrd1 = data->cs_h2formation_ncrd1[threadID];
    double *h2formation_ncrd2 = data->cs_h2formation_ncrd2[threadID];
    double *h2formation_ncrn = data->cs_h2formation_ncrn[threadID];
    double *reHeII1_reHeII1 = data->cs_reHeII1_reHeII1[threadID];
    double *reHeII2_reHeII2 = data->cs_reHeII2_reHeII2[threadID];
    double *reHeIII_reHeIII = data->cs_reHeIII_reHeIII[threadID];
    double *reHII_reHII = data->cs_reHII_reHII[threadID];
    double h2_optical_depth_approx;    
    
    double cie_optical_depth_approx;
    

    double z;
    double T;

    double mh = 1.66054e-24;
    double mdensity, inv_mdensity;

   
    for ( i = 0; i < nstrip; i++ ){
        
        T            = data->Ts[threadID][i];
        z            = data->current_z;
        mdensity     = data->mdensity[threadID][i];
        inv_mdensity = data->inv_mdensity[threadID][i];
        h2_optical_depth_approx = data->h2_optical_depth_approx[threadID][i];
        
        cie_optical_depth_approx = data->cie_optical_depth_approx[threadID][i];
        

        j = i * nchem;
        H2_1 = y_arr[j];
        j++;
        H2_2 = y_arr[j];
        j++;
        H_1 = y_arr[j];
        j++;
        H_2 = y_arr[j];
        j++;
        H_m0 = y_arr[j];
        j++;
        He_1 = y_arr[j];
        j++;
        He_2 = y_arr[j];
        j++;
        He_3 = y_arr[j];
        j++;
        de = y_arr[j];
        j++;
        ge = y_arr[j];
        j++;
    
        j = i * nchem;
        // Species: H2_1
        ydot_ptr[j] = k08[i]*H_1*H_m0 + k10[i]*H2_2*H_1 - k11[i]*H2_1*H_2 - k12[i]*H2_1*de - k13[i]*H2_1*H_1 + k19[i]*H2_2*H_m0 + k21[i]*H2_1*H_1*H_1 + k22[i]*H_1*H_1*H_1 - k23[i]*H2_1*H2_1;
        #ifdef SCALE_INPUT
        ydot_ptr[j] *= inv_scale[j];
        #endif
        
        //fprintf(stderr, "H2_1: %0.5g\n", scale[j]);
        //fprintf(stderr, "ydot = %0.5g \n", ydot_ptr[j]*scale[j] );
        j++;
        
        // Species: H2_2
        ydot_ptr[j] = k09[i]*H_1*H_2 - k10[i]*H2_2*H_1 + k11[i]*H2_1*H_2 + k17[i]*H_2*H_m0 - k18[i]*H2_2*de - k19[i]*H2_2*H_m0;
        #ifdef SCALE_INPUT
        ydot_ptr[j] *= inv_scale[j];
        #endif
        
        //fprintf(stderr, "H2_2: %0.5g\n", scale[j]);
        //fprintf(stderr, "ydot = %0.5g \n", ydot_ptr[j]*scale[j] );
        j++;
        
        // Species: H_1
        ydot_ptr[j] = -k01[i]*H_1*de + k02[i]*H_2*de - k07[i]*H_1*de - k08[i]*H_1*H_m0 - k09[i]*H_1*H_2 - k10[i]*H2_2*H_1 + k11[i]*H2_1*H_2 + 2*k12[i]*H2_1*de + 2*k13[i]*H2_1*H_1 + k14[i]*H_m0*de + k15[i]*H_1*H_m0 + 2*k16[i]*H_2*H_m0 + 2*k18[i]*H2_2*de + k19[i]*H2_2*H_m0 - 2*k21[i]*H2_1*H_1*H_1 - 2*k22[i]*H_1*H_1*H_1 + 2*k23[i]*H2_1*H2_1;
        #ifdef SCALE_INPUT
        ydot_ptr[j] *= inv_scale[j];
        #endif
        
        //fprintf(stderr, "H_1: %0.5g\n", scale[j]);
        //fprintf(stderr, "ydot = %0.5g \n", ydot_ptr[j]*scale[j] );
        j++;
        
        // Species: H_2
        ydot_ptr[j] = k01[i]*H_1*de - k02[i]*H_2*de - k09[i]*H_1*H_2 + k10[i]*H2_2*H_1 - k11[i]*H2_1*H_2 - k16[i]*H_2*H_m0 - k17[i]*H_2*H_m0;
        #ifdef SCALE_INPUT
        ydot_ptr[j] *= inv_scale[j];
        #endif
        
        //fprintf(stderr, "H_2: %0.5g\n", scale[j]);
        //fprintf(stderr, "ydot = %0.5g \n", ydot_ptr[j]*scale[j] );
        j++;
        
        // Species: H_m0
        ydot_ptr[j] = k07[i]*H_1*de - k08[i]*H_1*H_m0 - k14[i]*H_m0*de - k15[i]*H_1*H_m0 - k16[i]*H_2*H_m0 - k17[i]*H_2*H_m0 - k19[i]*H2_2*H_m0;
        #ifdef SCALE_INPUT
        ydot_ptr[j] *= inv_scale[j];
        #endif
        
        //fprintf(stderr, "H_m0: %0.5g\n", scale[j]);
        //fprintf(stderr, "ydot = %0.5g \n", ydot_ptr[j]*scale[j] );
        j++;
        
        // Species: He_1
        ydot_ptr[j] = -k03[i]*He_1*de + k04[i]*He_2*de;
        #ifdef SCALE_INPUT
        ydot_ptr[j] *= inv_scale[j];
        #endif
        
        //fprintf(stderr, "He_1: %0.5g\n", scale[j]);
        //fprintf(stderr, "ydot = %0.5g \n", ydot_ptr[j]*scale[j] );
        j++;
        
        // Species: He_2
        ydot_ptr[j] = k03[i]*He_1*de - k04[i]*He_2*de - k05[i]*He_2*de + k06[i]*He_3*de;
        #ifdef SCALE_INPUT
        ydot_ptr[j] *= inv_scale[j];
        #endif
        
        //fprintf(stderr, "He_2: %0.5g\n", scale[j]);
        //fprintf(stderr, "ydot = %0.5g \n", ydot_ptr[j]*scale[j] );
        j++;
        
        // Species: He_3
        ydot_ptr[j] = k05[i]*He_2*de - k06[i]*He_3*de;
        #ifdef SCALE_INPUT
        ydot_ptr[j] *= inv_scale[j];
        #endif
        
        //fprintf(stderr, "He_3: %0.5g\n", scale[j]);
        //fprintf(stderr, "ydot = %0.5g \n", ydot_ptr[j]*scale[j] );
        j++;
        
        // Species: de
        ydot_ptr[j] = k01[i]*H_1*de - k02[i]*H_2*de + k03[i]*He_1*de - k04[i]*He_2*de + k05[i]*He_2*de - k06[i]*He_3*de - k07[i]*H_1*de + k08[i]*H_1*H_m0 + k14[i]*H_m0*de + k15[i]*H_1*H_m0 + k17[i]*H_2*H_m0 - k18[i]*H2_2*de;
        #ifdef SCALE_INPUT
        ydot_ptr[j] *= inv_scale[j];
        #endif
        
        //fprintf(stderr, "de: %0.5g\n", scale[j]);
        //fprintf(stderr, "ydot = %0.5g \n", ydot_ptr[j]*scale[j] );
        j++;
        
        // Species: ge
        ydot_ptr[j] = -brem_brem[i]*cie_optical_depth_approx*de*(H_2 + He_2 + 4.0*He_3) - ceHI_ceHI[i]*H_1*cie_optical_depth_approx*de - ceHeII_ceHeII[i]*He_2*cie_optical_depth_approx*de - ceHeI_ceHeI[i]*He_2*cie_optical_depth_approx*pow(de, 2) - ciHI_ciHI[i]*H_1*cie_optical_depth_approx*de - ciHeII_ciHeII[i]*He_2*cie_optical_depth_approx*de - ciHeIS_ciHeIS[i]*He_2*cie_optical_depth_approx*pow(de, 2) - ciHeI_ciHeI[i]*He_1*cie_optical_depth_approx*de - 2.0158800000000001*cie_cooling_cieco[i]*H2_1*cie_optical_depth_approx*mdensity - compton_comp_[i]*cie_optical_depth_approx*de*pow(z + 1.0, 4)*(T - 2.73*z - 2.73) - gloverabel08_h2lte[i]*H2_1*cie_optical_depth_approx*h2_optical_depth_approx/(gloverabel08_h2lte[i]/(gloverabel08_gaH2[i]*H2_1 + gloverabel08_gaHI[i]*H_1 + gloverabel08_gaHe[i]*He_1 + gloverabel08_gaHp[i]*H_2 + gloverabel08_gael[i]*de) + 1.0) - reHII_reHII[i]*H_2*cie_optical_depth_approx*de - reHeII1_reHeII1[i]*He_2*cie_optical_depth_approx*de - reHeII2_reHeII2[i]*He_2*cie_optical_depth_approx*de - reHeIII_reHeIII[i]*He_3*cie_optical_depth_approx*de + 0.5*1.0/(h2formation_ncrn[i]/(h2formation_ncrd1[i]*H_1 + h2formation_ncrd2[i]*H2_1) + 1.0)*(-h2formation_h2mcool[i]*H2_1*H_1 + h2formation_h2mheat[i]*pow(H_1, 3));
        #ifdef SCALE_INPUT
        ydot_ptr[j] *= inv_scale[j];
        #endif
        
        ydot_ptr[j] *= inv_mdensity;
        
        //fprintf(stderr, "ge: %0.5g\n", scale[j]);
        //fprintf(stderr, "ydot = %0.5g \n", ydot_ptr[j]*scale[j] );
        j++;
        
    
    //fprintf(stderr, "----------------\n");
    }
    return 0;
    }



int calculate_jacobian_primordial( realtype t,
                                        N_Vector y, N_Vector fy,
                                        SUNMatrix J, void *user_data,
                                        N_Vector tmp1, N_Vector tmp2,
                                        N_Vector tmp3)
{
    /* We iterate over all of the rates */
    /* Calcuate temperature first */
    

    primordial_data *data = (primordial_data*)user_data; 
    
    #ifdef _OPENMP
    int threadID = omp_get_thread_num();
    #else
    int threadID = 0;
    #endif

    int nchem = 10;
    int nstrip = data->nstrip;
    int i, j;

    /* change N_Vector back to an array */
    double *yvec_ptr = N_VGetArrayPointer(y);
    double y_arr[ 10 * nstrip ];
 
    #ifdef SCALE_INPUT
    double *scale     = data->scale[threadID];
    double *inv_scale = data->inv_scale[threadID];   
    double scale1, inv_scale2;
    for ( i = 0; i < nstrip*nchem; i++ ){
        y_arr[i] = yvec_ptr[i]*scale[i];
    }
    #else
    for ( i = 0; i < nstrip*nchem; i++ ){
        y_arr[i] = yvec_ptr[i];
    }
    #endif


    /*
    int flag;
    flag = primordial_calculate_temperature(data, y_arr, nstrip, nchem );
    if (flag > 0){
        // check if the temperature failed to converged
        return -1;    
    }
    primordial_interpolate_rates(data, nstrip);
    */

    // primordial_calculate_temperature(data, y_arr, nstrip, nchem);
    // primordial_interpolate_rates(data, nstrip);

    /* Now We set up some temporaries */
    double *Tge = data->dTs_ge[threadID];
    double *k01 = data->rs_k01[threadID];
    double *rk01= data->drs_k01[threadID];
    double *k02 = data->rs_k02[threadID];
    double *rk02= data->drs_k02[threadID];
    double *k03 = data->rs_k03[threadID];
    double *rk03= data->drs_k03[threadID];
    double *k04 = data->rs_k04[threadID];
    double *rk04= data->drs_k04[threadID];
    double *k05 = data->rs_k05[threadID];
    double *rk05= data->drs_k05[threadID];
    double *k06 = data->rs_k06[threadID];
    double *rk06= data->drs_k06[threadID];
    double *k07 = data->rs_k07[threadID];
    double *rk07= data->drs_k07[threadID];
    double *k08 = data->rs_k08[threadID];
    double *rk08= data->drs_k08[threadID];
    double *k09 = data->rs_k09[threadID];
    double *rk09= data->drs_k09[threadID];
    double *k10 = data->rs_k10[threadID];
    double *rk10= data->drs_k10[threadID];
    double *k11 = data->rs_k11[threadID];
    double *rk11= data->drs_k11[threadID];
    double *k12 = data->rs_k12[threadID];
    double *rk12= data->drs_k12[threadID];
    double *k13 = data->rs_k13[threadID];
    double *rk13= data->drs_k13[threadID];
    double *k14 = data->rs_k14[threadID];
    double *rk14= data->drs_k14[threadID];
    double *k15 = data->rs_k15[threadID];
    double *rk15= data->drs_k15[threadID];
    double *k16 = data->rs_k16[threadID];
    double *rk16= data->drs_k16[threadID];
    double *k17 = data->rs_k17[threadID];
    double *rk17= data->drs_k17[threadID];
    double *k18 = data->rs_k18[threadID];
    double *rk18= data->drs_k18[threadID];
    double *k19 = data->rs_k19[threadID];
    double *rk19= data->drs_k19[threadID];
    double *k21 = data->rs_k21[threadID];
    double *rk21= data->drs_k21[threadID];
    double *k22 = data->rs_k22[threadID];
    double *rk22= data->drs_k22[threadID];
    double *k23 = data->rs_k23[threadID];
    double *rk23= data->drs_k23[threadID];
    double *brem_brem = data->cs_brem_brem[threadID];
    double *rbrem_brem = data->dcs_brem_brem[threadID];
    double *ceHeI_ceHeI = data->cs_ceHeI_ceHeI[threadID];
    double *rceHeI_ceHeI = data->dcs_ceHeI_ceHeI[threadID];
    double *ceHeII_ceHeII = data->cs_ceHeII_ceHeII[threadID];
    double *rceHeII_ceHeII = data->dcs_ceHeII_ceHeII[threadID];
    double *ceHI_ceHI = data->cs_ceHI_ceHI[threadID];
    double *rceHI_ceHI = data->dcs_ceHI_ceHI[threadID];
    double *cie_cooling_cieco = data->cs_cie_cooling_cieco[threadID];
    double *rcie_cooling_cieco = data->dcs_cie_cooling_cieco[threadID];
    double *ciHeI_ciHeI = data->cs_ciHeI_ciHeI[threadID];
    double *rciHeI_ciHeI = data->dcs_ciHeI_ciHeI[threadID];
    double *ciHeII_ciHeII = data->cs_ciHeII_ciHeII[threadID];
    double *rciHeII_ciHeII = data->dcs_ciHeII_ciHeII[threadID];
    double *ciHeIS_ciHeIS = data->cs_ciHeIS_ciHeIS[threadID];
    double *rciHeIS_ciHeIS = data->dcs_ciHeIS_ciHeIS[threadID];
    double *ciHI_ciHI = data->cs_ciHI_ciHI[threadID];
    double *rciHI_ciHI = data->dcs_ciHI_ciHI[threadID];
    double *compton_comp_ = data->cs_compton_comp_[threadID];
    double *rcompton_comp_ = data->dcs_compton_comp_[threadID];
    double *gammah_gammah = data->cs_gammah_gammah[threadID];
    double *rgammah_gammah = data->dcs_gammah_gammah[threadID];
    double *gloverabel08_gael = data->cs_gloverabel08_gael[threadID];
    double *rgloverabel08_gael = data->dcs_gloverabel08_gael[threadID];
    double *gloverabel08_gaH2 = data->cs_gloverabel08_gaH2[threadID];
    double *rgloverabel08_gaH2 = data->dcs_gloverabel08_gaH2[threadID];
    double *gloverabel08_gaHe = data->cs_gloverabel08_gaHe[threadID];
    double *rgloverabel08_gaHe = data->dcs_gloverabel08_gaHe[threadID];
    double *gloverabel08_gaHI = data->cs_gloverabel08_gaHI[threadID];
    double *rgloverabel08_gaHI = data->dcs_gloverabel08_gaHI[threadID];
    double *gloverabel08_gaHp = data->cs_gloverabel08_gaHp[threadID];
    double *rgloverabel08_gaHp = data->dcs_gloverabel08_gaHp[threadID];
    double *gloverabel08_gphdl = data->cs_gloverabel08_gphdl[threadID];
    double *rgloverabel08_gphdl = data->dcs_gloverabel08_gphdl[threadID];
    double *gloverabel08_gpldl = data->cs_gloverabel08_gpldl[threadID];
    double *rgloverabel08_gpldl = data->dcs_gloverabel08_gpldl[threadID];
    double *gloverabel08_h2lte = data->cs_gloverabel08_h2lte[threadID];
    double *rgloverabel08_h2lte = data->dcs_gloverabel08_h2lte[threadID];
    double *h2formation_h2mcool = data->cs_h2formation_h2mcool[threadID];
    double *rh2formation_h2mcool = data->dcs_h2formation_h2mcool[threadID];
    double *h2formation_h2mheat = data->cs_h2formation_h2mheat[threadID];
    double *rh2formation_h2mheat = data->dcs_h2formation_h2mheat[threadID];
    double *h2formation_ncrd1 = data->cs_h2formation_ncrd1[threadID];
    double *rh2formation_ncrd1 = data->dcs_h2formation_ncrd1[threadID];
    double *h2formation_ncrd2 = data->cs_h2formation_ncrd2[threadID];
    double *rh2formation_ncrd2 = data->dcs_h2formation_ncrd2[threadID];
    double *h2formation_ncrn = data->cs_h2formation_ncrn[threadID];
    double *rh2formation_ncrn = data->dcs_h2formation_ncrn[threadID];
    double *reHeII1_reHeII1 = data->cs_reHeII1_reHeII1[threadID];
    double *rreHeII1_reHeII1 = data->dcs_reHeII1_reHeII1[threadID];
    double *reHeII2_reHeII2 = data->cs_reHeII2_reHeII2[threadID];
    double *rreHeII2_reHeII2 = data->dcs_reHeII2_reHeII2[threadID];
    double *reHeIII_reHeIII = data->cs_reHeIII_reHeIII[threadID];
    double *rreHeIII_reHeIII = data->dcs_reHeIII_reHeIII[threadID];
    double *reHII_reHII = data->cs_reHII_reHII[threadID];
    double *rreHII_reHII = data->dcs_reHII_reHII[threadID];
    double H2_1;
    double H2_2;
    double H_1;
    double H_2;
    double H_m0;
    double He_1;
    double He_2;
    double He_3;
    double de;
    double ge;
    double z;
    double T;

    double mh = 1.66054e-24;
    double mdensity, inv_mdensity;
    
    double h2_optical_depth_approx;  
    
    
    double cie_optical_depth_approx; 
    
  

    j = 0;
    mdensity = 0.0;
    z = data->current_z;

    for ( i = 0; i < nstrip; i++ ){
        j = i * nchem;
        H2_1 = y_arr[j];
        j++;
        H2_2 = y_arr[j];
        j++;
        H_1 = y_arr[j];
        j++;
        H_2 = y_arr[j];
        j++;
        H_m0 = y_arr[j];
        j++;
        He_1 = y_arr[j];
        j++;
        He_2 = y_arr[j];
        j++;
        He_3 = y_arr[j];
        j++;
        de = y_arr[j];
        j++;
        ge = y_arr[j];
        j++;
	
        mdensity = data->mdensity[threadID][i];
        inv_mdensity = 1.0 / mdensity;
        
        h2_optical_depth_approx  = data->h2_optical_depth_approx[threadID][i];  
        
        
        cie_optical_depth_approx = data->cie_optical_depth_approx[threadID][i]; 
        

        j = i * nchem;
        //
        // Species: H2_1
        //
        
        
        // H2_1 by H2_1
        
        SM_ELEMENT_D( J, j + 0, j + 0 ) = -k11[i]*H_2 - k12[i]*de - k13[i]*H_1 + k21[i]*pow(H_1, 2) - 2*k23[i]*H2_1;
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 0];
        scale1     = scale    [ j + 0];
        SM_ELEMENT_D( J, j + 0, j + 0) *= inv_scale2*scale1;
        #endif
        
        // H2_1 by H2_2
        
        SM_ELEMENT_D( J, j + 0, j + 1 ) = k10[i]*H_1 + k19[i]*H_m0;
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 0];
        scale1     = scale    [ j + 1];
        SM_ELEMENT_D( J, j + 0, j + 1) *= inv_scale2*scale1;
        #endif
        
        // H2_1 by H_1
        
        SM_ELEMENT_D( J, j + 0, j + 2 ) = k08[i]*H_m0 + k10[i]*H2_2 - k13[i]*H2_1 + 2*k21[i]*H2_1*H_1 + 3*k22[i]*pow(H_1, 2);
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 0];
        scale1     = scale    [ j + 2];
        SM_ELEMENT_D( J, j + 0, j + 2) *= inv_scale2*scale1;
        #endif
        
        // H2_1 by H_2
        
        SM_ELEMENT_D( J, j + 0, j + 3 ) = -k11[i]*H2_1;
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 0];
        scale1     = scale    [ j + 3];
        SM_ELEMENT_D( J, j + 0, j + 3) *= inv_scale2*scale1;
        #endif
        
        // H2_1 by H_m0
        
        SM_ELEMENT_D( J, j + 0, j + 4 ) = k08[i]*H_1 + k19[i]*H2_2;
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 0];
        scale1     = scale    [ j + 4];
        SM_ELEMENT_D( J, j + 0, j + 4) *= inv_scale2*scale1;
        #endif
        
        // H2_1 by He_1
        
        // because the Jacobian is initialized to zeros by default
        // SM_ELEMENT_D( J, j + 0, j + 5 ) = 0;
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 0];
        scale1     = scale    [ j + 5];
        SM_ELEMENT_D( J, j + 0, j + 5) *= inv_scale2*scale1;
        #endif
        
        // H2_1 by He_2
        
        // because the Jacobian is initialized to zeros by default
        // SM_ELEMENT_D( J, j + 0, j + 6 ) = 0;
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 0];
        scale1     = scale    [ j + 6];
        SM_ELEMENT_D( J, j + 0, j + 6) *= inv_scale2*scale1;
        #endif
        
        // H2_1 by He_3
        
        // because the Jacobian is initialized to zeros by default
        // SM_ELEMENT_D( J, j + 0, j + 7 ) = 0;
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 0];
        scale1     = scale    [ j + 7];
        SM_ELEMENT_D( J, j + 0, j + 7) *= inv_scale2*scale1;
        #endif
        
        // H2_1 by de
        
        SM_ELEMENT_D( J, j + 0, j + 8 ) = -k12[i]*H2_1;
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 0];
        scale1     = scale    [ j + 8];
        SM_ELEMENT_D( J, j + 0, j + 8) *= inv_scale2*scale1;
        #endif
        
        // H2_1 by ge
        
        SM_ELEMENT_D( J, j + 0, j + 9 ) = -pow(H2_1, 2)*rk23[i] + H2_1*pow(H_1, 2)*rk21[i] - H2_1*H_1*rk13[i] - H2_1*H_2*rk11[i] - H2_1*de*rk12[i] + H2_2*H_1*rk10[i] + H2_2*H_m0*rk19[i] + pow(H_1, 3)*rk22[i] + H_1*H_m0*rk08[i];
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 0];
        scale1     = scale    [ j + 9];
        SM_ELEMENT_D( J, j + 0, j + 9) *= inv_scale2*scale1;
        #endif
        SM_ELEMENT_D( J, j + 0, j + 9 ) *= Tge[i];
        //
        // Species: H2_2
        //
        
        
        // H2_2 by H2_1
        
        SM_ELEMENT_D( J, j + 1, j + 0 ) = k11[i]*H_2;
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 1];
        scale1     = scale    [ j + 0];
        SM_ELEMENT_D( J, j + 1, j + 0) *= inv_scale2*scale1;
        #endif
        
        // H2_2 by H2_2
        
        SM_ELEMENT_D( J, j + 1, j + 1 ) = -k10[i]*H_1 - k18[i]*de - k19[i]*H_m0;
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 1];
        scale1     = scale    [ j + 1];
        SM_ELEMENT_D( J, j + 1, j + 1) *= inv_scale2*scale1;
        #endif
        
        // H2_2 by H_1
        
        SM_ELEMENT_D( J, j + 1, j + 2 ) = k09[i]*H_2 - k10[i]*H2_2;
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 1];
        scale1     = scale    [ j + 2];
        SM_ELEMENT_D( J, j + 1, j + 2) *= inv_scale2*scale1;
        #endif
        
        // H2_2 by H_2
        
        SM_ELEMENT_D( J, j + 1, j + 3 ) = k09[i]*H_1 + k11[i]*H2_1 + k17[i]*H_m0;
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 1];
        scale1     = scale    [ j + 3];
        SM_ELEMENT_D( J, j + 1, j + 3) *= inv_scale2*scale1;
        #endif
        
        // H2_2 by H_m0
        
        SM_ELEMENT_D( J, j + 1, j + 4 ) = k17[i]*H_2 - k19[i]*H2_2;
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 1];
        scale1     = scale    [ j + 4];
        SM_ELEMENT_D( J, j + 1, j + 4) *= inv_scale2*scale1;
        #endif
        
        // H2_2 by He_1
        
        // because the Jacobian is initialized to zeros by default
        // SM_ELEMENT_D( J, j + 1, j + 5 ) = 0;
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 1];
        scale1     = scale    [ j + 5];
        SM_ELEMENT_D( J, j + 1, j + 5) *= inv_scale2*scale1;
        #endif
        
        // H2_2 by He_2
        
        // because the Jacobian is initialized to zeros by default
        // SM_ELEMENT_D( J, j + 1, j + 6 ) = 0;
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 1];
        scale1     = scale    [ j + 6];
        SM_ELEMENT_D( J, j + 1, j + 6) *= inv_scale2*scale1;
        #endif
        
        // H2_2 by He_3
        
        // because the Jacobian is initialized to zeros by default
        // SM_ELEMENT_D( J, j + 1, j + 7 ) = 0;
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 1];
        scale1     = scale    [ j + 7];
        SM_ELEMENT_D( J, j + 1, j + 7) *= inv_scale2*scale1;
        #endif
        
        // H2_2 by de
        
        SM_ELEMENT_D( J, j + 1, j + 8 ) = -k18[i]*H2_2;
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 1];
        scale1     = scale    [ j + 8];
        SM_ELEMENT_D( J, j + 1, j + 8) *= inv_scale2*scale1;
        #endif
        
        // H2_2 by ge
        
        SM_ELEMENT_D( J, j + 1, j + 9 ) = H2_1*H_2*rk11[i] - H2_2*H_1*rk10[i] - H2_2*H_m0*rk19[i] - H2_2*de*rk18[i] + H_1*H_2*rk09[i] + H_2*H_m0*rk17[i];
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 1];
        scale1     = scale    [ j + 9];
        SM_ELEMENT_D( J, j + 1, j + 9) *= inv_scale2*scale1;
        #endif
        SM_ELEMENT_D( J, j + 1, j + 9 ) *= Tge[i];
        //
        // Species: H_1
        //
        
        
        // H_1 by H2_1
        
        SM_ELEMENT_D( J, j + 2, j + 0 ) = k11[i]*H_2 + 2*k12[i]*de + 2*k13[i]*H_1 - 2*k21[i]*pow(H_1, 2) + 4*k23[i]*H2_1;
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 2];
        scale1     = scale    [ j + 0];
        SM_ELEMENT_D( J, j + 2, j + 0) *= inv_scale2*scale1;
        #endif
        
        // H_1 by H2_2
        
        SM_ELEMENT_D( J, j + 2, j + 1 ) = -k10[i]*H_1 + 2*k18[i]*de + k19[i]*H_m0;
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 2];
        scale1     = scale    [ j + 1];
        SM_ELEMENT_D( J, j + 2, j + 1) *= inv_scale2*scale1;
        #endif
        
        // H_1 by H_1
        
        SM_ELEMENT_D( J, j + 2, j + 2 ) = -k01[i]*de - k07[i]*de - k08[i]*H_m0 - k09[i]*H_2 - k10[i]*H2_2 + 2*k13[i]*H2_1 + k15[i]*H_m0 - 4*k21[i]*H2_1*H_1 - 6*k22[i]*pow(H_1, 2);
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 2];
        scale1     = scale    [ j + 2];
        SM_ELEMENT_D( J, j + 2, j + 2) *= inv_scale2*scale1;
        #endif
        
        // H_1 by H_2
        
        SM_ELEMENT_D( J, j + 2, j + 3 ) = k02[i]*de - k09[i]*H_1 + k11[i]*H2_1 + 2*k16[i]*H_m0;
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 2];
        scale1     = scale    [ j + 3];
        SM_ELEMENT_D( J, j + 2, j + 3) *= inv_scale2*scale1;
        #endif
        
        // H_1 by H_m0
        
        SM_ELEMENT_D( J, j + 2, j + 4 ) = -k08[i]*H_1 + k14[i]*de + k15[i]*H_1 + 2*k16[i]*H_2 + k19[i]*H2_2;
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 2];
        scale1     = scale    [ j + 4];
        SM_ELEMENT_D( J, j + 2, j + 4) *= inv_scale2*scale1;
        #endif
        
        // H_1 by He_1
        
        // because the Jacobian is initialized to zeros by default
        // SM_ELEMENT_D( J, j + 2, j + 5 ) = 0;
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 2];
        scale1     = scale    [ j + 5];
        SM_ELEMENT_D( J, j + 2, j + 5) *= inv_scale2*scale1;
        #endif
        
        // H_1 by He_2
        
        // because the Jacobian is initialized to zeros by default
        // SM_ELEMENT_D( J, j + 2, j + 6 ) = 0;
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 2];
        scale1     = scale    [ j + 6];
        SM_ELEMENT_D( J, j + 2, j + 6) *= inv_scale2*scale1;
        #endif
        
        // H_1 by He_3
        
        // because the Jacobian is initialized to zeros by default
        // SM_ELEMENT_D( J, j + 2, j + 7 ) = 0;
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 2];
        scale1     = scale    [ j + 7];
        SM_ELEMENT_D( J, j + 2, j + 7) *= inv_scale2*scale1;
        #endif
        
        // H_1 by de
        
        SM_ELEMENT_D( J, j + 2, j + 8 ) = -k01[i]*H_1 + k02[i]*H_2 - k07[i]*H_1 + 2*k12[i]*H2_1 + k14[i]*H_m0 + 2*k18[i]*H2_2;
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 2];
        scale1     = scale    [ j + 8];
        SM_ELEMENT_D( J, j + 2, j + 8) *= inv_scale2*scale1;
        #endif
        
        // H_1 by ge
        
        SM_ELEMENT_D( J, j + 2, j + 9 ) = 2*pow(H2_1, 2)*rk23[i] - 2*H2_1*pow(H_1, 2)*rk21[i] + 2*H2_1*H_1*rk13[i] + H2_1*H_2*rk11[i] + 2*H2_1*de*rk12[i] - H2_2*H_1*rk10[i] + H2_2*H_m0*rk19[i] + 2*H2_2*de*rk18[i] - 2*pow(H_1, 3)*rk22[i] - H_1*H_2*rk09[i] - H_1*H_m0*rk08[i] + H_1*H_m0*rk15[i] - H_1*de*rk01[i] - H_1*de*rk07[i] + 2*H_2*H_m0*rk16[i] + H_2*de*rk02[i] + H_m0*de*rk14[i];
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 2];
        scale1     = scale    [ j + 9];
        SM_ELEMENT_D( J, j + 2, j + 9) *= inv_scale2*scale1;
        #endif
        SM_ELEMENT_D( J, j + 2, j + 9 ) *= Tge[i];
        //
        // Species: H_2
        //
        
        
        // H_2 by H2_1
        
        SM_ELEMENT_D( J, j + 3, j + 0 ) = -k11[i]*H_2;
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 3];
        scale1     = scale    [ j + 0];
        SM_ELEMENT_D( J, j + 3, j + 0) *= inv_scale2*scale1;
        #endif
        
        // H_2 by H2_2
        
        SM_ELEMENT_D( J, j + 3, j + 1 ) = k10[i]*H_1;
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 3];
        scale1     = scale    [ j + 1];
        SM_ELEMENT_D( J, j + 3, j + 1) *= inv_scale2*scale1;
        #endif
        
        // H_2 by H_1
        
        SM_ELEMENT_D( J, j + 3, j + 2 ) = k01[i]*de - k09[i]*H_2 + k10[i]*H2_2;
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 3];
        scale1     = scale    [ j + 2];
        SM_ELEMENT_D( J, j + 3, j + 2) *= inv_scale2*scale1;
        #endif
        
        // H_2 by H_2
        
        SM_ELEMENT_D( J, j + 3, j + 3 ) = -k02[i]*de - k09[i]*H_1 - k11[i]*H2_1 - k16[i]*H_m0 - k17[i]*H_m0;
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 3];
        scale1     = scale    [ j + 3];
        SM_ELEMENT_D( J, j + 3, j + 3) *= inv_scale2*scale1;
        #endif
        
        // H_2 by H_m0
        
        SM_ELEMENT_D( J, j + 3, j + 4 ) = -k16[i]*H_2 - k17[i]*H_2;
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 3];
        scale1     = scale    [ j + 4];
        SM_ELEMENT_D( J, j + 3, j + 4) *= inv_scale2*scale1;
        #endif
        
        // H_2 by He_1
        
        // because the Jacobian is initialized to zeros by default
        // SM_ELEMENT_D( J, j + 3, j + 5 ) = 0;
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 3];
        scale1     = scale    [ j + 5];
        SM_ELEMENT_D( J, j + 3, j + 5) *= inv_scale2*scale1;
        #endif
        
        // H_2 by He_2
        
        // because the Jacobian is initialized to zeros by default
        // SM_ELEMENT_D( J, j + 3, j + 6 ) = 0;
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 3];
        scale1     = scale    [ j + 6];
        SM_ELEMENT_D( J, j + 3, j + 6) *= inv_scale2*scale1;
        #endif
        
        // H_2 by He_3
        
        // because the Jacobian is initialized to zeros by default
        // SM_ELEMENT_D( J, j + 3, j + 7 ) = 0;
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 3];
        scale1     = scale    [ j + 7];
        SM_ELEMENT_D( J, j + 3, j + 7) *= inv_scale2*scale1;
        #endif
        
        // H_2 by de
        
        SM_ELEMENT_D( J, j + 3, j + 8 ) = k01[i]*H_1 - k02[i]*H_2;
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 3];
        scale1     = scale    [ j + 8];
        SM_ELEMENT_D( J, j + 3, j + 8) *= inv_scale2*scale1;
        #endif
        
        // H_2 by ge
        
        SM_ELEMENT_D( J, j + 3, j + 9 ) = -H2_1*H_2*rk11[i] + H2_2*H_1*rk10[i] - H_1*H_2*rk09[i] + H_1*de*rk01[i] - H_2*H_m0*rk16[i] - H_2*H_m0*rk17[i] - H_2*de*rk02[i];
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 3];
        scale1     = scale    [ j + 9];
        SM_ELEMENT_D( J, j + 3, j + 9) *= inv_scale2*scale1;
        #endif
        SM_ELEMENT_D( J, j + 3, j + 9 ) *= Tge[i];
        //
        // Species: H_m0
        //
        
        
        // H_m0 by H2_1
        
        // because the Jacobian is initialized to zeros by default
        // SM_ELEMENT_D( J, j + 4, j + 0 ) = 0;
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 4];
        scale1     = scale    [ j + 0];
        SM_ELEMENT_D( J, j + 4, j + 0) *= inv_scale2*scale1;
        #endif
        
        // H_m0 by H2_2
        
        SM_ELEMENT_D( J, j + 4, j + 1 ) = -k19[i]*H_m0;
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 4];
        scale1     = scale    [ j + 1];
        SM_ELEMENT_D( J, j + 4, j + 1) *= inv_scale2*scale1;
        #endif
        
        // H_m0 by H_1
        
        SM_ELEMENT_D( J, j + 4, j + 2 ) = k07[i]*de - k08[i]*H_m0 - k15[i]*H_m0;
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 4];
        scale1     = scale    [ j + 2];
        SM_ELEMENT_D( J, j + 4, j + 2) *= inv_scale2*scale1;
        #endif
        
        // H_m0 by H_2
        
        SM_ELEMENT_D( J, j + 4, j + 3 ) = -k16[i]*H_m0 - k17[i]*H_m0;
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 4];
        scale1     = scale    [ j + 3];
        SM_ELEMENT_D( J, j + 4, j + 3) *= inv_scale2*scale1;
        #endif
        
        // H_m0 by H_m0
        
        SM_ELEMENT_D( J, j + 4, j + 4 ) = -k08[i]*H_1 - k14[i]*de - k15[i]*H_1 - k16[i]*H_2 - k17[i]*H_2 - k19[i]*H2_2;
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 4];
        scale1     = scale    [ j + 4];
        SM_ELEMENT_D( J, j + 4, j + 4) *= inv_scale2*scale1;
        #endif
        
        // H_m0 by He_1
        
        // because the Jacobian is initialized to zeros by default
        // SM_ELEMENT_D( J, j + 4, j + 5 ) = 0;
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 4];
        scale1     = scale    [ j + 5];
        SM_ELEMENT_D( J, j + 4, j + 5) *= inv_scale2*scale1;
        #endif
        
        // H_m0 by He_2
        
        // because the Jacobian is initialized to zeros by default
        // SM_ELEMENT_D( J, j + 4, j + 6 ) = 0;
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 4];
        scale1     = scale    [ j + 6];
        SM_ELEMENT_D( J, j + 4, j + 6) *= inv_scale2*scale1;
        #endif
        
        // H_m0 by He_3
        
        // because the Jacobian is initialized to zeros by default
        // SM_ELEMENT_D( J, j + 4, j + 7 ) = 0;
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 4];
        scale1     = scale    [ j + 7];
        SM_ELEMENT_D( J, j + 4, j + 7) *= inv_scale2*scale1;
        #endif
        
        // H_m0 by de
        
        SM_ELEMENT_D( J, j + 4, j + 8 ) = k07[i]*H_1 - k14[i]*H_m0;
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 4];
        scale1     = scale    [ j + 8];
        SM_ELEMENT_D( J, j + 4, j + 8) *= inv_scale2*scale1;
        #endif
        
        // H_m0 by ge
        
        SM_ELEMENT_D( J, j + 4, j + 9 ) = -H2_2*H_m0*rk19[i] - H_1*H_m0*rk08[i] - H_1*H_m0*rk15[i] + H_1*de*rk07[i] - H_2*H_m0*rk16[i] - H_2*H_m0*rk17[i] - H_m0*de*rk14[i];
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 4];
        scale1     = scale    [ j + 9];
        SM_ELEMENT_D( J, j + 4, j + 9) *= inv_scale2*scale1;
        #endif
        SM_ELEMENT_D( J, j + 4, j + 9 ) *= Tge[i];
        //
        // Species: He_1
        //
        
        
        // He_1 by H2_1
        
        // because the Jacobian is initialized to zeros by default
        // SM_ELEMENT_D( J, j + 5, j + 0 ) = 0;
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 5];
        scale1     = scale    [ j + 0];
        SM_ELEMENT_D( J, j + 5, j + 0) *= inv_scale2*scale1;
        #endif
        
        // He_1 by H2_2
        
        // because the Jacobian is initialized to zeros by default
        // SM_ELEMENT_D( J, j + 5, j + 1 ) = 0;
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 5];
        scale1     = scale    [ j + 1];
        SM_ELEMENT_D( J, j + 5, j + 1) *= inv_scale2*scale1;
        #endif
        
        // He_1 by H_1
        
        // because the Jacobian is initialized to zeros by default
        // SM_ELEMENT_D( J, j + 5, j + 2 ) = 0;
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 5];
        scale1     = scale    [ j + 2];
        SM_ELEMENT_D( J, j + 5, j + 2) *= inv_scale2*scale1;
        #endif
        
        // He_1 by H_2
        
        // because the Jacobian is initialized to zeros by default
        // SM_ELEMENT_D( J, j + 5, j + 3 ) = 0;
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 5];
        scale1     = scale    [ j + 3];
        SM_ELEMENT_D( J, j + 5, j + 3) *= inv_scale2*scale1;
        #endif
        
        // He_1 by H_m0
        
        // because the Jacobian is initialized to zeros by default
        // SM_ELEMENT_D( J, j + 5, j + 4 ) = 0;
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 5];
        scale1     = scale    [ j + 4];
        SM_ELEMENT_D( J, j + 5, j + 4) *= inv_scale2*scale1;
        #endif
        
        // He_1 by He_1
        
        SM_ELEMENT_D( J, j + 5, j + 5 ) = -k03[i]*de;
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 5];
        scale1     = scale    [ j + 5];
        SM_ELEMENT_D( J, j + 5, j + 5) *= inv_scale2*scale1;
        #endif
        
        // He_1 by He_2
        
        SM_ELEMENT_D( J, j + 5, j + 6 ) = k04[i]*de;
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 5];
        scale1     = scale    [ j + 6];
        SM_ELEMENT_D( J, j + 5, j + 6) *= inv_scale2*scale1;
        #endif
        
        // He_1 by He_3
        
        // because the Jacobian is initialized to zeros by default
        // SM_ELEMENT_D( J, j + 5, j + 7 ) = 0;
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 5];
        scale1     = scale    [ j + 7];
        SM_ELEMENT_D( J, j + 5, j + 7) *= inv_scale2*scale1;
        #endif
        
        // He_1 by de
        
        SM_ELEMENT_D( J, j + 5, j + 8 ) = -k03[i]*He_1 + k04[i]*He_2;
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 5];
        scale1     = scale    [ j + 8];
        SM_ELEMENT_D( J, j + 5, j + 8) *= inv_scale2*scale1;
        #endif
        
        // He_1 by ge
        
        SM_ELEMENT_D( J, j + 5, j + 9 ) = -He_1*de*rk03[i] + He_2*de*rk04[i];
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 5];
        scale1     = scale    [ j + 9];
        SM_ELEMENT_D( J, j + 5, j + 9) *= inv_scale2*scale1;
        #endif
        SM_ELEMENT_D( J, j + 5, j + 9 ) *= Tge[i];
        //
        // Species: He_2
        //
        
        
        // He_2 by H2_1
        
        // because the Jacobian is initialized to zeros by default
        // SM_ELEMENT_D( J, j + 6, j + 0 ) = 0;
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 6];
        scale1     = scale    [ j + 0];
        SM_ELEMENT_D( J, j + 6, j + 0) *= inv_scale2*scale1;
        #endif
        
        // He_2 by H2_2
        
        // because the Jacobian is initialized to zeros by default
        // SM_ELEMENT_D( J, j + 6, j + 1 ) = 0;
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 6];
        scale1     = scale    [ j + 1];
        SM_ELEMENT_D( J, j + 6, j + 1) *= inv_scale2*scale1;
        #endif
        
        // He_2 by H_1
        
        // because the Jacobian is initialized to zeros by default
        // SM_ELEMENT_D( J, j + 6, j + 2 ) = 0;
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 6];
        scale1     = scale    [ j + 2];
        SM_ELEMENT_D( J, j + 6, j + 2) *= inv_scale2*scale1;
        #endif
        
        // He_2 by H_2
        
        // because the Jacobian is initialized to zeros by default
        // SM_ELEMENT_D( J, j + 6, j + 3 ) = 0;
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 6];
        scale1     = scale    [ j + 3];
        SM_ELEMENT_D( J, j + 6, j + 3) *= inv_scale2*scale1;
        #endif
        
        // He_2 by H_m0
        
        // because the Jacobian is initialized to zeros by default
        // SM_ELEMENT_D( J, j + 6, j + 4 ) = 0;
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 6];
        scale1     = scale    [ j + 4];
        SM_ELEMENT_D( J, j + 6, j + 4) *= inv_scale2*scale1;
        #endif
        
        // He_2 by He_1
        
        SM_ELEMENT_D( J, j + 6, j + 5 ) = k03[i]*de;
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 6];
        scale1     = scale    [ j + 5];
        SM_ELEMENT_D( J, j + 6, j + 5) *= inv_scale2*scale1;
        #endif
        
        // He_2 by He_2
        
        SM_ELEMENT_D( J, j + 6, j + 6 ) = -k04[i]*de - k05[i]*de;
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 6];
        scale1     = scale    [ j + 6];
        SM_ELEMENT_D( J, j + 6, j + 6) *= inv_scale2*scale1;
        #endif
        
        // He_2 by He_3
        
        SM_ELEMENT_D( J, j + 6, j + 7 ) = k06[i]*de;
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 6];
        scale1     = scale    [ j + 7];
        SM_ELEMENT_D( J, j + 6, j + 7) *= inv_scale2*scale1;
        #endif
        
        // He_2 by de
        
        SM_ELEMENT_D( J, j + 6, j + 8 ) = k03[i]*He_1 - k04[i]*He_2 - k05[i]*He_2 + k06[i]*He_3;
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 6];
        scale1     = scale    [ j + 8];
        SM_ELEMENT_D( J, j + 6, j + 8) *= inv_scale2*scale1;
        #endif
        
        // He_2 by ge
        
        SM_ELEMENT_D( J, j + 6, j + 9 ) = He_1*de*rk03[i] - He_2*de*rk04[i] - He_2*de*rk05[i] + He_3*de*rk06[i];
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 6];
        scale1     = scale    [ j + 9];
        SM_ELEMENT_D( J, j + 6, j + 9) *= inv_scale2*scale1;
        #endif
        SM_ELEMENT_D( J, j + 6, j + 9 ) *= Tge[i];
        //
        // Species: He_3
        //
        
        
        // He_3 by H2_1
        
        // because the Jacobian is initialized to zeros by default
        // SM_ELEMENT_D( J, j + 7, j + 0 ) = 0;
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 7];
        scale1     = scale    [ j + 0];
        SM_ELEMENT_D( J, j + 7, j + 0) *= inv_scale2*scale1;
        #endif
        
        // He_3 by H2_2
        
        // because the Jacobian is initialized to zeros by default
        // SM_ELEMENT_D( J, j + 7, j + 1 ) = 0;
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 7];
        scale1     = scale    [ j + 1];
        SM_ELEMENT_D( J, j + 7, j + 1) *= inv_scale2*scale1;
        #endif
        
        // He_3 by H_1
        
        // because the Jacobian is initialized to zeros by default
        // SM_ELEMENT_D( J, j + 7, j + 2 ) = 0;
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 7];
        scale1     = scale    [ j + 2];
        SM_ELEMENT_D( J, j + 7, j + 2) *= inv_scale2*scale1;
        #endif
        
        // He_3 by H_2
        
        // because the Jacobian is initialized to zeros by default
        // SM_ELEMENT_D( J, j + 7, j + 3 ) = 0;
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 7];
        scale1     = scale    [ j + 3];
        SM_ELEMENT_D( J, j + 7, j + 3) *= inv_scale2*scale1;
        #endif
        
        // He_3 by H_m0
        
        // because the Jacobian is initialized to zeros by default
        // SM_ELEMENT_D( J, j + 7, j + 4 ) = 0;
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 7];
        scale1     = scale    [ j + 4];
        SM_ELEMENT_D( J, j + 7, j + 4) *= inv_scale2*scale1;
        #endif
        
        // He_3 by He_1
        
        // because the Jacobian is initialized to zeros by default
        // SM_ELEMENT_D( J, j + 7, j + 5 ) = 0;
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 7];
        scale1     = scale    [ j + 5];
        SM_ELEMENT_D( J, j + 7, j + 5) *= inv_scale2*scale1;
        #endif
        
        // He_3 by He_2
        
        SM_ELEMENT_D( J, j + 7, j + 6 ) = k05[i]*de;
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 7];
        scale1     = scale    [ j + 6];
        SM_ELEMENT_D( J, j + 7, j + 6) *= inv_scale2*scale1;
        #endif
        
        // He_3 by He_3
        
        SM_ELEMENT_D( J, j + 7, j + 7 ) = -k06[i]*de;
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 7];
        scale1     = scale    [ j + 7];
        SM_ELEMENT_D( J, j + 7, j + 7) *= inv_scale2*scale1;
        #endif
        
        // He_3 by de
        
        SM_ELEMENT_D( J, j + 7, j + 8 ) = k05[i]*He_2 - k06[i]*He_3;
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 7];
        scale1     = scale    [ j + 8];
        SM_ELEMENT_D( J, j + 7, j + 8) *= inv_scale2*scale1;
        #endif
        
        // He_3 by ge
        
        SM_ELEMENT_D( J, j + 7, j + 9 ) = He_2*de*rk05[i] - He_3*de*rk06[i];
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 7];
        scale1     = scale    [ j + 9];
        SM_ELEMENT_D( J, j + 7, j + 9) *= inv_scale2*scale1;
        #endif
        SM_ELEMENT_D( J, j + 7, j + 9 ) *= Tge[i];
        //
        // Species: de
        //
        
        
        // de by H2_1
        
        // because the Jacobian is initialized to zeros by default
        // SM_ELEMENT_D( J, j + 8, j + 0 ) = 0;
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 8];
        scale1     = scale    [ j + 0];
        SM_ELEMENT_D( J, j + 8, j + 0) *= inv_scale2*scale1;
        #endif
        
        // de by H2_2
        
        SM_ELEMENT_D( J, j + 8, j + 1 ) = -k18[i]*de;
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 8];
        scale1     = scale    [ j + 1];
        SM_ELEMENT_D( J, j + 8, j + 1) *= inv_scale2*scale1;
        #endif
        
        // de by H_1
        
        SM_ELEMENT_D( J, j + 8, j + 2 ) = k01[i]*de - k07[i]*de + k08[i]*H_m0 + k15[i]*H_m0;
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 8];
        scale1     = scale    [ j + 2];
        SM_ELEMENT_D( J, j + 8, j + 2) *= inv_scale2*scale1;
        #endif
        
        // de by H_2
        
        SM_ELEMENT_D( J, j + 8, j + 3 ) = -k02[i]*de + k17[i]*H_m0;
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 8];
        scale1     = scale    [ j + 3];
        SM_ELEMENT_D( J, j + 8, j + 3) *= inv_scale2*scale1;
        #endif
        
        // de by H_m0
        
        SM_ELEMENT_D( J, j + 8, j + 4 ) = k08[i]*H_1 + k14[i]*de + k15[i]*H_1 + k17[i]*H_2;
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 8];
        scale1     = scale    [ j + 4];
        SM_ELEMENT_D( J, j + 8, j + 4) *= inv_scale2*scale1;
        #endif
        
        // de by He_1
        
        SM_ELEMENT_D( J, j + 8, j + 5 ) = k03[i]*de;
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 8];
        scale1     = scale    [ j + 5];
        SM_ELEMENT_D( J, j + 8, j + 5) *= inv_scale2*scale1;
        #endif
        
        // de by He_2
        
        SM_ELEMENT_D( J, j + 8, j + 6 ) = -k04[i]*de + k05[i]*de;
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 8];
        scale1     = scale    [ j + 6];
        SM_ELEMENT_D( J, j + 8, j + 6) *= inv_scale2*scale1;
        #endif
        
        // de by He_3
        
        SM_ELEMENT_D( J, j + 8, j + 7 ) = -k06[i]*de;
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 8];
        scale1     = scale    [ j + 7];
        SM_ELEMENT_D( J, j + 8, j + 7) *= inv_scale2*scale1;
        #endif
        
        // de by de
        
        SM_ELEMENT_D( J, j + 8, j + 8 ) = k01[i]*H_1 - k02[i]*H_2 + k03[i]*He_1 - k04[i]*He_2 + k05[i]*He_2 - k06[i]*He_3 - k07[i]*H_1 + k14[i]*H_m0 - k18[i]*H2_2;
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 8];
        scale1     = scale    [ j + 8];
        SM_ELEMENT_D( J, j + 8, j + 8) *= inv_scale2*scale1;
        #endif
        
        // de by ge
        
        SM_ELEMENT_D( J, j + 8, j + 9 ) = -H2_2*de*rk18[i] + H_1*H_m0*rk08[i] + H_1*H_m0*rk15[i] + H_1*de*rk01[i] - H_1*de*rk07[i] + H_2*H_m0*rk17[i] - H_2*de*rk02[i] + H_m0*de*rk14[i] + He_1*de*rk03[i] - He_2*de*rk04[i] + He_2*de*rk05[i] - He_3*de*rk06[i];
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 8];
        scale1     = scale    [ j + 9];
        SM_ELEMENT_D( J, j + 8, j + 9) *= inv_scale2*scale1;
        #endif
        SM_ELEMENT_D( J, j + 8, j + 9 ) *= Tge[i];
        //
        // Species: ge
        //
        
        
        // ge by H2_1
        
        SM_ELEMENT_D( J, j + 9, j + 0 ) = -2.0158800000000001*cie_cooling_cieco[i]*mdensity - gloverabel08_gaH2[i]*pow(gloverabel08_h2lte[i], 2)*H2_1*h2_optical_depth_approx/(pow(gloverabel08_h2lte[i]/(gloverabel08_gaH2[i]*H2_1 + gloverabel08_gaHI[i]*H_1 + gloverabel08_gaHe[i]*He_1 + gloverabel08_gaHp[i]*H_2 + gloverabel08_gael[i]*de) + 1.0, 2)*pow(gloverabel08_gaH2[i]*H2_1 + gloverabel08_gaHI[i]*H_1 + gloverabel08_gaHe[i]*He_1 + gloverabel08_gaHp[i]*H_2 + gloverabel08_gael[i]*de, 2)) - gloverabel08_h2lte[i]*h2_optical_depth_approx/(gloverabel08_h2lte[i]/(gloverabel08_gaH2[i]*H2_1 + gloverabel08_gaHI[i]*H_1 + gloverabel08_gaHe[i]*He_1 + gloverabel08_gaHp[i]*H_2 + gloverabel08_gael[i]*de) + 1.0) - 0.5*h2formation_h2mcool[i]*H_1*1.0/(h2formation_ncrn[i]/(h2formation_ncrd1[i]*H_1 + h2formation_ncrd2[i]*H2_1) + 1.0) + 0.5*h2formation_ncrd2[i]*h2formation_ncrn[i]*pow(h2formation_ncrn[i]/(h2formation_ncrd1[i]*H_1 + h2formation_ncrd2[i]*H2_1) + 1.0, -2.0)*(-h2formation_h2mcool[i]*H2_1*H_1 + h2formation_h2mheat[i]*pow(H_1, 3))/pow(h2formation_ncrd1[i]*H_1 + h2formation_ncrd2[i]*H2_1, 2);
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 9];
        scale1     = scale    [ j + 0];
        SM_ELEMENT_D( J, j + 9, j + 0) *= inv_scale2*scale1;
        #endif
        SM_ELEMENT_D( J, j + 9, j + 0 ) *= inv_mdensity;
        
        // ge by H2_2
        
        // because the Jacobian is initialized to zeros by default
        // SM_ELEMENT_D( J, j + 9, j + 1 ) = 0;
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 9];
        scale1     = scale    [ j + 1];
        SM_ELEMENT_D( J, j + 9, j + 1) *= inv_scale2*scale1;
        #endif
        SM_ELEMENT_D( J, j + 9, j + 1 ) *= inv_mdensity;
        
        // ge by H_1
        
        SM_ELEMENT_D( J, j + 9, j + 2 ) = -ceHI_ceHI[i]*de - ciHI_ciHI[i]*de - gloverabel08_gaHI[i]*pow(gloverabel08_h2lte[i], 2)*H2_1*h2_optical_depth_approx/(pow(gloverabel08_h2lte[i]/(gloverabel08_gaH2[i]*H2_1 + gloverabel08_gaHI[i]*H_1 + gloverabel08_gaHe[i]*He_1 + gloverabel08_gaHp[i]*H_2 + gloverabel08_gael[i]*de) + 1.0, 2)*pow(gloverabel08_gaH2[i]*H2_1 + gloverabel08_gaHI[i]*H_1 + gloverabel08_gaHe[i]*He_1 + gloverabel08_gaHp[i]*H_2 + gloverabel08_gael[i]*de, 2)) + 0.5*h2formation_ncrd1[i]*h2formation_ncrn[i]*pow(h2formation_ncrn[i]/(h2formation_ncrd1[i]*H_1 + h2formation_ncrd2[i]*H2_1) + 1.0, -2.0)*(-h2formation_h2mcool[i]*H2_1*H_1 + h2formation_h2mheat[i]*pow(H_1, 3))/pow(h2formation_ncrd1[i]*H_1 + h2formation_ncrd2[i]*H2_1, 2) + 0.5*(-h2formation_h2mcool[i]*H2_1 + 3*h2formation_h2mheat[i]*pow(H_1, 2))*1.0/(h2formation_ncrn[i]/(h2formation_ncrd1[i]*H_1 + h2formation_ncrd2[i]*H2_1) + 1.0);
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 9];
        scale1     = scale    [ j + 2];
        SM_ELEMENT_D( J, j + 9, j + 2) *= inv_scale2*scale1;
        #endif
        SM_ELEMENT_D( J, j + 9, j + 2 ) *= inv_mdensity;
        
        // ge by H_2
        
        SM_ELEMENT_D( J, j + 9, j + 3 ) = -brem_brem[i]*de - gloverabel08_gaHp[i]*pow(gloverabel08_h2lte[i], 2)*H2_1*h2_optical_depth_approx/(pow(gloverabel08_h2lte[i]/(gloverabel08_gaH2[i]*H2_1 + gloverabel08_gaHI[i]*H_1 + gloverabel08_gaHe[i]*He_1 + gloverabel08_gaHp[i]*H_2 + gloverabel08_gael[i]*de) + 1.0, 2)*pow(gloverabel08_gaH2[i]*H2_1 + gloverabel08_gaHI[i]*H_1 + gloverabel08_gaHe[i]*He_1 + gloverabel08_gaHp[i]*H_2 + gloverabel08_gael[i]*de, 2)) - reHII_reHII[i]*de;
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 9];
        scale1     = scale    [ j + 3];
        SM_ELEMENT_D( J, j + 9, j + 3) *= inv_scale2*scale1;
        #endif
        SM_ELEMENT_D( J, j + 9, j + 3 ) *= inv_mdensity;
        
        // ge by H_m0
        
        // because the Jacobian is initialized to zeros by default
        // SM_ELEMENT_D( J, j + 9, j + 4 ) = 0;
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 9];
        scale1     = scale    [ j + 4];
        SM_ELEMENT_D( J, j + 9, j + 4) *= inv_scale2*scale1;
        #endif
        SM_ELEMENT_D( J, j + 9, j + 4 ) *= inv_mdensity;
        
        // ge by He_1
        
        SM_ELEMENT_D( J, j + 9, j + 5 ) = -ciHeI_ciHeI[i]*de - gloverabel08_gaHe[i]*pow(gloverabel08_h2lte[i], 2)*H2_1*h2_optical_depth_approx/(pow(gloverabel08_h2lte[i]/(gloverabel08_gaH2[i]*H2_1 + gloverabel08_gaHI[i]*H_1 + gloverabel08_gaHe[i]*He_1 + gloverabel08_gaHp[i]*H_2 + gloverabel08_gael[i]*de) + 1.0, 2)*pow(gloverabel08_gaH2[i]*H2_1 + gloverabel08_gaHI[i]*H_1 + gloverabel08_gaHe[i]*He_1 + gloverabel08_gaHp[i]*H_2 + gloverabel08_gael[i]*de, 2));
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 9];
        scale1     = scale    [ j + 5];
        SM_ELEMENT_D( J, j + 9, j + 5) *= inv_scale2*scale1;
        #endif
        SM_ELEMENT_D( J, j + 9, j + 5 ) *= inv_mdensity;
        
        // ge by He_2
        
        SM_ELEMENT_D( J, j + 9, j + 6 ) = -brem_brem[i]*de - ceHeII_ceHeII[i]*de - ceHeI_ceHeI[i]*pow(de, 2) - ciHeII_ciHeII[i]*de - ciHeIS_ciHeIS[i]*pow(de, 2) - reHeII1_reHeII1[i]*de - reHeII2_reHeII2[i]*de;
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 9];
        scale1     = scale    [ j + 6];
        SM_ELEMENT_D( J, j + 9, j + 6) *= inv_scale2*scale1;
        #endif
        SM_ELEMENT_D( J, j + 9, j + 6 ) *= inv_mdensity;
        
        // ge by He_3
        
        SM_ELEMENT_D( J, j + 9, j + 7 ) = -4.0*brem_brem[i]*de - reHeIII_reHeIII[i]*de;
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 9];
        scale1     = scale    [ j + 7];
        SM_ELEMENT_D( J, j + 9, j + 7) *= inv_scale2*scale1;
        #endif
        SM_ELEMENT_D( J, j + 9, j + 7 ) *= inv_mdensity;
        
        // ge by de
        
        SM_ELEMENT_D( J, j + 9, j + 8 ) = brem_brem[i]*(-H_2 - He_2 - 4.0*He_3) - ceHI_ceHI[i]*H_1 - ceHeII_ceHeII[i]*He_2 - 2*ceHeI_ceHeI[i]*He_2*de - ciHI_ciHI[i]*H_1 - ciHeII_ciHeII[i]*He_2 - 2*ciHeIS_ciHeIS[i]*He_2*de - ciHeI_ciHeI[i]*He_1 - compton_comp_[i]*pow(z + 1.0, 4)*(T - 2.73*z - 2.73) - gloverabel08_gael[i]*pow(gloverabel08_h2lte[i], 2)*H2_1*h2_optical_depth_approx/(pow(gloverabel08_h2lte[i]/(gloverabel08_gaH2[i]*H2_1 + gloverabel08_gaHI[i]*H_1 + gloverabel08_gaHe[i]*He_1 + gloverabel08_gaHp[i]*H_2 + gloverabel08_gael[i]*de) + 1.0, 2)*pow(gloverabel08_gaH2[i]*H2_1 + gloverabel08_gaHI[i]*H_1 + gloverabel08_gaHe[i]*He_1 + gloverabel08_gaHp[i]*H_2 + gloverabel08_gael[i]*de, 2)) - reHII_reHII[i]*H_2 - reHeII1_reHeII1[i]*He_2 - reHeII2_reHeII2[i]*He_2 - reHeIII_reHeIII[i]*He_3;
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 9];
        scale1     = scale    [ j + 8];
        SM_ELEMENT_D( J, j + 9, j + 8) *= inv_scale2*scale1;
        #endif
        SM_ELEMENT_D( J, j + 9, j + 8 ) *= inv_mdensity;
        
        // ge by ge
        
        SM_ELEMENT_D( J, j + 9, j + 9 ) = -gloverabel08_h2lte[i]*H2_1*h2_optical_depth_approx*(-gloverabel08_h2lte[i]*(-H2_1*rgloverabel08_gaH2[i] - H_1*rgloverabel08_gaHI[i] - H_2*rgloverabel08_gaHp[i] - He_1*rgloverabel08_gaHe[i] - de*rgloverabel08_gael[i])/pow(gloverabel08_gaH2[i]*H2_1 + gloverabel08_gaHI[i]*H_1 + gloverabel08_gaHe[i]*He_1 + gloverabel08_gaHp[i]*H_2 + gloverabel08_gael[i]*de, 2) - rgloverabel08_h2lte[i]/(gloverabel08_gaH2[i]*H2_1 + gloverabel08_gaHI[i]*H_1 + gloverabel08_gaHe[i]*He_1 + gloverabel08_gaHp[i]*H_2 + gloverabel08_gael[i]*de))/pow(gloverabel08_h2lte[i]/(gloverabel08_gaH2[i]*H2_1 + gloverabel08_gaHI[i]*H_1 + gloverabel08_gaHe[i]*He_1 + gloverabel08_gaHp[i]*H_2 + gloverabel08_gael[i]*de) + 1.0, 2) - H2_1*h2_optical_depth_approx*rgloverabel08_h2lte[i]/(gloverabel08_h2lte[i]/(gloverabel08_gaH2[i]*H2_1 + gloverabel08_gaHI[i]*H_1 + gloverabel08_gaHe[i]*He_1 + gloverabel08_gaHp[i]*H_2 + gloverabel08_gael[i]*de) + 1.0) - 2.0158800000000001*H2_1*mdensity*rcie_cooling_cieco[i] - H_1*de*rceHI_ceHI[i] - H_1*de*rciHI_ciHI[i] - H_2*de*rreHII_reHII[i] - He_1*de*rciHeI_ciHeI[i] - He_2*pow(de, 2)*rceHeI_ceHeI[i] - He_2*pow(de, 2)*rciHeIS_ciHeIS[i] - He_2*de*rceHeII_ceHeII[i] - He_2*de*rciHeII_ciHeII[i] - He_2*de*rreHeII1_reHeII1[i] - He_2*de*rreHeII2_reHeII2[i] - He_3*de*rreHeIII_reHeIII[i] - de*rbrem_brem[i]*(H_2 + He_2 + 4.0*He_3) - de*rcompton_comp_[i]*pow(z + 1.0, 4)*(T - 2.73*z - 2.73) + 0.5*pow(h2formation_ncrn[i]/(h2formation_ncrd1[i]*H_1 + h2formation_ncrd2[i]*H2_1) + 1.0, -2.0)*(-h2formation_h2mcool[i]*H2_1*H_1 + h2formation_h2mheat[i]*pow(H_1, 3))*(-1.0*h2formation_ncrn[i]*(-H2_1*rh2formation_ncrd2[i] - H_1*rh2formation_ncrd1[i])/pow(h2formation_ncrd1[i]*H_1 + h2formation_ncrd2[i]*H2_1, 2) - 1.0*rh2formation_ncrn[i]/(h2formation_ncrd1[i]*H_1 + h2formation_ncrd2[i]*H2_1)) + 0.5*1.0/(h2formation_ncrn[i]/(h2formation_ncrd1[i]*H_1 + h2formation_ncrd2[i]*H2_1) + 1.0)*(-H2_1*H_1*rh2formation_h2mcool[i] + pow(H_1, 3)*rh2formation_h2mheat[i]);
        
        

        #ifdef SCALE_INPUT
        inv_scale2 = inv_scale[ j + 9];
        scale1     = scale    [ j + 9];
        SM_ELEMENT_D( J, j + 9, j + 9) *= inv_scale2*scale1;
        #endif
        SM_ELEMENT_D( J, j + 9, j + 9 ) *= inv_mdensity;
        SM_ELEMENT_D( J, j + 9, j + 9 ) *= Tge[i];
    }
    return 0;
}



#ifdef CVKLU
int calculate_sparse_jacobian_primordial( realtype t,
                                        N_Vector y, N_Vector fy,
                                        SUNMatrix J, void *user_data,
                                        N_Vector tmp1, N_Vector tmp2,
                                        N_Vector tmp3)
{
    /* We iterate over all of the rates */
    /* Calcuate temperature first */
    

    primordial_data *data = (primordial_data*)user_data; 
    
    int nchem = 10;
    int nstrip = data->nstrip;
    int i, j;
    int NSPARSE = 64;
    
    #ifdef _OPENMP
    int threadID = omp_get_thread_num();
    #else
    int threadID = 0;
    #endif
    
    /* change N_Vector back to an array */
    double y_arr[ 10 * nstrip ];
    double *scale     = data->scale[threadID];
    double *inv_scale = data->inv_scale[threadID];
    double h2_optical_depth_approx;
    double cie_optical_depth_approx;

    //TODO: Here we assumed during the evaluation of jacobian
    // temperature is approximately constant, 
    // i.e. close enough to the point evaluation of f(y)
    // such that the rates and temperature need not be interpolated or evalulated 
    // again during the jacobian evaluation.
    // We havent really fully explored the effect of this `assumption`...
    // But it definitely boost the performance 
    
    /*
    int flag;
    flag = primordial_calculate_temperature(data, y_arr , nstrip, nchem );
    if (flag > 0){
        // check if the temperature failed to converged
        return -1;    
    }
    primordial_interpolate_rates(data, nstrip);
    */

    // primordial_calculate_temperature(data, y_arr, nstrip, nchem);
    // primordial_interpolate_rates(data, nstrip);

    /* Now We set up some temporaries */
    // CSR is what we choose
    sunindextype *rowptrs = SUNSparseMatrix_IndexPointers(J);
    sunindextype *colvals = SUNSparseMatrix_IndexValues(J);
    realtype *matrix_data = SUNSparseMatrix_Data(J);
    
    SUNMatZero(J);
   
    double *Tge = data->dTs_ge[threadID];
    double *k01 = data->rs_k01[threadID];
    double *rk01= data->drs_k01[threadID];
    double *k02 = data->rs_k02[threadID];
    double *rk02= data->drs_k02[threadID];
    double *k03 = data->rs_k03[threadID];
    double *rk03= data->drs_k03[threadID];
    double *k04 = data->rs_k04[threadID];
    double *rk04= data->drs_k04[threadID];
    double *k05 = data->rs_k05[threadID];
    double *rk05= data->drs_k05[threadID];
    double *k06 = data->rs_k06[threadID];
    double *rk06= data->drs_k06[threadID];
    double *k07 = data->rs_k07[threadID];
    double *rk07= data->drs_k07[threadID];
    double *k08 = data->rs_k08[threadID];
    double *rk08= data->drs_k08[threadID];
    double *k09 = data->rs_k09[threadID];
    double *rk09= data->drs_k09[threadID];
    double *k10 = data->rs_k10[threadID];
    double *rk10= data->drs_k10[threadID];
    double *k11 = data->rs_k11[threadID];
    double *rk11= data->drs_k11[threadID];
    double *k12 = data->rs_k12[threadID];
    double *rk12= data->drs_k12[threadID];
    double *k13 = data->rs_k13[threadID];
    double *rk13= data->drs_k13[threadID];
    double *k14 = data->rs_k14[threadID];
    double *rk14= data->drs_k14[threadID];
    double *k15 = data->rs_k15[threadID];
    double *rk15= data->drs_k15[threadID];
    double *k16 = data->rs_k16[threadID];
    double *rk16= data->drs_k16[threadID];
    double *k17 = data->rs_k17[threadID];
    double *rk17= data->drs_k17[threadID];
    double *k18 = data->rs_k18[threadID];
    double *rk18= data->drs_k18[threadID];
    double *k19 = data->rs_k19[threadID];
    double *rk19= data->drs_k19[threadID];
    double *k21 = data->rs_k21[threadID];
    double *rk21= data->drs_k21[threadID];
    double *k22 = data->rs_k22[threadID];
    double *rk22= data->drs_k22[threadID];
    double *k23 = data->rs_k23[threadID];
    double *rk23= data->drs_k23[threadID];
    double *brem_brem = data->cs_brem_brem[threadID];
    double *rbrem_brem = data->dcs_brem_brem[threadID];
    double *ceHeI_ceHeI = data->cs_ceHeI_ceHeI[threadID];
    double *rceHeI_ceHeI = data->dcs_ceHeI_ceHeI[threadID];
    double *ceHeII_ceHeII = data->cs_ceHeII_ceHeII[threadID];
    double *rceHeII_ceHeII = data->dcs_ceHeII_ceHeII[threadID];
    double *ceHI_ceHI = data->cs_ceHI_ceHI[threadID];
    double *rceHI_ceHI = data->dcs_ceHI_ceHI[threadID];
    double *cie_cooling_cieco = data->cs_cie_cooling_cieco[threadID];
    double *rcie_cooling_cieco = data->dcs_cie_cooling_cieco[threadID];
    double *ciHeI_ciHeI = data->cs_ciHeI_ciHeI[threadID];
    double *rciHeI_ciHeI = data->dcs_ciHeI_ciHeI[threadID];
    double *ciHeII_ciHeII = data->cs_ciHeII_ciHeII[threadID];
    double *rciHeII_ciHeII = data->dcs_ciHeII_ciHeII[threadID];
    double *ciHeIS_ciHeIS = data->cs_ciHeIS_ciHeIS[threadID];
    double *rciHeIS_ciHeIS = data->dcs_ciHeIS_ciHeIS[threadID];
    double *ciHI_ciHI = data->cs_ciHI_ciHI[threadID];
    double *rciHI_ciHI = data->dcs_ciHI_ciHI[threadID];
    double *compton_comp_ = data->cs_compton_comp_[threadID];
    double *rcompton_comp_ = data->dcs_compton_comp_[threadID];
    double *gammah_gammah = data->cs_gammah_gammah[threadID];
    double *rgammah_gammah = data->dcs_gammah_gammah[threadID];
    double *gloverabel08_gael = data->cs_gloverabel08_gael[threadID];
    double *rgloverabel08_gael = data->dcs_gloverabel08_gael[threadID];
    double *gloverabel08_gaH2 = data->cs_gloverabel08_gaH2[threadID];
    double *rgloverabel08_gaH2 = data->dcs_gloverabel08_gaH2[threadID];
    double *gloverabel08_gaHe = data->cs_gloverabel08_gaHe[threadID];
    double *rgloverabel08_gaHe = data->dcs_gloverabel08_gaHe[threadID];
    double *gloverabel08_gaHI = data->cs_gloverabel08_gaHI[threadID];
    double *rgloverabel08_gaHI = data->dcs_gloverabel08_gaHI[threadID];
    double *gloverabel08_gaHp = data->cs_gloverabel08_gaHp[threadID];
    double *rgloverabel08_gaHp = data->dcs_gloverabel08_gaHp[threadID];
    double *gloverabel08_gphdl = data->cs_gloverabel08_gphdl[threadID];
    double *rgloverabel08_gphdl = data->dcs_gloverabel08_gphdl[threadID];
    double *gloverabel08_gpldl = data->cs_gloverabel08_gpldl[threadID];
    double *rgloverabel08_gpldl = data->dcs_gloverabel08_gpldl[threadID];
    double *gloverabel08_h2lte = data->cs_gloverabel08_h2lte[threadID];
    double *rgloverabel08_h2lte = data->dcs_gloverabel08_h2lte[threadID];
    double *h2formation_h2mcool = data->cs_h2formation_h2mcool[threadID];
    double *rh2formation_h2mcool = data->dcs_h2formation_h2mcool[threadID];
    double *h2formation_h2mheat = data->cs_h2formation_h2mheat[threadID];
    double *rh2formation_h2mheat = data->dcs_h2formation_h2mheat[threadID];
    double *h2formation_ncrd1 = data->cs_h2formation_ncrd1[threadID];
    double *rh2formation_ncrd1 = data->dcs_h2formation_ncrd1[threadID];
    double *h2formation_ncrd2 = data->cs_h2formation_ncrd2[threadID];
    double *rh2formation_ncrd2 = data->dcs_h2formation_ncrd2[threadID];
    double *h2formation_ncrn = data->cs_h2formation_ncrn[threadID];
    double *rh2formation_ncrn = data->dcs_h2formation_ncrn[threadID];
    double *reHeII1_reHeII1 = data->cs_reHeII1_reHeII1[threadID];
    double *rreHeII1_reHeII1 = data->dcs_reHeII1_reHeII1[threadID];
    double *reHeII2_reHeII2 = data->cs_reHeII2_reHeII2[threadID];
    double *rreHeII2_reHeII2 = data->dcs_reHeII2_reHeII2[threadID];
    double *reHeIII_reHeIII = data->cs_reHeIII_reHeIII[threadID];
    double *rreHeIII_reHeIII = data->dcs_reHeIII_reHeIII[threadID];
    double *reHII_reHII = data->cs_reHII_reHII[threadID];
    double *rreHII_reHII = data->dcs_reHII_reHII[threadID];
    double H2_1;
    double H2_2;
    double H_1;
    double H_2;
    double H_m0;
    double He_1;
    double He_2;
    double He_3;
    double de;
    double ge;
    double z;
    double T;

    double mh = 1.66054e-24;
    double mdensity, inv_mdensity;
    
    double scale2, inv_scale1;

    j = 0;
    mdensity = 0.0;
    z = data->current_z;
   
    int k = 0;
    
    double *yvec_ptr = N_VGetArrayPointer(y);

    for ( i = 0; i < nstrip; i++ ){

        #ifdef SCALE_INPUT
        j = i * nchem;
        H2_1 = yvec_ptr[j]*scale[j];
        j++;
        
        H2_2 = yvec_ptr[j]*scale[j];
        j++;
        
        H_1 = yvec_ptr[j]*scale[j];
        j++;
        
        H_2 = yvec_ptr[j]*scale[j];
        j++;
        
        H_m0 = yvec_ptr[j]*scale[j];
        j++;
        
        He_1 = yvec_ptr[j]*scale[j];
        j++;
        
        He_2 = yvec_ptr[j]*scale[j];
        j++;
        
        He_3 = yvec_ptr[j]*scale[j];
        j++;
        
        de = yvec_ptr[j]*scale[j];
        j++;
        
        ge = yvec_ptr[j]*scale[j];
        j++;
        
        #else
        j = i * nchem;
        H2_1 = yvec_ptr[j];
        j++;
        
        H2_2 = yvec_ptr[j];
        j++;
        
        H_1 = yvec_ptr[j];
        j++;
        
        H_2 = yvec_ptr[j];
        j++;
        
        H_m0 = yvec_ptr[j];
        j++;
        
        He_1 = yvec_ptr[j];
        j++;
        
        He_2 = yvec_ptr[j];
        j++;
        
        He_3 = yvec_ptr[j];
        j++;
        
        de = yvec_ptr[j];
        j++;
        
        ge = yvec_ptr[j];
        j++;
        
        #endif

        mdensity = data->mdensity[threadID][i];
        inv_mdensity = 1.0 / mdensity; 
        
        h2_optical_depth_approx = data->h2_optical_depth_approx[threadID][i];  
        
        
        
        cie_optical_depth_approx = data->cie_optical_depth_approx[threadID][i];
        

        j = i * NSPARSE;
        // H2_1 by H2_1
        colvals[j + 0] = i * nchem + 0 ;
        matrix_data[ j + 0 ] = -k11[i]*H_2 - k12[i]*de - k13[i]*H_1 + k21[i]*pow(H_1, 2) - 2*k23[i]*H2_1;

        
        // H2_1 by H2_2
        colvals[j + 1] = i * nchem + 1 ;
        matrix_data[ j + 1 ] = k10[i]*H_1 + k19[i]*H_m0;

        
        // H2_1 by H_1
        colvals[j + 2] = i * nchem + 2 ;
        matrix_data[ j + 2 ] = k08[i]*H_m0 + k10[i]*H2_2 - k13[i]*H2_1 + 2*k21[i]*H2_1*H_1 + 3*k22[i]*pow(H_1, 2);

        
        // H2_1 by H_2
        colvals[j + 3] = i * nchem + 3 ;
        matrix_data[ j + 3 ] = -k11[i]*H2_1;

        
        // H2_1 by H_m0
        colvals[j + 4] = i * nchem + 4 ;
        matrix_data[ j + 4 ] = k08[i]*H_1 + k19[i]*H2_2;

        
        // H2_1 by de
        colvals[j + 5] = i * nchem + 8 ;
        matrix_data[ j + 5 ] = -k12[i]*H2_1;

        
        // H2_1 by ge
        colvals[j + 6] = i * nchem + 9 ;
        matrix_data[ j + 6 ] = -pow(H2_1, 2)*rk23[i] + H2_1*pow(H_1, 2)*rk21[i] - H2_1*H_1*rk13[i] - H2_1*H_2*rk11[i] - H2_1*de*rk12[i] + H2_2*H_1*rk10[i] + H2_2*H_m0*rk19[i] + pow(H_1, 3)*rk22[i] + H_1*H_m0*rk08[i];

        
        matrix_data[ j + 6] *= Tge[i];
        // H2_2 by H2_1
        colvals[j + 7] = i * nchem + 0 ;
        matrix_data[ j + 7 ] = k11[i]*H_2;

        
        // H2_2 by H2_2
        colvals[j + 8] = i * nchem + 1 ;
        matrix_data[ j + 8 ] = -k10[i]*H_1 - k18[i]*de - k19[i]*H_m0;

        
        // H2_2 by H_1
        colvals[j + 9] = i * nchem + 2 ;
        matrix_data[ j + 9 ] = k09[i]*H_2 - k10[i]*H2_2;

        
        // H2_2 by H_2
        colvals[j + 10] = i * nchem + 3 ;
        matrix_data[ j + 10 ] = k09[i]*H_1 + k11[i]*H2_1 + k17[i]*H_m0;

        
        // H2_2 by H_m0
        colvals[j + 11] = i * nchem + 4 ;
        matrix_data[ j + 11 ] = k17[i]*H_2 - k19[i]*H2_2;

        
        // H2_2 by de
        colvals[j + 12] = i * nchem + 8 ;
        matrix_data[ j + 12 ] = -k18[i]*H2_2;

        
        // H2_2 by ge
        colvals[j + 13] = i * nchem + 9 ;
        matrix_data[ j + 13 ] = H2_1*H_2*rk11[i] - H2_2*H_1*rk10[i] - H2_2*H_m0*rk19[i] - H2_2*de*rk18[i] + H_1*H_2*rk09[i] + H_2*H_m0*rk17[i];

        
        matrix_data[ j + 13] *= Tge[i];
        // H_1 by H2_1
        colvals[j + 14] = i * nchem + 0 ;
        matrix_data[ j + 14 ] = k11[i]*H_2 + 2*k12[i]*de + 2*k13[i]*H_1 - 2*k21[i]*pow(H_1, 2) + 4*k23[i]*H2_1;

        
        // H_1 by H2_2
        colvals[j + 15] = i * nchem + 1 ;
        matrix_data[ j + 15 ] = -k10[i]*H_1 + 2*k18[i]*de + k19[i]*H_m0;

        
        // H_1 by H_1
        colvals[j + 16] = i * nchem + 2 ;
        matrix_data[ j + 16 ] = -k01[i]*de - k07[i]*de - k08[i]*H_m0 - k09[i]*H_2 - k10[i]*H2_2 + 2*k13[i]*H2_1 + k15[i]*H_m0 - 4*k21[i]*H2_1*H_1 - 6*k22[i]*pow(H_1, 2);

        
        // H_1 by H_2
        colvals[j + 17] = i * nchem + 3 ;
        matrix_data[ j + 17 ] = k02[i]*de - k09[i]*H_1 + k11[i]*H2_1 + 2*k16[i]*H_m0;

        
        // H_1 by H_m0
        colvals[j + 18] = i * nchem + 4 ;
        matrix_data[ j + 18 ] = -k08[i]*H_1 + k14[i]*de + k15[i]*H_1 + 2*k16[i]*H_2 + k19[i]*H2_2;

        
        // H_1 by de
        colvals[j + 19] = i * nchem + 8 ;
        matrix_data[ j + 19 ] = -k01[i]*H_1 + k02[i]*H_2 - k07[i]*H_1 + 2*k12[i]*H2_1 + k14[i]*H_m0 + 2*k18[i]*H2_2;

        
        // H_1 by ge
        colvals[j + 20] = i * nchem + 9 ;
        matrix_data[ j + 20 ] = 2*pow(H2_1, 2)*rk23[i] - 2*H2_1*pow(H_1, 2)*rk21[i] + 2*H2_1*H_1*rk13[i] + H2_1*H_2*rk11[i] + 2*H2_1*de*rk12[i] - H2_2*H_1*rk10[i] + H2_2*H_m0*rk19[i] + 2*H2_2*de*rk18[i] - 2*pow(H_1, 3)*rk22[i] - H_1*H_2*rk09[i] - H_1*H_m0*rk08[i] + H_1*H_m0*rk15[i] - H_1*de*rk01[i] - H_1*de*rk07[i] + 2*H_2*H_m0*rk16[i] + H_2*de*rk02[i] + H_m0*de*rk14[i];

        
        matrix_data[ j + 20] *= Tge[i];
        // H_2 by H2_1
        colvals[j + 21] = i * nchem + 0 ;
        matrix_data[ j + 21 ] = -k11[i]*H_2;

        
        // H_2 by H2_2
        colvals[j + 22] = i * nchem + 1 ;
        matrix_data[ j + 22 ] = k10[i]*H_1;

        
        // H_2 by H_1
        colvals[j + 23] = i * nchem + 2 ;
        matrix_data[ j + 23 ] = k01[i]*de - k09[i]*H_2 + k10[i]*H2_2;

        
        // H_2 by H_2
        colvals[j + 24] = i * nchem + 3 ;
        matrix_data[ j + 24 ] = -k02[i]*de - k09[i]*H_1 - k11[i]*H2_1 - k16[i]*H_m0 - k17[i]*H_m0;

        
        // H_2 by H_m0
        colvals[j + 25] = i * nchem + 4 ;
        matrix_data[ j + 25 ] = -k16[i]*H_2 - k17[i]*H_2;

        
        // H_2 by de
        colvals[j + 26] = i * nchem + 8 ;
        matrix_data[ j + 26 ] = k01[i]*H_1 - k02[i]*H_2;

        
        // H_2 by ge
        colvals[j + 27] = i * nchem + 9 ;
        matrix_data[ j + 27 ] = -H2_1*H_2*rk11[i] + H2_2*H_1*rk10[i] - H_1*H_2*rk09[i] + H_1*de*rk01[i] - H_2*H_m0*rk16[i] - H_2*H_m0*rk17[i] - H_2*de*rk02[i];

        
        matrix_data[ j + 27] *= Tge[i];
        // H_m0 by H2_2
        colvals[j + 28] = i * nchem + 1 ;
        matrix_data[ j + 28 ] = -k19[i]*H_m0;

        
        // H_m0 by H_1
        colvals[j + 29] = i * nchem + 2 ;
        matrix_data[ j + 29 ] = k07[i]*de - k08[i]*H_m0 - k15[i]*H_m0;

        
        // H_m0 by H_2
        colvals[j + 30] = i * nchem + 3 ;
        matrix_data[ j + 30 ] = -k16[i]*H_m0 - k17[i]*H_m0;

        
        // H_m0 by H_m0
        colvals[j + 31] = i * nchem + 4 ;
        matrix_data[ j + 31 ] = -k08[i]*H_1 - k14[i]*de - k15[i]*H_1 - k16[i]*H_2 - k17[i]*H_2 - k19[i]*H2_2;

        
        // H_m0 by de
        colvals[j + 32] = i * nchem + 8 ;
        matrix_data[ j + 32 ] = k07[i]*H_1 - k14[i]*H_m0;

        
        // H_m0 by ge
        colvals[j + 33] = i * nchem + 9 ;
        matrix_data[ j + 33 ] = -H2_2*H_m0*rk19[i] - H_1*H_m0*rk08[i] - H_1*H_m0*rk15[i] + H_1*de*rk07[i] - H_2*H_m0*rk16[i] - H_2*H_m0*rk17[i] - H_m0*de*rk14[i];

        
        matrix_data[ j + 33] *= Tge[i];
        // He_1 by He_1
        colvals[j + 34] = i * nchem + 5 ;
        matrix_data[ j + 34 ] = -k03[i]*de;

        
        // He_1 by He_2
        colvals[j + 35] = i * nchem + 6 ;
        matrix_data[ j + 35 ] = k04[i]*de;

        
        // He_1 by de
        colvals[j + 36] = i * nchem + 8 ;
        matrix_data[ j + 36 ] = -k03[i]*He_1 + k04[i]*He_2;

        
        // He_1 by ge
        colvals[j + 37] = i * nchem + 9 ;
        matrix_data[ j + 37 ] = -He_1*de*rk03[i] + He_2*de*rk04[i];

        
        matrix_data[ j + 37] *= Tge[i];
        // He_2 by He_1
        colvals[j + 38] = i * nchem + 5 ;
        matrix_data[ j + 38 ] = k03[i]*de;

        
        // He_2 by He_2
        colvals[j + 39] = i * nchem + 6 ;
        matrix_data[ j + 39 ] = -k04[i]*de - k05[i]*de;

        
        // He_2 by He_3
        colvals[j + 40] = i * nchem + 7 ;
        matrix_data[ j + 40 ] = k06[i]*de;

        
        // He_2 by de
        colvals[j + 41] = i * nchem + 8 ;
        matrix_data[ j + 41 ] = k03[i]*He_1 - k04[i]*He_2 - k05[i]*He_2 + k06[i]*He_3;

        
        // He_2 by ge
        colvals[j + 42] = i * nchem + 9 ;
        matrix_data[ j + 42 ] = He_1*de*rk03[i] - He_2*de*rk04[i] - He_2*de*rk05[i] + He_3*de*rk06[i];

        
        matrix_data[ j + 42] *= Tge[i];
        // He_3 by He_2
        colvals[j + 43] = i * nchem + 6 ;
        matrix_data[ j + 43 ] = k05[i]*de;

        
        // He_3 by He_3
        colvals[j + 44] = i * nchem + 7 ;
        matrix_data[ j + 44 ] = -k06[i]*de;

        
        // He_3 by de
        colvals[j + 45] = i * nchem + 8 ;
        matrix_data[ j + 45 ] = k05[i]*He_2 - k06[i]*He_3;

        
        // He_3 by ge
        colvals[j + 46] = i * nchem + 9 ;
        matrix_data[ j + 46 ] = He_2*de*rk05[i] - He_3*de*rk06[i];

        
        matrix_data[ j + 46] *= Tge[i];
        // de by H2_2
        colvals[j + 47] = i * nchem + 1 ;
        matrix_data[ j + 47 ] = -k18[i]*de;

        
        // de by H_1
        colvals[j + 48] = i * nchem + 2 ;
        matrix_data[ j + 48 ] = k01[i]*de - k07[i]*de + k08[i]*H_m0 + k15[i]*H_m0;

        
        // de by H_2
        colvals[j + 49] = i * nchem + 3 ;
        matrix_data[ j + 49 ] = -k02[i]*de + k17[i]*H_m0;

        
        // de by H_m0
        colvals[j + 50] = i * nchem + 4 ;
        matrix_data[ j + 50 ] = k08[i]*H_1 + k14[i]*de + k15[i]*H_1 + k17[i]*H_2;

        
        // de by He_1
        colvals[j + 51] = i * nchem + 5 ;
        matrix_data[ j + 51 ] = k03[i]*de;

        
        // de by He_2
        colvals[j + 52] = i * nchem + 6 ;
        matrix_data[ j + 52 ] = -k04[i]*de + k05[i]*de;

        
        // de by He_3
        colvals[j + 53] = i * nchem + 7 ;
        matrix_data[ j + 53 ] = -k06[i]*de;

        
        // de by de
        colvals[j + 54] = i * nchem + 8 ;
        matrix_data[ j + 54 ] = k01[i]*H_1 - k02[i]*H_2 + k03[i]*He_1 - k04[i]*He_2 + k05[i]*He_2 - k06[i]*He_3 - k07[i]*H_1 + k14[i]*H_m0 - k18[i]*H2_2;

        
        // de by ge
        colvals[j + 55] = i * nchem + 9 ;
        matrix_data[ j + 55 ] = -H2_2*de*rk18[i] + H_1*H_m0*rk08[i] + H_1*H_m0*rk15[i] + H_1*de*rk01[i] - H_1*de*rk07[i] + H_2*H_m0*rk17[i] - H_2*de*rk02[i] + H_m0*de*rk14[i] + He_1*de*rk03[i] - He_2*de*rk04[i] + He_2*de*rk05[i] - He_3*de*rk06[i];

        
        matrix_data[ j + 55] *= Tge[i];
        // ge by H2_1
        colvals[j + 56] = i * nchem + 0 ;
        matrix_data[ j + 56 ] = -2.0158800000000001*cie_cooling_cieco[i]*mdensity - gloverabel08_gaH2[i]*pow(gloverabel08_h2lte[i], 2)*H2_1*h2_optical_depth_approx/(pow(gloverabel08_h2lte[i]/(gloverabel08_gaH2[i]*H2_1 + gloverabel08_gaHI[i]*H_1 + gloverabel08_gaHe[i]*He_1 + gloverabel08_gaHp[i]*H_2 + gloverabel08_gael[i]*de) + 1.0, 2)*pow(gloverabel08_gaH2[i]*H2_1 + gloverabel08_gaHI[i]*H_1 + gloverabel08_gaHe[i]*He_1 + gloverabel08_gaHp[i]*H_2 + gloverabel08_gael[i]*de, 2)) - gloverabel08_h2lte[i]*h2_optical_depth_approx/(gloverabel08_h2lte[i]/(gloverabel08_gaH2[i]*H2_1 + gloverabel08_gaHI[i]*H_1 + gloverabel08_gaHe[i]*He_1 + gloverabel08_gaHp[i]*H_2 + gloverabel08_gael[i]*de) + 1.0) - 0.5*h2formation_h2mcool[i]*H_1*1.0/(h2formation_ncrn[i]/(h2formation_ncrd1[i]*H_1 + h2formation_ncrd2[i]*H2_1) + 1.0) + 0.5*h2formation_ncrd2[i]*h2formation_ncrn[i]*pow(h2formation_ncrn[i]/(h2formation_ncrd1[i]*H_1 + h2formation_ncrd2[i]*H2_1) + 1.0, -2.0)*(-h2formation_h2mcool[i]*H2_1*H_1 + h2formation_h2mheat[i]*pow(H_1, 3))/pow(h2formation_ncrd1[i]*H_1 + h2formation_ncrd2[i]*H2_1, 2);

        
        matrix_data[j + 56] *= inv_mdensity;
        // ge by H_1
        colvals[j + 57] = i * nchem + 2 ;
        matrix_data[ j + 57 ] = -ceHI_ceHI[i]*de - ciHI_ciHI[i]*de - gloverabel08_gaHI[i]*pow(gloverabel08_h2lte[i], 2)*H2_1*h2_optical_depth_approx/(pow(gloverabel08_h2lte[i]/(gloverabel08_gaH2[i]*H2_1 + gloverabel08_gaHI[i]*H_1 + gloverabel08_gaHe[i]*He_1 + gloverabel08_gaHp[i]*H_2 + gloverabel08_gael[i]*de) + 1.0, 2)*pow(gloverabel08_gaH2[i]*H2_1 + gloverabel08_gaHI[i]*H_1 + gloverabel08_gaHe[i]*He_1 + gloverabel08_gaHp[i]*H_2 + gloverabel08_gael[i]*de, 2)) + 0.5*h2formation_ncrd1[i]*h2formation_ncrn[i]*pow(h2formation_ncrn[i]/(h2formation_ncrd1[i]*H_1 + h2formation_ncrd2[i]*H2_1) + 1.0, -2.0)*(-h2formation_h2mcool[i]*H2_1*H_1 + h2formation_h2mheat[i]*pow(H_1, 3))/pow(h2formation_ncrd1[i]*H_1 + h2formation_ncrd2[i]*H2_1, 2) + 0.5*(-h2formation_h2mcool[i]*H2_1 + 3*h2formation_h2mheat[i]*pow(H_1, 2))*1.0/(h2formation_ncrn[i]/(h2formation_ncrd1[i]*H_1 + h2formation_ncrd2[i]*H2_1) + 1.0);

        
        matrix_data[j + 57] *= inv_mdensity;
        // ge by H_2
        colvals[j + 58] = i * nchem + 3 ;
        matrix_data[ j + 58 ] = -brem_brem[i]*de - gloverabel08_gaHp[i]*pow(gloverabel08_h2lte[i], 2)*H2_1*h2_optical_depth_approx/(pow(gloverabel08_h2lte[i]/(gloverabel08_gaH2[i]*H2_1 + gloverabel08_gaHI[i]*H_1 + gloverabel08_gaHe[i]*He_1 + gloverabel08_gaHp[i]*H_2 + gloverabel08_gael[i]*de) + 1.0, 2)*pow(gloverabel08_gaH2[i]*H2_1 + gloverabel08_gaHI[i]*H_1 + gloverabel08_gaHe[i]*He_1 + gloverabel08_gaHp[i]*H_2 + gloverabel08_gael[i]*de, 2)) - reHII_reHII[i]*de;

        
        matrix_data[j + 58] *= inv_mdensity;
        // ge by He_1
        colvals[j + 59] = i * nchem + 5 ;
        matrix_data[ j + 59 ] = -ciHeI_ciHeI[i]*de - gloverabel08_gaHe[i]*pow(gloverabel08_h2lte[i], 2)*H2_1*h2_optical_depth_approx/(pow(gloverabel08_h2lte[i]/(gloverabel08_gaH2[i]*H2_1 + gloverabel08_gaHI[i]*H_1 + gloverabel08_gaHe[i]*He_1 + gloverabel08_gaHp[i]*H_2 + gloverabel08_gael[i]*de) + 1.0, 2)*pow(gloverabel08_gaH2[i]*H2_1 + gloverabel08_gaHI[i]*H_1 + gloverabel08_gaHe[i]*He_1 + gloverabel08_gaHp[i]*H_2 + gloverabel08_gael[i]*de, 2));

        
        matrix_data[j + 59] *= inv_mdensity;
        // ge by He_2
        colvals[j + 60] = i * nchem + 6 ;
        matrix_data[ j + 60 ] = -brem_brem[i]*de - ceHeII_ceHeII[i]*de - ceHeI_ceHeI[i]*pow(de, 2) - ciHeII_ciHeII[i]*de - ciHeIS_ciHeIS[i]*pow(de, 2) - reHeII1_reHeII1[i]*de - reHeII2_reHeII2[i]*de;

        
        matrix_data[j + 60] *= inv_mdensity;
        // ge by He_3
        colvals[j + 61] = i * nchem + 7 ;
        matrix_data[ j + 61 ] = -4.0*brem_brem[i]*de - reHeIII_reHeIII[i]*de;

        
        matrix_data[j + 61] *= inv_mdensity;
        // ge by de
        colvals[j + 62] = i * nchem + 8 ;
        matrix_data[ j + 62 ] = brem_brem[i]*(-H_2 - He_2 - 4.0*He_3) - ceHI_ceHI[i]*H_1 - ceHeII_ceHeII[i]*He_2 - 2*ceHeI_ceHeI[i]*He_2*de - ciHI_ciHI[i]*H_1 - ciHeII_ciHeII[i]*He_2 - 2*ciHeIS_ciHeIS[i]*He_2*de - ciHeI_ciHeI[i]*He_1 - compton_comp_[i]*pow(z + 1.0, 4)*(T - 2.73*z - 2.73) - gloverabel08_gael[i]*pow(gloverabel08_h2lte[i], 2)*H2_1*h2_optical_depth_approx/(pow(gloverabel08_h2lte[i]/(gloverabel08_gaH2[i]*H2_1 + gloverabel08_gaHI[i]*H_1 + gloverabel08_gaHe[i]*He_1 + gloverabel08_gaHp[i]*H_2 + gloverabel08_gael[i]*de) + 1.0, 2)*pow(gloverabel08_gaH2[i]*H2_1 + gloverabel08_gaHI[i]*H_1 + gloverabel08_gaHe[i]*He_1 + gloverabel08_gaHp[i]*H_2 + gloverabel08_gael[i]*de, 2)) - reHII_reHII[i]*H_2 - reHeII1_reHeII1[i]*He_2 - reHeII2_reHeII2[i]*He_2 - reHeIII_reHeIII[i]*He_3;

        
        matrix_data[j + 62] *= inv_mdensity;
        // ge by ge
        colvals[j + 63] = i * nchem + 9 ;
        matrix_data[ j + 63 ] = -gloverabel08_h2lte[i]*H2_1*h2_optical_depth_approx*(-gloverabel08_h2lte[i]*(-H2_1*rgloverabel08_gaH2[i] - H_1*rgloverabel08_gaHI[i] - H_2*rgloverabel08_gaHp[i] - He_1*rgloverabel08_gaHe[i] - de*rgloverabel08_gael[i])/pow(gloverabel08_gaH2[i]*H2_1 + gloverabel08_gaHI[i]*H_1 + gloverabel08_gaHe[i]*He_1 + gloverabel08_gaHp[i]*H_2 + gloverabel08_gael[i]*de, 2) - rgloverabel08_h2lte[i]/(gloverabel08_gaH2[i]*H2_1 + gloverabel08_gaHI[i]*H_1 + gloverabel08_gaHe[i]*He_1 + gloverabel08_gaHp[i]*H_2 + gloverabel08_gael[i]*de))/pow(gloverabel08_h2lte[i]/(gloverabel08_gaH2[i]*H2_1 + gloverabel08_gaHI[i]*H_1 + gloverabel08_gaHe[i]*He_1 + gloverabel08_gaHp[i]*H_2 + gloverabel08_gael[i]*de) + 1.0, 2) - H2_1*h2_optical_depth_approx*rgloverabel08_h2lte[i]/(gloverabel08_h2lte[i]/(gloverabel08_gaH2[i]*H2_1 + gloverabel08_gaHI[i]*H_1 + gloverabel08_gaHe[i]*He_1 + gloverabel08_gaHp[i]*H_2 + gloverabel08_gael[i]*de) + 1.0) - 2.0158800000000001*H2_1*mdensity*rcie_cooling_cieco[i] - H_1*de*rceHI_ceHI[i] - H_1*de*rciHI_ciHI[i] - H_2*de*rreHII_reHII[i] - He_1*de*rciHeI_ciHeI[i] - He_2*pow(de, 2)*rceHeI_ceHeI[i] - He_2*pow(de, 2)*rciHeIS_ciHeIS[i] - He_2*de*rceHeII_ceHeII[i] - He_2*de*rciHeII_ciHeII[i] - He_2*de*rreHeII1_reHeII1[i] - He_2*de*rreHeII2_reHeII2[i] - He_3*de*rreHeIII_reHeIII[i] - de*rbrem_brem[i]*(H_2 + He_2 + 4.0*He_3) - de*rcompton_comp_[i]*pow(z + 1.0, 4)*(T - 2.73*z - 2.73) + 0.5*pow(h2formation_ncrn[i]/(h2formation_ncrd1[i]*H_1 + h2formation_ncrd2[i]*H2_1) + 1.0, -2.0)*(-h2formation_h2mcool[i]*H2_1*H_1 + h2formation_h2mheat[i]*pow(H_1, 3))*(-1.0*h2formation_ncrn[i]*(-H2_1*rh2formation_ncrd2[i] - H_1*rh2formation_ncrd1[i])/pow(h2formation_ncrd1[i]*H_1 + h2formation_ncrd2[i]*H2_1, 2) - 1.0*rh2formation_ncrn[i]/(h2formation_ncrd1[i]*H_1 + h2formation_ncrd2[i]*H2_1)) + 0.5*1.0/(h2formation_ncrn[i]/(h2formation_ncrd1[i]*H_1 + h2formation_ncrd2[i]*H2_1) + 1.0)*(-H2_1*H_1*rh2formation_h2mcool[i] + pow(H_1, 3)*rh2formation_h2mheat[i]);

        
        matrix_data[j + 63] *= inv_mdensity;
        matrix_data[ j + 63] *= Tge[i];
        rowptrs[ i * nchem +  0] = i * NSPARSE + 0;
        rowptrs[ i * nchem +  1] = i * NSPARSE + 7;
        rowptrs[ i * nchem +  2] = i * NSPARSE + 14;
        rowptrs[ i * nchem +  3] = i * NSPARSE + 21;
        rowptrs[ i * nchem +  4] = i * NSPARSE + 28;
        rowptrs[ i * nchem +  5] = i * NSPARSE + 34;
        rowptrs[ i * nchem +  6] = i * NSPARSE + 38;
        rowptrs[ i * nchem +  7] = i * NSPARSE + 43;
        rowptrs[ i * nchem +  8] = i * NSPARSE + 47;
        rowptrs[ i * nchem +  9] = i * NSPARSE + 56;
       
        #ifdef SCALE_INPUT
        j = i * nchem;
        inv_scale1 = inv_scale[ j + 0 ];
        scale2     = scale    [ j + 0 ];
        matrix_data[ i * NSPARSE + 0]  *= inv_scale1*scale2;
        inv_scale1 = inv_scale[ j + 0 ];
        scale2     = scale    [ j + 1 ];
        matrix_data[ i * NSPARSE + 1]  *= inv_scale1*scale2;
        inv_scale1 = inv_scale[ j + 0 ];
        scale2     = scale    [ j + 2 ];
        matrix_data[ i * NSPARSE + 2]  *= inv_scale1*scale2;
        inv_scale1 = inv_scale[ j + 0 ];
        scale2     = scale    [ j + 3 ];
        matrix_data[ i * NSPARSE + 3]  *= inv_scale1*scale2;
        inv_scale1 = inv_scale[ j + 0 ];
        scale2     = scale    [ j + 4 ];
        matrix_data[ i * NSPARSE + 4]  *= inv_scale1*scale2;
        inv_scale1 = inv_scale[ j + 0 ];
        scale2     = scale    [ j + 8 ];
        matrix_data[ i * NSPARSE + 5]  *= inv_scale1*scale2;
        inv_scale1 = inv_scale[ j + 0 ];
        scale2     = scale    [ j + 9 ];
        matrix_data[ i * NSPARSE + 6]  *= inv_scale1*scale2;
        inv_scale1 = inv_scale[ j + 1 ];
        scale2     = scale    [ j + 0 ];
        matrix_data[ i * NSPARSE + 7]  *= inv_scale1*scale2;
        inv_scale1 = inv_scale[ j + 1 ];
        scale2     = scale    [ j + 1 ];
        matrix_data[ i * NSPARSE + 8]  *= inv_scale1*scale2;
        inv_scale1 = inv_scale[ j + 1 ];
        scale2     = scale    [ j + 2 ];
        matrix_data[ i * NSPARSE + 9]  *= inv_scale1*scale2;
        inv_scale1 = inv_scale[ j + 1 ];
        scale2     = scale    [ j + 3 ];
        matrix_data[ i * NSPARSE + 10]  *= inv_scale1*scale2;
        inv_scale1 = inv_scale[ j + 1 ];
        scale2     = scale    [ j + 4 ];
        matrix_data[ i * NSPARSE + 11]  *= inv_scale1*scale2;
        inv_scale1 = inv_scale[ j + 1 ];
        scale2     = scale    [ j + 8 ];
        matrix_data[ i * NSPARSE + 12]  *= inv_scale1*scale2;
        inv_scale1 = inv_scale[ j + 1 ];
        scale2     = scale    [ j + 9 ];
        matrix_data[ i * NSPARSE + 13]  *= inv_scale1*scale2;
        inv_scale1 = inv_scale[ j + 2 ];
        scale2     = scale    [ j + 0 ];
        matrix_data[ i * NSPARSE + 14]  *= inv_scale1*scale2;
        inv_scale1 = inv_scale[ j + 2 ];
        scale2     = scale    [ j + 1 ];
        matrix_data[ i * NSPARSE + 15]  *= inv_scale1*scale2;
        inv_scale1 = inv_scale[ j + 2 ];
        scale2     = scale    [ j + 2 ];
        matrix_data[ i * NSPARSE + 16]  *= inv_scale1*scale2;
        inv_scale1 = inv_scale[ j + 2 ];
        scale2     = scale    [ j + 3 ];
        matrix_data[ i * NSPARSE + 17]  *= inv_scale1*scale2;
        inv_scale1 = inv_scale[ j + 2 ];
        scale2     = scale    [ j + 4 ];
        matrix_data[ i * NSPARSE + 18]  *= inv_scale1*scale2;
        inv_scale1 = inv_scale[ j + 2 ];
        scale2     = scale    [ j + 8 ];
        matrix_data[ i * NSPARSE + 19]  *= inv_scale1*scale2;
        inv_scale1 = inv_scale[ j + 2 ];
        scale2     = scale    [ j + 9 ];
        matrix_data[ i * NSPARSE + 20]  *= inv_scale1*scale2;
        inv_scale1 = inv_scale[ j + 3 ];
        scale2     = scale    [ j + 0 ];
        matrix_data[ i * NSPARSE + 21]  *= inv_scale1*scale2;
        inv_scale1 = inv_scale[ j + 3 ];
        scale2     = scale    [ j + 1 ];
        matrix_data[ i * NSPARSE + 22]  *= inv_scale1*scale2;
        inv_scale1 = inv_scale[ j + 3 ];
        scale2     = scale    [ j + 2 ];
        matrix_data[ i * NSPARSE + 23]  *= inv_scale1*scale2;
        inv_scale1 = inv_scale[ j + 3 ];
        scale2     = scale    [ j + 3 ];
        matrix_data[ i * NSPARSE + 24]  *= inv_scale1*scale2;
        inv_scale1 = inv_scale[ j + 3 ];
        scale2     = scale    [ j + 4 ];
        matrix_data[ i * NSPARSE + 25]  *= inv_scale1*scale2;
        inv_scale1 = inv_scale[ j + 3 ];
        scale2     = scale    [ j + 8 ];
        matrix_data[ i * NSPARSE + 26]  *= inv_scale1*scale2;
        inv_scale1 = inv_scale[ j + 3 ];
        scale2     = scale    [ j + 9 ];
        matrix_data[ i * NSPARSE + 27]  *= inv_scale1*scale2;
        inv_scale1 = inv_scale[ j + 4 ];
        scale2     = scale    [ j + 1 ];
        matrix_data[ i * NSPARSE + 28]  *= inv_scale1*scale2;
        inv_scale1 = inv_scale[ j + 4 ];
        scale2     = scale    [ j + 2 ];
        matrix_data[ i * NSPARSE + 29]  *= inv_scale1*scale2;
        inv_scale1 = inv_scale[ j + 4 ];
        scale2     = scale    [ j + 3 ];
        matrix_data[ i * NSPARSE + 30]  *= inv_scale1*scale2;
        inv_scale1 = inv_scale[ j + 4 ];
        scale2     = scale    [ j + 4 ];
        matrix_data[ i * NSPARSE + 31]  *= inv_scale1*scale2;
        inv_scale1 = inv_scale[ j + 4 ];
        scale2     = scale    [ j + 8 ];
        matrix_data[ i * NSPARSE + 32]  *= inv_scale1*scale2;
        inv_scale1 = inv_scale[ j + 4 ];
        scale2     = scale    [ j + 9 ];
        matrix_data[ i * NSPARSE + 33]  *= inv_scale1*scale2;
        inv_scale1 = inv_scale[ j + 5 ];
        scale2     = scale    [ j + 5 ];
        matrix_data[ i * NSPARSE + 34]  *= inv_scale1*scale2;
        inv_scale1 = inv_scale[ j + 5 ];
        scale2     = scale    [ j + 6 ];
        matrix_data[ i * NSPARSE + 35]  *= inv_scale1*scale2;
        inv_scale1 = inv_scale[ j + 5 ];
        scale2     = scale    [ j + 8 ];
        matrix_data[ i * NSPARSE + 36]  *= inv_scale1*scale2;
        inv_scale1 = inv_scale[ j + 5 ];
        scale2     = scale    [ j + 9 ];
        matrix_data[ i * NSPARSE + 37]  *= inv_scale1*scale2;
        inv_scale1 = inv_scale[ j + 6 ];
        scale2     = scale    [ j + 5 ];
        matrix_data[ i * NSPARSE + 38]  *= inv_scale1*scale2;
        inv_scale1 = inv_scale[ j + 6 ];
        scale2     = scale    [ j + 6 ];
        matrix_data[ i * NSPARSE + 39]  *= inv_scale1*scale2;
        inv_scale1 = inv_scale[ j + 6 ];
        scale2     = scale    [ j + 7 ];
        matrix_data[ i * NSPARSE + 40]  *= inv_scale1*scale2;
        inv_scale1 = inv_scale[ j + 6 ];
        scale2     = scale    [ j + 8 ];
        matrix_data[ i * NSPARSE + 41]  *= inv_scale1*scale2;
        inv_scale1 = inv_scale[ j + 6 ];
        scale2     = scale    [ j + 9 ];
        matrix_data[ i * NSPARSE + 42]  *= inv_scale1*scale2;
        inv_scale1 = inv_scale[ j + 7 ];
        scale2     = scale    [ j + 6 ];
        matrix_data[ i * NSPARSE + 43]  *= inv_scale1*scale2;
        inv_scale1 = inv_scale[ j + 7 ];
        scale2     = scale    [ j + 7 ];
        matrix_data[ i * NSPARSE + 44]  *= inv_scale1*scale2;
        inv_scale1 = inv_scale[ j + 7 ];
        scale2     = scale    [ j + 8 ];
        matrix_data[ i * NSPARSE + 45]  *= inv_scale1*scale2;
        inv_scale1 = inv_scale[ j + 7 ];
        scale2     = scale    [ j + 9 ];
        matrix_data[ i * NSPARSE + 46]  *= inv_scale1*scale2;
        inv_scale1 = inv_scale[ j + 8 ];
        scale2     = scale    [ j + 1 ];
        matrix_data[ i * NSPARSE + 47]  *= inv_scale1*scale2;
        inv_scale1 = inv_scale[ j + 8 ];
        scale2     = scale    [ j + 2 ];
        matrix_data[ i * NSPARSE + 48]  *= inv_scale1*scale2;
        inv_scale1 = inv_scale[ j + 8 ];
        scale2     = scale    [ j + 3 ];
        matrix_data[ i * NSPARSE + 49]  *= inv_scale1*scale2;
        inv_scale1 = inv_scale[ j + 8 ];
        scale2     = scale    [ j + 4 ];
        matrix_data[ i * NSPARSE + 50]  *= inv_scale1*scale2;
        inv_scale1 = inv_scale[ j + 8 ];
        scale2     = scale    [ j + 5 ];
        matrix_data[ i * NSPARSE + 51]  *= inv_scale1*scale2;
        inv_scale1 = inv_scale[ j + 8 ];
        scale2     = scale    [ j + 6 ];
        matrix_data[ i * NSPARSE + 52]  *= inv_scale1*scale2;
        inv_scale1 = inv_scale[ j + 8 ];
        scale2     = scale    [ j + 7 ];
        matrix_data[ i * NSPARSE + 53]  *= inv_scale1*scale2;
        inv_scale1 = inv_scale[ j + 8 ];
        scale2     = scale    [ j + 8 ];
        matrix_data[ i * NSPARSE + 54]  *= inv_scale1*scale2;
        inv_scale1 = inv_scale[ j + 8 ];
        scale2     = scale    [ j + 9 ];
        matrix_data[ i * NSPARSE + 55]  *= inv_scale1*scale2;
        inv_scale1 = inv_scale[ j + 9 ];
        scale2     = scale    [ j + 0 ];
        matrix_data[ i * NSPARSE + 56]  *= inv_scale1*scale2;
        inv_scale1 = inv_scale[ j + 9 ];
        scale2     = scale    [ j + 2 ];
        matrix_data[ i * NSPARSE + 57]  *= inv_scale1*scale2;
        inv_scale1 = inv_scale[ j + 9 ];
        scale2     = scale    [ j + 3 ];
        matrix_data[ i * NSPARSE + 58]  *= inv_scale1*scale2;
        inv_scale1 = inv_scale[ j + 9 ];
        scale2     = scale    [ j + 5 ];
        matrix_data[ i * NSPARSE + 59]  *= inv_scale1*scale2;
        inv_scale1 = inv_scale[ j + 9 ];
        scale2     = scale    [ j + 6 ];
        matrix_data[ i * NSPARSE + 60]  *= inv_scale1*scale2;
        inv_scale1 = inv_scale[ j + 9 ];
        scale2     = scale    [ j + 7 ];
        matrix_data[ i * NSPARSE + 61]  *= inv_scale1*scale2;
        inv_scale1 = inv_scale[ j + 9 ];
        scale2     = scale    [ j + 8 ];
        matrix_data[ i * NSPARSE + 62]  *= inv_scale1*scale2;
        inv_scale1 = inv_scale[ j + 9 ];
        scale2     = scale    [ j + 9 ];
        matrix_data[ i * NSPARSE + 63]  *= inv_scale1*scale2;
        #endif

    }

    rowptrs[ i * nchem ] = i * NSPARSE ;
    return 0;
}

#endif




void setting_up_extra_variables( primordial_data * data, double * input, int nstrip ){
    //-------------------------------------------------------------------------    
    // Function: setting_up_extra_variables
    // Desciption: calculating variables that are independent on the state with time
    //             to avoid repeated evaluation. Examples here are h2 optical depth
    //             and cie_optical depth. Functions that depends only on density
    //-------------------------------------------------------------------------    
    #ifdef _OPENMP
    int threadID = omp_get_thread_num();
    #else
    int threadID = 0;
    #endif

    int i, j;
    double mh = 1.66054e-24;
    double mdensity;
    // TODO: maybe this should be filled out by Dengo as well
    for ( i = 0; i < nstrip; i++){
        data->mdensity[threadID][i] = 0;
        j = i * 10;
        // species: H2_1
        data->mdensity[threadID][i] += input[j] * 2.0;
        j++;
        // species: H2_2
        data->mdensity[threadID][i] += input[j] * 2.0;
        j++;
        // species: H_1
        data->mdensity[threadID][i] += input[j] * 1.00794;
        j++;
        // species: H_2
        data->mdensity[threadID][i] += input[j] * 1.00794;
        j++;
        // species: H_m0
        data->mdensity[threadID][i] += input[j] * 1.00794;
        j++;
        // species: He_1
        data->mdensity[threadID][i] += input[j] * 4.002602;
        j++;
        // species: He_2
        data->mdensity[threadID][i] += input[j] * 4.002602;
        j++;
        // species: He_3
        data->mdensity[threadID][i] += input[j] * 4.002602;
        j++;
        j++;
        j++;
        // TODO: update temperature and rates to obtain abundances of equilibrium states
	// for more detail: go to the _calculate_temperature 
        data->mdensity[threadID][i] *= mh;
        data->inv_mdensity[threadID][i] = 1.0 / data->mdensity[threadID][i];
    }
    double tau;
    for ( i = 0; i < nstrip; i++){
        
        mdensity = data->mdensity[threadID][i];
        tau      = pow( (mdensity / 3.3e-8 ), 2.8);
        tau      = fmax( tau, 1.0e-5 );
        data->cie_optical_depth_approx[threadID][i] = fmin( 1.0, (1.0 - exp(-tau) ) / tau );
    }
    for ( i = 0; i < nstrip; i++ ){
        mdensity = data->mdensity[threadID][i];
        data->h2_optical_depth_approx[threadID][i] = fmin( 1.0, pow( (mdensity / (1.34e-14) )  , -0.45) );
    }
}



///////////////////////////////////////////////////////////////////////////////
/////////////////// Sturcture Incoming Data ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// Reshape/ Flatten data to match the data shape/ type 
// required by the solver in dengo


int flatten_dengo_field_data(code_units *units, dengo_field_data *field_data, double *input){

    //-----------------------------------------------------
    // Function     :   flatten_dengo_field_data 
    // Parameter    :   
    //                  code_units: units from the incoming field_data 
    //                  field_data: dengo_field_data class that contains pointer to species array 
    //                  input     : 1D array that flattens the field_data, 
    //                              i.e. 
    //                              s = species, s0 = 0th species, with d dimensions
    //                              [s0, s1, s2..., , s0_1, s1_1, ... s0_d, s1_d, ...] 
    //-----------------------------------------------------
    //
    unsigned long d, dims, i, j;
    dims = field_data->ncells;

    // code unit in terms of erg/g
    double UNIT_E_per_M = units->velocity_units * units->velocity_units;
    double m_amu = 1.66053904e-24;
    double dom = units->density_units/m_amu;
    int N = 10;

    #pragma omp parallel for private (i, j ,d) num_threads(NTHREADS) schedule(static,1)
    for ( d = 0; d< dims; d++  ){
        j = d*N;
        // this should be the normalized 
        // by the input units later
        // atol = input * rtol;
        // which again should be set by dengo
        // input in *mass density per amu* 
        // and energy in the units of (erg / g)
        input[j]  = field_data->H2_1_density[d] ;
        input[j] *= dom;
        j++;
        input[j]  = field_data->H2_2_density[d] ;
        input[j] *= dom;
        j++;
        input[j]  = field_data->H_1_density[d] ;
        input[j] *= dom;
        j++;
        input[j]  = field_data->H_2_density[d] ;
        input[j] *= dom;
        j++;
        input[j]  = field_data->H_m0_density[d] ;
        input[j] *= dom;
        j++;
        input[j]  = field_data->He_1_density[d] ;
        input[j] *= dom;
        j++;
        input[j]  = field_data->He_2_density[d] ;
        input[j] *= dom;
        j++;
        input[j]  = field_data->He_3_density[d] ;
        input[j] *= dom;
        j++;
        input[j]  = field_data->de_density[d] ;
        input[j] *= dom;
        j++;
        input[j]  = field_data->ge_density[d] ;
        input[j] *= UNIT_E_per_M;
        j++;
    }
}



int reshape_to_dengo_field_data( code_units* units, dengo_field_data *field_data, double* input ){
    //------------------------------------------------------------------------------------
    // Function   :     reshape_to_dengo_field_data
    // Description:     reshape the 1d output array from solver to a dengo_field_data object  
    //                  and covert them to code units 
    //                  i.e. ge_density in erg /g
    //                       H_1_density in g / cm^-3 / amu (mass density per amu)
    // Parameter  :     code_units
    //                  dengo_field_data
    //                  input
    //------------------------------------------------------------------------------------

    unsigned long int i, j, d, dims;
    int N = 10;
    dims = field_data->ncells; // total number of strips to be evaluated
    

    // code unit in terms of erg/g
    double UNIT_E_per_M = units->velocity_units * units->velocity_units;
    double m_amu = 1.66053904e-24;
    double dom = units->density_units/ m_amu;
   
    #pragma omp parallel for private (i, j ,d) num_threads(NTHREADS) schedule(static, 1)
    for ( d = 0; d< dims; d++  ){
        j = d*N;
        field_data->H2_1_density[d] = input[j];
        field_data->H2_1_density[d] /= dom;
        j++;
        field_data->H2_2_density[d] = input[j];
        field_data->H2_2_density[d] /= dom;
        j++;
        field_data->H_1_density[d] = input[j];
        field_data->H_1_density[d] /= dom;
        j++;
        field_data->H_2_density[d] = input[j];
        field_data->H_2_density[d] /= dom;
        j++;
        field_data->H_m0_density[d] = input[j];
        field_data->H_m0_density[d] /= dom;
        j++;
        field_data->He_1_density[d] = input[j];
        field_data->He_1_density[d] /= dom;
        j++;
        field_data->He_2_density[d] = input[j];
        field_data->He_2_density[d] /= dom;
        j++;
        field_data->He_3_density[d] = input[j];
        field_data->He_3_density[d] /= dom;
        j++;
        field_data->de_density[d] = input[j];
        field_data->de_density[d] /= dom;
        j++;
        field_data->ge_density[d] = input[j];
        field_data->ge_density[d] /= UNIT_E_per_M;
        j++;
    }
   
    return 0;
}



// and a enzo version
//

int flatten_dengo_field_data_enzo(code_units *units, dengo_field_data *field_data, double *input){

    //-----------------------------------------------------
    // Function     :   flatten_dengo_field_data_enzo
    // Description  :   To read in data from Enzo pointers 
    // Parameter    :   
    //                  code_units: units from the incoming field_data 
    //                  field_data: dengo_field_data class that contains pointer to species array 
    //                  input     : 1D array that flattens the field_data, 
    //                              i.e. 
    //                              s = species, s0 = 0th species, with d dimensions
    //                              [s0, s1, s2..., , s0_1, s1_1, ... s0_d, s1_d, ...] 
    //                              abundances in units of mass density / m_amu
    //                              m_amu is in atomic mass units                           
    //-----------------------------------------------------
    //
    
    int is, ie, js, je, ks, ke;
    int i, j, k, N;
    int ni, nj, nk, idim, jdim, kdim;
    unsigned long dims, ccount, c, idx;

    N = 10;
    // avoid ghost zones
    is = field_data->grid_start[0];
    ie = field_data->grid_end[0];

    js = field_data->grid_start[1];
    je = field_data->grid_end[1];

    ks = field_data->grid_start[2];
    ke = field_data->grid_end[2];

    idim = field_data->grid_dimension[0];
    jdim = field_data->grid_dimension[1];
    kdim = field_data->grid_dimension[2];

    // number of cells that actually required calculations
    ni = ie - is + 1;
    nj = je - js + 1;
    nk = ke - ks + 1;
    dims = ni*nj*nk;
    field_data->ncells = dims;

    // code unit in terms of erg/g
    double UNIT_E_per_M = units->velocity_units * units->velocity_units;
    double m_amu = 1.66053904e-24;
    double dom = units->density_units/m_amu;


    ccount = 0;
    for (k = ks; k <= ke; k++){
    for (j = js; j <= je; j++){
    for (i = is; i <= ie; i++){
        c = ccount * N;
	idx = ((k* jdim + j)*idim + i);
        input[c]  = field_data->H2_1_density[idx];
        input[c] *= dom;
    c++;
        input[c]  = field_data->H2_2_density[idx];
        input[c] *= dom;
    c++;
        input[c]  = field_data->H_1_density[idx];
        input[c] *= dom;
    c++;
        input[c]  = field_data->H_2_density[idx];
        input[c] *= dom;
    c++;
        input[c]  = field_data->H_m0_density[idx];
        input[c] *= dom;
    c++;
        input[c]  = field_data->He_1_density[idx];
        input[c] *= dom;
    c++;
        input[c]  = field_data->He_2_density[idx];
        input[c] *= dom;
    c++;
        input[c]  = field_data->He_3_density[idx];
        input[c] *= dom;
    c++;
        input[c]  = field_data->de_density[idx];
        input[c] *= dom;
    c++;
        input[c]  = field_data->ge_density[idx];
        input[c] *= UNIT_E_per_M;
    c++;
    ccount += 1;

    }}}
}



int reshape_to_dengo_field_data_enzo( code_units* units, dengo_field_data *field_data, double* input, double *temp ){
    //------------------------------------------------------------------------------------
    // Function   :     reshape_to_dengo_field_data
    // Description:     reshape the 1d output array from solver to a dengo_field_data object  
    //                  and covert them to code units 
    //                  i.e. ge_density in erg /g
    //                       H_1_density in g / cm^-3 / amu (mass density per amu)
    // Parameter  :     code_units
    //                  dengo_field_data
    //                  input
    //------------------------------------------------------------------------------------

    unsigned long i, j, k, d; 
    unsigned long idx, ccount, c,dims;
    int is, ie, js, je, ks, ke;
    int ni, nj,nk;
    int N = 10;

    is = field_data->grid_start[0];
    ie = field_data->grid_end[0];

    js = field_data->grid_start[1];
    je = field_data->grid_end[1];

    ks = field_data->grid_start[2];
    ke = field_data->grid_end[2];

    int idim = field_data->grid_dimension[0];
    int jdim = field_data->grid_dimension[1];
    int kdim = field_data->grid_dimension[2];

    ni = ie - is + 1;
    nj = je - js + 1;
    nk = ke - ks + 1;
    dims = ni*nj*nk;
    field_data->ncells = dims;
    dims = field_data->ncells; // total number of strips to be evaluated

    // code unit in terms of erg/g
    double UNIT_E_per_M = units->velocity_units * units->velocity_units;
    double m_amu = 1.66053904e-24;
    double dom = units->density_units/ m_amu;
   
    // extra bits for calculating  conservation
    ////////////////////////////////////////////
    ccount = 0;
    for (k = ks; k <= ke; k++){
    for (j = js; j <= je; j++){
    for (i = is; i <= ie; i++){
        c = ccount * N;
	idx = ((k* jdim + j)*idim + i);
        field_data->H2_1_density[idx] = input[c];
        field_data->H2_1_density[idx] /= dom;
	c++;
        field_data->H2_2_density[idx] = input[c];
        field_data->H2_2_density[idx] /= dom;
	c++;
        field_data->H_1_density[idx] = input[c];
        field_data->H_1_density[idx] /= dom;
	c++;
        field_data->H_2_density[idx] = input[c];
        field_data->H_2_density[idx] /= dom;
	c++;
        field_data->H_m0_density[idx] = input[c];
        field_data->H_m0_density[idx] /= dom;
	c++;
        field_data->He_1_density[idx] = input[c];
        field_data->He_1_density[idx] /= dom;
	c++;
        field_data->He_2_density[idx] = input[c];
        field_data->He_2_density[idx] /= dom;
	c++;
        field_data->He_3_density[idx] = input[c];
        field_data->He_3_density[idx] /= dom;
	c++;
        field_data->de_density[idx] = input[c];
        field_data->de_density[idx] /= dom;
	c++;
        field_data->ge_density[idx] = input[c];
        field_data->ge_density[idx] /= UNIT_E_per_M;
	c++;
	ccount += 1;
    }}}

    return 0;
}




int read_init_data_to_dengo( dengo_field_data *field_data, char const *filename){

    // this reads the initial abundances of the data from
    // a hdf5 file, and initialize a field_data object

    hid_t file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);

    if (file_id < 0){
        fprintf(stderr, "failed to open %s so dying. \n", filename);
        return (1);
        }

    hsize_t dims;
    /* Check gas energy to get the number of cells */
    fprintf(stderr, "Getting dimensionality from ge:\n");
    herr_t status = H5LTget_dataset_info(file_id, "/ge", &dims, NULL, NULL);
    if(status == -1) {
        fprintf(stderr, "Error opening initial conditions file.\n");
        return 1;
    }
    fprintf(stderr, "  ncells = % 3i\n", (int) dims);
    
    field_data->ncells = (int) dims;
    int N = 10;
    double *atol, *rtol;
    atol = (double *) malloc(N * dims * sizeof(double));
    rtol = (double *) malloc(N * dims * sizeof(double));

    double *tics = (double *) malloc(dims * sizeof(double));
    double *ics = (double *) malloc(dims * N * sizeof(double));
    double *input = (double *) malloc(dims * N * sizeof(double));
    
    unsigned int i = 0, j;
    double *H2_1 = (double *) malloc(dims * sizeof(double));    
    field_data->H2_1_density = H2_1;
    
    double *H2_2 = (double *) malloc(dims * sizeof(double));    
    field_data->H2_2_density = H2_2;
    
    double *H_1 = (double *) malloc(dims * sizeof(double));    
    field_data->H_1_density = H_1;
    
    double *H_2 = (double *) malloc(dims * sizeof(double));    
    field_data->H_2_density = H_2;
    
    double *H_m0 = (double *) malloc(dims * sizeof(double));    
    field_data->H_m0_density = H_m0;
    
    double *He_1 = (double *) malloc(dims * sizeof(double));    
    field_data->He_1_density = He_1;
    
    double *He_2 = (double *) malloc(dims * sizeof(double));    
    field_data->He_2_density = He_2;
    
    double *He_3 = (double *) malloc(dims * sizeof(double));    
    field_data->He_3_density = He_3;
    
    double *de = (double *) malloc(dims * sizeof(double));    
    field_data->de_density = de;
    
    double *ge = (double *) malloc(dims * sizeof(double));    
    field_data->ge_density = ge;
    
    fprintf(stderr, "Reading I.C. for /H2_1\n");
    H5LTread_dataset_double(file_id, "/H2_1", tics);
    for ( j = 0; j < dims; j++ ){
        field_data->H2_1_density[j] = tics[j];
        if (j == 0){
            fprintf(stderr, "H2_1[0] = %0.3g \n", tics[j] );
        }
    }
    i++;
    
    fprintf(stderr, "Reading I.C. for /H2_2\n");
    H5LTread_dataset_double(file_id, "/H2_2", tics);
    for ( j = 0; j < dims; j++ ){
        field_data->H2_2_density[j] = tics[j];
        if (j == 0){
            fprintf(stderr, "H2_2[0] = %0.3g \n", tics[j] );
        }
    }
    i++;
    
    fprintf(stderr, "Reading I.C. for /H_1\n");
    H5LTread_dataset_double(file_id, "/H_1", tics);
    for ( j = 0; j < dims; j++ ){
        field_data->H_1_density[j] = tics[j];
        if (j == 0){
            fprintf(stderr, "H_1[0] = %0.3g \n", tics[j] );
        }
    }
    i++;
    
    fprintf(stderr, "Reading I.C. for /H_2\n");
    H5LTread_dataset_double(file_id, "/H_2", tics);
    for ( j = 0; j < dims; j++ ){
        field_data->H_2_density[j] = tics[j];
        if (j == 0){
            fprintf(stderr, "H_2[0] = %0.3g \n", tics[j] );
        }
    }
    i++;
    
    fprintf(stderr, "Reading I.C. for /H_m0\n");
    H5LTread_dataset_double(file_id, "/H_m0", tics);
    for ( j = 0; j < dims; j++ ){
        field_data->H_m0_density[j] = tics[j];
        if (j == 0){
            fprintf(stderr, "H_m0[0] = %0.3g \n", tics[j] );
        }
    }
    i++;
    
    fprintf(stderr, "Reading I.C. for /He_1\n");
    H5LTread_dataset_double(file_id, "/He_1", tics);
    for ( j = 0; j < dims; j++ ){
        field_data->He_1_density[j] = tics[j];
        if (j == 0){
            fprintf(stderr, "He_1[0] = %0.3g \n", tics[j] );
        }
    }
    i++;
    
    fprintf(stderr, "Reading I.C. for /He_2\n");
    H5LTread_dataset_double(file_id, "/He_2", tics);
    for ( j = 0; j < dims; j++ ){
        field_data->He_2_density[j] = tics[j];
        if (j == 0){
            fprintf(stderr, "He_2[0] = %0.3g \n", tics[j] );
        }
    }
    i++;
    
    fprintf(stderr, "Reading I.C. for /He_3\n");
    H5LTread_dataset_double(file_id, "/He_3", tics);
    for ( j = 0; j < dims; j++ ){
        field_data->He_3_density[j] = tics[j];
        if (j == 0){
            fprintf(stderr, "He_3[0] = %0.3g \n", tics[j] );
        }
    }
    i++;
    
    fprintf(stderr, "Reading I.C. for /de\n");
    H5LTread_dataset_double(file_id, "/de", tics);
    for ( j = 0; j < dims; j++ ){
        field_data->de_density[j] = tics[j];
        if (j == 0){
            fprintf(stderr, "de[0] = %0.3g \n", tics[j] );
        }
    }
    i++;
    
    fprintf(stderr, "Reading I.C. for /ge\n");
    H5LTread_dataset_double(file_id, "/ge", tics);
    for ( j = 0; j < dims; j++ ){
        field_data->ge_density[j] = tics[j];
        if (j == 0){
            fprintf(stderr, "ge[0] = %0.3g \n", tics[j] );
        }
    }
    i++;
    
    
    H5Fclose(file_id);
    free(input);
    free(ics);    
    return 1;
}

