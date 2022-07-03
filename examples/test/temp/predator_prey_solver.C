
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


#include "predator_prey_solver.h"


///////////////////////////////////////////////////////////////////////////////
/////////// Setup the reaction, cooling rate data table ///////////////////////
///////////////////////////////////////////////////////////////////////////////
predator_prey_data *predator_prey_setup_data( const char *FileLocation, int *NumberOfFields, char ***FieldNames)
{

    //-----------------------------------------------------
    // Function : predator_prey_setup_data
    // Description: Initialize a data object that stores the reaction/ cooling rate data 
    //-----------------------------------------------------

    int i, n;
    
    predator_prey_data *data = (predator_prey_data *) malloc(sizeof(predator_prey_data));
    
    // point the module to look for predator_prey_tables.h5
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
    
    predator_prey_read_rate_tables(data);
    //fprintf(stderr, "Successfully read in rate tables.\n");

    predator_prey_read_cooling_tables(data);
    //fprintf(stderr, "Successfully read in cooling rate tables.\n");
    
    predator_prey_read_gamma(data);
    //fprintf(stderr, "Successfully read in gamma tables. \n");

    if (FieldNames != NULL && NumberOfFields != NULL) {
        NumberOfFields[0] = 5;
        FieldNames[0] = new char*[5];
        i = 0;
        FieldNames[0][i++] = strdup("dead_predator");
        
        FieldNames[0][i++] = strdup("dead_prey");
        
        FieldNames[0][i++] = strdup("ge");
        
        FieldNames[0][i++] = strdup("predator");
        
        FieldNames[0][i++] = strdup("prey");
        
    }

    data->dengo_data_file = NULL;

    return data;

}


void predator_prey_read_rate_tables(predator_prey_data *data)
{
    const char * filedir;
    if (data->dengo_data_file != NULL){
        filedir =  data->dengo_data_file; 
    } else{
        filedir = "predator_prey_tables.h5";   
    }

    hid_t file_id = H5Fopen( filedir , H5F_ACC_RDONLY, H5P_DEFAULT);
    /* Allocate the correct number of rate tables */
    H5LTread_dataset_double(file_id, "/exp_growth_prey", data->r_exp_growth_prey);
    H5LTread_dataset_double(file_id, "/natural_death_predator", data->r_natural_death_predator);
    H5LTread_dataset_double(file_id, "/predation", data->r_predation);
    
    H5Fclose(file_id);
}


void predator_prey_read_cooling_tables(predator_prey_data *data)
{

    const char * filedir;
    if (data->dengo_data_file != NULL){
        filedir =  data->dengo_data_file; 
    } else{
        filedir = "predator_prey_tables.h5";   
    }
    hid_t file_id = H5Fopen( filedir , H5F_ACC_RDONLY, H5P_DEFAULT);
    /* Allocate the correct number of rate tables */

    H5Fclose(file_id);
}

void predator_prey_read_gamma(predator_prey_data *data)
{

    const char * filedir;
    if (data->dengo_data_file != NULL){
        filedir =  data->dengo_data_file; 
    } else{
        filedir = "predator_prey_tables.h5";   
    }
    
    hid_t file_id = H5Fopen( filedir , H5F_ACC_RDONLY, H5P_DEFAULT);
    /* Allocate the correct number of rate tables */

    H5Fclose(file_id);

}
 


/*
   This setup may be different than the user may anticipate, as a result
   of the lockstep timestep we use for a pencil beam through the grid.
   As such, it accepts the number of things to interpolate and makes
   assumptions about the sizes of the rates.
*/

/* This also requires no templating other than for the solver name...*/
void predator_prey_interpolate_rates(predator_prey_data *data,
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
        data->rs_exp_growth_prey[threadID][i] = data->r_exp_growth_prey[bin_id] +
            Tdef * (data->r_exp_growth_prey[bin_id+1] - data->r_exp_growth_prey[bin_id]);
        data->drs_exp_growth_prey[threadID][i] = (data->r_exp_growth_prey[bin_id+1] - data->r_exp_growth_prey[bin_id]);
        data->drs_exp_growth_prey[threadID][i] /= data->dT[threadID][i];
        data->drs_exp_growth_prey[threadID][i] *= data->invTs[threadID][i];
    }
    
    
    
    for ( i = 0; i < nstrip; i++ ){
        bin_id = data->bin_id[threadID][i];
        Tdef   = data->Tdef[threadID][i];
        data->rs_natural_death_predator[threadID][i] = data->r_natural_death_predator[bin_id] +
            Tdef * (data->r_natural_death_predator[bin_id+1] - data->r_natural_death_predator[bin_id]);
        data->drs_natural_death_predator[threadID][i] = (data->r_natural_death_predator[bin_id+1] - data->r_natural_death_predator[bin_id]);
        data->drs_natural_death_predator[threadID][i] /= data->dT[threadID][i];
        data->drs_natural_death_predator[threadID][i] *= data->invTs[threadID][i];
    }
    
    
    
    for ( i = 0; i < nstrip; i++ ){
        bin_id = data->bin_id[threadID][i];
        Tdef   = data->Tdef[threadID][i];
        data->rs_predation[threadID][i] = data->r_predation[bin_id] +
            Tdef * (data->r_predation[bin_id+1] - data->r_predation[bin_id]);
        data->drs_predation[threadID][i] = (data->r_predation[bin_id+1] - data->r_predation[bin_id]);
        data->drs_predation[threadID][i] /= data->dT[threadID][i];
        data->drs_predation[threadID][i] *= data->invTs[threadID][i];
    }
    
    
    
}
 


void predator_prey_interpolate_gamma(predator_prey_data *data,
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

    
       
    }



/////// MAIN

int predator_prey_main(int argc, char** argv)
{
    predator_prey_data *data = predator_prey_setup_data("predator_prey_tables.h5", NULL, NULL);

    /* Initial conditions */

    hid_t file_id = H5Fopen("predator_prey_initial_conditions.h5", H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file_id < 0) {fprintf(stderr, "Failed to open "
        "predator_prey_initial_conditions.h5 so dying.\n");
        return(1);}

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

    int N = 5;

    double *atol, *rtol;
    atol = (double *) alloca(N * dims * sizeof(double));
    rtol = (double *) alloca(N * dims * sizeof(double));

    double *tics = (double *) alloca(dims * sizeof(double));
    double *ics = (double *) alloca(dims * N * sizeof(double));
    double *input = (double *) alloca(dims * N * sizeof(double));
    
    unsigned int i = 0, j;
    
    fprintf(stderr, "Reading I.C. for /dead_predator\n");
    H5LTread_dataset_double(file_id, "/dead_predator", tics);
    for (j = 0; j < dims; j++) {
        ics[j * N + i] = tics[j]; 
        atol[j * N + i] = tics[j] * 1e-09;
        rtol[j * N + i] = 1e-09;
        if(j==0) {
            fprintf(stderr, "dead_predator[0] = %0.3g, atol => % 0.16g\n",
                    tics[j], atol[j]);
        }
    }
    i++;
    
    fprintf(stderr, "Reading I.C. for /dead_prey\n");
    H5LTread_dataset_double(file_id, "/dead_prey", tics);
    for (j = 0; j < dims; j++) {
        ics[j * N + i] = tics[j]; 
        atol[j * N + i] = tics[j] * 1e-09;
        rtol[j * N + i] = 1e-09;
        if(j==0) {
            fprintf(stderr, "dead_prey[0] = %0.3g, atol => % 0.16g\n",
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
    
    fprintf(stderr, "Reading I.C. for /predator\n");
    H5LTread_dataset_double(file_id, "/predator", tics);
    for (j = 0; j < dims; j++) {
        ics[j * N + i] = tics[j]; 
        atol[j * N + i] = tics[j] * 1e-09;
        rtol[j * N + i] = 1e-09;
        if(j==0) {
            fprintf(stderr, "predator[0] = %0.3g, atol => % 0.16g\n",
                    tics[j], atol[j]);
        }
    }
    i++;
    
    fprintf(stderr, "Reading I.C. for /prey\n");
    H5LTread_dataset_double(file_id, "/prey", tics);
    for (j = 0; j < dims; j++) {
        ics[j * N + i] = tics[j]; 
        atol[j * N + i] = tics[j] * 1e-09;
        rtol[j * N + i] = 1e-09;
        if(j==0) {
            fprintf(stderr, "prey[0] = %0.3g, atol => % 0.16g\n",
                    tics[j], atol[j]);
        }
    }
    i++;
    

    H5Fclose(file_id);

    double dtf = 31557000000000.0;
    double dt = -1.0;
    double z = -1.0;
    for (i = 0; i < dims * N; i++) input[i] = ics[i];
    double ttot;
    ttot = dengo_evolve_predator_prey(dtf, dt, z, input, rtol, atol, (long long) dims, data);

    /* Write results to HDF5 file */
    file_id = H5Fcreate("predator_prey_solution.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    hsize_t dimsarr[1];
    dimsarr[0] = dims;
    i = 0;
    
    double dead_predator[dims];
    for (j = 0; j < dims; j++) {
        dead_predator[j] = input[j * N + i]; 
    }
    fprintf(stderr, "Writing solution for /dead_predator\n");
    H5LTmake_dataset_double(file_id, "/dead_predator", 1, dimsarr, dead_predator);
    i++;
    
    double dead_prey[dims];
    for (j = 0; j < dims; j++) {
        dead_prey[j] = input[j * N + i]; 
    }
    fprintf(stderr, "Writing solution for /dead_prey\n");
    H5LTmake_dataset_double(file_id, "/dead_prey", 1, dimsarr, dead_prey);
    i++;
    
    double ge[dims];
    for (j = 0; j < dims; j++) {
        ge[j] = input[j * N + i]; 
    }
    fprintf(stderr, "Writing solution for /ge\n");
    H5LTmake_dataset_double(file_id, "/ge", 1, dimsarr, ge);
    i++;
    
    double predator[dims];
    for (j = 0; j < dims; j++) {
        predator[j] = input[j * N + i]; 
    }
    fprintf(stderr, "Writing solution for /predator\n");
    H5LTmake_dataset_double(file_id, "/predator", 1, dimsarr, predator);
    i++;
    
    double prey[dims];
    for (j = 0; j < dims; j++) {
        prey[j] = input[j * N + i]; 
    }
    fprintf(stderr, "Writing solution for /prey\n");
    H5LTmake_dataset_double(file_id, "/prey", 1, dimsarr, prey);
    i++;
    
    //TODO: Ts is linked to threadID
    double temperature[dims];
    for (j = 0; j < dims; j++) {
    	temperature[j] = data->Ts[0][j];
    }
    H5LTmake_dataset_double(file_id, "/T", 1, dimsarr, temperature);
    double time[1];
    time[0] = ttot;
    double timestep[1];
    timestep[0] = dt;
    H5LTset_attribute_double(file_id, "/", "time", time, 1); 
    H5LTset_attribute_double(file_id, "/", "timestep", timestep, 1);
    H5Fclose(file_id);
    
    return 0;
}
 



int dengo_evolve_predator_prey (double dtf, double dt, double z, double *input,
            double *rtol, double *atol, long long dims, predator_prey_data *data) {
    int i, j;
    hid_t file_id;
    /* fprintf(stderr, "  ncells = % 3i\n", (int) dims); */

    int N = 5;
    for (i = 0; i<dims; i++) {
      j = i * N;
      
      input[j] /= 1.0; // dead_predator;
      atol[j] /= 1.0; // dead_predator;
      j++;
      
      input[j] /= 1.0; // dead_prey;
      atol[j] /= 1.0; // dead_prey;
      j++;
      
      j++;
      
      input[j] /= 1.0; // predator;
      atol[j] /= 1.0; // predator;
      j++;
      
      input[j] /= 1.0; // prey;
      atol[j] /= 1.0; // prey;
      j++;
    }
    ensure_electron_consistency(input, dims, N);
    double floor_value = data->floor_value;
    for (j = 0; j< dims*N; j++) input[j] = fmax(input[j], floor_value);

    rhs_f f = calculate_rhs_predator_prey;
    jac_f jf = calculate_jacobian_predator_prey;

    unsigned long MAX_ITERATION = 1e4;
    if (dt < 0) dt = dtf / MAX_ITERATION;
    (z < 0) ? data->current_z = 0: data->current_z = z;
    int niter = 0;
    int siter = 0;
    double ttot = 0;
    double *scale = (double *) malloc(dims * N * sizeof(double));
    double *prev = (double *) malloc(dims * N * sizeof(double));
    for (i = 0; i < dims * N; i++) scale[i] = input[i];
    for (i = 0; i < dims * N; i++) prev[i] = input[i];
    double *u0 = (double *) malloc(N*dims*sizeof(double));
    double *s  = (double *) malloc(N*sizeof(double));
    double *gu = (double *) malloc(N*dims*sizeof(double));
    double *Ju = (double *) malloc(N*N*dims*sizeof(double));

    // setting up extra temporary variables
    setting_up_extra_variables(data, scale, dims);

    while (ttot < dtf) {
        int rv = BE_chem_solve(f, jf, input, dt, rtol, atol, dims, N,
                               scale, (void *) data, u0, s, gu, Ju);
        #ifdef DENGO_DEBUG
        fprintf(stderr, "Return value [%d]: %i.  %0.5g / %0.5g = %0.5g (%0.5g)\n",
                niter, rv, ttot, dtf, ttot/dtf, dt);
        fprintf(stderr, "Value[0] = %0.5g %0.5g\n",
                input[0], prev[0]);
        #endif
        for (i = 0; i < dims * N; i++) {
            if (input[i] < 0) {
                rv = 1;
                break;
            }
        }
        if (rv == 0) {
	    if (siter == MAX_ITERATION) break;
	    siter++;
            if (siter % 10000 == 0) {
                fprintf(stderr, "Successful Iteration[%d]: (%0.4g) %0.16g / %0.16g\n",
                        siter, dt, ttot, dtf);
            }
            ttot += dt;
	    dt = DMIN(dt * 1.1, dtf - ttot);
	    
	    for (i = 0; i < dims * N; i++) prev[i] = input[i];
            for (i = 0; i < dims * N; i++) {     
                if (input[i] < floor_value) {
                  input[i] = floor_value;
                }
            }
        } else {
            dt /= 2.0;
            for (i = 0; i < dims * N; i++) input[i] = prev[i];
            if (dt/dtf < 1e-15)  {
                fprintf(stderr, "Dying! dt/dtf = %0.5g\n", dt/dtf);
                break;
            }
        }
        niter++;
    }
    #ifdef DENGO_DEBUG
    fprintf(stderr, "End: %0.5g / %0.5g (%0.5g)\n",
       ttot, dtf, dtf-ttot);
    #endif
    for (i = 0; i<dims; i++) {
      j = i * N;
      input[j] *= 1.0; // dead_predator;
      atol[j] *= 1.0; // dead_predator
      j++;
      input[j] *= 1.0; // dead_prey;
      atol[j] *= 1.0; // dead_prey
      j++;
      j++;
      input[j] *= 1.0; // predator;
      atol[j] *= 1.0; // predator
      j++;
      input[j] *= 1.0; // prey;
      atol[j] *= 1.0; // prey
      j++;
    }
    free(scale);
    free(prev);
    free(u0);
    free(s);
    free(gu);
    free(Ju);

    if (ttot < dtf) return 1;
    return 0;
}
 

///////////////////////////////////////////////////////////////////////////////
//////////// Evaluate Temperature /////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////


int predator_prey_calculate_temperature(predator_prey_data *data,
                        double *input, int nstrip, int nchem)
{
    int i, j;
    double density, T, Tnew;
    double kb = 1.3806504e-16; // Boltzmann constant [erg/K]
    double mh = 1.66054e-24;
    double gamma = 5.e0/3.e0;
    double _gamma_m1 = 1.0 / (gamma - 1);

    double dge_dT;

    /* Calculate total density */
    double dead_predator;
    double dead_prey;
    double ge;
    double predator;
    double prey;
    
    i = 0;
    #ifdef _OPENMP
    int threadID = omp_get_thread_num();
    #else 
    int threadID = 0;
    #endif

    for ( i = 0; i < nstrip; i++ ){
        j = i * nchem;
        dead_predator = input[j];
        j++;
        dead_prey = input[j];
        j++;
        ge = input[j];
        j++;
        predator = input[j];
        j++;
        prey = input[j];
        j++;
    
        /*
        */
	
        // TODO: pull the rates from predator_prey_data
        // these species usually contribute negligbly to the number density (?)
	// it is a little tricky here,
	// since these species depends on the temperature 
	// and the abundance of the rest of the species
	// BUT, without their abundance, temperature CANNOT be evaluated....
	// FOR NOW, a not entirely correct physically, 
	// BUT a not-too-bad surrogate is:
	// assume these species has negligible abundance....

        density = 1.0*dead_predator + 1.0*dead_prey + 1.0*predator + 1.0*prey;
        
        // Requires iteration on the convergence of temperature
        // since the gammaH2 is not fixed 
        data->Ts[threadID][i] = density*ge*mh/(kb*(dead_predator/(gamma - 1.0) + dead_prey/(gamma - 1.0) + predator/(gamma - 1.0) + prey/(gamma - 1.0)));
        

        if (data->Ts[threadID][i] < data->bounds[0]) {
            data->Ts[threadID][i] = data->bounds[0];
        } else if (data->Ts[threadID][i] > data->bounds[1]) {
            data->Ts[threadID][i] = data->bounds[1];
        }
        data->logTs[threadID][i] = log(data->Ts[threadID][i]);
        data->invTs[threadID][i] = 1.0 / data->Ts[threadID][i];

        dge_dT = kb*(dead_predator/(gamma - 1.0) + dead_prey/(gamma - 1.0) + predator/(gamma - 1.0) + prey/(gamma - 1.0))/(density*mh);
        data->dTs_ge[threadID][i] = 1.0 / dge_dT;
    } // for i in nstrip loop
    return 0;
         
}
 




int calculate_rhs_predator_prey(double *input, double *rhs, int nstrip,
                  int nchem, void *sdata)
{
    /* We iterate over all of the rates */
    /* Calculate temperature first */
    predator_prey_data *data = (predator_prey_data*)sdata;
    int i, j;
    predator_prey_calculate_temperature(data, input, nstrip, nchem);

    predator_prey_interpolate_rates(data, nstrip);
    
    #ifdef _OPENMP
    int threadID = omp_get_thread_num();
    #else 
    int threadID = 0;
    #endif

    /* Now we set up some temporaries */
    double *exp_growth_prey = data->rs_exp_growth_prey[threadID];
    double *natural_death_predator = data->rs_natural_death_predator[threadID];
    double *predation = data->rs_predation[threadID];
    double dead_predator;
    double dead_prey;
    double ge;
    double predator;
    double prey;
    double z;
    double T;


    double mh = 1.67e-24;
    double mdensity, inv_mdensity;
    for (i = 0; i<nstrip; i++) {
        j = i * nchem;
        mdensity = 0.0;
        T = data->Ts[threadID][i];
        z = data->current_z;
        dead_predator = input[j];
        
        mdensity += dead_predator;
        
        if (dead_predator < 0.0) {
            /* fprintf(stderr, "RNegative[%d][dead_predator] = % 0.16g [%d]\n",
               i, dead_predator, j); */
            return 1;
          dead_predator = 1e-20;
        }
        j++;
    	
        dead_prey = input[j];
        
        mdensity += dead_prey;
        
        if (dead_prey < 0.0) {
            /* fprintf(stderr, "RNegative[%d][dead_prey] = % 0.16g [%d]\n",
               i, dead_prey, j); */
            return 1;
          dead_prey = 1e-20;
        }
        j++;
    	
        ge = input[j];
        
        if (ge < 0.0) {
            /* fprintf(stderr, "RNegative[%d][ge] = % 0.16g [%d]\n",
               i, ge, j); */
            return 1;
          ge = 1e-20;
        }
        j++;
    	
        predator = input[j];
        
        mdensity += predator;
        
        if (predator < 0.0) {
            /* fprintf(stderr, "RNegative[%d][predator] = % 0.16g [%d]\n",
               i, predator, j); */
            return 1;
          predator = 1e-20;
        }
        j++;
    	
        prey = input[j];
        
        mdensity += prey;
        
        if (prey < 0.0) {
            /* fprintf(stderr, "RNegative[%d][prey] = % 0.16g [%d]\n",
               i, prey, j); */
            return 1;
          prey = 1e-20;
        }
        j++;
    	


/*
*/
        mdensity *= mh;
        j = i * nchem;
        // 
        // Species: dead_predator
        // 
        rhs[j] = natural_death_predator[i]*predator;
        
        j++;
    
        // 
        // Species: dead_prey
        // 
        rhs[j] = predation[i]*predator*prey;
        
        j++;
    
        // 
        // Species: ge
        // 
        rhs[j] = 0;
        
	rhs[j] /= mdensity;
        
        j++;
    
        // 
        // Species: predator
        // 
        rhs[j] = -natural_death_predator[i]*predator + 0.75*predation[i]*predator*prey;
        
        j++;
    
        // 
        // Species: prey
        // 
        rhs[j] = exp_growth_prey[i]*prey - predation[i]*predator*prey;
        
        j++;
    
    }  

    for (i = 0; i < nchem*nstrip; i++){
        if (rhs[i] != rhs[i]){
            printf("FAILED rhs[%d] = %0.5g; T = %0.5g; z = %0.5g\n", i, rhs[i], T, z);
            return 1;
        }
    }

    return 0;
}




int calculate_jacobian_predator_prey(double *input, double *Joutput,
        int nstrip, int nchem, void *sdata)
{
    /* We iterate over all of the rates */
    /* Calculate temperature first */
    predator_prey_data *data = (predator_prey_data*)sdata;
    
    #ifdef _OPENMP
    int threadID = omp_get_thread_num();
    #else 
    int threadID = 0;
    #endif

    int i, j;
    /* Now we set up some temporaries */
    double *Tge = data->dTs_ge[threadID];
    double *exp_growth_prey = data->rs_exp_growth_prey[threadID];
    double *rexp_growth_prey = data->drs_exp_growth_prey[threadID];
    double *natural_death_predator = data->rs_natural_death_predator[threadID];
    double *rnatural_death_predator = data->drs_natural_death_predator[threadID];
    double *predation = data->rs_predation[threadID];
    double *rpredation = data->drs_predation[threadID];
    double dead_predator;
    double dead_prey;
    double ge;
    double predator;
    double prey;
    double z;
    double T;

    double mh = 1.67e-24;
    double mdensity, inv_mdensity;
    
    

    for (i = 0; i<nstrip; i++) {
        j = i * nchem;
        mdensity = 0.0;
        T = data->Ts[threadID][i];
        z = data->current_z;
	dead_predator = input[j];
        
        mdensity += dead_predator*1.0;
	



        if (dead_predator < 0.0) {
            fprintf(stderr, "JNegative[%d][dead_predator] = % 0.16g [%d]\n",
                    i, dead_predator, j);
            /*dead_predator = 0.0;*/
            dead_predator = 1e-20;
            return 1;
        }
        j++;
        
	dead_prey = input[j];
        
        mdensity += dead_prey*1.0;
	



        if (dead_prey < 0.0) {
            fprintf(stderr, "JNegative[%d][dead_prey] = % 0.16g [%d]\n",
                    i, dead_prey, j);
            /*dead_prey = 0.0;*/
            dead_prey = 1e-20;
            return 1;
        }
        j++;
        
	ge = input[j];
        



        if (ge < 0.0) {
            fprintf(stderr, "JNegative[%d][ge] = % 0.16g [%d]\n",
                    i, ge, j);
            /*ge = 0.0;*/
            ge = 1e-20;
            return 1;
        }
        j++;
        
	predator = input[j];
        
        mdensity += predator*1.0;
	



        if (predator < 0.0) {
            fprintf(stderr, "JNegative[%d][predator] = % 0.16g [%d]\n",
                    i, predator, j);
            /*predator = 0.0;*/
            predator = 1e-20;
            return 1;
        }
        j++;
        
	prey = input[j];
        
        mdensity += prey*1.0;
	



        if (prey < 0.0) {
            fprintf(stderr, "JNegative[%d][prey] = % 0.16g [%d]\n",
                    i, prey, j);
            /*prey = 0.0;*/
            prey = 1e-20;
            return 1;
        }
        j++;
        
        mdensity *= mh;

        j = i * nchem * nchem;
        // 
        // Species: dead_predator
        //
            // dead_predator by dead_predator
            Joutput[j] = 0;
	    
	    
            j++;
            // dead_prey by dead_predator
            Joutput[j] = 0;
	    
	    
            j++;
            // ge by dead_predator
            Joutput[j] = 0;
	    
	    Joutput[j] /= mdensity;
	    
	    
            j++;
            // predator by dead_predator
            Joutput[j] = 0;
	    
	    
            j++;
            // prey by dead_predator
            Joutput[j] = 0;
	    
	    
            j++;
    
        // 
        // Species: dead_prey
        //
            // dead_predator by dead_prey
            Joutput[j] = 0;
	    
	    
            j++;
            // dead_prey by dead_prey
            Joutput[j] = 0;
	    
	    
            j++;
            // ge by dead_prey
            Joutput[j] = 0;
	    
	    Joutput[j] /= mdensity;
	    
	    
            j++;
            // predator by dead_prey
            Joutput[j] = 0;
	    
	    
            j++;
            // prey by dead_prey
            Joutput[j] = 0;
	    
	    
            j++;
    
        // 
        // Species: ge
        //
            // dead_predator by ge
            Joutput[j] = predator*rnatural_death_predator[i];
	    
	    
            Joutput[j] *= Tge[i];
            
            j++;
            // dead_prey by ge
            Joutput[j] = predator*prey*rpredation[i];
	    
	    
            Joutput[j] *= Tge[i];
            
            j++;
            // ge by ge
            Joutput[j] = 0;
	    
	    Joutput[j] /= mdensity;
	    
	    
            Joutput[j] *= Tge[i];
            
            j++;
            // predator by ge
            Joutput[j] = 0.75*predator*prey*rpredation[i] - predator*rnatural_death_predator[i];
	    
	    
            Joutput[j] *= Tge[i];
            
            j++;
            // prey by ge
            Joutput[j] = -predator*prey*rpredation[i] + prey*rexp_growth_prey[i];
	    
	    
            Joutput[j] *= Tge[i];
            
            j++;
    
        // 
        // Species: predator
        //
            // dead_predator by predator
            Joutput[j] = natural_death_predator[i];
	    
	    
            j++;
            // dead_prey by predator
            Joutput[j] = predation[i]*prey;
	    
	    
            j++;
            // ge by predator
            Joutput[j] = 0;
	    
	    Joutput[j] /= mdensity;
	    
	    
            j++;
            // predator by predator
            Joutput[j] = -natural_death_predator[i] + 0.75*predation[i]*prey;
	    
	    
            j++;
            // prey by predator
            Joutput[j] = -predation[i]*prey;
	    
	    
            j++;
    
        // 
        // Species: prey
        //
            // dead_predator by prey
            Joutput[j] = 0;
	    
	    
            j++;
            // dead_prey by prey
            Joutput[j] = predation[i]*predator;
	    
	    
            j++;
            // ge by prey
            Joutput[j] = 0;
	    
	    Joutput[j] /= mdensity;
	    
	    
            j++;
            // predator by prey
            Joutput[j] = 0.75*predation[i]*predator;
	    
	    
            j++;
            // prey by prey
            Joutput[j] = exp_growth_prey[i] - predation[i]*predator;
	    
	    
            j++;
    
    }

    return 0;
    
}




void ensure_electron_consistency(double *input, long long nstrip, int nchem)
{
    int i, j;

    /* Now we set up some temporaries */
    double dead_predator;
    double dead_prey;
    double ge;
    double predator;
    double prey;
    double total_e = 0.0;
    int e_indx;

    for (i = 0; i<nstrip; i++) {
        total_e = 0.0;
        j = i * nchem;
        dead_predator = input[j];
        
        total_e += dead_predator * 0.0;
        j++;
        
        dead_prey = input[j];
        
        total_e += dead_prey * 0.0;
        j++;
        
        ge = input[j];
        
        j++;
        
        predator = input[j];
        
        total_e += predator * 0.0;
        j++;
        
        prey = input[j];
        
        total_e += prey * 0.0;
        j++;
        
        input[e_indx] = total_e;
    }  
}




void setting_up_extra_variables( predator_prey_data * data, double * input, int nstrip ){
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
        j = i * 5;
        // species: dead_predator
        data->mdensity[threadID][i] += input[j] * 1.0;
        j++;
        // species: dead_prey
        data->mdensity[threadID][i] += input[j] * 1.0;
        j++;
        j++;
        // species: predator
        data->mdensity[threadID][i] += input[j] * 1.0;
        j++;
        // species: prey
        data->mdensity[threadID][i] += input[j] * 1.0;
        j++;
        // TODO: update temperature and rates to obtain abundances of equilibrium states
	// for more detail: go to the _calculate_temperature 
        data->mdensity[threadID][i] *= mh;
        data->inv_mdensity[threadID][i] = 1.0 / data->mdensity[threadID][i];
    }
}





// Enzo Routine

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
    
    unsigned long is, ie, js, je, ks, ke;
    unsigned long i, j, k, N;
    unsigned long ni, nj, nk, idim, jdim, kdim;
    unsigned long dims, ccount, c, idx;

    N = 5;
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
        input[c]  = field_data->dead_predator_density[idx];
        input[c] *= dom;
    c++;
        input[c]  = field_data->dead_prey_density[idx];
        input[c] *= dom;
    c++;
        input[c]  = field_data->ge_density[idx];
        input[c] *= UNIT_E_per_M;
    c++;
        input[c]  = field_data->predator_density[idx];
        input[c] *= dom;
    c++;
        input[c]  = field_data->prey_density[idx];
        input[c] *= dom;
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
    unsigned long is, ie, js, je, ks, ke;
    unsigned long ni, nj,nk;
    int N = 5;

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
        field_data->dead_predator_density[idx] = input[c];
        field_data->dead_predator_density[idx] /= dom;
	c++;
        field_data->dead_prey_density[idx] = input[c];
        field_data->dead_prey_density[idx] /= dom;
	c++;
        field_data->ge_density[idx] = input[c];
        field_data->ge_density[idx] /= UNIT_E_per_M;
	c++;
        field_data->predator_density[idx] = input[c];
        field_data->predator_density[idx] /= dom;
	c++;
        field_data->prey_density[idx] = input[c];
        field_data->prey_density[idx] /= dom;
	c++;
	ccount += 1;
    }}}

    return 0;
}



int predator_prey_solve_chemistry_enzo( code_units *units, dengo_field_data *field_data, double dt ){
    //-------------------------------------------------------------------------
    // Function: predator_prey_solve_chemistry_enzo
    // Description: takes the same input as predator_prey_solve_chemistry
    //              BUT field_data needs to be specified with things like
    //              grid_start, grid_end, grid_dimension
    //              such that dengo can avoid evaluating ghost zones
    //-------------------------------------------------------------------------

    // to use this, reltol and floor_value must be specified through 
    // dengo_field_data 
    unsigned long int i, j, k, d, dims;
    int N = 5;
    
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


    double *input = (double *) malloc(dims * N * sizeof(double));
    double *temp  = (double *) malloc(dims * (0+1) * sizeof(double) );
    double *atol  = (double *) malloc(dims * N * sizeof(double));
    double *rtol  = (double *) malloc(dims * N * sizeof(double));
    
    
    // fills the `input` array with field data
    // note that input are still some kind of mass density, rho/ m_amu  
    flatten_dengo_field_data_enzo(units, field_data, input);
    
    for (int j = 0; j< N*dims; j++) atol[j] = field_data->reltol*input[j];
    for (int j = 0; j< N*dims; j++) rtol[j] = field_data->reltol;

    const char * FileLocation = field_data->dengo_data_file;

    // chop up dims into chunks of "MAX_NCELLS"
    // and omp parallel it
    unsigned long ntimes = dims / MAX_NCELLS;
    unsigned long nstrip_res = dims % MAX_NCELLS;

    // reading in rate data
    predator_prey_data *data = predator_prey_setup_data( field_data->dengo_data_file, NULL, NULL);
    // specify the relative tolerance and floor value
    data->reltol       = field_data->reltol;
    double floor_value = field_data->floor_value;

    // specifiy the redshift
    int flag;
    double a = units->a_value * units->a_units;
    double z = 0.0; // 1./a - 1.;
    double dtf;

    // convert code time to seconds 
    dt *= units->time_units;
    dtf = dt;

    // update the rate table location
    data->dengo_data_file = field_data->dengo_data_file; 
    data->floor_value = floor_value;

    int flag_array[ntimes+1];
    #pragma omp parallel for private(dt, flag)
    for (i = 0; i < ntimes; i++)
    {
        #ifdef _OPENMP
        int threadID = omp_get_thread_num();
        #else
        int threadID = 0;
        #endif
        
        dt = dtf;
        // evolve chemistry in dengo
        flag = dengo_evolve_predator_prey(dtf, dt, z, &input[i*N*MAX_NCELLS], rtol, &atol[i*N*MAX_NCELLS], MAX_NCELLS, data);
        #ifdef DENGO_DEBUG
        printf("tID(%d) flag = %d; %lu/ %lu\n", threadID, flag, i, ntimes);
        #endif
        flag_array[i] = flag;
    }

    if (nstrip_res > 0){
        int threadID = 0;

        dt = dtf;
        // evolve chemistry in dengo
        flag = dengo_evolve_predator_prey(dtf, dt, z, &input[ntimes*N*MAX_NCELLS], rtol, &atol[ntimes*N*MAX_NCELLS], nstrip_res, data);
        #ifdef DENGO_DEBUG
        printf("tID(%d) flag = %d; %lu/ %lu\n", threadID, flag, i, nstrip_res);
        #endif
        flag_array[ntimes] = flag;
    } else{
        flag_array[ntimes] = 0;
    }

    // fill results in `input` back to field data
    // in appropriate units
    int f = 0;
    for (i = 0; i < ntimes; i++)
        f = DMAX(f, flag_array[i]);

    if (f == 0)
        reshape_to_dengo_field_data_enzo(units, field_data, input, temp);
        
    // free all pointers
    free(input);
    free(temp);
    free(data);
    free(atol);
    free(rtol);

    // in Enzo; 0 == FAIL
    //          1 == SUCCESS
    int Enzo_Fail    = 0;
    int Enzo_Success = 1;
    if (f == 0) return Enzo_Success;
    return Enzo_Fail;
}



int predator_prey_calculate_cooling_timescale( double *cooling_time, double *input, int nstrip, predator_prey_data *data){
    
    #ifdef _OPENMP
    int threadID = omp_get_thread_num(); 
    #else
    int threadID = 0;
    #endif

    unsigned long i, j, dims;
    int flag;
    int nchem = 5;
    /* Now we set up some temporaries */
    // make sure the input are in number density
    for (i = 0; i < nstrip; i++) {
        j = i * nchem;
        input[j] /= 1.0; // dead_predator
        j++;
        input[j] /= 1.0; // dead_prey
        j++;
        j++;
        input[j] /= 1.0; // predator
        j++;
        input[j] /= 1.0; // prey
        j++;
    }
    
    // calculate temperature and cooling rate
    setting_up_extra_variables(data, input, nstrip );
    flag = predator_prey_calculate_temperature(data,  input , nstrip, nchem );
    if (flag > 0){
        // check if the temperature failed to converged
        return -1;    
    }
    predator_prey_interpolate_rates(data, nstrip);
    double predator;
    double dead_prey;
    double ge;
    double prey;
    double dead_predator;

    // this might be redundant
    // should only select relavant rates
    double *exp_growth_prey = data->rs_exp_growth_prey[threadID];
    double *natural_death_predator = data->rs_natural_death_predator[threadID];
    double *predation = data->rs_predation[threadID];
    
    
    
    
    double z;
    double T;
    double mdensity, inv_mdensity, dge_dt;

    for ( i = 0; i < nstrip; i++ ){
        
        T            = data->Ts[threadID][i];
        z            = data->current_z;
        mdensity     = data->mdensity[threadID][i];
        inv_mdensity = data->inv_mdensity[threadID][i];
        
        

        j = i * nchem;
        dead_predator = input[j];
        j++;
        dead_prey = input[j];
        j++;
        ge = input[j];
        j++;
        predator = input[j];
        j++;
        prey = input[j];
        j++;
   	
        // obtain a quasis-equilibrium estimate

        //
        // Species: ge
        //
        dge_dt = 0;
        dge_dt *= inv_mdensity;
        cooling_time[i] = fabs( ge / dge_dt);
    
    //fprintf(stderr, "----------------\n");
    }
}



int dengo_estimate_cooling_time_enzo( code_units* units, dengo_field_data *field_data ){
    

    unsigned long int i, j, k, d, dims;
    int nchem = 5;
    
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
    predator_prey_data *data = predator_prey_setup_data( FileLocation, NULL, NULL);

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
        predator_prey_calculate_cooling_timescale( cooling_time_batch, input_batch, nstrip, data);
    }

    if (nstrip_res > 0){
        input_batch        = &input[ntimes* nstrip * nchem];
        cooling_time_batch = &field_data->CoolingTime[ntimes*nstrip]; 
        predator_prey_calculate_cooling_timescale( cooling_time_batch, input_batch, nstrip_res, data );    
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
