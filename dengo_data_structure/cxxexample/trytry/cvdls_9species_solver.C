
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


#include "cvdls_9species_solver.h"

cvdls_9species_data *cvdls_9species_setup_data(
    int *NumberOfFields, char ***FieldNames)
{
    int i;

    cvdls_9species_data *data = (cvdls_9species_data *) malloc(sizeof(cvdls_9species_data));
    
    /* allocate space for the scale related pieces */
    for (i = 0; i< 10 ; i++){
    data->scale[i] = 1.0;
    }
    
    /*initialize temperature so it wont crash*/
    data->Ts[0] = 1000.0;

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
    
    cvdls_9species_read_rate_tables(data);
    fprintf(stderr, "Successfully read in rate tables.\n");

    cvdls_9species_read_cooling_tables(data);
    fprintf(stderr, "Successfully read in cooling rate tables.\n");
    
    cvdls_9species_read_gamma(data);
    fprintf(stderr, "Successfully read in gamma tables. \n");

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
    return data;

}


int cvdls_9species_main(int argc, char** argv)
{
    cvdls_9species_data *data = cvdls_9species_setup_data(NULL, NULL);

    /* Initial conditions */
    
    char const *filename;
    if (argc > 1){
        filename = argv[1];
    } else{
        filename = "cvdls_9species_initial_conditions.h5";
    }
    
    hid_t file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    
    if (file_id < 0) {fprintf(stderr, "Failed to open %s so dying.\n", filename);
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

    int N = 10;

    double *atol, *rtol;
    atol = (double *) alloca(N * dims * sizeof(double));
    rtol = (double *) alloca(N * dims * sizeof(double));

    double *tics = (double *) alloca(dims * sizeof(double));
    double *ics = (double *) alloca(dims * N * sizeof(double));
    double *input = (double *) alloca(dims * N * sizeof(double));
    
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
                    tics[j*N + i], atol[j*N + i]);
        }
    }
    i++;
    

    H5Fclose(file_id);

    double dtf = 3.1557e+10;
    double dt = -1.0;
    double z = -1.0;
    for (i = 0; i < dims * N; i++) input[i] = ics[i];
    double ttot;
    ttot = dengo_evolve_cvdls_9species(dtf, dt, z, input, rtol, atol, dims, data);

    /* Write results to HDF5 file */
    file_id = H5Fcreate("cvdls_9species_solution.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
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
    
    double temperature[dims];
    for (j = 0; j < dims; j++) {
    	temperature[j] = data->Ts[j];
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
 



double dengo_evolve_cvdls_9species (double dtf, double &dt, double z, double *input,
            double *rtol, double *atol, long long dims, cvdls_9species_data *data) {
    int i, j;
    hid_t file_id;
    fprintf(stderr, "  ncells = % 3i\n", (int) dims);

    int N = 10;
    for (i = 0; i<dims; i++) {
      j = i * N;
        
          input[j] /= 2.01588 ;
          atol[j] /= 2.01588  ;
        
        j++;
      
        
          input[j] /= 2.01588 ;
          atol[j] /= 2.01588  ;
        
        j++;
      
        
          input[j] /= 1.00794 ;
          atol[j] /= 1.00794  ;
        
        j++;
      
        
          input[j] /= 1.00794 ;
          atol[j] /= 1.00794  ;
        
        j++;
      
        
          input[j] /= 1.00794 ;
          atol[j] /= 1.00794  ;
        
        j++;
      
        
          input[j] /= 4.002602 ;
          atol[j] /= 4.002602  ;
        
        j++;
      
        
          input[j] /= 4.002602 ;
          atol[j] /= 4.002602  ;
        
        j++;
      
        
          input[j] /= 4.002602 ;
          atol[j] /= 4.002602  ;
        
        j++;
      
        
          input[j] /= 1.0 ;
          atol[j] /= 1.0  ;
        
        j++;
      
        
        j++;
      
    }
    ensure_electron_consistency(input, dims, N);

    rhs_f f = calculate_rhs_cvdls_9species;
    jac_f jf = calculate_jacobian_cvdls_9species;
    if (dt < 0) dt = dtf / 1.0e3;
    data->current_z = z;
    int niter = 0;
    int siter = 0;
    double ttot = 0;
    double *scale = (double *) alloca(dims * N * sizeof(double));
    double *prev = (double *) alloca(dims * N * sizeof(double));
    for (i = 0; i < dims * N; i++) scale[i] = input[i];
    for (i = 0; i < dims * N; i++) prev[i] = input[i];
    double *u0 = (double *) alloca(N*dims*sizeof(double));
    double *s  = (double *) alloca(N*sizeof(double));
    double *gu = (double *) alloca(N*dims*sizeof(double));
    double *Ju = (double *) alloca(N*N*dims*sizeof(double));
    double floor_value = 1e-25;
            

        /* f and jf are function for evaluating rhs and jac for the solver
        * input: {double array}
        * rtol : {double} scalar (relative tolerance)
        * atol : {double array} vector (absolue tolerance)
        * NSPECIES: {int}
        * 
        */
        // int rv = 1;
        
        //fprintf(stderr, "ttot: %0.5g\n", ttot);
        // fprintf(stderr, "---------dt = %0.5g-------\n", dt); 
        
        double *input_N = (double *) alloca(N * sizeof(double)); 
        double *atol_N  = (double *) alloca(N * sizeof(double));
        double *rtol_N  = (double *) alloca(N * sizeof(double));
        double *prev_N  = (double *) alloca(N * sizeof(double));

        
        for (int d = 0; d < dims; d++){

            // fprintf(stderr, "nth strip: %d", d);

            // copy array which can be passed to the solver
            for (i = 0; i < N; i++){ 
                input_N[i] = input[d*N + i];
                atol_N[ i]  = atol[ d*N + i];
                rtol_N[ i]  = rtol[ d*N + i];
                prev_N[ i]  = prev[ d*N + i];
                
                // this is being passed around 
                // passively by the "dengo_rate_data" 
                // will have to fix it for openmp
                data->scale[i] = input[d*N + i];
            }
            // initialize a dt for the solver
            dt = dtf / 1.0e3;
            ttot = 0.0;
            siter = 0;
            //fprintf(stderr, "%d th strip calculation \n", d);
            while (ttot < dtf) {
                int rv = cvodes_main_solver( f, jf, input_N, rtol_N ,atol_N, NSPECIES, data, &dt);
                // fprintf(stderr, "%d th strip: %d iterations, time: %0.5g\n", d, siter, ttot );

                for (i = 0; i < N; i++) {
                    if (input[i] < 0) {
                        rv = 1;
                        break;
                    }
                }

                ttot += dt;
	            dt = DMIN(dt * 1.1, dtf - ttot);
                if (rv == 0) {
                // fprintf(stderr, "succesful integratoin\n");

	                if (siter == 50000) break;
	                siter++;
                    if (siter % 10000 == 0) {
                        fprintf(stderr, "Successful Iteration[%d]: (%0.4g) %0.16g / %0.16g\n", 
                                siter, dt, ttot, dtf);
                    }

	    
	            for (i = 0; i < N; i++) prev_N[i] = input_N[i];
                // for (i = 0; i < dims * N; i++) {     
                //    if (input[i] < floor_value) {
                //      input[i] = floor_value;
                //    }
                //}
                } else {
                    // fprintf(stderr, "failed\n");
                    dt /= 2.0;
                    // for (i = 0; i < N; i++) input_N[i] = prev_N[i];
                    // if (dt/dtf < 1e-30)  {
                    //     fprintf(stderr, "Dying! dt/dtf = %0.5g\n", dt/dtf);
                    //     break;
                    //}
                }
                siter++;
            }
            //fprintf(stderr, "%d strip finished at t = %0.5g \n",d, ttot);
            // should copy the results back to input[i] from input_N
            for (i = 0; i < N; i++){ 
                input[d*N +i] = input_N[i] ;
                // fprintf(stderr, "%d: %0.5g\n",i, input_N[i]);
            }
        }
    /* fprintf(stderr, "End: %0.5g / %0.5g (%0.5g)\n",
       ttot, dtf, dtf-ttot); */
    for (i = 0; i<dims; i++) {
      j = i * N;
        
          input[j] *= 2.01588 ;
          atol[j] *= 2.01588 ;
        
        j++;
      
        
          input[j] *= 2.01588;
          atol[j] *= 2.01588 ;
        
        j++;
      
        
          input[j] *= 1.00794;
          atol[j] *= 1.00794 ;
        
        j++;
      
        
          input[j] *= 1.00794;
          atol[j] *= 1.00794 ;
        
        j++;
      
        
          input[j] *= 1.00794;
          atol[j] *= 1.00794 ;
        
        j++;
      
        
          input[j] *= 4.002602;
          atol[j] *= 4.002602 ;
        
        j++;
      
        
          input[j] *= 4.002602;
          atol[j] *= 4.002602 ;
        
        j++;
      
        
          input[j] *= 4.002602;
          atol[j] *= 4.002602 ;
        
        j++;
      
        
          input[j] *= 1.0 ;
          atol[j] *= 1.0  ;
        
        j++;
      
        
        j++;
      
    }
    return ttot;
}
 


void cvdls_9species_read_rate_tables(cvdls_9species_data *data)
{
    hid_t file_id = H5Fopen("cvdls_9species_tables.h5", H5F_ACC_RDONLY, H5P_DEFAULT);
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

void cvdls_9species_read_cooling_tables(cvdls_9species_data *data)
{

    hid_t file_id = H5Fopen("cvdls_9species_tables.h5", H5F_ACC_RDONLY, H5P_DEFAULT);
    /* Allocate the correct number of rate tables */

    H5Fclose(file_id);
}

void cvdls_9species_read_gamma(cvdls_9species_data *data)
{

    hid_t file_id = H5Fopen("cvdls_9species_tables.h5", H5F_ACC_RDONLY, H5P_DEFAULT);
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
 


void cvdls_9species_calculate_temperature(cvdls_9species_data *data,
                        double *input, int nstrip, int nchem)
{
    int i, j;
    double density;
    double kb = 1.3806504e-16; // Boltzmann constant [erg/K]
    double mh = 1.67e-24;
    double gamma = 5.e0/3.e0;
    
    
    double gammaH2 = 7.e0/5.e0; // Should be a function of temperature
    	   	     		// this is a temporary solution
    double T,x, expx, Tnew;
    

    
    
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
    
    /* define scale */
    double scale;

    for (i = 0; i<nstrip; i++) {
        j = i * nchem;
        scale = data->scale[j];
        H2_1 = input[j]*scale;
        j++;
    
        scale = data->scale[j];
        H2_2 = input[j]*scale;
        j++;
    
        scale = data->scale[j];
        H_1 = input[j]*scale;
        j++;
    
        scale = data->scale[j];
        H_2 = input[j]*scale;
        j++;
    
        scale = data->scale[j];
        H_m0 = input[j]*scale;
        j++;
    
        scale = data->scale[j];
        He_1 = input[j]*scale;
        j++;
    
        scale = data->scale[j];
        He_2 = input[j]*scale;
        j++;
    
        scale = data->scale[j];
        He_3 = input[j]*scale;
        j++;
    
        scale = data->scale[j];
        de = input[j]*scale;
        j++;
    
        scale = data->scale[j];
        ge = input[j]*scale;
        j++;
    
        density = 2.01588*H2_1 + 2.01588*H2_2 + 1.00794*H_1 + 1.00794*H_2 + 1.00794*H_m0 + 4.002602*He_1 + 4.002602*He_2 + 4.002602*He_3;
        
    
        
        
        // Initiate the "guess" temperature
        T = data->Ts[i];
        Tnew = T + 1.0;
        double dge_dT;
        double dge;

        

        double gammaH2_1;
        double dgammaH2_1_dT;
        

        double gammaH2_2;
        double dgammaH2_2_dT;
        
        
        
        while ( abs(T - Tnew) > 0.1 ){
        // We do Newton's Iteration to calculate the temperature
        // Since gammaH2 is dependent on the temperature too!

        T = data->Ts[i];

        cvdls_9species_interpolate_gamma(data, i);
        
        gammaH2_1 = data->gammaH2_1[i];
        dgammaH2_1_dT = data->dgammaH2_1_dT[i];
        // fprintf(stderr, ":gammaH2_1 %0.5g , dgammaH2_1_dT: %.5g \n", gammaH2_1, dgammaH2_1_dT  );
        
        gammaH2_2 = data->gammaH2_2[i];
        dgammaH2_2_dT = data->dgammaH2_2_dT[i];
        // fprintf(stderr, ":gammaH2_2 %0.5g , dgammaH2_2_dT: %.5g \n", gammaH2_2, dgammaH2_2_dT);
        
       
        
        // update gammaH2

        

        // The derivatives of  sum (nkT/(gamma - 1)/mh/density) - ge
        // This is the function we want to minimize
        // which should only be dependent on the first part
        dge_dT = T*kb*(-H2_1*dgammaH2_1_dT/pow(gammaH2_1 - 1.0, 2) - H2_2*dgammaH2_2_dT/pow(gammaH2_2 - 1.0, 2))/(density*mh) + kb*(H2_1/(gammaH2_1 - 1.0) + H2_2/(gammaH2_2 - 1.0) + H_1/(gamma - 1.0) + H_2/(gamma - 1.0) + H_m0/(gamma - 1.0) + He_1/(gamma - 1.0) + He_2/(gamma - 1.0) + He_3/(gamma - 1.0) + de/(gamma - 1.0))/(density*mh);
        
        //This is the change in ge for each iteration
        dge = T*kb*(H2_1/(gammaH2_1 - 1.0) + H2_2/(gammaH2_2 - 1.0) + H_1/(gamma - 1.0) + H_2/(gamma - 1.0) + H_m0/(gamma - 1.0) + He_1/(gamma - 1.0) + He_2/(gamma - 1.0) + He_3/(gamma - 1.0) + de/(gamma - 1.0))/(density*mh) - ge;

        Tnew = T - dge/dge_dT;
        data->Ts[i] = Tnew;

        // fprintf(stderr, "T: %0.5g ; Tnew: %0.5g; dge_dT: %.5g, dge: %.5g, ge: %.5g \n", T,Tnew, dge_dT, dge, ge);
        }
        // fprintf(stderr,"---------------------\n");
        data->Ts[i] = Tnew;


        // fprintf(stderr,"T : %0.5g, density : %0.5g, d_gammaH2: %0.5g \n", Tnew, density, gammaH2 - 7./5.);


        

        if (data->Ts[i] < data->bounds[0]) {
            data->Ts[i] = data->bounds[0];
        } else if (data->Ts[i] > data->bounds[1]) {
            data->Ts[i] = data->bounds[1];
        }
        data->logTs[i] = log(data->Ts[i]);
        data->invTs[i] = 1.0 / data->Ts[i];
	    data->dTs_ge[i] = 
        density*mh/(kb*(H2_1/(gammaH2 - 1.0) + H2_2/(gammaH2 - 1.0) + H_1/(gamma - 1.0) + H_2/(gamma - 1.0) + H_m0/(gamma - 1.0) + He_1/(gamma - 1.0) + He_2/(gamma - 1.0) + He_3/(gamma - 1.0) + de/(gamma - 1.0)));
        /*fprintf(stderr, "T[%d] = % 0.16g, density = % 0.16g\n",
                i, data->Ts[i], density);*/
    }
         
}
 


/*
   This setup may be different than the user may anticipate, as a result
   of the lockstep timestep we use for a pencil beam through the grid.
   As such, it accepts the number of things to interpolate and makes
   assumptions about the sizes of the rates.
*/

/* This also requires no templating other than for the solver name...*/
void cvdls_9species_interpolate_rates(cvdls_9species_data *data,
                    int nstrip)
{
    int i, bin_id, zbin_id;
    double lb, t1, t2;
    double lbz, z1, z2;
    int no_photo = 0;
    lb = log(data->bounds[0]);
    lbz = log(data->z_bounds[0] + 1.0);
    /*fprintf(stderr, "lb = % 0.16g, ub = % 0.16g\n", lb, ub);*/
    for (i = 0; i < nstrip; i++) {
        data->bin_id[i] = bin_id = (int) (data->idbin * (data->logTs[i] - lb));
        if (data->bin_id[i] <= 0) {
            data->bin_id[i] = 0;
        } else if (data->bin_id[i] >= data->nbins) {
            data->bin_id[i] = data->nbins - 1;
        }
        t1 = (lb + (bin_id    ) * data->dbin);
        t2 = (lb + (bin_id + 1) * data->dbin);
        data->Tdef[i] = (data->logTs[i] - t1)/(t2 - t1);
        data->dT[i] = (t2 - t1);
        /*fprintf(stderr, "INTERP: %d, bin_id = %d, dT = % 0.16g, T = % 0.16g, logT = % 0.16g\n",
                i, data->bin_id[i], data->dT[i], data->Ts[i],
                data->logTs[i]);*/
    }
    
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
    
    for (i = 0; i < nstrip; i++) {
        bin_id = data->bin_id[i];
        data->rs_k01[i] = data->r_k01[bin_id] +
            data->Tdef[i] * (data->r_k01[bin_id+1] - data->r_k01[bin_id]);
        data->drs_k01[i] = (data->r_k01[bin_id+1] - data->r_k01[bin_id]);
        data->drs_k01[i] /= data->dT[i];
	data->drs_k01[i] *= data->invTs[i];
    }
    
    for (i = 0; i < nstrip; i++) {
        bin_id = data->bin_id[i];
        data->rs_k02[i] = data->r_k02[bin_id] +
            data->Tdef[i] * (data->r_k02[bin_id+1] - data->r_k02[bin_id]);
        data->drs_k02[i] = (data->r_k02[bin_id+1] - data->r_k02[bin_id]);
        data->drs_k02[i] /= data->dT[i];
	data->drs_k02[i] *= data->invTs[i];
    }
    
    for (i = 0; i < nstrip; i++) {
        bin_id = data->bin_id[i];
        data->rs_k03[i] = data->r_k03[bin_id] +
            data->Tdef[i] * (data->r_k03[bin_id+1] - data->r_k03[bin_id]);
        data->drs_k03[i] = (data->r_k03[bin_id+1] - data->r_k03[bin_id]);
        data->drs_k03[i] /= data->dT[i];
	data->drs_k03[i] *= data->invTs[i];
    }
    
    for (i = 0; i < nstrip; i++) {
        bin_id = data->bin_id[i];
        data->rs_k04[i] = data->r_k04[bin_id] +
            data->Tdef[i] * (data->r_k04[bin_id+1] - data->r_k04[bin_id]);
        data->drs_k04[i] = (data->r_k04[bin_id+1] - data->r_k04[bin_id]);
        data->drs_k04[i] /= data->dT[i];
	data->drs_k04[i] *= data->invTs[i];
    }
    
    for (i = 0; i < nstrip; i++) {
        bin_id = data->bin_id[i];
        data->rs_k05[i] = data->r_k05[bin_id] +
            data->Tdef[i] * (data->r_k05[bin_id+1] - data->r_k05[bin_id]);
        data->drs_k05[i] = (data->r_k05[bin_id+1] - data->r_k05[bin_id]);
        data->drs_k05[i] /= data->dT[i];
	data->drs_k05[i] *= data->invTs[i];
    }
    
    for (i = 0; i < nstrip; i++) {
        bin_id = data->bin_id[i];
        data->rs_k06[i] = data->r_k06[bin_id] +
            data->Tdef[i] * (data->r_k06[bin_id+1] - data->r_k06[bin_id]);
        data->drs_k06[i] = (data->r_k06[bin_id+1] - data->r_k06[bin_id]);
        data->drs_k06[i] /= data->dT[i];
	data->drs_k06[i] *= data->invTs[i];
    }
    
    for (i = 0; i < nstrip; i++) {
        bin_id = data->bin_id[i];
        data->rs_k07[i] = data->r_k07[bin_id] +
            data->Tdef[i] * (data->r_k07[bin_id+1] - data->r_k07[bin_id]);
        data->drs_k07[i] = (data->r_k07[bin_id+1] - data->r_k07[bin_id]);
        data->drs_k07[i] /= data->dT[i];
	data->drs_k07[i] *= data->invTs[i];
    }
    
    for (i = 0; i < nstrip; i++) {
        bin_id = data->bin_id[i];
        data->rs_k08[i] = data->r_k08[bin_id] +
            data->Tdef[i] * (data->r_k08[bin_id+1] - data->r_k08[bin_id]);
        data->drs_k08[i] = (data->r_k08[bin_id+1] - data->r_k08[bin_id]);
        data->drs_k08[i] /= data->dT[i];
	data->drs_k08[i] *= data->invTs[i];
    }
    
    for (i = 0; i < nstrip; i++) {
        bin_id = data->bin_id[i];
        data->rs_k09[i] = data->r_k09[bin_id] +
            data->Tdef[i] * (data->r_k09[bin_id+1] - data->r_k09[bin_id]);
        data->drs_k09[i] = (data->r_k09[bin_id+1] - data->r_k09[bin_id]);
        data->drs_k09[i] /= data->dT[i];
	data->drs_k09[i] *= data->invTs[i];
    }
    
    for (i = 0; i < nstrip; i++) {
        bin_id = data->bin_id[i];
        data->rs_k10[i] = data->r_k10[bin_id] +
            data->Tdef[i] * (data->r_k10[bin_id+1] - data->r_k10[bin_id]);
        data->drs_k10[i] = (data->r_k10[bin_id+1] - data->r_k10[bin_id]);
        data->drs_k10[i] /= data->dT[i];
	data->drs_k10[i] *= data->invTs[i];
    }
    
    for (i = 0; i < nstrip; i++) {
        bin_id = data->bin_id[i];
        data->rs_k11[i] = data->r_k11[bin_id] +
            data->Tdef[i] * (data->r_k11[bin_id+1] - data->r_k11[bin_id]);
        data->drs_k11[i] = (data->r_k11[bin_id+1] - data->r_k11[bin_id]);
        data->drs_k11[i] /= data->dT[i];
	data->drs_k11[i] *= data->invTs[i];
    }
    
    for (i = 0; i < nstrip; i++) {
        bin_id = data->bin_id[i];
        data->rs_k12[i] = data->r_k12[bin_id] +
            data->Tdef[i] * (data->r_k12[bin_id+1] - data->r_k12[bin_id]);
        data->drs_k12[i] = (data->r_k12[bin_id+1] - data->r_k12[bin_id]);
        data->drs_k12[i] /= data->dT[i];
	data->drs_k12[i] *= data->invTs[i];
    }
    
    for (i = 0; i < nstrip; i++) {
        bin_id = data->bin_id[i];
        data->rs_k13[i] = data->r_k13[bin_id] +
            data->Tdef[i] * (data->r_k13[bin_id+1] - data->r_k13[bin_id]);
        data->drs_k13[i] = (data->r_k13[bin_id+1] - data->r_k13[bin_id]);
        data->drs_k13[i] /= data->dT[i];
	data->drs_k13[i] *= data->invTs[i];
    }
    
    for (i = 0; i < nstrip; i++) {
        bin_id = data->bin_id[i];
        data->rs_k14[i] = data->r_k14[bin_id] +
            data->Tdef[i] * (data->r_k14[bin_id+1] - data->r_k14[bin_id]);
        data->drs_k14[i] = (data->r_k14[bin_id+1] - data->r_k14[bin_id]);
        data->drs_k14[i] /= data->dT[i];
	data->drs_k14[i] *= data->invTs[i];
    }
    
    for (i = 0; i < nstrip; i++) {
        bin_id = data->bin_id[i];
        data->rs_k15[i] = data->r_k15[bin_id] +
            data->Tdef[i] * (data->r_k15[bin_id+1] - data->r_k15[bin_id]);
        data->drs_k15[i] = (data->r_k15[bin_id+1] - data->r_k15[bin_id]);
        data->drs_k15[i] /= data->dT[i];
	data->drs_k15[i] *= data->invTs[i];
    }
    
    for (i = 0; i < nstrip; i++) {
        bin_id = data->bin_id[i];
        data->rs_k16[i] = data->r_k16[bin_id] +
            data->Tdef[i] * (data->r_k16[bin_id+1] - data->r_k16[bin_id]);
        data->drs_k16[i] = (data->r_k16[bin_id+1] - data->r_k16[bin_id]);
        data->drs_k16[i] /= data->dT[i];
	data->drs_k16[i] *= data->invTs[i];
    }
    
    for (i = 0; i < nstrip; i++) {
        bin_id = data->bin_id[i];
        data->rs_k17[i] = data->r_k17[bin_id] +
            data->Tdef[i] * (data->r_k17[bin_id+1] - data->r_k17[bin_id]);
        data->drs_k17[i] = (data->r_k17[bin_id+1] - data->r_k17[bin_id]);
        data->drs_k17[i] /= data->dT[i];
	data->drs_k17[i] *= data->invTs[i];
    }
    
    for (i = 0; i < nstrip; i++) {
        bin_id = data->bin_id[i];
        data->rs_k18[i] = data->r_k18[bin_id] +
            data->Tdef[i] * (data->r_k18[bin_id+1] - data->r_k18[bin_id]);
        data->drs_k18[i] = (data->r_k18[bin_id+1] - data->r_k18[bin_id]);
        data->drs_k18[i] /= data->dT[i];
	data->drs_k18[i] *= data->invTs[i];
    }
    
    for (i = 0; i < nstrip; i++) {
        bin_id = data->bin_id[i];
        data->rs_k19[i] = data->r_k19[bin_id] +
            data->Tdef[i] * (data->r_k19[bin_id+1] - data->r_k19[bin_id]);
        data->drs_k19[i] = (data->r_k19[bin_id+1] - data->r_k19[bin_id]);
        data->drs_k19[i] /= data->dT[i];
	data->drs_k19[i] *= data->invTs[i];
    }
    
    for (i = 0; i < nstrip; i++) {
        bin_id = data->bin_id[i];
        data->rs_k21[i] = data->r_k21[bin_id] +
            data->Tdef[i] * (data->r_k21[bin_id+1] - data->r_k21[bin_id]);
        data->drs_k21[i] = (data->r_k21[bin_id+1] - data->r_k21[bin_id]);
        data->drs_k21[i] /= data->dT[i];
	data->drs_k21[i] *= data->invTs[i];
    }
    
    for (i = 0; i < nstrip; i++) {
        bin_id = data->bin_id[i];
        data->rs_k22[i] = data->r_k22[bin_id] +
            data->Tdef[i] * (data->r_k22[bin_id+1] - data->r_k22[bin_id]);
        data->drs_k22[i] = (data->r_k22[bin_id+1] - data->r_k22[bin_id]);
        data->drs_k22[i] /= data->dT[i];
	data->drs_k22[i] *= data->invTs[i];
    }
    

}
 


void cvdls_9species_interpolate_gamma(cvdls_9species_data *data,
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
    
    data->bin_id[i] = bin_id = (int) (data->idbin * (data->logTs[i] - lb));
    if (data->bin_id[i] <= 0) {
        data->bin_id[i] = 0;
    } else if (data->bin_id[i] >= data->nbins) {
        data->bin_id[i] = data->nbins - 1;
    }
    t1 = (lb + (bin_id    ) * data->dbin);
    t2 = (lb + (bin_id + 1) * data->dbin);
    data->Tdef[i] = (data->logTs[i] - t1)/(t2 - t1);
    data->dT[i] = (t2 - t1);

    
    
    bin_id = data->bin_id[i];
    data->gammaH2_2[i] = data->g_gammaH2_2[bin_id] +
        data->Tdef[i] * (data->g_gammaH2_2[bin_id+1] - data->g_gammaH2_2[bin_id]);

    data->dgammaH2_2_dT[i] = data->g_dgammaH2_2_dT[bin_id] +
        data->Tdef[i] * (data->g_dgammaH2_2_dT[bin_id+1] 
        - data->g_dgammaH2_2_dT[bin_id]);
    
    
    bin_id = data->bin_id[i];
    data->gammaH2_1[i] = data->g_gammaH2_1[bin_id] +
        data->Tdef[i] * (data->g_gammaH2_1[bin_id+1] - data->g_gammaH2_1[bin_id]);

    data->dgammaH2_1_dT[i] = data->g_dgammaH2_1_dT[bin_id] +
        data->Tdef[i] * (data->g_dgammaH2_1_dT[bin_id+1] 
        - data->g_dgammaH2_1_dT[bin_id]);
    
       
    }






void ensure_electron_consistency(double *input, int nstrip, int nchem)
{
    int i, j;

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
    
    double scale;

    for (i = 0; i<nstrip; i++) {
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




void temperature_from_mass_density(double *input, int nstrip,
                                   int nchem, double *strip_temperature)
{
    int i, j;
    double density;
    double kb = 1.3806504e-16; // Boltzmann constant [erg/K]
    double mh = 1.67e-24;
    double gamma = 5.e0/3.e0;
    
    double gammaH2 = 7.e0/5.e0; // Should be a function of temperature
    	   	     		// this is a temporary solution
    
    double T =  1000.0; // THIS IS TEMPORARY!!! DELTETE!!

    
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
    
    double scale;

    for (i = 0; i<nstrip; i++) {
        j = i * nchem;
        H2_1 = input[j];
        
        H2_1 /= 2.01588 * mh;
        
        /*fprintf(stderr, "H2_1[%d] = % 0.16g\n",
                i, H2_1);*/
        j++;
    
        H2_2 = input[j];
        
        H2_2 /= 2.01588 * mh;
        
        /*fprintf(stderr, "H2_2[%d] = % 0.16g\n",
                i, H2_2);*/
        j++;
    
        H_1 = input[j];
        
        H_1 /= 1.00794 * mh;
        
        /*fprintf(stderr, "H_1[%d] = % 0.16g\n",
                i, H_1);*/
        j++;
    
        H_2 = input[j];
        
        H_2 /= 1.00794 * mh;
        
        /*fprintf(stderr, "H_2[%d] = % 0.16g\n",
                i, H_2);*/
        j++;
    
        H_m0 = input[j];
        
        H_m0 /= 1.00794 * mh;
        
        /*fprintf(stderr, "H_m0[%d] = % 0.16g\n",
                i, H_m0);*/
        j++;
    
        He_1 = input[j];
        
        He_1 /= 4.002602 * mh;
        
        /*fprintf(stderr, "He_1[%d] = % 0.16g\n",
                i, He_1);*/
        j++;
    
        He_2 = input[j];
        
        He_2 /= 4.002602 * mh;
        
        /*fprintf(stderr, "He_2[%d] = % 0.16g\n",
                i, He_2);*/
        j++;
    
        He_3 = input[j];
        
        He_3 /= 4.002602 * mh;
        
        /*fprintf(stderr, "He_3[%d] = % 0.16g\n",
                i, He_3);*/
        j++;
    
        de = input[j];
        
        de /= 1.0 * mh;
        
        /*fprintf(stderr, "de[%d] = % 0.16g\n",
                i, de);*/
        j++;
    
        ge = input[j];
        
        /*fprintf(stderr, "ge[%d] = % 0.16g\n",
                i, ge);*/
        j++;
    
        density = 2.01588*H2_1 + 2.01588*H2_2 + 1.00794*H_1 + 1.00794*H_2 + 1.00794*H_m0 + 4.002602*He_1 + 4.002602*He_2 + 4.002602*He_3;
        strip_temperature[i] = density*ge*mh/(kb*(H2_1/(gammaH2 - 1.0) + H2_2/(gammaH2 - 1.0) + H_1/(gamma - 1.0) + H_2/(gamma - 1.0) + H_m0/(gamma - 1.0) + He_1/(gamma - 1.0) + He_2/(gamma - 1.0) + He_3/(gamma - 1.0) + de/(gamma - 1.0)));
        if (strip_temperature[i] < 1.0)
            strip_temperature[i] = 1.0;
    }
         
}
 


int calculate_jacobian_cvdls_9species
                                        ( realtype t,
                                        N_Vector y, N_Vector fy,
                                        SUNMatrix J, void *user_data,
                                        N_Vector tmp1, N_Vector tmp2,
                                        N_Vector tmp3)
{
    /* We iterate over all of the rates */
    /* Calcuate temperature first */
    
    int nstrip = 1;
    int nchem = 10;

    cvdls_9species_data *data = (cvdls_9species_data*)user_data; 
    

    int i, j;
    j = 0;
    /* change N_Vector back to an array */
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

    cvdls_9species_calculate_temperature(data, y_arr, nstrip, nchem);

    cvdls_9species_interpolate_rates(data, nstrip);

    /* Now We set up some temporaries */

    double *Tge = data->dTs_ge;
    double *k01 = data->rs_k01;
    double *rk01 = data->drs_k01;
    double *k02 = data->rs_k02;
    double *rk02 = data->drs_k02;
    double *k03 = data->rs_k03;
    double *rk03 = data->drs_k03;
    double *k04 = data->rs_k04;
    double *rk04 = data->drs_k04;
    double *k05 = data->rs_k05;
    double *rk05 = data->drs_k05;
    double *k06 = data->rs_k06;
    double *rk06 = data->drs_k06;
    double *k07 = data->rs_k07;
    double *rk07 = data->drs_k07;
    double *k08 = data->rs_k08;
    double *rk08 = data->drs_k08;
    double *k09 = data->rs_k09;
    double *rk09 = data->drs_k09;
    double *k10 = data->rs_k10;
    double *rk10 = data->drs_k10;
    double *k11 = data->rs_k11;
    double *rk11 = data->drs_k11;
    double *k12 = data->rs_k12;
    double *rk12 = data->drs_k12;
    double *k13 = data->rs_k13;
    double *rk13 = data->drs_k13;
    double *k14 = data->rs_k14;
    double *rk14 = data->drs_k14;
    double *k15 = data->rs_k15;
    double *rk15 = data->drs_k15;
    double *k16 = data->rs_k16;
    double *rk16 = data->drs_k16;
    double *k17 = data->rs_k17;
    double *rk17 = data->drs_k17;
    double *k18 = data->rs_k18;
    double *rk18 = data->drs_k18;
    double *k19 = data->rs_k19;
    double *rk19 = data->drs_k19;
    double *k21 = data->rs_k21;
    double *rk21 = data->drs_k21;
    double *k22 = data->rs_k22;
    double *rk22 = data->drs_k22;
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

    double mh = 1.67e-24;
    double mdensity;
    
    int jj;
    jj = 0;
    double scale, scale1, scale2;

    for (i = 0; i<nstrip; i++) {
        j = i*nchem;
        mdensity = 0.0;
        z = data->current_z;
        scale = data->scale[j];
        H2_1 = Ith( y, 1  )*scale;
        //fprintf(stderr,"from jac: H2_1 = %.3g\n, Ith(y1) = %.3g \n" , H2_1, Ith(y,1 ) );
        
        mdensity += H2_1 * 2.01588;
        

        j++;
        
        scale = data->scale[j];
        H2_2 = Ith( y, 2  )*scale;
        //fprintf(stderr,"from jac: H2_2 = %.3g\n, Ith(y2) = %.3g \n" , H2_2, Ith(y,2 ) );
        
        mdensity += H2_2 * 2.01588;
        

        j++;
        
        scale = data->scale[j];
        H_1 = Ith( y, 3  )*scale;
        //fprintf(stderr,"from jac: H_1 = %.3g\n, Ith(y3) = %.3g \n" , H_1, Ith(y,3 ) );
        
        mdensity += H_1 * 1.00794;
        

        j++;
        
        scale = data->scale[j];
        H_2 = Ith( y, 4  )*scale;
        //fprintf(stderr,"from jac: H_2 = %.3g\n, Ith(y4) = %.3g \n" , H_2, Ith(y,4 ) );
        
        mdensity += H_2 * 1.00794;
        

        j++;
        
        scale = data->scale[j];
        H_m0 = Ith( y, 5  )*scale;
        //fprintf(stderr,"from jac: H_m0 = %.3g\n, Ith(y5) = %.3g \n" , H_m0, Ith(y,5 ) );
        
        mdensity += H_m0 * 1.00794;
        

        j++;
        
        scale = data->scale[j];
        He_1 = Ith( y, 6  )*scale;
        //fprintf(stderr,"from jac: He_1 = %.3g\n, Ith(y6) = %.3g \n" , He_1, Ith(y,6 ) );
        
        mdensity += He_1 * 4.002602;
        

        j++;
        
        scale = data->scale[j];
        He_2 = Ith( y, 7  )*scale;
        //fprintf(stderr,"from jac: He_2 = %.3g\n, Ith(y7) = %.3g \n" , He_2, Ith(y,7 ) );
        
        mdensity += He_2 * 4.002602;
        

        j++;
        
        scale = data->scale[j];
        He_3 = Ith( y, 8  )*scale;
        //fprintf(stderr,"from jac: He_3 = %.3g\n, Ith(y8) = %.3g \n" , He_3, Ith(y,8 ) );
        
        mdensity += He_3 * 4.002602;
        

        j++;
        
        scale = data->scale[j];
        de = Ith( y, 9  )*scale;
        //fprintf(stderr,"from jac: de = %.3g\n, Ith(y9) = %.3g \n" , de, Ith(y,9 ) );
        

        j++;
        
        scale = data->scale[j];
        ge = Ith( y, 10  )*scale;
        //fprintf(stderr,"from jac: ge = %.3g\n, Ith(y10) = %.3g \n" , ge, Ith(y,10 ) );
        

        j++;
        
        double nH; 
        nH = H_1 + H_2 + 2.0*H2_1 + 2.0*H2_2;
        


        mdensity *= mh;

        j = i * nchem * nchem;
        //
        // Species: H2_1
        //
        
        // H2_1 by H2_1
        
        IJth(J, 1, 1 ) = -k11[i]*H_2 - k12[i]*de - k13[i]*H_1 + k21[i]*pow(H_1, 2);
        
        scale2 = data->scale[1 - 1];
        scale1 = data->scale[1 - 1];
        IJth(J, 1, 1) /= scale2/scale1;

        

        
        // H2_2 by H2_1
        
        IJth(J, 2, 1 ) = k11[i]*H_2;
        
        scale2 = data->scale[2 - 1];
        scale1 = data->scale[1 - 1];
        IJth(J, 2, 1) /= scale2/scale1;

        

        
        // H_1 by H2_1
        
        IJth(J, 3, 1 ) = k11[i]*H_2 + 2*k12[i]*de + 2*k13[i]*H_1 - 2*k21[i]*pow(H_1, 2);
        
        scale2 = data->scale[3 - 1];
        scale1 = data->scale[1 - 1];
        IJth(J, 3, 1) /= scale2/scale1;

        

        
        // H_2 by H2_1
        
        IJth(J, 4, 1 ) = -k11[i]*H_2;
        
        scale2 = data->scale[4 - 1];
        scale1 = data->scale[1 - 1];
        IJth(J, 4, 1) /= scale2/scale1;

        

        
        // H_m0 by H2_1
        
        IJth(J, 5, 1 ) = 0;
        
        scale2 = data->scale[5 - 1];
        scale1 = data->scale[1 - 1];
        IJth(J, 5, 1) /= scale2/scale1;

        

        
        // He_1 by H2_1
        
        IJth(J, 6, 1 ) = 0;
        
        scale2 = data->scale[6 - 1];
        scale1 = data->scale[1 - 1];
        IJth(J, 6, 1) /= scale2/scale1;

        

        
        // He_2 by H2_1
        
        IJth(J, 7, 1 ) = 0;
        
        scale2 = data->scale[7 - 1];
        scale1 = data->scale[1 - 1];
        IJth(J, 7, 1) /= scale2/scale1;

        

        
        // He_3 by H2_1
        
        IJth(J, 8, 1 ) = 0;
        
        scale2 = data->scale[8 - 1];
        scale1 = data->scale[1 - 1];
        IJth(J, 8, 1) /= scale2/scale1;

        

        
        // de by H2_1
        
        IJth(J, 9, 1 ) = 0;
        
        scale2 = data->scale[9 - 1];
        scale1 = data->scale[1 - 1];
        IJth(J, 9, 1) /= scale2/scale1;

        

        
        // ge by H2_1
        
        IJth(J, 10, 1 ) = 0;
        
        scale2 = data->scale[10 - 1];
        scale1 = data->scale[1 - 1];
        IJth(J, 10, 1) /= scale2/scale1;

        
        IJth(J, 10, 1 ) /= mdensity;
        

        
        
        //
        // Species: H2_2
        //
        
        // H2_1 by H2_2
        
        IJth(J, 1, 2 ) = k10[i]*H_1 + k19[i]*H_m0;
        
        scale2 = data->scale[1 - 1];
        scale1 = data->scale[2 - 1];
        IJth(J, 1, 2) /= scale2/scale1;

        

        
        // H2_2 by H2_2
        
        IJth(J, 2, 2 ) = -k10[i]*H_1 - k18[i]*de - k19[i]*H_m0;
        
        scale2 = data->scale[2 - 1];
        scale1 = data->scale[2 - 1];
        IJth(J, 2, 2) /= scale2/scale1;

        

        
        // H_1 by H2_2
        
        IJth(J, 3, 2 ) = -k10[i]*H_1 + 2*k18[i]*de + k19[i]*H_m0;
        
        scale2 = data->scale[3 - 1];
        scale1 = data->scale[2 - 1];
        IJth(J, 3, 2) /= scale2/scale1;

        

        
        // H_2 by H2_2
        
        IJth(J, 4, 2 ) = k10[i]*H_1;
        
        scale2 = data->scale[4 - 1];
        scale1 = data->scale[2 - 1];
        IJth(J, 4, 2) /= scale2/scale1;

        

        
        // H_m0 by H2_2
        
        IJth(J, 5, 2 ) = -k19[i]*H_m0;
        
        scale2 = data->scale[5 - 1];
        scale1 = data->scale[2 - 1];
        IJth(J, 5, 2) /= scale2/scale1;

        

        
        // He_1 by H2_2
        
        IJth(J, 6, 2 ) = 0;
        
        scale2 = data->scale[6 - 1];
        scale1 = data->scale[2 - 1];
        IJth(J, 6, 2) /= scale2/scale1;

        

        
        // He_2 by H2_2
        
        IJth(J, 7, 2 ) = 0;
        
        scale2 = data->scale[7 - 1];
        scale1 = data->scale[2 - 1];
        IJth(J, 7, 2) /= scale2/scale1;

        

        
        // He_3 by H2_2
        
        IJth(J, 8, 2 ) = 0;
        
        scale2 = data->scale[8 - 1];
        scale1 = data->scale[2 - 1];
        IJth(J, 8, 2) /= scale2/scale1;

        

        
        // de by H2_2
        
        IJth(J, 9, 2 ) = -k18[i]*de;
        
        scale2 = data->scale[9 - 1];
        scale1 = data->scale[2 - 1];
        IJth(J, 9, 2) /= scale2/scale1;

        

        
        // ge by H2_2
        
        IJth(J, 10, 2 ) = 0;
        
        scale2 = data->scale[10 - 1];
        scale1 = data->scale[2 - 1];
        IJth(J, 10, 2) /= scale2/scale1;

        
        IJth(J, 10, 2 ) /= mdensity;
        

        
        
        //
        // Species: H_1
        //
        
        // H2_1 by H_1
        
        IJth(J, 1, 3 ) = k08[i]*H_m0 + k10[i]*H2_2 - k13[i]*H2_1 + 2*k21[i]*H2_1*H_1 + 3*k22[i]*pow(H_1, 2);
        
        scale2 = data->scale[1 - 1];
        scale1 = data->scale[3 - 1];
        IJth(J, 1, 3) /= scale2/scale1;

        

        
        // H2_2 by H_1
        
        IJth(J, 2, 3 ) = k09[i]*H_2 - k10[i]*H2_2;
        
        scale2 = data->scale[2 - 1];
        scale1 = data->scale[3 - 1];
        IJth(J, 2, 3) /= scale2/scale1;

        

        
        // H_1 by H_1
        
        IJth(J, 3, 3 ) = -k01[i]*de - k07[i]*de - k08[i]*H_m0 - k09[i]*H_2 - k10[i]*H2_2 + 2*k13[i]*H2_1 + k15[i]*H_m0 - 4*k21[i]*H2_1*H_1 - 6*k22[i]*pow(H_1, 2);
        
        scale2 = data->scale[3 - 1];
        scale1 = data->scale[3 - 1];
        IJth(J, 3, 3) /= scale2/scale1;

        

        
        // H_2 by H_1
        
        IJth(J, 4, 3 ) = k01[i]*de - k09[i]*H_2 + k10[i]*H2_2;
        
        scale2 = data->scale[4 - 1];
        scale1 = data->scale[3 - 1];
        IJth(J, 4, 3) /= scale2/scale1;

        

        
        // H_m0 by H_1
        
        IJth(J, 5, 3 ) = k07[i]*de - k08[i]*H_m0 - k15[i]*H_m0;
        
        scale2 = data->scale[5 - 1];
        scale1 = data->scale[3 - 1];
        IJth(J, 5, 3) /= scale2/scale1;

        

        
        // He_1 by H_1
        
        IJth(J, 6, 3 ) = 0;
        
        scale2 = data->scale[6 - 1];
        scale1 = data->scale[3 - 1];
        IJth(J, 6, 3) /= scale2/scale1;

        

        
        // He_2 by H_1
        
        IJth(J, 7, 3 ) = 0;
        
        scale2 = data->scale[7 - 1];
        scale1 = data->scale[3 - 1];
        IJth(J, 7, 3) /= scale2/scale1;

        

        
        // He_3 by H_1
        
        IJth(J, 8, 3 ) = 0;
        
        scale2 = data->scale[8 - 1];
        scale1 = data->scale[3 - 1];
        IJth(J, 8, 3) /= scale2/scale1;

        

        
        // de by H_1
        
        IJth(J, 9, 3 ) = k01[i]*de - k07[i]*de + k08[i]*H_m0 + k15[i]*H_m0;
        
        scale2 = data->scale[9 - 1];
        scale1 = data->scale[3 - 1];
        IJth(J, 9, 3) /= scale2/scale1;

        

        
        // ge by H_1
        
        IJth(J, 10, 3 ) = 0;
        
        scale2 = data->scale[10 - 1];
        scale1 = data->scale[3 - 1];
        IJth(J, 10, 3) /= scale2/scale1;

        
        IJth(J, 10, 3 ) /= mdensity;
        

        
        
        //
        // Species: H_2
        //
        
        // H2_1 by H_2
        
        IJth(J, 1, 4 ) = -k11[i]*H2_1;
        
        scale2 = data->scale[1 - 1];
        scale1 = data->scale[4 - 1];
        IJth(J, 1, 4) /= scale2/scale1;

        

        
        // H2_2 by H_2
        
        IJth(J, 2, 4 ) = k09[i]*H_1 + k11[i]*H2_1 + k17[i]*H_m0;
        
        scale2 = data->scale[2 - 1];
        scale1 = data->scale[4 - 1];
        IJth(J, 2, 4) /= scale2/scale1;

        

        
        // H_1 by H_2
        
        IJth(J, 3, 4 ) = k02[i]*de - k09[i]*H_1 + k11[i]*H2_1 + 2*k16[i]*H_m0;
        
        scale2 = data->scale[3 - 1];
        scale1 = data->scale[4 - 1];
        IJth(J, 3, 4) /= scale2/scale1;

        

        
        // H_2 by H_2
        
        IJth(J, 4, 4 ) = -k02[i]*de - k09[i]*H_1 - k11[i]*H2_1 - k16[i]*H_m0 - k17[i]*H_m0;
        
        scale2 = data->scale[4 - 1];
        scale1 = data->scale[4 - 1];
        IJth(J, 4, 4) /= scale2/scale1;

        

        
        // H_m0 by H_2
        
        IJth(J, 5, 4 ) = -k16[i]*H_m0 - k17[i]*H_m0;
        
        scale2 = data->scale[5 - 1];
        scale1 = data->scale[4 - 1];
        IJth(J, 5, 4) /= scale2/scale1;

        

        
        // He_1 by H_2
        
        IJth(J, 6, 4 ) = 0;
        
        scale2 = data->scale[6 - 1];
        scale1 = data->scale[4 - 1];
        IJth(J, 6, 4) /= scale2/scale1;

        

        
        // He_2 by H_2
        
        IJth(J, 7, 4 ) = 0;
        
        scale2 = data->scale[7 - 1];
        scale1 = data->scale[4 - 1];
        IJth(J, 7, 4) /= scale2/scale1;

        

        
        // He_3 by H_2
        
        IJth(J, 8, 4 ) = 0;
        
        scale2 = data->scale[8 - 1];
        scale1 = data->scale[4 - 1];
        IJth(J, 8, 4) /= scale2/scale1;

        

        
        // de by H_2
        
        IJth(J, 9, 4 ) = -k02[i]*de + k17[i]*H_m0;
        
        scale2 = data->scale[9 - 1];
        scale1 = data->scale[4 - 1];
        IJth(J, 9, 4) /= scale2/scale1;

        

        
        // ge by H_2
        
        IJth(J, 10, 4 ) = 0;
        
        scale2 = data->scale[10 - 1];
        scale1 = data->scale[4 - 1];
        IJth(J, 10, 4) /= scale2/scale1;

        
        IJth(J, 10, 4 ) /= mdensity;
        

        
        
        //
        // Species: H_m0
        //
        
        // H2_1 by H_m0
        
        IJth(J, 1, 5 ) = k08[i]*H_1 + k19[i]*H2_2;
        
        scale2 = data->scale[1 - 1];
        scale1 = data->scale[5 - 1];
        IJth(J, 1, 5) /= scale2/scale1;

        

        
        // H2_2 by H_m0
        
        IJth(J, 2, 5 ) = k17[i]*H_2 - k19[i]*H2_2;
        
        scale2 = data->scale[2 - 1];
        scale1 = data->scale[5 - 1];
        IJth(J, 2, 5) /= scale2/scale1;

        

        
        // H_1 by H_m0
        
        IJth(J, 3, 5 ) = -k08[i]*H_1 + k14[i]*de + k15[i]*H_1 + 2*k16[i]*H_2 + k19[i]*H2_2;
        
        scale2 = data->scale[3 - 1];
        scale1 = data->scale[5 - 1];
        IJth(J, 3, 5) /= scale2/scale1;

        

        
        // H_2 by H_m0
        
        IJth(J, 4, 5 ) = -k16[i]*H_2 - k17[i]*H_2;
        
        scale2 = data->scale[4 - 1];
        scale1 = data->scale[5 - 1];
        IJth(J, 4, 5) /= scale2/scale1;

        

        
        // H_m0 by H_m0
        
        IJth(J, 5, 5 ) = -k08[i]*H_1 - k14[i]*de - k15[i]*H_1 - k16[i]*H_2 - k17[i]*H_2 - k19[i]*H2_2;
        
        scale2 = data->scale[5 - 1];
        scale1 = data->scale[5 - 1];
        IJth(J, 5, 5) /= scale2/scale1;

        

        
        // He_1 by H_m0
        
        IJth(J, 6, 5 ) = 0;
        
        scale2 = data->scale[6 - 1];
        scale1 = data->scale[5 - 1];
        IJth(J, 6, 5) /= scale2/scale1;

        

        
        // He_2 by H_m0
        
        IJth(J, 7, 5 ) = 0;
        
        scale2 = data->scale[7 - 1];
        scale1 = data->scale[5 - 1];
        IJth(J, 7, 5) /= scale2/scale1;

        

        
        // He_3 by H_m0
        
        IJth(J, 8, 5 ) = 0;
        
        scale2 = data->scale[8 - 1];
        scale1 = data->scale[5 - 1];
        IJth(J, 8, 5) /= scale2/scale1;

        

        
        // de by H_m0
        
        IJth(J, 9, 5 ) = k08[i]*H_1 + k14[i]*de + k15[i]*H_1 + k17[i]*H_2;
        
        scale2 = data->scale[9 - 1];
        scale1 = data->scale[5 - 1];
        IJth(J, 9, 5) /= scale2/scale1;

        

        
        // ge by H_m0
        
        IJth(J, 10, 5 ) = 0;
        
        scale2 = data->scale[10 - 1];
        scale1 = data->scale[5 - 1];
        IJth(J, 10, 5) /= scale2/scale1;

        
        IJth(J, 10, 5 ) /= mdensity;
        

        
        
        //
        // Species: He_1
        //
        
        // H2_1 by He_1
        
        IJth(J, 1, 6 ) = 0;
        
        scale2 = data->scale[1 - 1];
        scale1 = data->scale[6 - 1];
        IJth(J, 1, 6) /= scale2/scale1;

        

        
        // H2_2 by He_1
        
        IJth(J, 2, 6 ) = 0;
        
        scale2 = data->scale[2 - 1];
        scale1 = data->scale[6 - 1];
        IJth(J, 2, 6) /= scale2/scale1;

        

        
        // H_1 by He_1
        
        IJth(J, 3, 6 ) = 0;
        
        scale2 = data->scale[3 - 1];
        scale1 = data->scale[6 - 1];
        IJth(J, 3, 6) /= scale2/scale1;

        

        
        // H_2 by He_1
        
        IJth(J, 4, 6 ) = 0;
        
        scale2 = data->scale[4 - 1];
        scale1 = data->scale[6 - 1];
        IJth(J, 4, 6) /= scale2/scale1;

        

        
        // H_m0 by He_1
        
        IJth(J, 5, 6 ) = 0;
        
        scale2 = data->scale[5 - 1];
        scale1 = data->scale[6 - 1];
        IJth(J, 5, 6) /= scale2/scale1;

        

        
        // He_1 by He_1
        
        IJth(J, 6, 6 ) = -k03[i]*de;
        
        scale2 = data->scale[6 - 1];
        scale1 = data->scale[6 - 1];
        IJth(J, 6, 6) /= scale2/scale1;

        

        
        // He_2 by He_1
        
        IJth(J, 7, 6 ) = k03[i]*de;
        
        scale2 = data->scale[7 - 1];
        scale1 = data->scale[6 - 1];
        IJth(J, 7, 6) /= scale2/scale1;

        

        
        // He_3 by He_1
        
        IJth(J, 8, 6 ) = 0;
        
        scale2 = data->scale[8 - 1];
        scale1 = data->scale[6 - 1];
        IJth(J, 8, 6) /= scale2/scale1;

        

        
        // de by He_1
        
        IJth(J, 9, 6 ) = k03[i]*de;
        
        scale2 = data->scale[9 - 1];
        scale1 = data->scale[6 - 1];
        IJth(J, 9, 6) /= scale2/scale1;

        

        
        // ge by He_1
        
        IJth(J, 10, 6 ) = 0;
        
        scale2 = data->scale[10 - 1];
        scale1 = data->scale[6 - 1];
        IJth(J, 10, 6) /= scale2/scale1;

        
        IJth(J, 10, 6 ) /= mdensity;
        

        
        
        //
        // Species: He_2
        //
        
        // H2_1 by He_2
        
        IJth(J, 1, 7 ) = 0;
        
        scale2 = data->scale[1 - 1];
        scale1 = data->scale[7 - 1];
        IJth(J, 1, 7) /= scale2/scale1;

        

        
        // H2_2 by He_2
        
        IJth(J, 2, 7 ) = 0;
        
        scale2 = data->scale[2 - 1];
        scale1 = data->scale[7 - 1];
        IJth(J, 2, 7) /= scale2/scale1;

        

        
        // H_1 by He_2
        
        IJth(J, 3, 7 ) = 0;
        
        scale2 = data->scale[3 - 1];
        scale1 = data->scale[7 - 1];
        IJth(J, 3, 7) /= scale2/scale1;

        

        
        // H_2 by He_2
        
        IJth(J, 4, 7 ) = 0;
        
        scale2 = data->scale[4 - 1];
        scale1 = data->scale[7 - 1];
        IJth(J, 4, 7) /= scale2/scale1;

        

        
        // H_m0 by He_2
        
        IJth(J, 5, 7 ) = 0;
        
        scale2 = data->scale[5 - 1];
        scale1 = data->scale[7 - 1];
        IJth(J, 5, 7) /= scale2/scale1;

        

        
        // He_1 by He_2
        
        IJth(J, 6, 7 ) = k04[i]*de;
        
        scale2 = data->scale[6 - 1];
        scale1 = data->scale[7 - 1];
        IJth(J, 6, 7) /= scale2/scale1;

        

        
        // He_2 by He_2
        
        IJth(J, 7, 7 ) = -k04[i]*de - k05[i]*de;
        
        scale2 = data->scale[7 - 1];
        scale1 = data->scale[7 - 1];
        IJth(J, 7, 7) /= scale2/scale1;

        

        
        // He_3 by He_2
        
        IJth(J, 8, 7 ) = k05[i]*de;
        
        scale2 = data->scale[8 - 1];
        scale1 = data->scale[7 - 1];
        IJth(J, 8, 7) /= scale2/scale1;

        

        
        // de by He_2
        
        IJth(J, 9, 7 ) = -k04[i]*de + k05[i]*de;
        
        scale2 = data->scale[9 - 1];
        scale1 = data->scale[7 - 1];
        IJth(J, 9, 7) /= scale2/scale1;

        

        
        // ge by He_2
        
        IJth(J, 10, 7 ) = 0;
        
        scale2 = data->scale[10 - 1];
        scale1 = data->scale[7 - 1];
        IJth(J, 10, 7) /= scale2/scale1;

        
        IJth(J, 10, 7 ) /= mdensity;
        

        
        
        //
        // Species: He_3
        //
        
        // H2_1 by He_3
        
        IJth(J, 1, 8 ) = 0;
        
        scale2 = data->scale[1 - 1];
        scale1 = data->scale[8 - 1];
        IJth(J, 1, 8) /= scale2/scale1;

        

        
        // H2_2 by He_3
        
        IJth(J, 2, 8 ) = 0;
        
        scale2 = data->scale[2 - 1];
        scale1 = data->scale[8 - 1];
        IJth(J, 2, 8) /= scale2/scale1;

        

        
        // H_1 by He_3
        
        IJth(J, 3, 8 ) = 0;
        
        scale2 = data->scale[3 - 1];
        scale1 = data->scale[8 - 1];
        IJth(J, 3, 8) /= scale2/scale1;

        

        
        // H_2 by He_3
        
        IJth(J, 4, 8 ) = 0;
        
        scale2 = data->scale[4 - 1];
        scale1 = data->scale[8 - 1];
        IJth(J, 4, 8) /= scale2/scale1;

        

        
        // H_m0 by He_3
        
        IJth(J, 5, 8 ) = 0;
        
        scale2 = data->scale[5 - 1];
        scale1 = data->scale[8 - 1];
        IJth(J, 5, 8) /= scale2/scale1;

        

        
        // He_1 by He_3
        
        IJth(J, 6, 8 ) = 0;
        
        scale2 = data->scale[6 - 1];
        scale1 = data->scale[8 - 1];
        IJth(J, 6, 8) /= scale2/scale1;

        

        
        // He_2 by He_3
        
        IJth(J, 7, 8 ) = k06[i]*de;
        
        scale2 = data->scale[7 - 1];
        scale1 = data->scale[8 - 1];
        IJth(J, 7, 8) /= scale2/scale1;

        

        
        // He_3 by He_3
        
        IJth(J, 8, 8 ) = -k06[i]*de;
        
        scale2 = data->scale[8 - 1];
        scale1 = data->scale[8 - 1];
        IJth(J, 8, 8) /= scale2/scale1;

        

        
        // de by He_3
        
        IJth(J, 9, 8 ) = -k06[i]*de;
        
        scale2 = data->scale[9 - 1];
        scale1 = data->scale[8 - 1];
        IJth(J, 9, 8) /= scale2/scale1;

        

        
        // ge by He_3
        
        IJth(J, 10, 8 ) = 0;
        
        scale2 = data->scale[10 - 1];
        scale1 = data->scale[8 - 1];
        IJth(J, 10, 8) /= scale2/scale1;

        
        IJth(J, 10, 8 ) /= mdensity;
        

        
        
        //
        // Species: de
        //
        
        // H2_1 by de
        
        IJth(J, 1, 9 ) = -k12[i]*H2_1;
        
        scale2 = data->scale[1 - 1];
        scale1 = data->scale[9 - 1];
        IJth(J, 1, 9) /= scale2/scale1;

        

        
        // H2_2 by de
        
        IJth(J, 2, 9 ) = -k18[i]*H2_2;
        
        scale2 = data->scale[2 - 1];
        scale1 = data->scale[9 - 1];
        IJth(J, 2, 9) /= scale2/scale1;

        

        
        // H_1 by de
        
        IJth(J, 3, 9 ) = -k01[i]*H_1 + k02[i]*H_2 - k07[i]*H_1 + 2*k12[i]*H2_1 + k14[i]*H_m0 + 2*k18[i]*H2_2;
        
        scale2 = data->scale[3 - 1];
        scale1 = data->scale[9 - 1];
        IJth(J, 3, 9) /= scale2/scale1;

        

        
        // H_2 by de
        
        IJth(J, 4, 9 ) = k01[i]*H_1 - k02[i]*H_2;
        
        scale2 = data->scale[4 - 1];
        scale1 = data->scale[9 - 1];
        IJth(J, 4, 9) /= scale2/scale1;

        

        
        // H_m0 by de
        
        IJth(J, 5, 9 ) = k07[i]*H_1 - k14[i]*H_m0;
        
        scale2 = data->scale[5 - 1];
        scale1 = data->scale[9 - 1];
        IJth(J, 5, 9) /= scale2/scale1;

        

        
        // He_1 by de
        
        IJth(J, 6, 9 ) = -k03[i]*He_1 + k04[i]*He_2;
        
        scale2 = data->scale[6 - 1];
        scale1 = data->scale[9 - 1];
        IJth(J, 6, 9) /= scale2/scale1;

        

        
        // He_2 by de
        
        IJth(J, 7, 9 ) = k03[i]*He_1 - k04[i]*He_2 - k05[i]*He_2 + k06[i]*He_3;
        
        scale2 = data->scale[7 - 1];
        scale1 = data->scale[9 - 1];
        IJth(J, 7, 9) /= scale2/scale1;

        

        
        // He_3 by de
        
        IJth(J, 8, 9 ) = k05[i]*He_2 - k06[i]*He_3;
        
        scale2 = data->scale[8 - 1];
        scale1 = data->scale[9 - 1];
        IJth(J, 8, 9) /= scale2/scale1;

        

        
        // de by de
        
        IJth(J, 9, 9 ) = k01[i]*H_1 - k02[i]*H_2 + k03[i]*He_1 - k04[i]*He_2 + k05[i]*He_2 - k06[i]*He_3 - k07[i]*H_1 + k14[i]*H_m0 - k18[i]*H2_2;
        
        scale2 = data->scale[9 - 1];
        scale1 = data->scale[9 - 1];
        IJth(J, 9, 9) /= scale2/scale1;

        

        
        // ge by de
        
        IJth(J, 10, 9 ) = 0;
        
        scale2 = data->scale[10 - 1];
        scale1 = data->scale[9 - 1];
        IJth(J, 10, 9) /= scale2/scale1;

        
        IJth(J, 10, 9 ) /= mdensity;
        

        
        
        //
        // Species: ge
        //
        
        // H2_1 by ge
        
        IJth(J, 1, 10 ) = 0;
        
        scale2 = data->scale[1 - 1];
        scale1 = data->scale[10 - 1];
        IJth(J, 1, 10) /= scale2/scale1;

        

        
        IJth(J, 1, 10 ) *= Tge[i];
        
        // H2_2 by ge
        
        IJth(J, 2, 10 ) = 0;
        
        scale2 = data->scale[2 - 1];
        scale1 = data->scale[10 - 1];
        IJth(J, 2, 10) /= scale2/scale1;

        

        
        IJth(J, 2, 10 ) *= Tge[i];
        
        // H_1 by ge
        
        IJth(J, 3, 10 ) = 0;
        
        scale2 = data->scale[3 - 1];
        scale1 = data->scale[10 - 1];
        IJth(J, 3, 10) /= scale2/scale1;

        

        
        IJth(J, 3, 10 ) *= Tge[i];
        
        // H_2 by ge
        
        IJth(J, 4, 10 ) = 0;
        
        scale2 = data->scale[4 - 1];
        scale1 = data->scale[10 - 1];
        IJth(J, 4, 10) /= scale2/scale1;

        

        
        IJth(J, 4, 10 ) *= Tge[i];
        
        // H_m0 by ge
        
        IJth(J, 5, 10 ) = 0;
        
        scale2 = data->scale[5 - 1];
        scale1 = data->scale[10 - 1];
        IJth(J, 5, 10) /= scale2/scale1;

        

        
        IJth(J, 5, 10 ) *= Tge[i];
        
        // He_1 by ge
        
        IJth(J, 6, 10 ) = 0;
        
        scale2 = data->scale[6 - 1];
        scale1 = data->scale[10 - 1];
        IJth(J, 6, 10) /= scale2/scale1;

        

        
        IJth(J, 6, 10 ) *= Tge[i];
        
        // He_2 by ge
        
        IJth(J, 7, 10 ) = 0;
        
        scale2 = data->scale[7 - 1];
        scale1 = data->scale[10 - 1];
        IJth(J, 7, 10) /= scale2/scale1;

        

        
        IJth(J, 7, 10 ) *= Tge[i];
        
        // He_3 by ge
        
        IJth(J, 8, 10 ) = 0;
        
        scale2 = data->scale[8 - 1];
        scale1 = data->scale[10 - 1];
        IJth(J, 8, 10) /= scale2/scale1;

        

        
        IJth(J, 8, 10 ) *= Tge[i];
        
        // de by ge
        
        IJth(J, 9, 10 ) = 0;
        
        scale2 = data->scale[9 - 1];
        scale1 = data->scale[10 - 1];
        IJth(J, 9, 10) /= scale2/scale1;

        

        
        IJth(J, 9, 10 ) *= Tge[i];
        
        // ge by ge
        
        IJth(J, 10, 10 ) = 0;
        
        scale2 = data->scale[10 - 1];
        scale1 = data->scale[10 - 1];
        IJth(J, 10, 10) /= scale2/scale1;

        
        IJth(J, 10, 10 ) /= mdensity;
        

        
        IJth(J, 10, 10 ) *= Tge[i];
        
        
    }
    return 0;
}






int calculate_rhs_cvdls_9species(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
    cvdls_9species_data *data = (cvdls_9species_data* ) user_data;
    int i, j;

    int nchem = 10;
    int nstrip = 1;
    
    /* change N_Vector back to an array */
    double y_arr[ 10 ];
    /* the variable is ALREADY scaled in "calculate temperature" */
    y_arr[0] = Ith(y , 1);
    // fprintf(stderr, "scale: %.3g \n", data->scale[0]);
    /* the variable is ALREADY scaled in "calculate temperature" */
    y_arr[1] = Ith(y , 2);
    // fprintf(stderr, "scale: %.3g \n", data->scale[1]);
    /* the variable is ALREADY scaled in "calculate temperature" */
    y_arr[2] = Ith(y , 3);
    // fprintf(stderr, "scale: %.3g \n", data->scale[2]);
    /* the variable is ALREADY scaled in "calculate temperature" */
    y_arr[3] = Ith(y , 4);
    // fprintf(stderr, "scale: %.3g \n", data->scale[3]);
    /* the variable is ALREADY scaled in "calculate temperature" */
    y_arr[4] = Ith(y , 5);
    // fprintf(stderr, "scale: %.3g \n", data->scale[4]);
    /* the variable is ALREADY scaled in "calculate temperature" */
    y_arr[5] = Ith(y , 6);
    // fprintf(stderr, "scale: %.3g \n", data->scale[5]);
    /* the variable is ALREADY scaled in "calculate temperature" */
    y_arr[6] = Ith(y , 7);
    // fprintf(stderr, "scale: %.3g \n", data->scale[6]);
    /* the variable is ALREADY scaled in "calculate temperature" */
    y_arr[7] = Ith(y , 8);
    // fprintf(stderr, "scale: %.3g \n", data->scale[7]);
    /* the variable is ALREADY scaled in "calculate temperature" */
    y_arr[8] = Ith(y , 9);
    // fprintf(stderr, "scale: %.3g \n", data->scale[8]);
    /* the variable is ALREADY scaled in "calculate temperature" */
    y_arr[9] = Ith(y , 10);
    // fprintf(stderr, "scale: %.3g \n", data->scale[9]);

    cvdls_9species_calculate_temperature(data, y_arr , nstrip, nchem );
    cvdls_9species_interpolate_rates(data, nstrip);


    /* Now we set up some temporaries */
    double *k01 = data->rs_k01;
    double *k02 = data->rs_k02;
    double *k03 = data->rs_k03;
    double *k04 = data->rs_k04;
    double *k05 = data->rs_k05;
    double *k06 = data->rs_k06;
    double *k07 = data->rs_k07;
    double *k08 = data->rs_k08;
    double *k09 = data->rs_k09;
    double *k10 = data->rs_k10;
    double *k11 = data->rs_k11;
    double *k12 = data->rs_k12;
    double *k13 = data->rs_k13;
    double *k14 = data->rs_k14;
    double *k15 = data->rs_k15;
    double *k16 = data->rs_k16;
    double *k17 = data->rs_k17;
    double *k18 = data->rs_k18;
    double *k19 = data->rs_k19;
    double *k21 = data->rs_k21;
    double *k22 = data->rs_k22;
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

    double mh = 1.67e-24;
    double mdensity;
    
    // i = nstrip;
    i = 0;
    T = data->Ts[i];
    z = data->current_z;
   
    // fprintf(stderr, "T from rhs: %0.5g\n", T);

    double scale;
    int jj =0;
    scale = data->scale[jj];
    H2_1 = Ith( y,1 )*scale;
    jj++;
    
    mdensity += H2_1 * 2.01588;
    


    
    scale = data->scale[jj];
    H2_2 = Ith( y,2 )*scale;
    jj++;
    
    mdensity += H2_2 * 2.01588;
    


    
    scale = data->scale[jj];
    H_1 = Ith( y,3 )*scale;
    jj++;
    
    mdensity += H_1 * 1.00794;
    


    
    scale = data->scale[jj];
    H_2 = Ith( y,4 )*scale;
    jj++;
    
    mdensity += H_2 * 1.00794;
    


    
    scale = data->scale[jj];
    H_m0 = Ith( y,5 )*scale;
    jj++;
    
    mdensity += H_m0 * 1.00794;
    


    
    scale = data->scale[jj];
    He_1 = Ith( y,6 )*scale;
    jj++;
    
    mdensity += He_1 * 4.002602;
    


    
    scale = data->scale[jj];
    He_2 = Ith( y,7 )*scale;
    jj++;
    
    mdensity += He_2 * 4.002602;
    


    
    scale = data->scale[jj];
    He_3 = Ith( y,8 )*scale;
    jj++;
    
    mdensity += He_3 * 4.002602;
    


    
    scale = data->scale[jj];
    de = Ith( y,9 )*scale;
    jj++;
    


    
    scale = data->scale[jj];
    ge = Ith( y,10 )*scale;
    jj++;
    


    
    double nH; 
    nH = H_1 + H_2 + 2.0*H2_1 + 2.0*H2_2;
    
    
    mdensity *= mh;
    //
    // Species: H2_1
    //
    Ith(ydot, 1) = k08[i]*H_1*H_m0 + k10[i]*H2_2*H_1 - k11[i]*H2_1*H_2 - k12[i]*H2_1*de - k13[i]*H2_1*H_1 + k19[i]*H2_2*H_m0 + k21[i]*H2_1*pow(H_1, 2) + k22[i]*pow(H_1, 3);
 
    scale = data->scale[1 - 1];
    Ith(ydot, 1) /= scale;

    
    
    //
    // Species: H2_2
    //
    Ith(ydot, 2) = k09[i]*H_1*H_2 - k10[i]*H2_2*H_1 + k11[i]*H2_1*H_2 + k17[i]*H_2*H_m0 - k18[i]*H2_2*de - k19[i]*H2_2*H_m0;
 
    scale = data->scale[2 - 1];
    Ith(ydot, 2) /= scale;

    
    
    //
    // Species: H_1
    //
    Ith(ydot, 3) = -k01[i]*H_1*de + k02[i]*H_2*de - k07[i]*H_1*de - k08[i]*H_1*H_m0 - k09[i]*H_1*H_2 - k10[i]*H2_2*H_1 + k11[i]*H2_1*H_2 + 2*k12[i]*H2_1*de + 2*k13[i]*H2_1*H_1 + k14[i]*H_m0*de + k15[i]*H_1*H_m0 + 2*k16[i]*H_2*H_m0 + 2*k18[i]*H2_2*de + k19[i]*H2_2*H_m0 - 2*k21[i]*H2_1*pow(H_1, 2) - 2*k22[i]*pow(H_1, 3);
 
    scale = data->scale[3 - 1];
    Ith(ydot, 3) /= scale;

    
    
    //
    // Species: H_2
    //
    Ith(ydot, 4) = k01[i]*H_1*de - k02[i]*H_2*de - k09[i]*H_1*H_2 + k10[i]*H2_2*H_1 - k11[i]*H2_1*H_2 - k16[i]*H_2*H_m0 - k17[i]*H_2*H_m0;
 
    scale = data->scale[4 - 1];
    Ith(ydot, 4) /= scale;

    
    
    //
    // Species: H_m0
    //
    Ith(ydot, 5) = k07[i]*H_1*de - k08[i]*H_1*H_m0 - k14[i]*H_m0*de - k15[i]*H_1*H_m0 - k16[i]*H_2*H_m0 - k17[i]*H_2*H_m0 - k19[i]*H2_2*H_m0;
 
    scale = data->scale[5 - 1];
    Ith(ydot, 5) /= scale;

    
    
    //
    // Species: He_1
    //
    Ith(ydot, 6) = -k03[i]*He_1*de + k04[i]*He_2*de;
 
    scale = data->scale[6 - 1];
    Ith(ydot, 6) /= scale;

    
    
    //
    // Species: He_2
    //
    Ith(ydot, 7) = k03[i]*He_1*de - k04[i]*He_2*de - k05[i]*He_2*de + k06[i]*He_3*de;
 
    scale = data->scale[7 - 1];
    Ith(ydot, 7) /= scale;

    
    
    //
    // Species: He_3
    //
    Ith(ydot, 8) = k05[i]*He_2*de - k06[i]*He_3*de;
 
    scale = data->scale[8 - 1];
    Ith(ydot, 8) /= scale;

    
    
    //
    // Species: de
    //
    Ith(ydot, 9) = k01[i]*H_1*de - k02[i]*H_2*de + k03[i]*He_1*de - k04[i]*He_2*de + k05[i]*He_2*de - k06[i]*He_3*de - k07[i]*H_1*de + k08[i]*H_1*H_m0 + k14[i]*H_m0*de + k15[i]*H_1*H_m0 + k17[i]*H_2*H_m0 - k18[i]*H2_2*de;
 
    scale = data->scale[9 - 1];
    Ith(ydot, 9) /= scale;

    
    
    //
    // Species: ge
    //
    Ith(ydot, 10) = 0;
 
    scale = data->scale[10 - 1];
    Ith(ydot, 10) /= scale;

    
    Ith(ydot, 10) /= mdensity;
    
    
    // fprintf(stderr, "k22: %0.5g, T: %0.5g \n", k22[i], T);
    return 0;
    }
