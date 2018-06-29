
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


#include "eq_H2_2_solver.h"

void *setup_cvode_solver( rhs_f f, jac_f Jac, int NEQ,
                    eq_H2_2_data *data , SUNLinearSolver LS, SUNMatrix A, N_Vector , double, N_Vector);
int cvode_solver( void *cvode_mem, double *output, int NEQ, double *dt, eq_H2_2_data *, N_Vector, double, N_Vector);

eq_H2_2_data *eq_H2_2_setup_data(
    int *NumberOfFields, char ***FieldNames)
{
    int i;

    eq_H2_2_data *data = (eq_H2_2_data *) malloc(sizeof(eq_H2_2_data));
    
    /* allocate space for the scale related pieces */
    for (i = 0; i< 10 ; i++){
        data->scale[i] = 1.0;
        data->inv_scale[i] = 1.0;
    }
    
    /*initialize temperature so it wont crash*/
    data->Ts[0] = 1000.0;
    data->logTs[0] = log(1000.0);

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
    
    eq_H2_2_read_rate_tables(data);
    fprintf(stderr, "Successfully read in rate tables.\n");

    eq_H2_2_read_cooling_tables(data);
    fprintf(stderr, "Successfully read in cooling rate tables.\n");
    
    eq_H2_2_read_gamma(data);
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


int eq_H2_2_main(int argc, char** argv)
{
    eq_H2_2_data *data = eq_H2_2_setup_data(NULL, NULL);

    /* Initial conditions */

    hid_t file_id = H5Fopen("eq_H2_2_initial_conditions.h5", H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file_id < 0) {fprintf(stderr, "Failed to open "
        "eq_H2_2_initial_conditions.h5 so dying.\n");
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
    

    H5Fclose(file_id);

    double dtf = 299204917.32712233;
    double dt = -1.0;
    double z = -1.0;
    for (i = 0; i < dims * N; i++) input[i] = ics[i];
    double ttot;
    ttot = dengo_evolve_eq_H2_2(dtf, dt, z, input, rtol, atol, dims, data, temp);

    /* Write results to HDF5 file */
    file_id = H5Fcreate("eq_H2_2_solution.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
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
    
    free(tics);
    free(ics);
    free(input);
    free(data);
    free(atol);
    free(rtol);
    return 0;
}
 



double dengo_evolve_eq_H2_2 (double dtf, double &dt, double z, double *input,
            double *rtol, double *atol, long long dims, eq_H2_2_data *data, double *temp_array ) {
    int i, j;
    hid_t file_id;
    /* fprintf(stderr, "  ncells = % 3i\n", (int) dims); */

    int N = 10;
    for (i = 0; i<dims; i++) {
      j = i * N;
        
          input[j] /= 2.01588 ;
          atol[j] /= 2.01588 ;
        
        j++;
      
        
          input[j] /= 2.01588 ;
          atol[j] /= 2.01588 ;
        
        j++;
      
        
          input[j] /= 1.00794 ;
          atol[j] /= 1.00794 ;
        
        j++;
      
        
          input[j] /= 1.00794 ;
          atol[j] /= 1.00794 ;
        
        j++;
      
        
          input[j] /= 1.00794 ;
          atol[j] /= 1.00794 ;
        
        j++;
      
        
          input[j] /= 4.002602 ;
          atol[j] /= 4.002602 ;
        
        j++;
      
        
          input[j] /= 4.002602 ;
          atol[j] /= 4.002602 ;
        
        j++;
      
        
          input[j] /= 4.002602 ;
          atol[j] /= 4.002602 ;
        
        j++;
      
        
          input[j] /= 1.0 ;
          atol[j] /= 1.0 ;
        
        j++;
      
        
        j++;
      
    }
    //ensure_electron_consistency(input, dims, N);

    rhs_f f = calculate_rhs_eq_H2_2;
    jac_f jf = calculate_jacobian_eq_H2_2;
    if (dt < 0) dt = dtf / 1e0;
    data->current_z = z;
    int niter = 0;
    int siter = 0;
    double ttot = 0;
    
    double *ttot_all = (double *) malloc( dims * sizeof(double) );

  
    
    double floor_value = 1e-25;

    // Initialize a CVODE object, memory spaces
    // and attach rhs, jac to them
    int flag;
    double reltol = 1.0e-3;
    void *cvode_mem;
    int MAX_ITERATION = 100; 
    double y[10];
    
    SUNLinearSolver LS;
    SUNMatrix A;
    N_Vector y_vec, abstol;
    
    y_vec  = NULL;   
    LS     = NULL;
    A      = NULL;
    abstol = NULL;
    
    y_vec  = N_VNew_Serial(N);
    abstol = N_VNew_Serial(N);


    for (i=0; i<N; i++) {
        y[i] = 1.0;
        NV_Ith_S(y_vec , i )   = 1.0;
        NV_Ith_S(abstol, i )   = 1.0;

    }

    A = SUNDenseMatrix(N, N);
    LS = SUNDenseLinearSolver(y_vec, A);
    cvode_mem = setup_cvode_solver( f, jf, N, data, LS, A, y_vec, reltol, abstol);
    
    double h_density  = 1.0e14;
    double he_density = 3.0e14;
    double intermediate_solution;
    
    double q_density = input[8];

    for (int d = 0; d < dims; d++){
        


        // copy array which can be passed to the solver
        for (i = 0; i < N; i++){ 
            // this is being passed around 
            // passively by the "dengo_rate_data" 
            // will have to fix it for openmp
            //
            //
            //
            //
            /*
                 data->scale[i]     = input[d*N + i];
                data->inv_scale[i] = 1.0 / data->scale[i];
                NV_Ith_S(y_vec , i )   = input[ d*N+i ] * data->inv_scale[i];
                NV_Ith_S(abstol, i )   = reltol ;
                */
            
            
            if ( i == 0  || i == 2){
                data->scale[i] = h_density;
                data->inv_scale[i] = 1.0 / data->scale[i];
            
            NV_Ith_S(y_vec , i )   = input[ d*N+i ] * data->inv_scale[i];
            NV_Ith_S(abstol, i )   = NV_Ith_S(y_vec , i ) * reltol * 1.0e-4;
            
            fprintf(stderr, "y_vec = %0.5g; abstol = %0.5g \n", NV_Ith_S(y_vec , i ), NV_Ith_S(abstol , i ));  

            }
            else{
                 data->scale[i]     = input[d*N + i];
                data->inv_scale[i] = 1.0 / data->scale[i];
                NV_Ith_S(y_vec , i )   = input[ d*N+i ] * data->inv_scale[i];
                NV_Ith_S(abstol, i )   = reltol ;
            }
        }
        /*
            if ( i == 5 ){
                data->scale[i] = he_density;
                data->inv_scale[i] = 1.0 / data->scale[i];
            
                NV_Ith_S(y_vec , i )   = input[ d*N+i ] * data->inv_scale[i];
                NV_Ith_S(abstol, i )   = reltol ;

            
            }
            if ( i == 1 || i == 3 || i == 4 || i == 6 || i == 7 || i == 8 ){
                data->scale[i] = q_density;
                data->inv_scale[i] = 1.0 / data->scale[i];
            
                NV_Ith_S(y_vec , i )   = input[ d*N+i ] * data->inv_scale[i];
                NV_Ith_S(abstol, i )   = NV_Ith_S(y_vec , i ) * reltol ;
                
                
                fprintf(stderr, "input[%d] : %0.5g \n", i, input[d*N + i]);
                fprintf(stderr, "y_vec [%d]: %0.5g \n", i, NV_Ith_S(y_vec , i));
                fprintf(stderr, "abstol[%d]: %0.5g \n", i, NV_Ith_S(abstol, i));
                
 
            }


            if ( i  > 8) {
                data->scale[i]     = input[d*N + i];
                data->inv_scale[i] = 1.0 / data->scale[i];
                NV_Ith_S(y_vec , i )   = input[ d*N+i ] * data->inv_scale[i];
                NV_Ith_S(abstol, i )   = reltol ;

           }
            } 

            // data->scale[i]     = input[d*N + i];
            // data->inv_scale[i] = 1.0 / data->scale[i];
           

            

        }*/
        
        // initialize a dt for the solver    
        dt = dtf;
        ttot = 0.0;
        siter = 0;
            
        while (ttot < dtf) { 
            fprintf(stderr, "%d th strip: %d iterations, time: %0.5g\n", d, siter, ttot );
            
            //
            
            // we are not returning the CVode flag
            // 0: success
            // 1: failure
            flag = cvode_solver( cvode_mem, y, N, &dt, data, y_vec, reltol, abstol);
            for (i = 0; i < N; i++) {
                if (y[i] < 0) {
                    flag = 1;
                    fprintf(stderr, "failed, negative \n ");
                    break;
                }
            }

            if (flag < 1){

                for (i = 0; i < N; i++){

                    
                    intermediate_solution = y[i] * data->scale[i];
fprintf(stderr, "y[%d]: %0.5g \n",i, y[i]);

                    data->scale[i]     = intermediate_solution;
                    data->inv_scale[i] = 1.0 / data->scale[i];
                    

                    // fprintf(stderr, "y_vec [%d]: %0.5g \n", i, NV_Ith_S(y_vec , i));
                    NV_Ith_S(y_vec , i )   = intermediate_solution * data->inv_scale[i];
                    NV_Ith_S(abstol, i )   = NV_Ith_S(y_vec, i) * reltol ;
                

                    //fprintf(stderr, "abstol[%d]: %0.5g \n", i, NV_Ith_S(abstol, i));
  
                
                }
                
                ttot += dt;

            } else{
                fprintf(stderr, "real fail\n");
                dt /= 2.0;
            }

	        dt = DMIN(dt * 1.1, dtf - ttot);

            if (siter == MAX_ITERATION) break;

        
            siter++;
        } // while loop for each strip
        
        fprintf(stderr, "%d the strip = %0.5g\n", d, ttot);
        temp_array[ d ] = data->Ts[0];
        ttot_all[d] = ttot;

        for (i = 0; i < N; i++){
            input[d*N +i] = data->scale[i] ;
            
        } // copy data back to the input array
    } // for d dims loop
    CVodeFree(&cvode_mem);
    SUNLinSolFree(LS);
    SUNMatDestroy(A);
    N_VDestroy(y_vec);
    N_VDestroy(abstol);

    for (i = 0; i < dims; i++) {
      j = i * N;
        
        input[j] *= 2.01588;
        atol[j] *= 2.01588;
        
        j++;
      
        
        input[j] *= 2.01588;
        atol[j] *= 2.01588;
        
        j++;
      
        
        input[j] *= 1.00794;
        atol[j] *= 1.00794;
        
        j++;
      
        
        input[j] *= 1.00794;
        atol[j] *= 1.00794;
        
        j++;
      
        
        input[j] *= 1.00794;
        atol[j] *= 1.00794;
        
        j++;
      
        
        input[j] *= 4.002602;
        atol[j] *= 4.002602;
        
        j++;
      
        
        input[j] *= 4.002602;
        atol[j] *= 4.002602;
        
        j++;
      
        
        input[j] *= 4.002602;
        atol[j] *= 4.002602;
        
        j++;
      
        
        input[j] *= 1.0;
        atol[j] *= 1.0;
        
        j++;
      
        
        j++;
      
    }

    double dt_final = dtf;
    
    for (int d = 0; d < dims; d++){
        if (ttot_all[d] < dt_final) dt_final = ttot_all[d];    
    }
    free(ttot_all);

    return dt_final;
}
 


void eq_H2_2_read_rate_tables(eq_H2_2_data *data)
{
    hid_t file_id = H5Fopen("eq_H2_2_tables.h5", H5F_ACC_RDONLY, H5P_DEFAULT);
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

void eq_H2_2_read_cooling_tables(eq_H2_2_data *data)
{

    hid_t file_id = H5Fopen("eq_H2_2_tables.h5", H5F_ACC_RDONLY, H5P_DEFAULT);
    /* Allocate the correct number of rate tables */
    H5LTread_dataset_double(file_id, "/brem_brem",
                            data->c_brem_brem);
    H5LTread_dataset_double(file_id, "/ceHeI_ceHeI",
                            data->c_ceHeI_ceHeI);
    H5LTread_dataset_double(file_id, "/ceHeII_ceHeII",
                            data->c_ceHeII_ceHeII);
    H5LTread_dataset_double(file_id, "/ceHI_ceHI",
                            data->c_ceHI_ceHI);
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

void eq_H2_2_read_gamma(eq_H2_2_data *data)
{

    hid_t file_id = H5Fopen("eq_H2_2_tables.h5", H5F_ACC_RDONLY, H5P_DEFAULT);
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
 


void eq_H2_2_calculate_temperature(eq_H2_2_data *data,
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
    double scale = 1.0;
    
    i = 0;
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
    
        density = 2.01588*H2_1 + 2.01588*H2_2 + 1.00794*H_1 + 1.00794*H_2 + 1.00794*H_m0 + 4.002602*He_1 + 4.002602*He_2 + 4.002602*He_3;
        
    
        
        
        // Initiate the "guess" temperature
        T = data->Ts[i];

        // fprintf(stderr, "Temp: %0.5g\n", T);

        Tnew = T + 1.0;
        double dge_dT;
        double dge;

        

        double gammaH2_1;
        double dgammaH2_1_dT;
        

        double gammaH2_2;
        double dgammaH2_2_dT;
        
        
        double Tdiff = 1.0; 
        while ( Tdiff > 0.1 ){
        // We do Newton's Iteration to calculate the temperature
        // Since gammaH2 is dependent on the temperature too!

        T = data->Ts[i];
        
        eq_H2_2_interpolate_gamma(data, i);
        
        gammaH2_1 = data->gammaH2_1[i];
        dgammaH2_1_dT = data->dgammaH2_1_dT[i];
        // fprintf(stderr, ":gammaH2_1 %0.5g , dgammaH2_1_dT: %.5g \n", gammaH2_1, dgammaH2_1_dT  );
        
        gammaH2_2 = data->gammaH2_2[i];
        dgammaH2_2_dT = data->dgammaH2_2_dT[i];
        // fprintf(stderr, ":gammaH2_2 %0.5g , dgammaH2_2_dT: %.5g \n", gammaH2_2, dgammaH2_2_dT  );
        
       
        
        // update gammaH2

        

        // The derivatives of  sum (nkT/(gamma - 1)/mh/density) - ge
        // This is the function we want to minimize
        // which should only be dependent on the first part
        dge_dT = T*kb*(-H2_1*dgammaH2_1_dT/pow(gammaH2_1 - 1.0, 2) - H2_2*dgammaH2_2_dT/pow(gammaH2_2 - 1.0, 2))/(density*mh) + kb*(H2_1/(gammaH2_1 - 1.0) + H2_2/(gammaH2_2 - 1.0) + H_1/(gamma - 1.0) + H_2/(gamma - 1.0) + H_m0/(gamma - 1.0) + He_1/(gamma - 1.0) + He_2/(gamma - 1.0) + He_3/(gamma - 1.0) + de/(gamma - 1.0))/(density*mh);
        
        //This is the change in ge for each iteration
        dge = T*kb*(H2_1/(gammaH2_1 - 1.0) + H2_2/(gammaH2_2 - 1.0) + H_1/(gamma - 1.0) + H_2/(gamma - 1.0) + H_m0/(gamma - 1.0) + He_1/(gamma - 1.0) + He_2/(gamma - 1.0) + He_3/(gamma - 1.0) + de/(gamma - 1.0))/(density*mh) - ge;

        Tnew = T - dge/dge_dT;
        data->Ts[i] = Tnew;
        
        Tdiff = fabs(T - Tnew);
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
        /*fprintf(stderr, "T[%d] = % 0.16g, density = % 0.16g\n",
                i, data->Ts[i], density);*/
         
}
 


/*
   This setup may be different than the user may anticipate, as a result
   of the lockstep timestep we use for a pencil beam through the grid.
   As such, it accepts the number of things to interpolate and makes
   assumptions about the sizes of the rates.
*/

/* This also requires no templating other than for the solver name...*/
void eq_H2_2_interpolate_rates(eq_H2_2_data *data,
                    int nstrip)
{
    int i, bin_id, zbin_id;
    double lb, t1, t2;
    double lbz, z1, z2;
    int no_photo = 0;
    lb = log(data->bounds[0]);
    
    i = 0;
    lbz = log(data->z_bounds[0] + 1.0);
    /*fprintf(stderr, "lb = % 0.16g, ub = % 0.16g\n", lb, ub);*/
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
    
        bin_id = data->bin_id[i];
        data->rs_k01[i] = data->r_k01[bin_id] +
            data->Tdef[i] * (data->r_k01[bin_id+1] - data->r_k01[bin_id]);
        data->drs_k01[i] = (data->r_k01[bin_id+1] - data->r_k01[bin_id]);
        data->drs_k01[i] /= data->dT[i];
	data->drs_k01[i] *= data->invTs[i];
    
        data->rs_k02[i] = data->r_k02[bin_id] +
            data->Tdef[i] * (data->r_k02[bin_id+1] - data->r_k02[bin_id]);
        data->drs_k02[i] = (data->r_k02[bin_id+1] - data->r_k02[bin_id]);
        data->drs_k02[i] /= data->dT[i];
	data->drs_k02[i] *= data->invTs[i];
    
        data->rs_k03[i] = data->r_k03[bin_id] +
            data->Tdef[i] * (data->r_k03[bin_id+1] - data->r_k03[bin_id]);
        data->drs_k03[i] = (data->r_k03[bin_id+1] - data->r_k03[bin_id]);
        data->drs_k03[i] /= data->dT[i];
	data->drs_k03[i] *= data->invTs[i];
    
        data->rs_k04[i] = data->r_k04[bin_id] +
            data->Tdef[i] * (data->r_k04[bin_id+1] - data->r_k04[bin_id]);
        data->drs_k04[i] = (data->r_k04[bin_id+1] - data->r_k04[bin_id]);
        data->drs_k04[i] /= data->dT[i];
	data->drs_k04[i] *= data->invTs[i];
    
        data->rs_k05[i] = data->r_k05[bin_id] +
            data->Tdef[i] * (data->r_k05[bin_id+1] - data->r_k05[bin_id]);
        data->drs_k05[i] = (data->r_k05[bin_id+1] - data->r_k05[bin_id]);
        data->drs_k05[i] /= data->dT[i];
	data->drs_k05[i] *= data->invTs[i];
    
        data->rs_k06[i] = data->r_k06[bin_id] +
            data->Tdef[i] * (data->r_k06[bin_id+1] - data->r_k06[bin_id]);
        data->drs_k06[i] = (data->r_k06[bin_id+1] - data->r_k06[bin_id]);
        data->drs_k06[i] /= data->dT[i];
	data->drs_k06[i] *= data->invTs[i];
    
        data->rs_k07[i] = data->r_k07[bin_id] +
            data->Tdef[i] * (data->r_k07[bin_id+1] - data->r_k07[bin_id]);
        data->drs_k07[i] = (data->r_k07[bin_id+1] - data->r_k07[bin_id]);
        data->drs_k07[i] /= data->dT[i];
	data->drs_k07[i] *= data->invTs[i];
    
        data->rs_k08[i] = data->r_k08[bin_id] +
            data->Tdef[i] * (data->r_k08[bin_id+1] - data->r_k08[bin_id]);
        data->drs_k08[i] = (data->r_k08[bin_id+1] - data->r_k08[bin_id]);
        data->drs_k08[i] /= data->dT[i];
	data->drs_k08[i] *= data->invTs[i];
    
        data->rs_k09[i] = data->r_k09[bin_id] +
            data->Tdef[i] * (data->r_k09[bin_id+1] - data->r_k09[bin_id]);
        data->drs_k09[i] = (data->r_k09[bin_id+1] - data->r_k09[bin_id]);
        data->drs_k09[i] /= data->dT[i];
	data->drs_k09[i] *= data->invTs[i];
    
        data->rs_k10[i] = data->r_k10[bin_id] +
            data->Tdef[i] * (data->r_k10[bin_id+1] - data->r_k10[bin_id]);
        data->drs_k10[i] = (data->r_k10[bin_id+1] - data->r_k10[bin_id]);
        data->drs_k10[i] /= data->dT[i];
	data->drs_k10[i] *= data->invTs[i];
    
        data->rs_k11[i] = data->r_k11[bin_id] +
            data->Tdef[i] * (data->r_k11[bin_id+1] - data->r_k11[bin_id]);
        data->drs_k11[i] = (data->r_k11[bin_id+1] - data->r_k11[bin_id]);
        data->drs_k11[i] /= data->dT[i];
	data->drs_k11[i] *= data->invTs[i];
    
        data->rs_k12[i] = data->r_k12[bin_id] +
            data->Tdef[i] * (data->r_k12[bin_id+1] - data->r_k12[bin_id]);
        data->drs_k12[i] = (data->r_k12[bin_id+1] - data->r_k12[bin_id]);
        data->drs_k12[i] /= data->dT[i];
	data->drs_k12[i] *= data->invTs[i];
    
        data->rs_k13[i] = data->r_k13[bin_id] +
            data->Tdef[i] * (data->r_k13[bin_id+1] - data->r_k13[bin_id]);
        data->drs_k13[i] = (data->r_k13[bin_id+1] - data->r_k13[bin_id]);
        data->drs_k13[i] /= data->dT[i];
	data->drs_k13[i] *= data->invTs[i];
    
        data->rs_k14[i] = data->r_k14[bin_id] +
            data->Tdef[i] * (data->r_k14[bin_id+1] - data->r_k14[bin_id]);
        data->drs_k14[i] = (data->r_k14[bin_id+1] - data->r_k14[bin_id]);
        data->drs_k14[i] /= data->dT[i];
	data->drs_k14[i] *= data->invTs[i];
    
        data->rs_k15[i] = data->r_k15[bin_id] +
            data->Tdef[i] * (data->r_k15[bin_id+1] - data->r_k15[bin_id]);
        data->drs_k15[i] = (data->r_k15[bin_id+1] - data->r_k15[bin_id]);
        data->drs_k15[i] /= data->dT[i];
	data->drs_k15[i] *= data->invTs[i];
    
        data->rs_k16[i] = data->r_k16[bin_id] +
            data->Tdef[i] * (data->r_k16[bin_id+1] - data->r_k16[bin_id]);
        data->drs_k16[i] = (data->r_k16[bin_id+1] - data->r_k16[bin_id]);
        data->drs_k16[i] /= data->dT[i];
	data->drs_k16[i] *= data->invTs[i];
    
        data->rs_k17[i] = data->r_k17[bin_id] +
            data->Tdef[i] * (data->r_k17[bin_id+1] - data->r_k17[bin_id]);
        data->drs_k17[i] = (data->r_k17[bin_id+1] - data->r_k17[bin_id]);
        data->drs_k17[i] /= data->dT[i];
	data->drs_k17[i] *= data->invTs[i];
    
        data->rs_k18[i] = data->r_k18[bin_id] +
            data->Tdef[i] * (data->r_k18[bin_id+1] - data->r_k18[bin_id]);
        data->drs_k18[i] = (data->r_k18[bin_id+1] - data->r_k18[bin_id]);
        data->drs_k18[i] /= data->dT[i];
	data->drs_k18[i] *= data->invTs[i];
    
        data->rs_k19[i] = data->r_k19[bin_id] +
            data->Tdef[i] * (data->r_k19[bin_id+1] - data->r_k19[bin_id]);
        data->drs_k19[i] = (data->r_k19[bin_id+1] - data->r_k19[bin_id]);
        data->drs_k19[i] /= data->dT[i];
	data->drs_k19[i] *= data->invTs[i];
    
        data->rs_k21[i] = data->r_k21[bin_id] +
            data->Tdef[i] * (data->r_k21[bin_id+1] - data->r_k21[bin_id]);
        data->drs_k21[i] = (data->r_k21[bin_id+1] - data->r_k21[bin_id]);
        data->drs_k21[i] /= data->dT[i];
	data->drs_k21[i] *= data->invTs[i];
    
        data->rs_k22[i] = data->r_k22[bin_id] +
            data->Tdef[i] * (data->r_k22[bin_id+1] - data->r_k22[bin_id]);
        data->drs_k22[i] = (data->r_k22[bin_id+1] - data->r_k22[bin_id]);
        data->drs_k22[i] /= data->dT[i];
	data->drs_k22[i] *= data->invTs[i];
    
        data->cs_brem_brem[i] = data->c_brem_brem[bin_id] +
            data->Tdef[i] * (data->c_brem_brem[bin_id+1] - data->c_brem_brem[bin_id]);
        data->dcs_brem_brem[i] = (data->c_brem_brem[bin_id+1] - data->c_brem_brem[bin_id]);;
        data->dcs_brem_brem[i] /= data->dT[i];
	data->dcs_brem_brem[i] *= data->invTs[i];
    
        data->cs_ceHeI_ceHeI[i] = data->c_ceHeI_ceHeI[bin_id] +
            data->Tdef[i] * (data->c_ceHeI_ceHeI[bin_id+1] - data->c_ceHeI_ceHeI[bin_id]);
        data->dcs_ceHeI_ceHeI[i] = (data->c_ceHeI_ceHeI[bin_id+1] - data->c_ceHeI_ceHeI[bin_id]);;
        data->dcs_ceHeI_ceHeI[i] /= data->dT[i];
	data->dcs_ceHeI_ceHeI[i] *= data->invTs[i];
    
        data->cs_ceHeII_ceHeII[i] = data->c_ceHeII_ceHeII[bin_id] +
            data->Tdef[i] * (data->c_ceHeII_ceHeII[bin_id+1] - data->c_ceHeII_ceHeII[bin_id]);
        data->dcs_ceHeII_ceHeII[i] = (data->c_ceHeII_ceHeII[bin_id+1] - data->c_ceHeII_ceHeII[bin_id]);;
        data->dcs_ceHeII_ceHeII[i] /= data->dT[i];
	data->dcs_ceHeII_ceHeII[i] *= data->invTs[i];
    
        data->cs_ceHI_ceHI[i] = data->c_ceHI_ceHI[bin_id] +
            data->Tdef[i] * (data->c_ceHI_ceHI[bin_id+1] - data->c_ceHI_ceHI[bin_id]);
        data->dcs_ceHI_ceHI[i] = (data->c_ceHI_ceHI[bin_id+1] - data->c_ceHI_ceHI[bin_id]);;
        data->dcs_ceHI_ceHI[i] /= data->dT[i];
	data->dcs_ceHI_ceHI[i] *= data->invTs[i];
    
        data->cs_ciHeI_ciHeI[i] = data->c_ciHeI_ciHeI[bin_id] +
            data->Tdef[i] * (data->c_ciHeI_ciHeI[bin_id+1] - data->c_ciHeI_ciHeI[bin_id]);
        data->dcs_ciHeI_ciHeI[i] = (data->c_ciHeI_ciHeI[bin_id+1] - data->c_ciHeI_ciHeI[bin_id]);;
        data->dcs_ciHeI_ciHeI[i] /= data->dT[i];
	data->dcs_ciHeI_ciHeI[i] *= data->invTs[i];
    
        data->cs_ciHeII_ciHeII[i] = data->c_ciHeII_ciHeII[bin_id] +
            data->Tdef[i] * (data->c_ciHeII_ciHeII[bin_id+1] - data->c_ciHeII_ciHeII[bin_id]);
        data->dcs_ciHeII_ciHeII[i] = (data->c_ciHeII_ciHeII[bin_id+1] - data->c_ciHeII_ciHeII[bin_id]);;
        data->dcs_ciHeII_ciHeII[i] /= data->dT[i];
	data->dcs_ciHeII_ciHeII[i] *= data->invTs[i];
    
        data->cs_ciHeIS_ciHeIS[i] = data->c_ciHeIS_ciHeIS[bin_id] +
            data->Tdef[i] * (data->c_ciHeIS_ciHeIS[bin_id+1] - data->c_ciHeIS_ciHeIS[bin_id]);
        data->dcs_ciHeIS_ciHeIS[i] = (data->c_ciHeIS_ciHeIS[bin_id+1] - data->c_ciHeIS_ciHeIS[bin_id]);;
        data->dcs_ciHeIS_ciHeIS[i] /= data->dT[i];
	data->dcs_ciHeIS_ciHeIS[i] *= data->invTs[i];
    
        data->cs_ciHI_ciHI[i] = data->c_ciHI_ciHI[bin_id] +
            data->Tdef[i] * (data->c_ciHI_ciHI[bin_id+1] - data->c_ciHI_ciHI[bin_id]);
        data->dcs_ciHI_ciHI[i] = (data->c_ciHI_ciHI[bin_id+1] - data->c_ciHI_ciHI[bin_id]);;
        data->dcs_ciHI_ciHI[i] /= data->dT[i];
	data->dcs_ciHI_ciHI[i] *= data->invTs[i];
    
        data->cs_compton_comp_[i] = data->c_compton_comp_[bin_id] +
            data->Tdef[i] * (data->c_compton_comp_[bin_id+1] - data->c_compton_comp_[bin_id]);
        data->dcs_compton_comp_[i] = (data->c_compton_comp_[bin_id+1] - data->c_compton_comp_[bin_id]);;
        data->dcs_compton_comp_[i] /= data->dT[i];
	data->dcs_compton_comp_[i] *= data->invTs[i];
    
        data->cs_gammah_gammah[i] = data->c_gammah_gammah[bin_id] +
            data->Tdef[i] * (data->c_gammah_gammah[bin_id+1] - data->c_gammah_gammah[bin_id]);
        data->dcs_gammah_gammah[i] = (data->c_gammah_gammah[bin_id+1] - data->c_gammah_gammah[bin_id]);;
        data->dcs_gammah_gammah[i] /= data->dT[i];
	data->dcs_gammah_gammah[i] *= data->invTs[i];
    
        data->cs_gloverabel08_gael[i] = data->c_gloverabel08_gael[bin_id] +
            data->Tdef[i] * (data->c_gloverabel08_gael[bin_id+1] - data->c_gloverabel08_gael[bin_id]);
        data->dcs_gloverabel08_gael[i] = (data->c_gloverabel08_gael[bin_id+1] - data->c_gloverabel08_gael[bin_id]);;
        data->dcs_gloverabel08_gael[i] /= data->dT[i];
	data->dcs_gloverabel08_gael[i] *= data->invTs[i];
    
        data->cs_gloverabel08_gaH2[i] = data->c_gloverabel08_gaH2[bin_id] +
            data->Tdef[i] * (data->c_gloverabel08_gaH2[bin_id+1] - data->c_gloverabel08_gaH2[bin_id]);
        data->dcs_gloverabel08_gaH2[i] = (data->c_gloverabel08_gaH2[bin_id+1] - data->c_gloverabel08_gaH2[bin_id]);;
        data->dcs_gloverabel08_gaH2[i] /= data->dT[i];
	data->dcs_gloverabel08_gaH2[i] *= data->invTs[i];
        
        data->cs_gloverabel08_gaHe[i] = data->c_gloverabel08_gaHe[bin_id] +
            data->Tdef[i] * (data->c_gloverabel08_gaHe[bin_id+1] - data->c_gloverabel08_gaHe[bin_id]);
        data->dcs_gloverabel08_gaHe[i] = (data->c_gloverabel08_gaHe[bin_id+1] - data->c_gloverabel08_gaHe[bin_id]);;
        data->dcs_gloverabel08_gaHe[i] /= data->dT[i];
	data->dcs_gloverabel08_gaHe[i] *= data->invTs[i];
        
        data->cs_gloverabel08_gaHI[i] = data->c_gloverabel08_gaHI[bin_id] +
            data->Tdef[i] * (data->c_gloverabel08_gaHI[bin_id+1] - data->c_gloverabel08_gaHI[bin_id]);
        data->dcs_gloverabel08_gaHI[i] = (data->c_gloverabel08_gaHI[bin_id+1] - data->c_gloverabel08_gaHI[bin_id]);;
        data->dcs_gloverabel08_gaHI[i] /= data->dT[i];
	data->dcs_gloverabel08_gaHI[i] *= data->invTs[i];
        
        data->cs_gloverabel08_gaHp[i] = data->c_gloverabel08_gaHp[bin_id] +
            data->Tdef[i] * (data->c_gloverabel08_gaHp[bin_id+1] - data->c_gloverabel08_gaHp[bin_id]);
        data->dcs_gloverabel08_gaHp[i] = (data->c_gloverabel08_gaHp[bin_id+1] - data->c_gloverabel08_gaHp[bin_id]);;
        data->dcs_gloverabel08_gaHp[i] /= data->dT[i];
	data->dcs_gloverabel08_gaHp[i] *= data->invTs[i];
        
        data->cs_gloverabel08_gphdl[i] = data->c_gloverabel08_gphdl[bin_id] +
            data->Tdef[i] * (data->c_gloverabel08_gphdl[bin_id+1] - data->c_gloverabel08_gphdl[bin_id]);
        data->dcs_gloverabel08_gphdl[i] = (data->c_gloverabel08_gphdl[bin_id+1] - data->c_gloverabel08_gphdl[bin_id]);;
        data->dcs_gloverabel08_gphdl[i] /= data->dT[i];
	data->dcs_gloverabel08_gphdl[i] *= data->invTs[i];
        
        data->cs_gloverabel08_gpldl[i] = data->c_gloverabel08_gpldl[bin_id] +
            data->Tdef[i] * (data->c_gloverabel08_gpldl[bin_id+1] - data->c_gloverabel08_gpldl[bin_id]);
        data->dcs_gloverabel08_gpldl[i] = (data->c_gloverabel08_gpldl[bin_id+1] - data->c_gloverabel08_gpldl[bin_id]);;
        data->dcs_gloverabel08_gpldl[i] /= data->dT[i];
	data->dcs_gloverabel08_gpldl[i] *= data->invTs[i];
        
        data->cs_gloverabel08_h2lte[i] = data->c_gloverabel08_h2lte[bin_id] +
            data->Tdef[i] * (data->c_gloverabel08_h2lte[bin_id+1] - data->c_gloverabel08_h2lte[bin_id]);
        data->dcs_gloverabel08_h2lte[i] = (data->c_gloverabel08_h2lte[bin_id+1] - data->c_gloverabel08_h2lte[bin_id]);;
        data->dcs_gloverabel08_h2lte[i] /= data->dT[i];
	data->dcs_gloverabel08_h2lte[i] *= data->invTs[i];
    
        data->cs_h2formation_h2mcool[i] = data->c_h2formation_h2mcool[bin_id] +
            data->Tdef[i] * (data->c_h2formation_h2mcool[bin_id+1] - data->c_h2formation_h2mcool[bin_id]);
        data->dcs_h2formation_h2mcool[i] = (data->c_h2formation_h2mcool[bin_id+1] - data->c_h2formation_h2mcool[bin_id]);;
        data->dcs_h2formation_h2mcool[i] /= data->dT[i];
	data->dcs_h2formation_h2mcool[i] *= data->invTs[i];
        
        data->cs_h2formation_h2mheat[i] = data->c_h2formation_h2mheat[bin_id] +
            data->Tdef[i] * (data->c_h2formation_h2mheat[bin_id+1] - data->c_h2formation_h2mheat[bin_id]);
        data->dcs_h2formation_h2mheat[i] = (data->c_h2formation_h2mheat[bin_id+1] - data->c_h2formation_h2mheat[bin_id]);;
        data->dcs_h2formation_h2mheat[i] /= data->dT[i];
	data->dcs_h2formation_h2mheat[i] *= data->invTs[i];
        
        data->cs_h2formation_ncrd1[i] = data->c_h2formation_ncrd1[bin_id] +
            data->Tdef[i] * (data->c_h2formation_ncrd1[bin_id+1] - data->c_h2formation_ncrd1[bin_id]);
        data->dcs_h2formation_ncrd1[i] = (data->c_h2formation_ncrd1[bin_id+1] - data->c_h2formation_ncrd1[bin_id]);;
        data->dcs_h2formation_ncrd1[i] /= data->dT[i];
	data->dcs_h2formation_ncrd1[i] *= data->invTs[i];
        
        data->cs_h2formation_ncrd2[i] = data->c_h2formation_ncrd2[bin_id] +
            data->Tdef[i] * (data->c_h2formation_ncrd2[bin_id+1] - data->c_h2formation_ncrd2[bin_id]);
        data->dcs_h2formation_ncrd2[i] = (data->c_h2formation_ncrd2[bin_id+1] - data->c_h2formation_ncrd2[bin_id]);;
        data->dcs_h2formation_ncrd2[i] /= data->dT[i];
	data->dcs_h2formation_ncrd2[i] *= data->invTs[i];
        
        data->cs_h2formation_ncrn[i] = data->c_h2formation_ncrn[bin_id] +
            data->Tdef[i] * (data->c_h2formation_ncrn[bin_id+1] - data->c_h2formation_ncrn[bin_id]);
        data->dcs_h2formation_ncrn[i] = (data->c_h2formation_ncrn[bin_id+1] - data->c_h2formation_ncrn[bin_id]);;
        data->dcs_h2formation_ncrn[i] /= data->dT[i];
	data->dcs_h2formation_ncrn[i] *= data->invTs[i];
    
        data->cs_reHeII1_reHeII1[i] = data->c_reHeII1_reHeII1[bin_id] +
            data->Tdef[i] * (data->c_reHeII1_reHeII1[bin_id+1] - data->c_reHeII1_reHeII1[bin_id]);
        data->dcs_reHeII1_reHeII1[i] = (data->c_reHeII1_reHeII1[bin_id+1] - data->c_reHeII1_reHeII1[bin_id]);;
        data->dcs_reHeII1_reHeII1[i] /= data->dT[i];
	data->dcs_reHeII1_reHeII1[i] *= data->invTs[i];
    
        data->cs_reHeII2_reHeII2[i] = data->c_reHeII2_reHeII2[bin_id] +
            data->Tdef[i] * (data->c_reHeII2_reHeII2[bin_id+1] - data->c_reHeII2_reHeII2[bin_id]);
        data->dcs_reHeII2_reHeII2[i] = (data->c_reHeII2_reHeII2[bin_id+1] - data->c_reHeII2_reHeII2[bin_id]);;
        data->dcs_reHeII2_reHeII2[i] /= data->dT[i];
	data->dcs_reHeII2_reHeII2[i] *= data->invTs[i];
    
        data->cs_reHeIII_reHeIII[i] = data->c_reHeIII_reHeIII[bin_id] +
            data->Tdef[i] * (data->c_reHeIII_reHeIII[bin_id+1] - data->c_reHeIII_reHeIII[bin_id]);
        data->dcs_reHeIII_reHeIII[i] = (data->c_reHeIII_reHeIII[bin_id+1] - data->c_reHeIII_reHeIII[bin_id]);;
        data->dcs_reHeIII_reHeIII[i] /= data->dT[i];
	data->dcs_reHeIII_reHeIII[i] *= data->invTs[i];
    
        data->cs_reHII_reHII[i] = data->c_reHII_reHII[bin_id] +
            data->Tdef[i] * (data->c_reHII_reHII[bin_id+1] - data->c_reHII_reHII[bin_id]);
        data->dcs_reHII_reHII[i] = (data->c_reHII_reHII[bin_id+1] - data->c_reHII_reHII[bin_id]);;
        data->dcs_reHII_reHII[i] /= data->dT[i];
	data->dcs_reHII_reHII[i] *= data->invTs[i];
    

}
 


void eq_H2_2_interpolate_gamma(eq_H2_2_data *data,
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
 


int calculate_jacobian_eq_H2_2( realtype t,
                                        N_Vector y, N_Vector fy,
                                        SUNMatrix J, void *user_data,
                                        N_Vector tmp1, N_Vector tmp2,
                                        N_Vector tmp3)
{
    /* We iterate over all of the rates */
    /* Calcuate temperature first */
    
    int nstrip = 1;
    int nchem = 10;

    eq_H2_2_data *data = (eq_H2_2_data*)user_data; 
    

    int i, j;
    j = 0;
    /* change N_Vector back to an array */
    double y_arr[ 10 ];
    y_arr[0] = NV_Ith_S(y , 0);
    y_arr[1] = NV_Ith_S(y , 1);
    y_arr[2] = NV_Ith_S(y , 2);
    y_arr[3] = NV_Ith_S(y , 3);
    y_arr[4] = NV_Ith_S(y , 4);
    y_arr[5] = NV_Ith_S(y , 5);
    y_arr[6] = NV_Ith_S(y , 6);
    y_arr[7] = NV_Ith_S(y , 7);
    y_arr[8] = NV_Ith_S(y , 8);
    y_arr[9] = NV_Ith_S(y , 9);

    // eq_H2_2_ealculate_temperature(data, y_arr, nstrip, nchem);
    // eq_H2_2_interpolate_rates(data, nstrip);

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
    double *brem_brem = data->cs_brem_brem;
    double *rbrem_brem = data->dcs_brem_brem;
    double *ceHeI_ceHeI = data->cs_ceHeI_ceHeI;
    double *rceHeI_ceHeI = data->dcs_ceHeI_ceHeI;
    double *ceHeII_ceHeII = data->cs_ceHeII_ceHeII;
    double *rceHeII_ceHeII = data->dcs_ceHeII_ceHeII;
    double *ceHI_ceHI = data->cs_ceHI_ceHI;
    double *rceHI_ceHI = data->dcs_ceHI_ceHI;
    double *ciHeI_ciHeI = data->cs_ciHeI_ciHeI;
    double *rciHeI_ciHeI = data->dcs_ciHeI_ciHeI;
    double *ciHeII_ciHeII = data->cs_ciHeII_ciHeII;
    double *rciHeII_ciHeII = data->dcs_ciHeII_ciHeII;
    double *ciHeIS_ciHeIS = data->cs_ciHeIS_ciHeIS;
    double *rciHeIS_ciHeIS = data->dcs_ciHeIS_ciHeIS;
    double *ciHI_ciHI = data->cs_ciHI_ciHI;
    double *rciHI_ciHI = data->dcs_ciHI_ciHI;
    double *compton_comp_ = data->cs_compton_comp_;
    double *rcompton_comp_ = data->dcs_compton_comp_;
    double *gammah_gammah = data->cs_gammah_gammah;
    double *rgammah_gammah = data->dcs_gammah_gammah;
    double *gloverabel08_gael = data->cs_gloverabel08_gael;
    double *rgloverabel08_gael = data->dcs_gloverabel08_gael;
    double *gloverabel08_gaH2 = data->cs_gloverabel08_gaH2;
    double *rgloverabel08_gaH2 = data->dcs_gloverabel08_gaH2;
    double *gloverabel08_gaHe = data->cs_gloverabel08_gaHe;
    double *rgloverabel08_gaHe = data->dcs_gloverabel08_gaHe;
    double *gloverabel08_gaHI = data->cs_gloverabel08_gaHI;
    double *rgloverabel08_gaHI = data->dcs_gloverabel08_gaHI;
    double *gloverabel08_gaHp = data->cs_gloverabel08_gaHp;
    double *rgloverabel08_gaHp = data->dcs_gloverabel08_gaHp;
    double *gloverabel08_gphdl = data->cs_gloverabel08_gphdl;
    double *rgloverabel08_gphdl = data->dcs_gloverabel08_gphdl;
    double *gloverabel08_gpldl = data->cs_gloverabel08_gpldl;
    double *rgloverabel08_gpldl = data->dcs_gloverabel08_gpldl;
    double *gloverabel08_h2lte = data->cs_gloverabel08_h2lte;
    double *rgloverabel08_h2lte = data->dcs_gloverabel08_h2lte;
    double *h2formation_h2mcool = data->cs_h2formation_h2mcool;
    double *rh2formation_h2mcool = data->dcs_h2formation_h2mcool;
    double *h2formation_h2mheat = data->cs_h2formation_h2mheat;
    double *rh2formation_h2mheat = data->dcs_h2formation_h2mheat;
    double *h2formation_ncrd1 = data->cs_h2formation_ncrd1;
    double *rh2formation_ncrd1 = data->dcs_h2formation_ncrd1;
    double *h2formation_ncrd2 = data->cs_h2formation_ncrd2;
    double *rh2formation_ncrd2 = data->dcs_h2formation_ncrd2;
    double *h2formation_ncrn = data->cs_h2formation_ncrn;
    double *rh2formation_ncrn = data->dcs_h2formation_ncrn;
    double *reHeII1_reHeII1 = data->cs_reHeII1_reHeII1;
    double *rreHeII1_reHeII1 = data->dcs_reHeII1_reHeII1;
    double *reHeII2_reHeII2 = data->cs_reHeII2_reHeII2;
    double *rreHeII2_reHeII2 = data->dcs_reHeII2_reHeII2;
    double *reHeIII_reHeIII = data->cs_reHeIII_reHeIII;
    double *rreHeIII_reHeIII = data->dcs_reHeIII_reHeIII;
    double *reHII_reHII = data->cs_reHII_reHII;
    double *rreHII_reHII = data->dcs_reHII_reHII;
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
    double mdensity, inv_mdensity;
    
    int jj;
    jj = 0;
    double scale, scale1, inv_scale2;

    j = 0;
    mdensity = 0.0;
    z = data->current_z;
    scale = data->scale[j];
    H2_1 = NV_Ith_S( y, 0  )*scale;
    
    mdensity += H2_1 * 2.01588;
    

    j++;
    
    scale = data->scale[j];
    H2_2 = NV_Ith_S( y, 1  )*scale;
    
    mdensity += H2_2 * 2.01588;
    

    j++;
    
    scale = data->scale[j];
    H_1 = NV_Ith_S( y, 2  )*scale;
    
    mdensity += H_1 * 1.00794;
    

    j++;
    
    scale = data->scale[j];
    H_2 = NV_Ith_S( y, 3  )*scale;
    
    mdensity += H_2 * 1.00794;
    

    j++;
    
    scale = data->scale[j];
    H_m0 = NV_Ith_S( y, 4  )*scale;
    
    mdensity += H_m0 * 1.00794;
    

    j++;
    
    scale = data->scale[j];
    He_1 = NV_Ith_S( y, 5  )*scale;
    
    mdensity += He_1 * 4.002602;
    

    j++;
    
    scale = data->scale[j];
    He_2 = NV_Ith_S( y, 6  )*scale;
    
    mdensity += He_2 * 4.002602;
    

    j++;
    
    scale = data->scale[j];
    He_3 = NV_Ith_S( y, 7  )*scale;
    
    mdensity += He_3 * 4.002602;
    

    j++;
    
    scale = data->scale[j];
    de = NV_Ith_S( y, 8  )*scale;
    

    j++;
    
    scale = data->scale[j];
    ge = NV_Ith_S( y, 9  )*scale;
    

    j++;
    

    mdensity *= mh;
    inv_mdensity = 1.0 / mdensity; 
    i = 0;
    j = i * nchem * nchem;
    //
    // Species: H2_1
    //
    
    // H2_1 by H2_1
    
    IJth(J, 1, 1 ) = -k11[i]*H_2 - k12[i]*de - k13[i]*H_1 + k21[i]*pow(H_1, 2);
    
    inv_scale2 = data->inv_scale[0];
    scale1     = data->scale[0];
    IJth(J, 1, 1) *= inv_scale2*scale1;

    

    
    // H2_2 by H2_1
    
    IJth(J, 2, 1 ) = k11[i]*H_2;
    
    inv_scale2 = data->inv_scale[1];
    scale1     = data->scale[0];
    IJth(J, 2, 1) *= inv_scale2*scale1;

    

    
    // H_1 by H2_1
    
    IJth(J, 3, 1 ) = k11[i]*H_2 + 2*k12[i]*de + 2*k13[i]*H_1 - 2*k21[i]*pow(H_1, 2);
    
    inv_scale2 = data->inv_scale[2];
    scale1     = data->scale[0];
    IJth(J, 3, 1) *= inv_scale2*scale1;

    

    
    // H_2 by H2_1
    
    IJth(J, 4, 1 ) = -k11[i]*H_2;
    
    inv_scale2 = data->inv_scale[3];
    scale1     = data->scale[0];
    IJth(J, 4, 1) *= inv_scale2*scale1;

    

    
    // H_m0 by H2_1
    
    IJth(J, 5, 1 ) = 0;
    
    inv_scale2 = data->inv_scale[4];
    scale1     = data->scale[0];
    IJth(J, 5, 1) *= inv_scale2*scale1;

    

    
    // He_1 by H2_1
    
    IJth(J, 6, 1 ) = 0;
    
    inv_scale2 = data->inv_scale[5];
    scale1     = data->scale[0];
    IJth(J, 6, 1) *= inv_scale2*scale1;

    

    
    // He_2 by H2_1
    
    IJth(J, 7, 1 ) = 0;
    
    inv_scale2 = data->inv_scale[6];
    scale1     = data->scale[0];
    IJth(J, 7, 1) *= inv_scale2*scale1;

    

    
    // He_3 by H2_1
    
    IJth(J, 8, 1 ) = 0;
    
    inv_scale2 = data->inv_scale[7];
    scale1     = data->scale[0];
    IJth(J, 8, 1) *= inv_scale2*scale1;

    

    
    // de by H2_1
    
    IJth(J, 9, 1 ) = 0;
    
    inv_scale2 = data->inv_scale[8];
    scale1     = data->scale[0];
    IJth(J, 9, 1) *= inv_scale2*scale1;

    

    
    // ge by H2_1
    
    IJth(J, 10, 1 ) = -H2_1*gloverabel08_gaH2[i]*pow(gloverabel08_h2lte[i], 2)/(pow(gloverabel08_h2lte[i]/(H2_1*gloverabel08_gaH2[i] + H_1*gloverabel08_gaHI[i] + H_2*gloverabel08_gaHp[i] + He_1*gloverabel08_gaHe[i] + de*gloverabel08_gael[i]) + 1.0, 2)*pow(H2_1*gloverabel08_gaH2[i] + H_1*gloverabel08_gaHI[i] + H_2*gloverabel08_gaHp[i] + He_1*gloverabel08_gaHe[i] + de*gloverabel08_gael[i], 2)) - 0.5*H_1*h2formation_h2mcool[i]*1.0/(h2formation_ncrn[i]/(H2_1*h2formation_ncrd2[i] + H_1*h2formation_ncrd1[i]) + 1.0) - gloverabel08_h2lte[i]/(gloverabel08_h2lte[i]/(H2_1*gloverabel08_gaH2[i] + H_1*gloverabel08_gaHI[i] + H_2*gloverabel08_gaHp[i] + He_1*gloverabel08_gaHe[i] + de*gloverabel08_gael[i]) + 1.0) + 0.5*h2formation_ncrd2[i]*h2formation_ncrn[i]*pow(h2formation_ncrn[i]/(H2_1*h2formation_ncrd2[i] + H_1*h2formation_ncrd1[i]) + 1.0, -2.0)*(-H2_1*H_1*h2formation_h2mcool[i] + pow(H_1, 3)*h2formation_h2mheat[i])/pow(H2_1*h2formation_ncrd2[i] + H_1*h2formation_ncrd1[i], 2);
    
    inv_scale2 = data->inv_scale[9];
    scale1     = data->scale[0];
    IJth(J, 10, 1) *= inv_scale2*scale1;

    
    IJth(J, 10, 1 ) *= inv_mdensity;
    

    
    
    //
    // Species: H2_2
    //
    
    // H2_1 by H2_2
    
    IJth(J, 1, 2 ) = k10[i]*H_1 + k19[i]*H_m0;
    
    inv_scale2 = data->inv_scale[0];
    scale1     = data->scale[1];
    IJth(J, 1, 2) *= inv_scale2*scale1;

    

    
    // H2_2 by H2_2
    
    IJth(J, 2, 2 ) = -k10[i]*H_1 - k18[i]*de - k19[i]*H_m0;
    
    inv_scale2 = data->inv_scale[1];
    scale1     = data->scale[1];
    IJth(J, 2, 2) *= inv_scale2*scale1;

    

    
    // H_1 by H2_2
    
    IJth(J, 3, 2 ) = -k10[i]*H_1 + 2*k18[i]*de + k19[i]*H_m0;
    
    inv_scale2 = data->inv_scale[2];
    scale1     = data->scale[1];
    IJth(J, 3, 2) *= inv_scale2*scale1;

    

    
    // H_2 by H2_2
    
    IJth(J, 4, 2 ) = k10[i]*H_1;
    
    inv_scale2 = data->inv_scale[3];
    scale1     = data->scale[1];
    IJth(J, 4, 2) *= inv_scale2*scale1;

    

    
    // H_m0 by H2_2
    
    IJth(J, 5, 2 ) = -k19[i]*H_m0;
    
    inv_scale2 = data->inv_scale[4];
    scale1     = data->scale[1];
    IJth(J, 5, 2) *= inv_scale2*scale1;

    

    
    // He_1 by H2_2
    
    IJth(J, 6, 2 ) = 0;
    
    inv_scale2 = data->inv_scale[5];
    scale1     = data->scale[1];
    IJth(J, 6, 2) *= inv_scale2*scale1;

    

    
    // He_2 by H2_2
    
    IJth(J, 7, 2 ) = 0;
    
    inv_scale2 = data->inv_scale[6];
    scale1     = data->scale[1];
    IJth(J, 7, 2) *= inv_scale2*scale1;

    

    
    // He_3 by H2_2
    
    IJth(J, 8, 2 ) = 0;
    
    inv_scale2 = data->inv_scale[7];
    scale1     = data->scale[1];
    IJth(J, 8, 2) *= inv_scale2*scale1;

    

    
    // de by H2_2
    
    IJth(J, 9, 2 ) = -k18[i]*de;
    
    inv_scale2 = data->inv_scale[8];
    scale1     = data->scale[1];
    IJth(J, 9, 2) *= inv_scale2*scale1;

    

    
    // ge by H2_2
    
    IJth(J, 10, 2 ) = 0;
    
    inv_scale2 = data->inv_scale[9];
    scale1     = data->scale[1];
    IJth(J, 10, 2) *= inv_scale2*scale1;

    
    IJth(J, 10, 2 ) *= inv_mdensity;
    

    
    
    //
    // Species: H_1
    //
    
    // H2_1 by H_1
    
    IJth(J, 1, 3 ) = k08[i]*H_m0 + k10[i]*H2_2 - k13[i]*H2_1 + 2*k21[i]*H2_1*H_1 + 3*k22[i]*pow(H_1, 2);
    
    inv_scale2 = data->inv_scale[0];
    scale1     = data->scale[2];
    IJth(J, 1, 3) *= inv_scale2*scale1;

    

    
    // H2_2 by H_1
    
    IJth(J, 2, 3 ) = k09[i]*H_2 - k10[i]*H2_2;
    
    inv_scale2 = data->inv_scale[1];
    scale1     = data->scale[2];
    IJth(J, 2, 3) *= inv_scale2*scale1;

    

    
    // H_1 by H_1
    
    IJth(J, 3, 3 ) = -k01[i]*de - k07[i]*de - k08[i]*H_m0 - k09[i]*H_2 - k10[i]*H2_2 + 2*k13[i]*H2_1 + k15[i]*H_m0 - 4*k21[i]*H2_1*H_1 - 6*k22[i]*pow(H_1, 2);
    
    inv_scale2 = data->inv_scale[2];
    scale1     = data->scale[2];
    IJth(J, 3, 3) *= inv_scale2*scale1;

    

    
    // H_2 by H_1
    
    IJth(J, 4, 3 ) = k01[i]*de - k09[i]*H_2 + k10[i]*H2_2;
    
    inv_scale2 = data->inv_scale[3];
    scale1     = data->scale[2];
    IJth(J, 4, 3) *= inv_scale2*scale1;

    

    
    // H_m0 by H_1
    
    IJth(J, 5, 3 ) = k07[i]*de - k08[i]*H_m0 - k15[i]*H_m0;
    
    inv_scale2 = data->inv_scale[4];
    scale1     = data->scale[2];
    IJth(J, 5, 3) *= inv_scale2*scale1;

    

    
    // He_1 by H_1
    
    IJth(J, 6, 3 ) = 0;
    
    inv_scale2 = data->inv_scale[5];
    scale1     = data->scale[2];
    IJth(J, 6, 3) *= inv_scale2*scale1;

    

    
    // He_2 by H_1
    
    IJth(J, 7, 3 ) = 0;
    
    inv_scale2 = data->inv_scale[6];
    scale1     = data->scale[2];
    IJth(J, 7, 3) *= inv_scale2*scale1;

    

    
    // He_3 by H_1
    
    IJth(J, 8, 3 ) = 0;
    
    inv_scale2 = data->inv_scale[7];
    scale1     = data->scale[2];
    IJth(J, 8, 3) *= inv_scale2*scale1;

    

    
    // de by H_1
    
    IJth(J, 9, 3 ) = k01[i]*de - k07[i]*de + k08[i]*H_m0 + k15[i]*H_m0;
    
    inv_scale2 = data->inv_scale[8];
    scale1     = data->scale[2];
    IJth(J, 9, 3) *= inv_scale2*scale1;

    

    
    // ge by H_1
    
    IJth(J, 10, 3 ) = -H2_1*gloverabel08_gaHI[i]*pow(gloverabel08_h2lte[i], 2)/(pow(gloverabel08_h2lte[i]/(H2_1*gloverabel08_gaH2[i] + H_1*gloverabel08_gaHI[i] + H_2*gloverabel08_gaHp[i] + He_1*gloverabel08_gaHe[i] + de*gloverabel08_gael[i]) + 1.0, 2)*pow(H2_1*gloverabel08_gaH2[i] + H_1*gloverabel08_gaHI[i] + H_2*gloverabel08_gaHp[i] + He_1*gloverabel08_gaHe[i] + de*gloverabel08_gael[i], 2)) - ceHI_ceHI[i]*de - ciHI_ciHI[i]*de + 0.5*h2formation_ncrd1[i]*h2formation_ncrn[i]*pow(h2formation_ncrn[i]/(H2_1*h2formation_ncrd2[i] + H_1*h2formation_ncrd1[i]) + 1.0, -2.0)*(-H2_1*H_1*h2formation_h2mcool[i] + pow(H_1, 3)*h2formation_h2mheat[i])/pow(H2_1*h2formation_ncrd2[i] + H_1*h2formation_ncrd1[i], 2) + 0.5*(-H2_1*h2formation_h2mcool[i] + 3*pow(H_1, 2)*h2formation_h2mheat[i])*1.0/(h2formation_ncrn[i]/(H2_1*h2formation_ncrd2[i] + H_1*h2formation_ncrd1[i]) + 1.0);
    
    inv_scale2 = data->inv_scale[9];
    scale1     = data->scale[2];
    IJth(J, 10, 3) *= inv_scale2*scale1;

    
    IJth(J, 10, 3 ) *= inv_mdensity;
    

    
    
    //
    // Species: H_2
    //
    
    // H2_1 by H_2
    
    IJth(J, 1, 4 ) = -k11[i]*H2_1;
    
    inv_scale2 = data->inv_scale[0];
    scale1     = data->scale[3];
    IJth(J, 1, 4) *= inv_scale2*scale1;

    

    
    // H2_2 by H_2
    
    IJth(J, 2, 4 ) = k09[i]*H_1 + k11[i]*H2_1 + k17[i]*H_m0;
    
    inv_scale2 = data->inv_scale[1];
    scale1     = data->scale[3];
    IJth(J, 2, 4) *= inv_scale2*scale1;

    

    
    // H_1 by H_2
    
    IJth(J, 3, 4 ) = k02[i]*de - k09[i]*H_1 + k11[i]*H2_1 + 2*k16[i]*H_m0;
    
    inv_scale2 = data->inv_scale[2];
    scale1     = data->scale[3];
    IJth(J, 3, 4) *= inv_scale2*scale1;

    

    
    // H_2 by H_2
    
    IJth(J, 4, 4 ) = -k02[i]*de - k09[i]*H_1 - k11[i]*H2_1 - k16[i]*H_m0 - k17[i]*H_m0;
    
    inv_scale2 = data->inv_scale[3];
    scale1     = data->scale[3];
    IJth(J, 4, 4) *= inv_scale2*scale1;

    

    
    // H_m0 by H_2
    
    IJth(J, 5, 4 ) = -k16[i]*H_m0 - k17[i]*H_m0;
    
    inv_scale2 = data->inv_scale[4];
    scale1     = data->scale[3];
    IJth(J, 5, 4) *= inv_scale2*scale1;

    

    
    // He_1 by H_2
    
    IJth(J, 6, 4 ) = 0;
    
    inv_scale2 = data->inv_scale[5];
    scale1     = data->scale[3];
    IJth(J, 6, 4) *= inv_scale2*scale1;

    

    
    // He_2 by H_2
    
    IJth(J, 7, 4 ) = 0;
    
    inv_scale2 = data->inv_scale[6];
    scale1     = data->scale[3];
    IJth(J, 7, 4) *= inv_scale2*scale1;

    

    
    // He_3 by H_2
    
    IJth(J, 8, 4 ) = 0;
    
    inv_scale2 = data->inv_scale[7];
    scale1     = data->scale[3];
    IJth(J, 8, 4) *= inv_scale2*scale1;

    

    
    // de by H_2
    
    IJth(J, 9, 4 ) = -k02[i]*de + k17[i]*H_m0;
    
    inv_scale2 = data->inv_scale[8];
    scale1     = data->scale[3];
    IJth(J, 9, 4) *= inv_scale2*scale1;

    

    
    // ge by H_2
    
    IJth(J, 10, 4 ) = -H2_1*gloverabel08_gaHp[i]*pow(gloverabel08_h2lte[i], 2)/(pow(gloverabel08_h2lte[i]/(H2_1*gloverabel08_gaH2[i] + H_1*gloverabel08_gaHI[i] + H_2*gloverabel08_gaHp[i] + He_1*gloverabel08_gaHe[i] + de*gloverabel08_gael[i]) + 1.0, 2)*pow(H2_1*gloverabel08_gaH2[i] + H_1*gloverabel08_gaHI[i] + H_2*gloverabel08_gaHp[i] + He_1*gloverabel08_gaHe[i] + de*gloverabel08_gael[i], 2)) - brem_brem[i]*de - de*reHII_reHII[i];
    
    inv_scale2 = data->inv_scale[9];
    scale1     = data->scale[3];
    IJth(J, 10, 4) *= inv_scale2*scale1;

    
    IJth(J, 10, 4 ) *= inv_mdensity;
    

    
    
    //
    // Species: H_m0
    //
    
    // H2_1 by H_m0
    
    IJth(J, 1, 5 ) = k08[i]*H_1 + k19[i]*H2_2;
    
    inv_scale2 = data->inv_scale[0];
    scale1     = data->scale[4];
    IJth(J, 1, 5) *= inv_scale2*scale1;

    

    
    // H2_2 by H_m0
    
    IJth(J, 2, 5 ) = k17[i]*H_2 - k19[i]*H2_2;
    
    inv_scale2 = data->inv_scale[1];
    scale1     = data->scale[4];
    IJth(J, 2, 5) *= inv_scale2*scale1;

    

    
    // H_1 by H_m0
    
    IJth(J, 3, 5 ) = -k08[i]*H_1 + k14[i]*de + k15[i]*H_1 + 2*k16[i]*H_2 + k19[i]*H2_2;
    
    inv_scale2 = data->inv_scale[2];
    scale1     = data->scale[4];
    IJth(J, 3, 5) *= inv_scale2*scale1;

    

    
    // H_2 by H_m0
    
    IJth(J, 4, 5 ) = -k16[i]*H_2 - k17[i]*H_2;
    
    inv_scale2 = data->inv_scale[3];
    scale1     = data->scale[4];
    IJth(J, 4, 5) *= inv_scale2*scale1;

    

    
    // H_m0 by H_m0
    
    IJth(J, 5, 5 ) = -k08[i]*H_1 - k14[i]*de - k15[i]*H_1 - k16[i]*H_2 - k17[i]*H_2 - k19[i]*H2_2;
    
    inv_scale2 = data->inv_scale[4];
    scale1     = data->scale[4];
    IJth(J, 5, 5) *= inv_scale2*scale1;

    

    
    // He_1 by H_m0
    
    IJth(J, 6, 5 ) = 0;
    
    inv_scale2 = data->inv_scale[5];
    scale1     = data->scale[4];
    IJth(J, 6, 5) *= inv_scale2*scale1;

    

    
    // He_2 by H_m0
    
    IJth(J, 7, 5 ) = 0;
    
    inv_scale2 = data->inv_scale[6];
    scale1     = data->scale[4];
    IJth(J, 7, 5) *= inv_scale2*scale1;

    

    
    // He_3 by H_m0
    
    IJth(J, 8, 5 ) = 0;
    
    inv_scale2 = data->inv_scale[7];
    scale1     = data->scale[4];
    IJth(J, 8, 5) *= inv_scale2*scale1;

    

    
    // de by H_m0
    
    IJth(J, 9, 5 ) = k08[i]*H_1 + k14[i]*de + k15[i]*H_1 + k17[i]*H_2;
    
    inv_scale2 = data->inv_scale[8];
    scale1     = data->scale[4];
    IJth(J, 9, 5) *= inv_scale2*scale1;

    

    
    // ge by H_m0
    
    IJth(J, 10, 5 ) = 0;
    
    inv_scale2 = data->inv_scale[9];
    scale1     = data->scale[4];
    IJth(J, 10, 5) *= inv_scale2*scale1;

    
    IJth(J, 10, 5 ) *= inv_mdensity;
    

    
    
    //
    // Species: He_1
    //
    
    // H2_1 by He_1
    
    IJth(J, 1, 6 ) = 0;
    
    inv_scale2 = data->inv_scale[0];
    scale1     = data->scale[5];
    IJth(J, 1, 6) *= inv_scale2*scale1;

    

    
    // H2_2 by He_1
    
    IJth(J, 2, 6 ) = 0;
    
    inv_scale2 = data->inv_scale[1];
    scale1     = data->scale[5];
    IJth(J, 2, 6) *= inv_scale2*scale1;

    

    
    // H_1 by He_1
    
    IJth(J, 3, 6 ) = 0;
    
    inv_scale2 = data->inv_scale[2];
    scale1     = data->scale[5];
    IJth(J, 3, 6) *= inv_scale2*scale1;

    

    
    // H_2 by He_1
    
    IJth(J, 4, 6 ) = 0;
    
    inv_scale2 = data->inv_scale[3];
    scale1     = data->scale[5];
    IJth(J, 4, 6) *= inv_scale2*scale1;

    

    
    // H_m0 by He_1
    
    IJth(J, 5, 6 ) = 0;
    
    inv_scale2 = data->inv_scale[4];
    scale1     = data->scale[5];
    IJth(J, 5, 6) *= inv_scale2*scale1;

    

    
    // He_1 by He_1
    
    IJth(J, 6, 6 ) = -k03[i]*de;
    
    inv_scale2 = data->inv_scale[5];
    scale1     = data->scale[5];
    IJth(J, 6, 6) *= inv_scale2*scale1;

    

    
    // He_2 by He_1
    
    IJth(J, 7, 6 ) = k03[i]*de;
    
    inv_scale2 = data->inv_scale[6];
    scale1     = data->scale[5];
    IJth(J, 7, 6) *= inv_scale2*scale1;

    

    
    // He_3 by He_1
    
    IJth(J, 8, 6 ) = 0;
    
    inv_scale2 = data->inv_scale[7];
    scale1     = data->scale[5];
    IJth(J, 8, 6) *= inv_scale2*scale1;

    

    
    // de by He_1
    
    IJth(J, 9, 6 ) = k03[i]*de;
    
    inv_scale2 = data->inv_scale[8];
    scale1     = data->scale[5];
    IJth(J, 9, 6) *= inv_scale2*scale1;

    

    
    // ge by He_1
    
    IJth(J, 10, 6 ) = -H2_1*gloverabel08_gaHe[i]*pow(gloverabel08_h2lte[i], 2)/(pow(gloverabel08_h2lte[i]/(H2_1*gloverabel08_gaH2[i] + H_1*gloverabel08_gaHI[i] + H_2*gloverabel08_gaHp[i] + He_1*gloverabel08_gaHe[i] + de*gloverabel08_gael[i]) + 1.0, 2)*pow(H2_1*gloverabel08_gaH2[i] + H_1*gloverabel08_gaHI[i] + H_2*gloverabel08_gaHp[i] + He_1*gloverabel08_gaHe[i] + de*gloverabel08_gael[i], 2)) - ciHeI_ciHeI[i]*de;
    
    inv_scale2 = data->inv_scale[9];
    scale1     = data->scale[5];
    IJth(J, 10, 6) *= inv_scale2*scale1;

    
    IJth(J, 10, 6 ) *= inv_mdensity;
    

    
    
    //
    // Species: He_2
    //
    
    // H2_1 by He_2
    
    IJth(J, 1, 7 ) = 0;
    
    inv_scale2 = data->inv_scale[0];
    scale1     = data->scale[6];
    IJth(J, 1, 7) *= inv_scale2*scale1;

    

    
    // H2_2 by He_2
    
    IJth(J, 2, 7 ) = 0;
    
    inv_scale2 = data->inv_scale[1];
    scale1     = data->scale[6];
    IJth(J, 2, 7) *= inv_scale2*scale1;

    

    
    // H_1 by He_2
    
    IJth(J, 3, 7 ) = 0;
    
    inv_scale2 = data->inv_scale[2];
    scale1     = data->scale[6];
    IJth(J, 3, 7) *= inv_scale2*scale1;

    

    
    // H_2 by He_2
    
    IJth(J, 4, 7 ) = 0;
    
    inv_scale2 = data->inv_scale[3];
    scale1     = data->scale[6];
    IJth(J, 4, 7) *= inv_scale2*scale1;

    

    
    // H_m0 by He_2
    
    IJth(J, 5, 7 ) = 0;
    
    inv_scale2 = data->inv_scale[4];
    scale1     = data->scale[6];
    IJth(J, 5, 7) *= inv_scale2*scale1;

    

    
    // He_1 by He_2
    
    IJth(J, 6, 7 ) = k04[i]*de;
    
    inv_scale2 = data->inv_scale[5];
    scale1     = data->scale[6];
    IJth(J, 6, 7) *= inv_scale2*scale1;

    

    
    // He_2 by He_2
    
    IJth(J, 7, 7 ) = -k04[i]*de - k05[i]*de;
    
    inv_scale2 = data->inv_scale[6];
    scale1     = data->scale[6];
    IJth(J, 7, 7) *= inv_scale2*scale1;

    

    
    // He_3 by He_2
    
    IJth(J, 8, 7 ) = k05[i]*de;
    
    inv_scale2 = data->inv_scale[7];
    scale1     = data->scale[6];
    IJth(J, 8, 7) *= inv_scale2*scale1;

    

    
    // de by He_2
    
    IJth(J, 9, 7 ) = -k04[i]*de + k05[i]*de;
    
    inv_scale2 = data->inv_scale[8];
    scale1     = data->scale[6];
    IJth(J, 9, 7) *= inv_scale2*scale1;

    

    
    // ge by He_2
    
    IJth(J, 10, 7 ) = -brem_brem[i]*de - ceHeII_ceHeII[i]*de - ceHeI_ceHeI[i]*pow(de, 2) - ciHeII_ciHeII[i]*de - ciHeIS_ciHeIS[i]*pow(de, 2) - de*reHeII1_reHeII1[i] - de*reHeII2_reHeII2[i];
    
    inv_scale2 = data->inv_scale[9];
    scale1     = data->scale[6];
    IJth(J, 10, 7) *= inv_scale2*scale1;

    
    IJth(J, 10, 7 ) *= inv_mdensity;
    

    
    
    //
    // Species: He_3
    //
    
    // H2_1 by He_3
    
    IJth(J, 1, 8 ) = 0;
    
    inv_scale2 = data->inv_scale[0];
    scale1     = data->scale[7];
    IJth(J, 1, 8) *= inv_scale2*scale1;

    

    
    // H2_2 by He_3
    
    IJth(J, 2, 8 ) = 0;
    
    inv_scale2 = data->inv_scale[1];
    scale1     = data->scale[7];
    IJth(J, 2, 8) *= inv_scale2*scale1;

    

    
    // H_1 by He_3
    
    IJth(J, 3, 8 ) = 0;
    
    inv_scale2 = data->inv_scale[2];
    scale1     = data->scale[7];
    IJth(J, 3, 8) *= inv_scale2*scale1;

    

    
    // H_2 by He_3
    
    IJth(J, 4, 8 ) = 0;
    
    inv_scale2 = data->inv_scale[3];
    scale1     = data->scale[7];
    IJth(J, 4, 8) *= inv_scale2*scale1;

    

    
    // H_m0 by He_3
    
    IJth(J, 5, 8 ) = 0;
    
    inv_scale2 = data->inv_scale[4];
    scale1     = data->scale[7];
    IJth(J, 5, 8) *= inv_scale2*scale1;

    

    
    // He_1 by He_3
    
    IJth(J, 6, 8 ) = 0;
    
    inv_scale2 = data->inv_scale[5];
    scale1     = data->scale[7];
    IJth(J, 6, 8) *= inv_scale2*scale1;

    

    
    // He_2 by He_3
    
    IJth(J, 7, 8 ) = k06[i]*de;
    
    inv_scale2 = data->inv_scale[6];
    scale1     = data->scale[7];
    IJth(J, 7, 8) *= inv_scale2*scale1;

    

    
    // He_3 by He_3
    
    IJth(J, 8, 8 ) = -k06[i]*de;
    
    inv_scale2 = data->inv_scale[7];
    scale1     = data->scale[7];
    IJth(J, 8, 8) *= inv_scale2*scale1;

    

    
    // de by He_3
    
    IJth(J, 9, 8 ) = -k06[i]*de;
    
    inv_scale2 = data->inv_scale[8];
    scale1     = data->scale[7];
    IJth(J, 9, 8) *= inv_scale2*scale1;

    

    
    // ge by He_3
    
    IJth(J, 10, 8 ) = -4.0*brem_brem[i]*de - de*reHeIII_reHeIII[i];
    
    inv_scale2 = data->inv_scale[9];
    scale1     = data->scale[7];
    IJth(J, 10, 8) *= inv_scale2*scale1;

    
    IJth(J, 10, 8 ) *= inv_mdensity;
    

    
    
    //
    // Species: de
    //
    
    // H2_1 by de
    
    IJth(J, 1, 9 ) = -k12[i]*H2_1;
    
    inv_scale2 = data->inv_scale[0];
    scale1     = data->scale[8];
    IJth(J, 1, 9) *= inv_scale2*scale1;

    

    
    // H2_2 by de
    
    IJth(J, 2, 9 ) = -k18[i]*H2_2;
    
    inv_scale2 = data->inv_scale[1];
    scale1     = data->scale[8];
    IJth(J, 2, 9) *= inv_scale2*scale1;

    

    
    // H_1 by de
    
    IJth(J, 3, 9 ) = -k01[i]*H_1 + k02[i]*H_2 - k07[i]*H_1 + 2*k12[i]*H2_1 + k14[i]*H_m0 + 2*k18[i]*H2_2;
    
    inv_scale2 = data->inv_scale[2];
    scale1     = data->scale[8];
    IJth(J, 3, 9) *= inv_scale2*scale1;

    

    
    // H_2 by de
    
    IJth(J, 4, 9 ) = k01[i]*H_1 - k02[i]*H_2;
    
    inv_scale2 = data->inv_scale[3];
    scale1     = data->scale[8];
    IJth(J, 4, 9) *= inv_scale2*scale1;

    

    
    // H_m0 by de
    
    IJth(J, 5, 9 ) = k07[i]*H_1 - k14[i]*H_m0;
    
    inv_scale2 = data->inv_scale[4];
    scale1     = data->scale[8];
    IJth(J, 5, 9) *= inv_scale2*scale1;

    

    
    // He_1 by de
    
    IJth(J, 6, 9 ) = -k03[i]*He_1 + k04[i]*He_2;
    
    inv_scale2 = data->inv_scale[5];
    scale1     = data->scale[8];
    IJth(J, 6, 9) *= inv_scale2*scale1;

    

    
    // He_2 by de
    
    IJth(J, 7, 9 ) = k03[i]*He_1 - k04[i]*He_2 - k05[i]*He_2 + k06[i]*He_3;
    
    inv_scale2 = data->inv_scale[6];
    scale1     = data->scale[8];
    IJth(J, 7, 9) *= inv_scale2*scale1;

    

    
    // He_3 by de
    
    IJth(J, 8, 9 ) = k05[i]*He_2 - k06[i]*He_3;
    
    inv_scale2 = data->inv_scale[7];
    scale1     = data->scale[8];
    IJth(J, 8, 9) *= inv_scale2*scale1;

    

    
    // de by de
    
    IJth(J, 9, 9 ) = k01[i]*H_1 - k02[i]*H_2 + k03[i]*He_1 - k04[i]*He_2 + k05[i]*He_2 - k06[i]*He_3 - k07[i]*H_1 + k14[i]*H_m0 - k18[i]*H2_2;
    
    inv_scale2 = data->inv_scale[8];
    scale1     = data->scale[8];
    IJth(J, 9, 9) *= inv_scale2*scale1;

    

    
    // ge by de
    
    IJth(J, 10, 9 ) = -H2_1*gloverabel08_gael[i]*pow(gloverabel08_h2lte[i], 2)/(pow(gloverabel08_h2lte[i]/(H2_1*gloverabel08_gaH2[i] + H_1*gloverabel08_gaHI[i] + H_2*gloverabel08_gaHp[i] + He_1*gloverabel08_gaHe[i] + de*gloverabel08_gael[i]) + 1.0, 2)*pow(H2_1*gloverabel08_gaH2[i] + H_1*gloverabel08_gaHI[i] + H_2*gloverabel08_gaHp[i] + He_1*gloverabel08_gaHe[i] + de*gloverabel08_gael[i], 2)) - H_1*ceHI_ceHI[i] - H_1*ciHI_ciHI[i] - H_2*reHII_reHII[i] - He_1*ciHeI_ciHeI[i] - He_2*ceHeII_ceHeII[i] - 2*He_2*ceHeI_ceHeI[i]*de - He_2*ciHeII_ciHeII[i] - 2*He_2*ciHeIS_ciHeIS[i]*de - He_2*reHeII1_reHeII1[i] - He_2*reHeII2_reHeII2[i] - He_3*reHeIII_reHeIII[i] - brem_brem[i]*(H_2 + He_2 + 4.0*He_3) - compton_comp_[i]*pow(z + 1.0, 4)*(T - 2.73*z - 2.73);
    
    inv_scale2 = data->inv_scale[9];
    scale1     = data->scale[8];
    IJth(J, 10, 9) *= inv_scale2*scale1;

    
    IJth(J, 10, 9 ) *= inv_mdensity;
    

    
    
    //
    // Species: ge
    //
    
    // H2_1 by ge
    
    IJth(J, 1, 10 ) = 0;
    
    inv_scale2 = data->inv_scale[0];
    scale1     = data->scale[9];
    IJth(J, 1, 10) *= inv_scale2*scale1;

    

    
    IJth(J, 1, 10 ) *= Tge[i];
    
    // H2_2 by ge
    
    IJth(J, 2, 10 ) = 0;
    
    inv_scale2 = data->inv_scale[1];
    scale1     = data->scale[9];
    IJth(J, 2, 10) *= inv_scale2*scale1;

    

    
    IJth(J, 2, 10 ) *= Tge[i];
    
    // H_1 by ge
    
    IJth(J, 3, 10 ) = 0;
    
    inv_scale2 = data->inv_scale[2];
    scale1     = data->scale[9];
    IJth(J, 3, 10) *= inv_scale2*scale1;

    

    
    IJth(J, 3, 10 ) *= Tge[i];
    
    // H_2 by ge
    
    IJth(J, 4, 10 ) = 0;
    
    inv_scale2 = data->inv_scale[3];
    scale1     = data->scale[9];
    IJth(J, 4, 10) *= inv_scale2*scale1;

    

    
    IJth(J, 4, 10 ) *= Tge[i];
    
    // H_m0 by ge
    
    IJth(J, 5, 10 ) = 0;
    
    inv_scale2 = data->inv_scale[4];
    scale1     = data->scale[9];
    IJth(J, 5, 10) *= inv_scale2*scale1;

    

    
    IJth(J, 5, 10 ) *= Tge[i];
    
    // He_1 by ge
    
    IJth(J, 6, 10 ) = 0;
    
    inv_scale2 = data->inv_scale[5];
    scale1     = data->scale[9];
    IJth(J, 6, 10) *= inv_scale2*scale1;

    

    
    IJth(J, 6, 10 ) *= Tge[i];
    
    // He_2 by ge
    
    IJth(J, 7, 10 ) = 0;
    
    inv_scale2 = data->inv_scale[6];
    scale1     = data->scale[9];
    IJth(J, 7, 10) *= inv_scale2*scale1;

    

    
    IJth(J, 7, 10 ) *= Tge[i];
    
    // He_3 by ge
    
    IJth(J, 8, 10 ) = 0;
    
    inv_scale2 = data->inv_scale[7];
    scale1     = data->scale[9];
    IJth(J, 8, 10) *= inv_scale2*scale1;

    

    
    IJth(J, 8, 10 ) *= Tge[i];
    
    // de by ge
    
    IJth(J, 9, 10 ) = 0;
    
    inv_scale2 = data->inv_scale[8];
    scale1     = data->scale[9];
    IJth(J, 9, 10) *= inv_scale2*scale1;

    

    
    IJth(J, 9, 10 ) *= Tge[i];
    
    // ge by ge
    
    IJth(J, 10, 10 ) = 0;
    
    inv_scale2 = data->inv_scale[9];
    scale1     = data->scale[9];
    IJth(J, 10, 10) *= inv_scale2*scale1;

    
    IJth(J, 10, 10 ) *= inv_mdensity;
    

    
    IJth(J, 10, 10 ) *= Tge[i];
    
    
    
    return 0;
}






int calculate_rhs_eq_H2_2(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
    eq_H2_2_data *data = (eq_H2_2_data* ) user_data;
    int i, j;

    int nchem = 10;
    int nstrip = 1;
    
    /* change N_Vector back to an array */
    double y_arr[ 10 ];
    /* the variable is ALREADY scaled in "calculate temperature" */
    
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
    double *brem_brem = data->cs_brem_brem;
    double *ceHeI_ceHeI = data->cs_ceHeI_ceHeI;
    double *ceHeII_ceHeII = data->cs_ceHeII_ceHeII;
    double *ceHI_ceHI = data->cs_ceHI_ceHI;
    double *ciHeI_ciHeI = data->cs_ciHeI_ciHeI;
    double *ciHeII_ciHeII = data->cs_ciHeII_ciHeII;
    double *ciHeIS_ciHeIS = data->cs_ciHeIS_ciHeIS;
    double *ciHI_ciHI = data->cs_ciHI_ciHI;
    double *compton_comp_ = data->cs_compton_comp_;
    double *gammah_gammah = data->cs_gammah_gammah;
    double *gloverabel08_gael = data->cs_gloverabel08_gael;
    double *gloverabel08_gaH2 = data->cs_gloverabel08_gaH2;
    double *gloverabel08_gaHe = data->cs_gloverabel08_gaHe;
    double *gloverabel08_gaHI = data->cs_gloverabel08_gaHI;
    double *gloverabel08_gaHp = data->cs_gloverabel08_gaHp;
    double *gloverabel08_gphdl = data->cs_gloverabel08_gphdl;
    double *gloverabel08_gpldl = data->cs_gloverabel08_gpldl;
    double *gloverabel08_h2lte = data->cs_gloverabel08_h2lte;
    double *h2formation_h2mcool = data->cs_h2formation_h2mcool;
    double *h2formation_h2mheat = data->cs_h2formation_h2mheat;
    double *h2formation_ncrd1 = data->cs_h2formation_ncrd1;
    double *h2formation_ncrd2 = data->cs_h2formation_ncrd2;
    double *h2formation_ncrn = data->cs_h2formation_ncrn;
    double *reHeII1_reHeII1 = data->cs_reHeII1_reHeII1;
    double *reHeII2_reHeII2 = data->cs_reHeII2_reHeII2;
    double *reHeIII_reHeIII = data->cs_reHeIII_reHeIII;
    double *reHII_reHII = data->cs_reHII_reHII;
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
    
    i = 0;
    T = data->Ts[i];
    z = data->current_z;
    
    double scale, inv_scale;
    int jj =0;
    scale = data->scale[jj];
    y_arr[jj] = H2_1 = NV_Ith_S( y,0 )*scale;
    jj++;
    
    mdensity += H2_1 * 2.01588;
    


    
    scale = data->scale[jj];
    y_arr[jj] = H2_2 = NV_Ith_S( y,1 )*scale;
    jj++;
    
    mdensity += H2_2 * 2.01588;
    


    
    scale = data->scale[jj];
    y_arr[jj] = H_1 = NV_Ith_S( y,2 )*scale;
    jj++;
    
    mdensity += H_1 * 1.00794;
    


    
    scale = data->scale[jj];
    y_arr[jj] = H_2 = NV_Ith_S( y,3 )*scale;
    jj++;
    
    mdensity += H_2 * 1.00794;
    


    
    scale = data->scale[jj];
    y_arr[jj] = H_m0 = NV_Ith_S( y,4 )*scale;
    jj++;
    
    mdensity += H_m0 * 1.00794;
    


    
    scale = data->scale[jj];
    y_arr[jj] = He_1 = NV_Ith_S( y,5 )*scale;
    jj++;
    
    mdensity += He_1 * 4.002602;
    


    
    scale = data->scale[jj];
    y_arr[jj] = He_2 = NV_Ith_S( y,6 )*scale;
    jj++;
    
    mdensity += He_2 * 4.002602;
    


    
    scale = data->scale[jj];
    y_arr[jj] = He_3 = NV_Ith_S( y,7 )*scale;
    jj++;
    
    mdensity += He_3 * 4.002602;
    


    
    scale = data->scale[jj];
    y_arr[jj] = de = NV_Ith_S( y,8 )*scale;
    jj++;
    


    
    scale = data->scale[jj];
    y_arr[jj] = ge = NV_Ith_S( y,9 )*scale;
    jj++;
    
    eq_H2_2_calculate_temperature(data, y_arr , nstrip, nchem );
    eq_H2_2_interpolate_rates(data, nstrip);


    
    
    mdensity *= mh;
    //
    // Species: H2_1
    //
    NV_Ith_S(ydot, 0) = k08[i]*H_1*H_m0 + k10[i]*H2_2*H_1 - k11[i]*H2_1*H_2 - k12[i]*H2_1*de - k13[i]*H2_1*H_1 + k19[i]*H2_2*H_m0 + k21[i]*H2_1*pow(H_1, 2) + k22[i]*pow(H_1, 3);
 
    inv_scale = data->inv_scale[0];
    NV_Ith_S(ydot, 0) *= inv_scale;

    
    
    //fprintf(stderr, "H2_1: %0.5g\n", scale);
    //fprintf(stderr, "ydot = %0.5g \n", NV_Ith_S(ydot, 0));
    
    //
    // Species: H2_2
    //
    NV_Ith_S(ydot, 1) = k09[i]*H_1*H_2 - k10[i]*H2_2*H_1 + k11[i]*H2_1*H_2 + k17[i]*H_2*H_m0 - k18[i]*H2_2*de - k19[i]*H2_2*H_m0;
 
    inv_scale = data->inv_scale[1];
    NV_Ith_S(ydot, 1) *= inv_scale;

    
    
    //fprintf(stderr, "H2_2: %0.5g\n", scale);
    //fprintf(stderr, "ydot = %0.5g \n", NV_Ith_S(ydot, 1));
    
    //
    // Species: H_1
    //
    NV_Ith_S(ydot, 2) = -k01[i]*H_1*de + k02[i]*H_2*de - k07[i]*H_1*de - k08[i]*H_1*H_m0 - k09[i]*H_1*H_2 - k10[i]*H2_2*H_1 + k11[i]*H2_1*H_2 + 2*k12[i]*H2_1*de + 2*k13[i]*H2_1*H_1 + k14[i]*H_m0*de + k15[i]*H_1*H_m0 + 2*k16[i]*H_2*H_m0 + 2*k18[i]*H2_2*de + k19[i]*H2_2*H_m0 - 2*k21[i]*H2_1*pow(H_1, 2) - 2*k22[i]*pow(H_1, 3);
 
    inv_scale = data->inv_scale[2];
    NV_Ith_S(ydot, 2) *= inv_scale;

    
    
    //fprintf(stderr, "H_1: %0.5g\n", scale);
    //fprintf(stderr, "ydot = %0.5g \n", NV_Ith_S(ydot, 2));
    
    //
    // Species: H_2
    //
    NV_Ith_S(ydot, 3) = k01[i]*H_1*de - k02[i]*H_2*de - k09[i]*H_1*H_2 + k10[i]*H2_2*H_1 - k11[i]*H2_1*H_2 - k16[i]*H_2*H_m0 - k17[i]*H_2*H_m0;
 
    inv_scale = data->inv_scale[3];
    NV_Ith_S(ydot, 3) *= inv_scale;

    
    
    //fprintf(stderr, "H_2: %0.5g\n", scale);
    //fprintf(stderr, "ydot = %0.5g \n", NV_Ith_S(ydot, 3));
    
    //
    // Species: H_m0
    //
    NV_Ith_S(ydot, 4) = k07[i]*H_1*de - k08[i]*H_1*H_m0 - k14[i]*H_m0*de - k15[i]*H_1*H_m0 - k16[i]*H_2*H_m0 - k17[i]*H_2*H_m0 - k19[i]*H2_2*H_m0;
 
    inv_scale = data->inv_scale[4];
    NV_Ith_S(ydot, 4) *= inv_scale;

    
    
    //fprintf(stderr, "H_m0: %0.5g\n", scale);
    //fprintf(stderr, "ydot = %0.5g \n", NV_Ith_S(ydot, 4));
    
    //
    // Species: He_1
    //
    NV_Ith_S(ydot, 5) = -k03[i]*He_1*de + k04[i]*He_2*de;
 
    inv_scale = data->inv_scale[5];
    NV_Ith_S(ydot, 5) *= inv_scale;

    
    
    //fprintf(stderr, "He_1: %0.5g\n", scale);
    //fprintf(stderr, "ydot = %0.5g \n", NV_Ith_S(ydot, 5));
    
    //
    // Species: He_2
    //
    NV_Ith_S(ydot, 6) = k03[i]*He_1*de - k04[i]*He_2*de - k05[i]*He_2*de + k06[i]*He_3*de;
 
    inv_scale = data->inv_scale[6];
    NV_Ith_S(ydot, 6) *= inv_scale;

    
    
    //fprintf(stderr, "He_2: %0.5g\n", scale);
    //fprintf(stderr, "ydot = %0.5g \n", NV_Ith_S(ydot, 6));
    
    //
    // Species: He_3
    //
    NV_Ith_S(ydot, 7) = k05[i]*He_2*de - k06[i]*He_3*de;
 
    inv_scale = data->inv_scale[7];
    NV_Ith_S(ydot, 7) *= inv_scale;

    
    
    //fprintf(stderr, "He_3: %0.5g\n", scale);
    //fprintf(stderr, "ydot = %0.5g \n", NV_Ith_S(ydot, 7));
    
    //
    // Species: de
    //
    NV_Ith_S(ydot, 8) = k01[i]*H_1*de - k02[i]*H_2*de + k03[i]*He_1*de - k04[i]*He_2*de + k05[i]*He_2*de - k06[i]*He_3*de - k07[i]*H_1*de + k08[i]*H_1*H_m0 + k14[i]*H_m0*de + k15[i]*H_1*H_m0 + k17[i]*H_2*H_m0 - k18[i]*H2_2*de;
 
    inv_scale = data->inv_scale[8];
    NV_Ith_S(ydot, 8) *= inv_scale;

    
    
    //fprintf(stderr, "de: %0.5g\n", scale);
    //fprintf(stderr, "ydot = %0.5g \n", NV_Ith_S(ydot, 8));
    
    //
    // Species: ge
    //
    NV_Ith_S(ydot, 9) = (-H2_1*gloverabel08_h2lte[i]/(gloverabel08_h2lte[i]/(H2_1*gloverabel08_gaH2[i] + H_1*gloverabel08_gaHI[i] + H_2*gloverabel08_gaHp[i] + He_1*gloverabel08_gaHe[i] + de*gloverabel08_gael[i]) + 1.0) - H_1*ceHI_ceHI[i]*de - H_1*ciHI_ciHI[i]*de - H_2*de*reHII_reHII[i] - He_1*ciHeI_ciHeI[i]*de - He_2*ceHeII_ceHeII[i]*de - He_2*ceHeI_ceHeI[i]*pow(de, 2) - He_2*ciHeII_ciHeII[i]*de - He_2*ciHeIS_ciHeIS[i]*pow(de, 2) - He_2*de*reHeII1_reHeII1[i] - He_2*de*reHeII2_reHeII2[i] - He_3*de*reHeIII_reHeIII[i] - brem_brem[i]*de*(H_2 + He_2 + 4.0*He_3) - compton_comp_[i]*de*pow(z + 1.0, 4)*(T - 2.73*z - 2.73) + 0.5*1.0/(h2formation_ncrn[i]/(H2_1*h2formation_ncrd2[i] + H_1*h2formation_ncrd1[i]) + 1.0)*(-H2_1*H_1*h2formation_h2mcool[i] + pow(H_1, 3)*h2formation_h2mheat[i]))*fmin(1.00000000000000, (1.0 - exp(-fmax(1.00000000000000e-5, mdensity)))/fmax(1.00000000000000e-5, mdensity));
 
    inv_scale = data->inv_scale[9];
    NV_Ith_S(ydot, 9) *= inv_scale;

    
    NV_Ith_S(ydot, 9) /= mdensity;
    
    
    //fprintf(stderr, "ge: %0.5g\n", scale);
    //fprintf(stderr, "ydot = %0.5g \n", NV_Ith_S(ydot, 9));
    
    
    //fprintf(stderr, "----------------\n");
    
    return 0;
    }




int eq_H2_2_solve_chemistry_dt( dengo_field_data *field_data, 
eq_H2_2_data *data, double dt ){
    
    int i, j;
    int N = 10;
    int dims = field_data->ncells; // total number of strips to be evaluated


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


    double *input = (double *) malloc(dims * N * sizeof(double));
    double *atol  = (double *) malloc(dims * N * sizeof(double));
    double *rtol  = (double *) malloc(dims * N * sizeof(double));
    double *temp  = (double *) malloc(dims * sizeof(double) );
    

    for ( int d = 0; d< dims; d++  ){
        j = d*N;
        // this should be the normalized 
        // by the input units later
        // atol = input * rtol;
        // which again should be set by dengo
        input[j] = field_data->H2_1_density[d];
        input[j] /= 2.01588;
        atol[j] = input[j] * 1.0e-9;
        j++;
        
        input[j] = field_data->H2_2_density[d];
        input[j] /= 2.01588;
        atol[j] = input[j] * 1.0e-9;
        j++;
        
        input[j] = field_data->H_1_density[d];
        input[j] /= 1.00794;
        atol[j] = input[j] * 1.0e-9;
        j++;
        
        input[j] = field_data->H_2_density[d];
        input[j] /= 1.00794;
        atol[j] = input[j] * 1.0e-9;
        j++;
        
        input[j] = field_data->H_m0_density[d];
        input[j] /= 1.00794;
        atol[j] = input[j] * 1.0e-9;
        j++;
        
        input[j] = field_data->He_1_density[d];
        input[j] /= 4.002602;
        atol[j] = input[j] * 1.0e-9;
        j++;
        
        input[j] = field_data->He_2_density[d];
        input[j] /= 4.002602;
        atol[j] = input[j] * 1.0e-9;
        j++;
        
        input[j] = field_data->He_3_density[d];
        input[j] /= 4.002602;
        atol[j] = input[j] * 1.0e-9;
        j++;
        
        input[j] = field_data->de_density[d];
        input[j] /= 1.0;
        atol[j] = input[j] * 1.0e-9;
        j++;
        
        input[j] = field_data->ge_density[d];
        input[j] /= 1.0;
        atol[j] = input[j] * 1.0e-9;
        j++;
        
    }

    double ttot;
    double z;

    ttot = dengo_evolve_eq_H2_2( dt, dt, z, input, rtol, atol, dims, data, temp );

    for ( int d = 0; d< dims; d++  ){
        j = d*N;
        // this should be the normalized 
        // by the input units later
        // atol = input * rtol;
        // which again should be set by dengo
        field_data->H2_1_density[d] = input[j] * 2.01588;
        j++;
        
        field_data->H2_2_density[d] = input[j] * 2.01588;
        j++;
        
        field_data->H_1_density[d] = input[j] * 1.00794;
        j++;
        
        field_data->H_2_density[d] = input[j] * 1.00794;
        j++;
        
        field_data->H_m0_density[d] = input[j] * 1.00794;
        j++;
        
        field_data->He_1_density[d] = input[j] * 4.002602;
        j++;
        
        field_data->He_2_density[d] = input[j] * 4.002602;
        j++;
        
        field_data->He_3_density[d] = input[j] * 4.002602;
        j++;
        
        field_data->de_density[d] = input[j] * 1.0;
        j++;
        
        field_data->ge_density[d] = input[j] * 1.0;
        j++;
        
    }


}






int read_init_data_to_dengo( dengo_field_data *field_data ){

    // this reads the initial abundances of the data from
    // a hdf5 file, and initialize a field_data object

    char const *filename;
    filename = "newIC.h5";

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
    free(atol);
    free(rtol);
    free(ics);
    
    return 1;

}
