
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


#include "sunny_solver.h"

sunny_data *sunny_setup_data(
    int *NumberOfFields, char ***FieldNames)
{
    int i;

    sunny_data *data = (sunny_data *) malloc(sizeof(sunny_data));

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
    
    sunny_read_rate_tables(data);
    fprintf(stderr, "Successfully read in rate tables.\n");

    sunny_read_cooling_tables(data);
    fprintf(stderr, "Successfully read in cooling rate tables.\n");

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


int sunny_main(int argc, char** argv)
{
    sunny_data *data = sunny_setup_data(NULL, NULL);

    /* Initial conditions */

    hid_t file_id = H5Fopen("sunny_initial_conditions.h5", H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file_id < 0) {fprintf(stderr, "Failed to open "
        "sunny_initial_conditions.h5 so dying.\n");
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
                    tics[j], atol[j]);
        }
    }
    i++;
    

    H5Fclose(file_id);

    double dtf = 3.1557e+13;
    double dt = -1.0;
    double z = -1.0;
    for (i = 0; i < dims * N; i++) input[i] = ics[i];
    double ttot;
    ttot = dengo_evolve_sunny(dtf, dt, z, input, rtol, atol, dims, data);

    /* Write results to HDF5 file */
    file_id = H5Fcreate("sunny_solution.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
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
 



double dengo_evolve_sunny (double dtf, double &dt, double z, double *input,
            double *rtol, double *atol, long long dims, sunny_data *data) {
    int i, j;
    hid_t file_id;
    /* fprintf(stderr, "  ncells = % 3i\n", (int) dims); */

    int N = 10;
    for (i = 0; i<dims; i++) {
      j = i * N;
        
          input[j] /= 2.01588 * 1.67e-24;
          atol[j] /= 2.01588 * 1.67e-24;
        
        j++;
      
        
          input[j] /= 2.01588 * 1.67e-24;
          atol[j] /= 2.01588 * 1.67e-24;
        
        j++;
      
        
          input[j] /= 1.00794 * 1.67e-24;
          atol[j] /= 1.00794 * 1.67e-24;
        
        j++;
      
        
          input[j] /= 1.00794 * 1.67e-24;
          atol[j] /= 1.00794 * 1.67e-24;
        
        j++;
      
        
          input[j] /= 1.00794 * 1.67e-24;
          atol[j] /= 1.00794 * 1.67e-24;
        
        j++;
      
        
          input[j] /= 4.002602 * 1.67e-24;
          atol[j] /= 4.002602 * 1.67e-24;
        
        j++;
      
        
          input[j] /= 4.002602 * 1.67e-24;
          atol[j] /= 4.002602 * 1.67e-24;
        
        j++;
      
        
          input[j] /= 4.002602 * 1.67e-24;
          atol[j] /= 4.002602 * 1.67e-24;
        
        j++;
      
        
          input[j] /= 1.0 * 1.67e-24;
          atol[j] /= 1.0 * 1.67e-24;
        
        j++;
      
        
        j++;
      
    }
    ensure_electron_consistency(input, dims, N);

    rhs_f f = calculate_rhs_sunny;
    jac_f jf = calculate_jacobian_sunny;
    if (dt < 0) dt = dtf / 1e5;
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
    while (ttot < dtf) {
        int rv = BE_chem_solve(f, jf, input, dt, rtol, atol, dims, N,
                               scale, (void *) data, u0, s, gu, Ju);
        /*
        fprintf(stderr, "Return value [%d]: %i.  %0.5g / %0.5g = %0.5g (%0.5g)\n",
                niter, rv, ttot, dtf, ttot/dtf, dt);
        fprintf(stderr, "Value[0] = %0.5g %0.5g\n",
                input[0], prev[0]);
        */
        for (i = 0; i < dims * N; i++) {
            if (input[i] < 0) {
                rv = 1;
                break;
            }
        }
        if (rv == 0) {
	    if (siter == 50000) break;
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
    /* fprintf(stderr, "End: %0.5g / %0.5g (%0.5g)\n",
       ttot, dtf, dtf-ttot); */
    for (i = 0; i<dims; i++) {
      j = i * N;
        
          input[j] *= 2.01588 * 1.67e-24;
          atol[j] *= 2.01588 * 1.67e-24;
        
        j++;
      
        
          input[j] *= 2.01588 * 1.67e-24;
          atol[j] *= 2.01588 * 1.67e-24;
        
        j++;
      
        
          input[j] *= 1.00794 * 1.67e-24;
          atol[j] *= 1.00794 * 1.67e-24;
        
        j++;
      
        
          input[j] *= 1.00794 * 1.67e-24;
          atol[j] *= 1.00794 * 1.67e-24;
        
        j++;
      
        
          input[j] *= 1.00794 * 1.67e-24;
          atol[j] *= 1.00794 * 1.67e-24;
        
        j++;
      
        
          input[j] *= 4.002602 * 1.67e-24;
          atol[j] *= 4.002602 * 1.67e-24;
        
        j++;
      
        
          input[j] *= 4.002602 * 1.67e-24;
          atol[j] *= 4.002602 * 1.67e-24;
        
        j++;
      
        
          input[j] *= 4.002602 * 1.67e-24;
          atol[j] *= 4.002602 * 1.67e-24;
        
        j++;
      
        
          input[j] *= 1.0 * 1.67e-24;
          atol[j] *= 1.0 * 1.67e-24;
        
        j++;
      
        
        j++;
      
    }
    return ttot;
}
 


void sunny_read_rate_tables(sunny_data *data)
{
    hid_t file_id = H5Fopen("sunny_tables.h5", H5F_ACC_RDONLY, H5P_DEFAULT);
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

void sunny_read_cooling_tables(sunny_data *data)
{

    hid_t file_id = H5Fopen("sunny_tables.h5", H5F_ACC_RDONLY, H5P_DEFAULT);
    /* Allocate the correct number of rate tables */

    H5Fclose(file_id);
}

 


void sunny_calculate_temperature(sunny_data *data,
                        double *input, int nstrip, int nchem)
{
    int i, j;
    double density;
    double kb = 1.3806504e-16; // Boltzmann constant [erg/K]
    double mh = 1.67e-24;
    double gamma = 5.e0/3.e0;
    
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

    for (i = 0; i<nstrip; i++) {
        j = i * nchem;
        H2_1 = input[j];
        /*fprintf(stderr, "H2_1[%d] = % 0.16g\n",
                i, H2_1);*/
        j++;
    
        H2_2 = input[j];
        /*fprintf(stderr, "H2_2[%d] = % 0.16g\n",
                i, H2_2);*/
        j++;
    
        H_1 = input[j];
        /*fprintf(stderr, "H_1[%d] = % 0.16g\n",
                i, H_1);*/
        j++;
    
        H_2 = input[j];
        /*fprintf(stderr, "H_2[%d] = % 0.16g\n",
                i, H_2);*/
        j++;
    
        H_m0 = input[j];
        /*fprintf(stderr, "H_m0[%d] = % 0.16g\n",
                i, H_m0);*/
        j++;
    
        He_1 = input[j];
        /*fprintf(stderr, "He_1[%d] = % 0.16g\n",
                i, He_1);*/
        j++;
    
        He_2 = input[j];
        /*fprintf(stderr, "He_2[%d] = % 0.16g\n",
                i, He_2);*/
        j++;
    
        He_3 = input[j];
        /*fprintf(stderr, "He_3[%d] = % 0.16g\n",
                i, He_3);*/
        j++;
    
        de = input[j];
        /*fprintf(stderr, "de[%d] = % 0.16g\n",
                i, de);*/
        j++;
    
        ge = input[j];
        /*fprintf(stderr, "ge[%d] = % 0.16g\n",
                i, ge);*/
        j++;
    
        density = 2.01588*H2_1 + 2.01588*H2_2 + 1.00794*H_1 + 1.00794*H_2 + 1.00794*H_m0 + 4.002602*He_1 + 4.002602*He_2 + 4.002602*He_3;
        data->Ts[i] = density*ge*mh/(kb*(H2_1/(gamma - 1.0) + H2_2/(gamma - 1.0) + H_1/(gamma - 1.0) + H_2/(gamma - 1.0) + H_m0/(gamma - 1.0) + He_1/(gamma - 1.0) + He_2/(gamma - 1.0) + He_3/(gamma - 1.0) + de/(gamma - 1.0)));
        if (data->Ts[i] < data->bounds[0]) {
            data->Ts[i] = data->bounds[0];
        } else if (data->Ts[i] > data->bounds[1]) {
            data->Ts[i] = data->bounds[1];
        }
        data->logTs[i] = log(data->Ts[i]);
        data->invTs[i] = 1.0 / data->Ts[i];
	data->dTs_ge[i] = 
        density*mh/(kb*(H2_1/(gamma - 1.0) + H2_2/(gamma - 1.0) + H_1/(gamma - 1.0) + H_2/(gamma - 1.0) + H_m0/(gamma - 1.0) + He_1/(gamma - 1.0) + He_2/(gamma - 1.0) + He_3/(gamma - 1.0) + de/(gamma - 1.0)));
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
void sunny_interpolate_rates(sunny_data *data,
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
    
    for (i = 0; i < nstrip; i++) {
        bin_id = data->bin_id[i];
        data->rs_k23[i] = data->r_k23[bin_id] +
            data->Tdef[i] * (data->r_k23[bin_id+1] - data->r_k23[bin_id]);
        data->drs_k23[i] = (data->r_k23[bin_id+1] - data->r_k23[bin_id]);
        data->drs_k23[i] /= data->dT[i];
	data->drs_k23[i] *= data->invTs[i];
    }
    

}
 



int calculate_rhs_sunny(double *input, double *rhs, int nstrip,
                  int nchem, void *sdata)
{
    /* We iterate over all of the rates */
    /* Calculate temperature first */
    sunny_data *data = (sunny_data*)sdata;
    int i, j;
    sunny_calculate_temperature(data, input, nstrip, nchem);

    sunny_interpolate_rates(data, nstrip);

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
    double *k23 = data->rs_k23;
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
    for (i = 0; i<nstrip; i++) {
        j = i * nchem;
        mdensity = 0.0;
        T = data->Ts[i];
        z = data->current_z;
        H2_1 = input[j];
        
        mdensity += H2_1;
        
        if (H2_1 < 0.0) {
            /* fprintf(stderr, "RNegative[%d][H2_1] = % 0.16g [%d]\n",
               i, H2_1, j); */
            return 1;
          H2_1 = 1e-20;
        }
        j++;
    
        H2_2 = input[j];
        
        mdensity += H2_2;
        
        if (H2_2 < 0.0) {
            /* fprintf(stderr, "RNegative[%d][H2_2] = % 0.16g [%d]\n",
               i, H2_2, j); */
            return 1;
          H2_2 = 1e-20;
        }
        j++;
    
        H_1 = input[j];
        
        mdensity += H_1;
        
        if (H_1 < 0.0) {
            /* fprintf(stderr, "RNegative[%d][H_1] = % 0.16g [%d]\n",
               i, H_1, j); */
            return 1;
          H_1 = 1e-20;
        }
        j++;
    
        H_2 = input[j];
        
        mdensity += H_2;
        
        if (H_2 < 0.0) {
            /* fprintf(stderr, "RNegative[%d][H_2] = % 0.16g [%d]\n",
               i, H_2, j); */
            return 1;
          H_2 = 1e-20;
        }
        j++;
    
        H_m0 = input[j];
        
        mdensity += H_m0;
        
        if (H_m0 < 0.0) {
            /* fprintf(stderr, "RNegative[%d][H_m0] = % 0.16g [%d]\n",
               i, H_m0, j); */
            return 1;
          H_m0 = 1e-20;
        }
        j++;
    
        He_1 = input[j];
        
        mdensity += He_1;
        
        if (He_1 < 0.0) {
            /* fprintf(stderr, "RNegative[%d][He_1] = % 0.16g [%d]\n",
               i, He_1, j); */
            return 1;
          He_1 = 1e-20;
        }
        j++;
    
        He_2 = input[j];
        
        mdensity += He_2;
        
        if (He_2 < 0.0) {
            /* fprintf(stderr, "RNegative[%d][He_2] = % 0.16g [%d]\n",
               i, He_2, j); */
            return 1;
          He_2 = 1e-20;
        }
        j++;
    
        He_3 = input[j];
        
        mdensity += He_3;
        
        if (He_3 < 0.0) {
            /* fprintf(stderr, "RNegative[%d][He_3] = % 0.16g [%d]\n",
               i, He_3, j); */
            return 1;
          He_3 = 1e-20;
        }
        j++;
    
        de = input[j];
        
        if (de < 0.0) {
            /* fprintf(stderr, "RNegative[%d][de] = % 0.16g [%d]\n",
               i, de, j); */
            return 1;
          de = 1e-20;
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
    
        mdensity *= mh;
        j = i * nchem;
        // 
        // Species: H2_1
        // 
        rhs[j] = k08[i]*H_1*H_m0 + k10[i]*H2_2*H_1 - k11[i]*H2_1*H_2 - k12[i]*H2_1*de - k13[i]*H2_1*H_1 + k19[i]*H2_2*H_m0 + k21[i]*H2_1*pow(H_1, 2) + k22[i]*pow(H_1, 3) - k23[i]*pow(H2_1, 2);
        
        j++;
    
        // 
        // Species: H2_2
        // 
        rhs[j] = k09[i]*H_1*H_2 - k10[i]*H2_2*H_1 + k11[i]*H2_1*H_2 + k17[i]*H_2*H_m0 - k18[i]*H2_2*de - k19[i]*H2_2*H_m0;
        
        j++;
    
        // 
        // Species: H_1
        // 
        rhs[j] = -k01[i]*H_1*de + k02[i]*H_2*de - k07[i]*H_1*de - k08[i]*H_1*H_m0 - k09[i]*H_1*H_2 - k10[i]*H2_2*H_1 + k11[i]*H2_1*H_2 + 2*k12[i]*H2_1*de + 2*k13[i]*H2_1*H_1 + k14[i]*H_m0*de + k15[i]*H_1*H_m0 + 2*k16[i]*H_2*H_m0 + 2*k18[i]*H2_2*de + k19[i]*H2_2*H_m0 - 2*k21[i]*H2_1*pow(H_1, 2) - 2*k22[i]*pow(H_1, 3) + 2*k23[i]*pow(H2_1, 2);
        
        j++;
    
        // 
        // Species: H_2
        // 
        rhs[j] = k01[i]*H_1*de - k02[i]*H_2*de - k09[i]*H_1*H_2 + k10[i]*H2_2*H_1 - k11[i]*H2_1*H_2 - k16[i]*H_2*H_m0 - k17[i]*H_2*H_m0;
        
        j++;
    
        // 
        // Species: H_m0
        // 
        rhs[j] = k07[i]*H_1*de - k08[i]*H_1*H_m0 - k14[i]*H_m0*de - k15[i]*H_1*H_m0 - k16[i]*H_2*H_m0 - k17[i]*H_2*H_m0 - k19[i]*H2_2*H_m0;
        
        j++;
    
        // 
        // Species: He_1
        // 
        rhs[j] = -k03[i]*He_1*de + k04[i]*He_2*de;
        
        j++;
    
        // 
        // Species: He_2
        // 
        rhs[j] = k03[i]*He_1*de - k04[i]*He_2*de - k05[i]*He_2*de + k06[i]*He_3*de;
        
        j++;
    
        // 
        // Species: He_3
        // 
        rhs[j] = k05[i]*He_2*de - k06[i]*He_3*de;
        
        j++;
    
        // 
        // Species: de
        // 
        rhs[j] = k01[i]*H_1*de - k02[i]*H_2*de + k03[i]*He_1*de - k04[i]*He_2*de + k05[i]*He_2*de - k06[i]*He_3*de - k07[i]*H_1*de + k08[i]*H_1*H_m0 + k14[i]*H_m0*de + k15[i]*H_1*H_m0 + k17[i]*H_2*H_m0 - k18[i]*H2_2*de;
        
        j++;
    
        // 
        // Species: ge
        // 
        rhs[j] = 0;
        
	rhs[j] /= mdensity;
        
        j++;
    
    }  
    return 0;
}




int calculate_jacobian_sunny(double *input, double *Joutput,
        int nstrip, int nchem, void *sdata)
{
    /* We iterate over all of the rates */
    /* Calculate temperature first */
    sunny_data *data = (sunny_data*)sdata;

    int i, j;
    /* Now we set up some temporaries */
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
    double *k23 = data->rs_k23;
    double *rk23 = data->drs_k23;
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
    for (i = 0; i<nstrip; i++) {
        j = i * nchem;
        mdensity = 0.0;
        T = data->Ts[i];
        z = data->current_z;
	H2_1 = input[j];
        
        mdensity += H2_1;
	
        if (H2_1 < 0.0) {
            fprintf(stderr, "JNegative[%d][H2_1] = % 0.16g [%d]\n",
                    i, H2_1, j);
            /*H2_1 = 0.0;*/
            H2_1 = 1e-20;
            return 1;
        }
        j++;
        
	H2_2 = input[j];
        
        mdensity += H2_2;
	
        if (H2_2 < 0.0) {
            fprintf(stderr, "JNegative[%d][H2_2] = % 0.16g [%d]\n",
                    i, H2_2, j);
            /*H2_2 = 0.0;*/
            H2_2 = 1e-20;
            return 1;
        }
        j++;
        
	H_1 = input[j];
        
        mdensity += H_1;
	
        if (H_1 < 0.0) {
            fprintf(stderr, "JNegative[%d][H_1] = % 0.16g [%d]\n",
                    i, H_1, j);
            /*H_1 = 0.0;*/
            H_1 = 1e-20;
            return 1;
        }
        j++;
        
	H_2 = input[j];
        
        mdensity += H_2;
	
        if (H_2 < 0.0) {
            fprintf(stderr, "JNegative[%d][H_2] = % 0.16g [%d]\n",
                    i, H_2, j);
            /*H_2 = 0.0;*/
            H_2 = 1e-20;
            return 1;
        }
        j++;
        
	H_m0 = input[j];
        
        mdensity += H_m0;
	
        if (H_m0 < 0.0) {
            fprintf(stderr, "JNegative[%d][H_m0] = % 0.16g [%d]\n",
                    i, H_m0, j);
            /*H_m0 = 0.0;*/
            H_m0 = 1e-20;
            return 1;
        }
        j++;
        
	He_1 = input[j];
        
        mdensity += He_1;
	
        if (He_1 < 0.0) {
            fprintf(stderr, "JNegative[%d][He_1] = % 0.16g [%d]\n",
                    i, He_1, j);
            /*He_1 = 0.0;*/
            He_1 = 1e-20;
            return 1;
        }
        j++;
        
	He_2 = input[j];
        
        mdensity += He_2;
	
        if (He_2 < 0.0) {
            fprintf(stderr, "JNegative[%d][He_2] = % 0.16g [%d]\n",
                    i, He_2, j);
            /*He_2 = 0.0;*/
            He_2 = 1e-20;
            return 1;
        }
        j++;
        
	He_3 = input[j];
        
        mdensity += He_3;
	
        if (He_3 < 0.0) {
            fprintf(stderr, "JNegative[%d][He_3] = % 0.16g [%d]\n",
                    i, He_3, j);
            /*He_3 = 0.0;*/
            He_3 = 1e-20;
            return 1;
        }
        j++;
        
	de = input[j];
        
        if (de < 0.0) {
            fprintf(stderr, "JNegative[%d][de] = % 0.16g [%d]\n",
                    i, de, j);
            /*de = 0.0;*/
            de = 1e-20;
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
        
        mdensity *= mh;
        
        j = i * nchem * nchem;
        // 
        // Species: H2_1
        //
            // H2_1 by H2_1
            Joutput[j] = -k11[i]*H_2 - k12[i]*de - k13[i]*H_1 + k21[i]*pow(H_1, 2) - 2*k23[i]*H2_1;
	    
	    
            j++;
            // H2_2 by H2_1
            Joutput[j] = k11[i]*H_2;
	    
	    
            j++;
            // H_1 by H2_1
            Joutput[j] = k11[i]*H_2 + 2*k12[i]*de + 2*k13[i]*H_1 - 2*k21[i]*pow(H_1, 2) + 4*k23[i]*H2_1;
	    
	    
            j++;
            // H_2 by H2_1
            Joutput[j] = -k11[i]*H_2;
	    
	    
            j++;
            // H_m0 by H2_1
            Joutput[j] = 0;
	    
	    
            j++;
            // He_1 by H2_1
            Joutput[j] = 0;
	    
	    
            j++;
            // He_2 by H2_1
            Joutput[j] = 0;
	    
	    
            j++;
            // He_3 by H2_1
            Joutput[j] = 0;
	    
	    
            j++;
            // de by H2_1
            Joutput[j] = 0;
	    
	    
            j++;
            // ge by H2_1
            Joutput[j] = 0;
	    
	    Joutput[j] /= mdensity;
	    
	    
            j++;
    
        // 
        // Species: H2_2
        //
            // H2_1 by H2_2
            Joutput[j] = k10[i]*H_1 + k19[i]*H_m0;
	    
	    
            j++;
            // H2_2 by H2_2
            Joutput[j] = -k10[i]*H_1 - k18[i]*de - k19[i]*H_m0;
	    
	    
            j++;
            // H_1 by H2_2
            Joutput[j] = -k10[i]*H_1 + 2*k18[i]*de + k19[i]*H_m0;
	    
	    
            j++;
            // H_2 by H2_2
            Joutput[j] = k10[i]*H_1;
	    
	    
            j++;
            // H_m0 by H2_2
            Joutput[j] = -k19[i]*H_m0;
	    
	    
            j++;
            // He_1 by H2_2
            Joutput[j] = 0;
	    
	    
            j++;
            // He_2 by H2_2
            Joutput[j] = 0;
	    
	    
            j++;
            // He_3 by H2_2
            Joutput[j] = 0;
	    
	    
            j++;
            // de by H2_2
            Joutput[j] = -k18[i]*de;
	    
	    
            j++;
            // ge by H2_2
            Joutput[j] = 0;
	    
	    Joutput[j] /= mdensity;
	    
	    
            j++;
    
        // 
        // Species: H_1
        //
            // H2_1 by H_1
            Joutput[j] = k08[i]*H_m0 + k10[i]*H2_2 - k13[i]*H2_1 + 2*k21[i]*H2_1*H_1 + 3*k22[i]*pow(H_1, 2);
	    
	    
            j++;
            // H2_2 by H_1
            Joutput[j] = k09[i]*H_2 - k10[i]*H2_2;
	    
	    
            j++;
            // H_1 by H_1
            Joutput[j] = -k01[i]*de - k07[i]*de - k08[i]*H_m0 - k09[i]*H_2 - k10[i]*H2_2 + 2*k13[i]*H2_1 + k15[i]*H_m0 - 4*k21[i]*H2_1*H_1 - 6*k22[i]*pow(H_1, 2);
	    
	    
            j++;
            // H_2 by H_1
            Joutput[j] = k01[i]*de - k09[i]*H_2 + k10[i]*H2_2;
	    
	    
            j++;
            // H_m0 by H_1
            Joutput[j] = k07[i]*de - k08[i]*H_m0 - k15[i]*H_m0;
	    
	    
            j++;
            // He_1 by H_1
            Joutput[j] = 0;
	    
	    
            j++;
            // He_2 by H_1
            Joutput[j] = 0;
	    
	    
            j++;
            // He_3 by H_1
            Joutput[j] = 0;
	    
	    
            j++;
            // de by H_1
            Joutput[j] = k01[i]*de - k07[i]*de + k08[i]*H_m0 + k15[i]*H_m0;
	    
	    
            j++;
            // ge by H_1
            Joutput[j] = 0;
	    
	    Joutput[j] /= mdensity;
	    
	    
            j++;
    
        // 
        // Species: H_2
        //
            // H2_1 by H_2
            Joutput[j] = -k11[i]*H2_1;
	    
	    
            j++;
            // H2_2 by H_2
            Joutput[j] = k09[i]*H_1 + k11[i]*H2_1 + k17[i]*H_m0;
	    
	    
            j++;
            // H_1 by H_2
            Joutput[j] = k02[i]*de - k09[i]*H_1 + k11[i]*H2_1 + 2*k16[i]*H_m0;
	    
	    
            j++;
            // H_2 by H_2
            Joutput[j] = -k02[i]*de - k09[i]*H_1 - k11[i]*H2_1 - k16[i]*H_m0 - k17[i]*H_m0;
	    
	    
            j++;
            // H_m0 by H_2
            Joutput[j] = -k16[i]*H_m0 - k17[i]*H_m0;
	    
	    
            j++;
            // He_1 by H_2
            Joutput[j] = 0;
	    
	    
            j++;
            // He_2 by H_2
            Joutput[j] = 0;
	    
	    
            j++;
            // He_3 by H_2
            Joutput[j] = 0;
	    
	    
            j++;
            // de by H_2
            Joutput[j] = -k02[i]*de + k17[i]*H_m0;
	    
	    
            j++;
            // ge by H_2
            Joutput[j] = 0;
	    
	    Joutput[j] /= mdensity;
	    
	    
            j++;
    
        // 
        // Species: H_m0
        //
            // H2_1 by H_m0
            Joutput[j] = k08[i]*H_1 + k19[i]*H2_2;
	    
	    
            j++;
            // H2_2 by H_m0
            Joutput[j] = k17[i]*H_2 - k19[i]*H2_2;
	    
	    
            j++;
            // H_1 by H_m0
            Joutput[j] = -k08[i]*H_1 + k14[i]*de + k15[i]*H_1 + 2*k16[i]*H_2 + k19[i]*H2_2;
	    
	    
            j++;
            // H_2 by H_m0
            Joutput[j] = -k16[i]*H_2 - k17[i]*H_2;
	    
	    
            j++;
            // H_m0 by H_m0
            Joutput[j] = -k08[i]*H_1 - k14[i]*de - k15[i]*H_1 - k16[i]*H_2 - k17[i]*H_2 - k19[i]*H2_2;
	    
	    
            j++;
            // He_1 by H_m0
            Joutput[j] = 0;
	    
	    
            j++;
            // He_2 by H_m0
            Joutput[j] = 0;
	    
	    
            j++;
            // He_3 by H_m0
            Joutput[j] = 0;
	    
	    
            j++;
            // de by H_m0
            Joutput[j] = k08[i]*H_1 + k14[i]*de + k15[i]*H_1 + k17[i]*H_2;
	    
	    
            j++;
            // ge by H_m0
            Joutput[j] = 0;
	    
	    Joutput[j] /= mdensity;
	    
	    
            j++;
    
        // 
        // Species: He_1
        //
            // H2_1 by He_1
            Joutput[j] = 0;
	    
	    
            j++;
            // H2_2 by He_1
            Joutput[j] = 0;
	    
	    
            j++;
            // H_1 by He_1
            Joutput[j] = 0;
	    
	    
            j++;
            // H_2 by He_1
            Joutput[j] = 0;
	    
	    
            j++;
            // H_m0 by He_1
            Joutput[j] = 0;
	    
	    
            j++;
            // He_1 by He_1
            Joutput[j] = -k03[i]*de;
	    
	    
            j++;
            // He_2 by He_1
            Joutput[j] = k03[i]*de;
	    
	    
            j++;
            // He_3 by He_1
            Joutput[j] = 0;
	    
	    
            j++;
            // de by He_1
            Joutput[j] = k03[i]*de;
	    
	    
            j++;
            // ge by He_1
            Joutput[j] = 0;
	    
	    Joutput[j] /= mdensity;
	    
	    
            j++;
    
        // 
        // Species: He_2
        //
            // H2_1 by He_2
            Joutput[j] = 0;
	    
	    
            j++;
            // H2_2 by He_2
            Joutput[j] = 0;
	    
	    
            j++;
            // H_1 by He_2
            Joutput[j] = 0;
	    
	    
            j++;
            // H_2 by He_2
            Joutput[j] = 0;
	    
	    
            j++;
            // H_m0 by He_2
            Joutput[j] = 0;
	    
	    
            j++;
            // He_1 by He_2
            Joutput[j] = k04[i]*de;
	    
	    
            j++;
            // He_2 by He_2
            Joutput[j] = -k04[i]*de - k05[i]*de;
	    
	    
            j++;
            // He_3 by He_2
            Joutput[j] = k05[i]*de;
	    
	    
            j++;
            // de by He_2
            Joutput[j] = -k04[i]*de + k05[i]*de;
	    
	    
            j++;
            // ge by He_2
            Joutput[j] = 0;
	    
	    Joutput[j] /= mdensity;
	    
	    
            j++;
    
        // 
        // Species: He_3
        //
            // H2_1 by He_3
            Joutput[j] = 0;
	    
	    
            j++;
            // H2_2 by He_3
            Joutput[j] = 0;
	    
	    
            j++;
            // H_1 by He_3
            Joutput[j] = 0;
	    
	    
            j++;
            // H_2 by He_3
            Joutput[j] = 0;
	    
	    
            j++;
            // H_m0 by He_3
            Joutput[j] = 0;
	    
	    
            j++;
            // He_1 by He_3
            Joutput[j] = 0;
	    
	    
            j++;
            // He_2 by He_3
            Joutput[j] = k06[i]*de;
	    
	    
            j++;
            // He_3 by He_3
            Joutput[j] = -k06[i]*de;
	    
	    
            j++;
            // de by He_3
            Joutput[j] = -k06[i]*de;
	    
	    
            j++;
            // ge by He_3
            Joutput[j] = 0;
	    
	    Joutput[j] /= mdensity;
	    
	    
            j++;
    
        // 
        // Species: de
        //
            // H2_1 by de
            Joutput[j] = -k12[i]*H2_1;
	    
	    
            j++;
            // H2_2 by de
            Joutput[j] = -k18[i]*H2_2;
	    
	    
            j++;
            // H_1 by de
            Joutput[j] = -k01[i]*H_1 + k02[i]*H_2 - k07[i]*H_1 + 2*k12[i]*H2_1 + k14[i]*H_m0 + 2*k18[i]*H2_2;
	    
	    
            j++;
            // H_2 by de
            Joutput[j] = k01[i]*H_1 - k02[i]*H_2;
	    
	    
            j++;
            // H_m0 by de
            Joutput[j] = k07[i]*H_1 - k14[i]*H_m0;
	    
	    
            j++;
            // He_1 by de
            Joutput[j] = -k03[i]*He_1 + k04[i]*He_2;
	    
	    
            j++;
            // He_2 by de
            Joutput[j] = k03[i]*He_1 - k04[i]*He_2 - k05[i]*He_2 + k06[i]*He_3;
	    
	    
            j++;
            // He_3 by de
            Joutput[j] = k05[i]*He_2 - k06[i]*He_3;
	    
	    
            j++;
            // de by de
            Joutput[j] = k01[i]*H_1 - k02[i]*H_2 + k03[i]*He_1 - k04[i]*He_2 + k05[i]*He_2 - k06[i]*He_3 - k07[i]*H_1 + k14[i]*H_m0 - k18[i]*H2_2;
	    
	    
            j++;
            // ge by de
            Joutput[j] = 0;
	    
	    Joutput[j] /= mdensity;
	    
	    
            j++;
    
        // 
        // Species: ge
        //
            // H2_1 by ge
            Joutput[j] = 0;
	    
	    
            Joutput[j] *= Tge[i];
            
            j++;
            // H2_2 by ge
            Joutput[j] = 0;
	    
	    
            Joutput[j] *= Tge[i];
            
            j++;
            // H_1 by ge
            Joutput[j] = 0;
	    
	    
            Joutput[j] *= Tge[i];
            
            j++;
            // H_2 by ge
            Joutput[j] = 0;
	    
	    
            Joutput[j] *= Tge[i];
            
            j++;
            // H_m0 by ge
            Joutput[j] = 0;
	    
	    
            Joutput[j] *= Tge[i];
            
            j++;
            // He_1 by ge
            Joutput[j] = 0;
	    
	    
            Joutput[j] *= Tge[i];
            
            j++;
            // He_2 by ge
            Joutput[j] = 0;
	    
	    
            Joutput[j] *= Tge[i];
            
            j++;
            // He_3 by ge
            Joutput[j] = 0;
	    
	    
            Joutput[j] *= Tge[i];
            
            j++;
            // de by ge
            Joutput[j] = 0;
	    
	    
            Joutput[j] *= Tge[i];
            
            j++;
            // ge by ge
            Joutput[j] = 0;
	    
	    Joutput[j] /= mdensity;
	    
	    
            Joutput[j] *= Tge[i];
            
            j++;
    
    }

    return 0;
    
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
        strip_temperature[i] = density*ge*mh/(kb*(H2_1/(gamma - 1.0) + H2_2/(gamma - 1.0) + H_1/(gamma - 1.0) + H_2/(gamma - 1.0) + H_m0/(gamma - 1.0) + He_1/(gamma - 1.0) + He_2/(gamma - 1.0) + He_3/(gamma - 1.0) + de/(gamma - 1.0)));
        if (strip_temperature[i] < 1.0)
            strip_temperature[i] = 1.0;
    }
         
}
 