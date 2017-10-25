
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


#include "sundials_solver.h"

sundials_data *sundials_setup_data(
    int *NumberOfFields, char ***FieldNames)
{
    int i;

    sundials_data *data = (sundials_data *) malloc(sizeof(sundials_data));

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
    
    sundials_read_rate_tables(data);
    fprintf(stderr, "Successfully read in rate tables.\n");

    sundials_read_cooling_tables(data);
    fprintf(stderr, "Successfully read in cooling rate tables.\n");

    if (FieldNames != NULL && NumberOfFields != NULL) {
        NumberOfFields[0] = 7;
        FieldNames[0] = new char*[7];
        i = 0;
        
        FieldNames[0][i++] = strdup("H2_1");
        
        FieldNames[0][i++] = strdup("H2_2");
        
        FieldNames[0][i++] = strdup("H_1");
        
        FieldNames[0][i++] = strdup("H_2");
        
        FieldNames[0][i++] = strdup("H_m0");
        
        FieldNames[0][i++] = strdup("de");
        
        FieldNames[0][i++] = strdup("ge");
        
    }
    return data;

}


int sundials_main(int argc, char** argv)
{
    sundials_data *data = sundials_setup_data(NULL, NULL);

    /* Initial conditions */

    hid_t file_id = H5Fopen("sundials_initial_conditions.h5", H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file_id < 0) {fprintf(stderr, "Failed to open "
        "sundials_initial_conditions.h5 so dying.\n");
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

    int N = 7;

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
    ttot = dengo_evolve_sundials(dtf, dt, z, input, rtol, atol, dims, data);

    /* Write results to HDF5 file */
    file_id = H5Fcreate("sundials_solution.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
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
 



double dengo_evolve_sundials (double dtf, double &dt, double z, double *input,
            double *rtol, double *atol, long long dims, sundials_data *data) {
    int i, j;
    hid_t file_id;
    /* fprintf(stderr, "  ncells = % 3i\n", (int) dims); */

    int N = 7;
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
      
        
          input[j] /= 1.0 * 1.67e-24;
          atol[j] /= 1.0 * 1.67e-24;
        
        j++;
      
        
        j++;
      
    }
    ensure_electron_consistency(input, dims, N);

    rhs_f f = calculate_rhs_sundials;
    jac_f jf = calculate_jacobian_sundials;
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
        

        /* f and jf are function for evaluating rhs and jac for the solver
        * input: {double array}
        * rtol : {double} scalar (relative tolerance)
        * atol : {double array} vector (absolue tolerance)
        * NSPECIES: {int}
        */
        /*int rv = cvodes_main_solver( f, jf, input, rtol ,  atol, NSPECIES, (void *) data, ttot , ttot + dt);*/
        
        int rv = 1;
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
      
        
          input[j] *= 1.0 * 1.67e-24;
          atol[j] *= 1.0 * 1.67e-24;
        
        j++;
      
        
        j++;
      
    }
    return ttot;
}
 


void sundials_read_rate_tables(sundials_data *data)
{
    hid_t file_id = H5Fopen("sundials_tables.h5", H5F_ACC_RDONLY, H5P_DEFAULT);
    /* Allocate the correct number of rate tables */
    H5LTread_dataset_double(file_id, "/k01", data->r_k01);
    H5LTread_dataset_double(file_id, "/k02", data->r_k02);
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

void sundials_read_cooling_tables(sundials_data *data)
{

    hid_t file_id = H5Fopen("sundials_tables.h5", H5F_ACC_RDONLY, H5P_DEFAULT);
    /* Allocate the correct number of rate tables */

    H5Fclose(file_id);
}

 


void sundials_calculate_temperature(sundials_data *data,
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
    
        de = input[j];
        /*fprintf(stderr, "de[%d] = % 0.16g\n",
                i, de);*/
        j++;
    
        ge = input[j];
        /*fprintf(stderr, "ge[%d] = % 0.16g\n",
                i, ge);*/
        j++;
    
        density = 2.01588*H2_1 + 2.01588*H2_2 + 1.00794*H_1 + 1.00794*H_2 + 1.00794*H_m0;
        data->Ts[i] = density*ge*mh/(kb*(H2_1/(gamma - 1.0) + H2_2/(gamma - 1.0) + H_1/(gamma - 1.0) + H_2/(gamma - 1.0) + H_m0/(gamma - 1.0) + de/(gamma - 1.0)));
        if (data->Ts[i] < data->bounds[0]) {
            data->Ts[i] = data->bounds[0];
        } else if (data->Ts[i] > data->bounds[1]) {
            data->Ts[i] = data->bounds[1];
        }
        data->logTs[i] = log(data->Ts[i]);
        data->invTs[i] = 1.0 / data->Ts[i];
	data->dTs_ge[i] = 
        density*mh/(kb*(H2_1/(gamma - 1.0) + H2_2/(gamma - 1.0) + H_1/(gamma - 1.0) + H_2/(gamma - 1.0) + H_m0/(gamma - 1.0) + de/(gamma - 1.0)));
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
void sundials_interpolate_rates(sundials_data *data,
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
 




void ensure_electron_consistency(double *input, int nstrip, int nchem)
{
    int i, j;

    /* Now we set up some temporaries */
    double H2_1;
    double H2_2;
    double H_1;
    double H_2;
    double H_m0;
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
    
        de = input[j];
        
        de /= 1.0 * mh;
        
        /*fprintf(stderr, "de[%d] = % 0.16g\n",
                i, de);*/
        j++;
    
        ge = input[j];
        
        /*fprintf(stderr, "ge[%d] = % 0.16g\n",
                i, ge);*/
        j++;
    
        density = 2.01588*H2_1 + 2.01588*H2_2 + 1.00794*H_1 + 1.00794*H_2 + 1.00794*H_m0;
        strip_temperature[i] = density*ge*mh/(kb*(H2_1/(gamma - 1.0) + H2_2/(gamma - 1.0) + H_1/(gamma - 1.0) + H_2/(gamma - 1.0) + H_m0/(gamma - 1.0) + de/(gamma - 1.0)));
        if (strip_temperature[i] < 1.0)
            strip_temperature[i] = 1.0;
    }
         
}
 


int calculate_jacobian_sundials
                                        ( long int N,  realtype t,
                                        N_Vector y, N_Vector fy,
                                        DlsMat J, void *user_data,
                                        N_Vector tmp1, N_Vector tmp2,
                                        N_Vector tmp3)
{
    /* We iterate over all of the rates */
    /* Calcuate temperature first */
    
    int nstrip = 1;
    int nchem = 7;

    sundials_data *data = (sundials_data*)user_data; 
    
    int i, j;
    /* Now We set up some temporaries */

    double *Tge = data->dTs_ge;
    double *k01 = data->rs_k01;
    double *rk01 = data->drs_k01;
    double *k02 = data->rs_k02;
    double *rk02 = data->drs_k02;
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
    double de;
    double ge;
    double z;
    double T;

    double mh = 1.67e-24;
    double mdensity;

    for (i = 0; i<nstrip; i++) {
        j = i*nchem;
        mdensity = 0.0;
        z = data->current_z;
        H2_1 = Ith( y, 0  );
        
        mdensity += H2_1;
        

        if (H2_1 < 0.0){
            fprintf(stderr, "JNegative[%d][H2_1] = %0.16g [%d]\n", i, H2_1, j);
            H2_1 = 1e-20;
            return 1;
        }
        j++;
        
        H2_2 = Ith( y, 1  );
        
        mdensity += H2_2;
        

        if (H2_2 < 0.0){
            fprintf(stderr, "JNegative[%d][H2_2] = %0.16g [%d]\n", i, H2_2, j);
            H2_2 = 1e-20;
            return 1;
        }
        j++;
        
        H_1 = Ith( y, 2  );
        
        mdensity += H_1;
        

        if (H_1 < 0.0){
            fprintf(stderr, "JNegative[%d][H_1] = %0.16g [%d]\n", i, H_1, j);
            H_1 = 1e-20;
            return 1;
        }
        j++;
        
        H_2 = Ith( y, 3  );
        
        mdensity += H_2;
        

        if (H_2 < 0.0){
            fprintf(stderr, "JNegative[%d][H_2] = %0.16g [%d]\n", i, H_2, j);
            H_2 = 1e-20;
            return 1;
        }
        j++;
        
        H_m0 = Ith( y, 4  );
        
        mdensity += H_m0;
        

        if (H_m0 < 0.0){
            fprintf(stderr, "JNegative[%d][H_m0] = %0.16g [%d]\n", i, H_m0, j);
            H_m0 = 1e-20;
            return 1;
        }
        j++;
        
        de = Ith( y, 5  );
        

        if (de < 0.0){
            fprintf(stderr, "JNegative[%d][de] = %0.16g [%d]\n", i, de, j);
            de = 1e-20;
            return 1;
        }
        j++;
        
        ge = Ith( y, 6  );
        

        if (ge < 0.0){
            fprintf(stderr, "JNegative[%d][ge] = %0.16g [%d]\n", i, ge, j);
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
            
            IJth(J, 1, 1 ) = -k11[i]*H_2 - k12[i]*de - k13[i]*H_1 + k21[i]*pow(H_1, 2) - 2*k23[i]*H2_1;
                

                
        // H2_2 by H2_1
            
            IJth(J, 2, 1 ) = k11[i]*H_2;
                

                
        // H_1 by H2_1
            
            IJth(J, 3, 1 ) = k11[i]*H_2 + 2*k12[i]*de + 2*k13[i]*H_1 - 2*k21[i]*pow(H_1, 2) + 4*k23[i]*H2_1;
                

                
        // H_2 by H2_1
            
            IJth(J, 4, 1 ) = -k11[i]*H_2;
                

                
        // H_m0 by H2_1
            
            IJth(J, 5, 1 ) = 0;
                

                
        // de by H2_1
            
            IJth(J, 6, 1 ) = 0;
                

                
        // ge by H2_1
            
            IJth(J, 7, 1 ) = 0;
                
                IJth(J, 7, 1 ) /= mdensity;
                

                
        
        //
        // Species: H2_2
        //
        
        // H2_1 by H2_2
            
            IJth(J, 1, 2 ) = k10[i]*H_1 + k19[i]*H_m0;
                

                
        // H2_2 by H2_2
            
            IJth(J, 2, 2 ) = -k10[i]*H_1 - k18[i]*de - k19[i]*H_m0;
                

                
        // H_1 by H2_2
            
            IJth(J, 3, 2 ) = -k10[i]*H_1 + 2*k18[i]*de + k19[i]*H_m0;
                

                
        // H_2 by H2_2
            
            IJth(J, 4, 2 ) = k10[i]*H_1;
                

                
        // H_m0 by H2_2
            
            IJth(J, 5, 2 ) = -k19[i]*H_m0;
                

                
        // de by H2_2
            
            IJth(J, 6, 2 ) = -k18[i]*de;
                

                
        // ge by H2_2
            
            IJth(J, 7, 2 ) = 0;
                
                IJth(J, 7, 2 ) /= mdensity;
                

                
        
        //
        // Species: H_1
        //
        
        // H2_1 by H_1
            
            IJth(J, 1, 3 ) = k08[i]*H_m0 + k10[i]*H2_2 - k13[i]*H2_1 + 2*k21[i]*H2_1*H_1 + 3*k22[i]*pow(H_1, 2);
                

                
        // H2_2 by H_1
            
            IJth(J, 2, 3 ) = k09[i]*H_2 - k10[i]*H2_2;
                

                
        // H_1 by H_1
            
            IJth(J, 3, 3 ) = -k01[i]*de - k07[i]*de - k08[i]*H_m0 - k09[i]*H_2 - k10[i]*H2_2 + 2*k13[i]*H2_1 + k15[i]*H_m0 - 4*k21[i]*H2_1*H_1 - 6*k22[i]*pow(H_1, 2);
                

                
        // H_2 by H_1
            
            IJth(J, 4, 3 ) = k01[i]*de - k09[i]*H_2 + k10[i]*H2_2;
                

                
        // H_m0 by H_1
            
            IJth(J, 5, 3 ) = k07[i]*de - k08[i]*H_m0 - k15[i]*H_m0;
                

                
        // de by H_1
            
            IJth(J, 6, 3 ) = k01[i]*de - k07[i]*de + k08[i]*H_m0 + k15[i]*H_m0;
                

                
        // ge by H_1
            
            IJth(J, 7, 3 ) = 0;
                
                IJth(J, 7, 3 ) /= mdensity;
                

                
        
        //
        // Species: H_2
        //
        
        // H2_1 by H_2
            
            IJth(J, 1, 4 ) = -k11[i]*H2_1;
                

                
        // H2_2 by H_2
            
            IJth(J, 2, 4 ) = k09[i]*H_1 + k11[i]*H2_1 + k17[i]*H_m0;
                

                
        // H_1 by H_2
            
            IJth(J, 3, 4 ) = k02[i]*de - k09[i]*H_1 + k11[i]*H2_1 + 2*k16[i]*H_m0;
                

                
        // H_2 by H_2
            
            IJth(J, 4, 4 ) = -k02[i]*de - k09[i]*H_1 - k11[i]*H2_1 - k16[i]*H_m0 - k17[i]*H_m0;
                

                
        // H_m0 by H_2
            
            IJth(J, 5, 4 ) = -k16[i]*H_m0 - k17[i]*H_m0;
                

                
        // de by H_2
            
            IJth(J, 6, 4 ) = -k02[i]*de + k17[i]*H_m0;
                

                
        // ge by H_2
            
            IJth(J, 7, 4 ) = 0;
                
                IJth(J, 7, 4 ) /= mdensity;
                

                
        
        //
        // Species: H_m0
        //
        
        // H2_1 by H_m0
            
            IJth(J, 1, 5 ) = k08[i]*H_1 + k19[i]*H2_2;
                

                
        // H2_2 by H_m0
            
            IJth(J, 2, 5 ) = k17[i]*H_2 - k19[i]*H2_2;
                

                
        // H_1 by H_m0
            
            IJth(J, 3, 5 ) = -k08[i]*H_1 + k14[i]*de + k15[i]*H_1 + 2*k16[i]*H_2 + k19[i]*H2_2;
                

                
        // H_2 by H_m0
            
            IJth(J, 4, 5 ) = -k16[i]*H_2 - k17[i]*H_2;
                

                
        // H_m0 by H_m0
            
            IJth(J, 5, 5 ) = -k08[i]*H_1 - k14[i]*de - k15[i]*H_1 - k16[i]*H_2 - k17[i]*H_2 - k19[i]*H2_2;
                

                
        // de by H_m0
            
            IJth(J, 6, 5 ) = k08[i]*H_1 + k14[i]*de + k15[i]*H_1 + k17[i]*H_2;
                

                
        // ge by H_m0
            
            IJth(J, 7, 5 ) = 0;
                
                IJth(J, 7, 5 ) /= mdensity;
                

                
        
        //
        // Species: de
        //
        
        // H2_1 by de
            
            IJth(J, 1, 6 ) = -k12[i]*H2_1;
                

                
        // H2_2 by de
            
            IJth(J, 2, 6 ) = -k18[i]*H2_2;
                

                
        // H_1 by de
            
            IJth(J, 3, 6 ) = -k01[i]*H_1 + k02[i]*H_2 - k07[i]*H_1 + 2*k12[i]*H2_1 + k14[i]*H_m0 + 2*k18[i]*H2_2;
                

                
        // H_2 by de
            
            IJth(J, 4, 6 ) = k01[i]*H_1 - k02[i]*H_2;
                

                
        // H_m0 by de
            
            IJth(J, 5, 6 ) = k07[i]*H_1 - k14[i]*H_m0;
                

                
        // de by de
            
            IJth(J, 6, 6 ) = k01[i]*H_1 - k02[i]*H_2 - k07[i]*H_1 + k14[i]*H_m0 - k18[i]*H2_2;
                

                
        // ge by de
            
            IJth(J, 7, 6 ) = 0;
                
                IJth(J, 7, 6 ) /= mdensity;
                

                
        
        //
        // Species: ge
        //
        
        // H2_1 by ge
            
            IJth(J, 1, 7 ) = 0;
                

                
                IJth(J, 1, 7 ) *= Tge[i];
                
        // H2_2 by ge
            
            IJth(J, 2, 7 ) = 0;
                

                
                IJth(J, 2, 7 ) *= Tge[i];
                
        // H_1 by ge
            
            IJth(J, 3, 7 ) = 0;
                

                
                IJth(J, 3, 7 ) *= Tge[i];
                
        // H_2 by ge
            
            IJth(J, 4, 7 ) = 0;
                

                
                IJth(J, 4, 7 ) *= Tge[i];
                
        // H_m0 by ge
            
            IJth(J, 5, 7 ) = 0;
                

                
                IJth(J, 5, 7 ) *= Tge[i];
                
        // de by ge
            
            IJth(J, 6, 7 ) = 0;
                

                
                IJth(J, 6, 7 ) *= Tge[i];
                
        // ge by ge
            
            IJth(J, 7, 7 ) = 0;
                
                IJth(J, 7, 7 ) /= mdensity;
                

                
                IJth(J, 7, 7 ) *= Tge[i];
                
        
    }
    return 0;
}






int calculate_rhs_sundials(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
    sundials_data *data = (sundials_data* ) user_data;
    int i, j;

    int nchem = 7;
    int nstrip = 1;
    
    /* change N_Vector back to an array */
    double y_arr[ 7 ];
    y_arr[0] = Ith(y , 0);
    y_arr[1] = Ith(y , 1);
    y_arr[2] = Ith(y , 2);
    y_arr[3] = Ith(y , 3);
    y_arr[4] = Ith(y , 4);
    y_arr[5] = Ith(y , 5);
    y_arr[6] = Ith(y , 6);

    sundials_calculate_temperature(data, y_arr , 1, nchem );
    sundials_interpolate_rates(data, nstrip);


    /* Now we set up some temporaries */
    double *k01 = data->rs_k01;
    double *k02 = data->rs_k02;
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
    double de;
    double ge;

    double z;
    double T;

    double mh = 1.67e-24;
    double mdensity;

    T = data->Ts[i];
    z = data->current_z;
    H2_1 = Ith( y,0 );
    
    mdensity += H2_1;
    

    if (H2_1 < 0.0 ){
        fprintf(stderr,"JNegative[%d][H2_1] = %0.16g [%d]\n", i, H2_1, j);
        H2_1 = 1e-20;

        return 1;

    }

    
    H2_2 = Ith( y,1 );
    
    mdensity += H2_2;
    

    if (H2_2 < 0.0 ){
        fprintf(stderr,"JNegative[%d][H2_2] = %0.16g [%d]\n", i, H2_2, j);
        H2_2 = 1e-20;

        return 1;

    }

    
    H_1 = Ith( y,2 );
    
    mdensity += H_1;
    

    if (H_1 < 0.0 ){
        fprintf(stderr,"JNegative[%d][H_1] = %0.16g [%d]\n", i, H_1, j);
        H_1 = 1e-20;

        return 1;

    }

    
    H_2 = Ith( y,3 );
    
    mdensity += H_2;
    

    if (H_2 < 0.0 ){
        fprintf(stderr,"JNegative[%d][H_2] = %0.16g [%d]\n", i, H_2, j);
        H_2 = 1e-20;

        return 1;

    }

    
    H_m0 = Ith( y,4 );
    
    mdensity += H_m0;
    

    if (H_m0 < 0.0 ){
        fprintf(stderr,"JNegative[%d][H_m0] = %0.16g [%d]\n", i, H_m0, j);
        H_m0 = 1e-20;

        return 1;

    }

    
    de = Ith( y,5 );
    

    if (de < 0.0 ){
        fprintf(stderr,"JNegative[%d][de] = %0.16g [%d]\n", i, de, j);
        de = 1e-20;

        return 1;

    }

    
    ge = Ith( y,6 );
    

    if (ge < 0.0 ){
        fprintf(stderr,"JNegative[%d][ge] = %0.16g [%d]\n", i, ge, j);
        ge = 1e-20;

        return 1;

    }

    

    mdensity *= mh;
    //
    // Species: H2_1
    //
    Ith(ydot, 1) = k08[i]*H_1*H_m0 + k10[i]*H2_2*H_1 - k11[i]*H2_1*H_2 - k12[i]*H2_1*de - k13[i]*H2_1*H_1 + k19[i]*H2_2*H_m0 + k21[i]*H2_1*pow(H_1, 2) + k22[i]*pow(H_1, 3) - k23[i]*pow(H2_1, 2);
    
    
    //
    // Species: H2_2
    //
    Ith(ydot, 2) = k09[i]*H_1*H_2 - k10[i]*H2_2*H_1 + k11[i]*H2_1*H_2 + k17[i]*H_2*H_m0 - k18[i]*H2_2*de - k19[i]*H2_2*H_m0;
    
    
    //
    // Species: H_1
    //
    Ith(ydot, 3) = -k01[i]*H_1*de + k02[i]*H_2*de - k07[i]*H_1*de - k08[i]*H_1*H_m0 - k09[i]*H_1*H_2 - k10[i]*H2_2*H_1 + k11[i]*H2_1*H_2 + 2*k12[i]*H2_1*de + 2*k13[i]*H2_1*H_1 + k14[i]*H_m0*de + k15[i]*H_1*H_m0 + 2*k16[i]*H_2*H_m0 + 2*k18[i]*H2_2*de + k19[i]*H2_2*H_m0 - 2*k21[i]*H2_1*pow(H_1, 2) - 2*k22[i]*pow(H_1, 3) + 2*k23[i]*pow(H2_1, 2);
    
    
    //
    // Species: H_2
    //
    Ith(ydot, 4) = k01[i]*H_1*de - k02[i]*H_2*de - k09[i]*H_1*H_2 + k10[i]*H2_2*H_1 - k11[i]*H2_1*H_2 - k16[i]*H_2*H_m0 - k17[i]*H_2*H_m0;
    
    
    //
    // Species: H_m0
    //
    Ith(ydot, 5) = k07[i]*H_1*de - k08[i]*H_1*H_m0 - k14[i]*H_m0*de - k15[i]*H_1*H_m0 - k16[i]*H_2*H_m0 - k17[i]*H_2*H_m0 - k19[i]*H2_2*H_m0;
    
    
    //
    // Species: de
    //
    Ith(ydot, 6) = k01[i]*H_1*de - k02[i]*H_2*de - k07[i]*H_1*de + k08[i]*H_1*H_m0 + k14[i]*H_m0*de + k15[i]*H_1*H_m0 + k17[i]*H_2*H_m0 - k18[i]*H2_2*de;
    
    
    //
    // Species: ge
    //
    Ith(ydot, 7) = 0;
    
    Ith(ydot, 7) /= mdensity;
    
    
    
    return 0;
}
