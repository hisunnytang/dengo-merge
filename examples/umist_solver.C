
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


#include "umist_solver.h"

umist_data *umist_setup_data(
    int *NumberOfFields, char ***FieldNames)
{
    int i;

    umist_data *data = (umist_data *) malloc(sizeof(umist_data));

    /* Temperature-related pieces */
    data->bounds[0] = 10.0;
    data->bounds[1] = 1000.0;
    data->nbins = 1024 - 1;
    data->dbin = (log(data->bounds[1]) - log(data->bounds[0])) / data->nbins;
    data->idbin = 1.0L / data->dbin;

    /* Redshift-related pieces */
    data->z_bounds[0] = 0.0;
    data->z_bounds[1] = 0.0;
    data->n_zbins = 0 - 1;
    data->d_zbin = (log(data->z_bounds[1] + 1.0) - log(data->z_bounds[0] + 1.0)) / data->n_zbins;
    data->id_zbin = 1.0L / data->d_zbin;
    
    umist_read_rate_tables(data);
    fprintf(stderr, "Successfully read in rate tables.\n");

    umist_read_cooling_tables(data);
    fprintf(stderr, "Successfully read in cooling rate tables.\n");

    if (FieldNames != NULL && NumberOfFields != NULL) {
        NumberOfFields[0] = 3;
        FieldNames[0] = new char*[3];
        i = 0;
        
        FieldNames[0][i++] = strdup("us_H_1");
        
        FieldNames[0][i++] = strdup("us_H2_1");
        
        FieldNames[0][i++] = strdup("ge");
        
    }
    return data;

}


int umist_main(int argc, char** argv)
{
    umist_data *data = umist_setup_data(NULL, NULL);

    /* Initial conditions */

    hid_t file_id = H5Fopen("umist_initial_conditions.h5", H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file_id < 0) {fprintf(stderr, "Failed to open "
        "umist_initial_conditions.h5 so dying.\n");
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

    int N = 3;

    double *atol, *rtol;
    atol = (double *) alloca(N * dims * sizeof(double));
    rtol = (double *) alloca(N * dims * sizeof(double));

    double *tics = (double *) alloca(dims * sizeof(double));
    double *ics = (double *) alloca(dims * N * sizeof(double));
    double *input = (double *) alloca(dims * N * sizeof(double));
    
    unsigned int i = 0, j;
    
    fprintf(stderr, "Reading I.C. for /us_H_1\n");
    H5LTread_dataset_double(file_id, "/us_H_1", tics);
    for (j = 0; j < dims; j++) {
        ics[j * N + i] = tics[j]; 
        atol[j * N + i] = tics[j] * 1e-09;
        rtol[j * N + i] = 1e-09;
        if(j==0) {
            fprintf(stderr, "us_H_1[0] = %0.3g, atol => % 0.16g\n",
                    tics[j], atol[j]);
        }
    }
    i++;
    
    fprintf(stderr, "Reading I.C. for /us_H2_1\n");
    H5LTread_dataset_double(file_id, "/us_H2_1", tics);
    for (j = 0; j < dims; j++) {
        ics[j * N + i] = tics[j]; 
        atol[j * N + i] = tics[j] * 1e-09;
        rtol[j * N + i] = 1e-09;
        if(j==0) {
            fprintf(stderr, "us_H2_1[0] = %0.3g, atol => % 0.16g\n",
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
    ttot = dengo_evolve_umist(dtf, dt, z, input, rtol, atol, dims, data);

    /* Write results to HDF5 file */
    file_id = H5Fcreate("umist_solution.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    hsize_t dimsarr[1];
    dimsarr[0] = dims;
    i = 0;
    
    double us_H_1[dims];
    for (j = 0; j < dims; j++) {
        us_H_1[j] = input[j * N + i]; 
    }
    fprintf(stderr, "Writing solution for /us_H_1\n");
    H5LTmake_dataset_double(file_id, "/us_H_1", 1, dimsarr, us_H_1);
    i++;
    
    double us_H2_1[dims];
    for (j = 0; j < dims; j++) {
        us_H2_1[j] = input[j * N + i]; 
    }
    fprintf(stderr, "Writing solution for /us_H2_1\n");
    H5LTmake_dataset_double(file_id, "/us_H2_1", 1, dimsarr, us_H2_1);
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
 



double dengo_evolve_umist (double dtf, double &dt, double z, double *input,
            double *rtol, double *atol, long long dims, umist_data *data) {
    int i, j;
    hid_t file_id;
    /* fprintf(stderr, "  ncells = % 3i\n", (int) dims); */

    int N = 3;
    for (i = 0; i<dims; i++) {
      j = i * N;
        
          input[j] /= 1 * 1.67e-24;
          atol[j] /= 1 * 1.67e-24;
        
        j++;
      
        
          input[j] /= 2 * 1.67e-24;
          atol[j] /= 2 * 1.67e-24;
        
        j++;
      
        
        j++;
      
    }
    ensure_electron_consistency(input, dims, N);

    rhs_f f = calculate_rhs_umist;
    jac_f jf = calculate_jacobian_umist;
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
        
          input[j] *= 1 * 1.67e-24;
          atol[j] *= 1 * 1.67e-24;
        
        j++;
      
        
          input[j] *= 2 * 1.67e-24;
          atol[j] *= 2 * 1.67e-24;
        
        j++;
      
        
        j++;
      
    }
    return ttot;
}
 


void umist_read_rate_tables(umist_data *data)
{
    hid_t file_id = H5Fopen("umist_tables.h5", H5F_ACC_RDONLY, H5P_DEFAULT);
    /* Allocate the correct number of rate tables */
    H5LTread_dataset_double(file_id, "/us_H2_1_plus_us_H2_1", data->r_us_H2_1_plus_us_H2_1);
    H5LTread_dataset_double(file_id, "/us_H_1_plus_us_H2_1", data->r_us_H_1_plus_us_H2_1);

    H5Fclose(file_id);
}

void umist_read_cooling_tables(umist_data *data)
{

    hid_t file_id = H5Fopen("umist_tables.h5", H5F_ACC_RDONLY, H5P_DEFAULT);
    /* Allocate the correct number of rate tables */

    H5Fclose(file_id);
}

 


void umist_calculate_temperature(umist_data *data,
                        double *input, int nstrip, int nchem)
{
    int i, j;
    double density;
    double kb = 1.3806504e-16; // Boltzmann constant [erg/K]
    double mh = 1.67e-24;
    double gamma = 5.e0/3.e0;
    
    /* Calculate total density */
    double us_H_1;
    double us_H2_1;
    double ge;

    for (i = 0; i<nstrip; i++) {
        j = i * nchem;
        us_H_1 = input[j];
        /*fprintf(stderr, "us_H_1[%d] = % 0.16g\n",
                i, us_H_1);*/
        j++;
    
        us_H2_1 = input[j];
        /*fprintf(stderr, "us_H2_1[%d] = % 0.16g\n",
                i, us_H2_1);*/
        j++;
    
        ge = input[j];
        /*fprintf(stderr, "ge[%d] = % 0.16g\n",
                i, ge);*/
        j++;
    
        density = 2*us_H2_1 + us_H_1;
        data->Ts[i] = density*ge*mh/(kb*(us_H2_1/(gamma - 1.0) + us_H_1/(gamma - 1.0)));
        if (data->Ts[i] < data->bounds[0]) {
            data->Ts[i] = data->bounds[0];
        } else if (data->Ts[i] > data->bounds[1]) {
            data->Ts[i] = data->bounds[1];
        }
        data->logTs[i] = log(data->Ts[i]);
        data->invTs[i] = 1.0 / data->Ts[i];
	data->dTs_ge[i] = 
        density*mh/(kb*(us_H2_1/(gamma - 1.0) + us_H_1/(gamma - 1.0)));
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
void umist_interpolate_rates(umist_data *data,
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
    
    if ((data->current_z > data->z_bounds[0]) && (data->current_z < data->z_bounds[1])) {
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
        data->rs_us_H2_1_plus_us_H2_1[i] = data->r_us_H2_1_plus_us_H2_1[bin_id] +
            data->Tdef[i] * (data->r_us_H2_1_plus_us_H2_1[bin_id+1] - data->r_us_H2_1_plus_us_H2_1[bin_id]);
        data->drs_us_H2_1_plus_us_H2_1[i] = (data->r_us_H2_1_plus_us_H2_1[bin_id+1] - data->r_us_H2_1_plus_us_H2_1[bin_id]);
        data->drs_us_H2_1_plus_us_H2_1[i] /= data->dT[i];
	data->drs_us_H2_1_plus_us_H2_1[i] *= data->invTs[i];
    }
    
    for (i = 0; i < nstrip; i++) {
        bin_id = data->bin_id[i];
        data->rs_us_H_1_plus_us_H2_1[i] = data->r_us_H_1_plus_us_H2_1[bin_id] +
            data->Tdef[i] * (data->r_us_H_1_plus_us_H2_1[bin_id+1] - data->r_us_H_1_plus_us_H2_1[bin_id]);
        data->drs_us_H_1_plus_us_H2_1[i] = (data->r_us_H_1_plus_us_H2_1[bin_id+1] - data->r_us_H_1_plus_us_H2_1[bin_id]);
        data->drs_us_H_1_plus_us_H2_1[i] /= data->dT[i];
	data->drs_us_H_1_plus_us_H2_1[i] *= data->invTs[i];
    }
    

}
 



int calculate_rhs_umist(double *input, double *rhs, int nstrip,
                  int nchem, void *sdata)
{
    /* We iterate over all of the rates */
    /* Calculate temperature first */
    umist_data *data = (umist_data*)sdata;
    int i, j;
    umist_calculate_temperature(data, input, nstrip, nchem);

    umist_interpolate_rates(data, nstrip);

    /* Now we set up some temporaries */
    double *us_H2_1_plus_us_H2_1 = data->rs_us_H2_1_plus_us_H2_1;
    double *us_H_1_plus_us_H2_1 = data->rs_us_H_1_plus_us_H2_1;
    double us_H_1;
    double us_H2_1;
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
        us_H_1 = input[j];
        
        mdensity += us_H_1;
        
        if (us_H_1 < 0.0) {
            /* fprintf(stderr, "RNegative[%d][us_H_1] = % 0.16g [%d]\n",
               i, us_H_1, j); */
            return 1;
          us_H_1 = 1e-20;
        }
        j++;
    
        us_H2_1 = input[j];
        
        mdensity += us_H2_1;
        
        if (us_H2_1 < 0.0) {
            /* fprintf(stderr, "RNegative[%d][us_H2_1] = % 0.16g [%d]\n",
               i, us_H2_1, j); */
            return 1;
          us_H2_1 = 1e-20;
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
        // Species: us_H_1
        // 
        rhs[j] = 2*us_H2_1_plus_us_H2_1[i]*pow(us_H2_1, 2) + 2*us_H_1_plus_us_H2_1[i]*us_H2_1*us_H_1;
        
        j++;
    
        // 
        // Species: us_H2_1
        // 
        rhs[j] = -us_H2_1_plus_us_H2_1[i]*pow(us_H2_1, 2) - us_H_1_plus_us_H2_1[i]*us_H2_1*us_H_1;
        
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




int calculate_jacobian_umist(double *input, double *Joutput,
        int nstrip, int nchem, void *sdata)
{
    /* We iterate over all of the rates */
    /* Calculate temperature first */
    umist_data *data = (umist_data*)sdata;

    int i, j;
    /* Now we set up some temporaries */
    double *Tge = data->dTs_ge;
    double *us_H2_1_plus_us_H2_1 = data->rs_us_H2_1_plus_us_H2_1;
    double *rus_H2_1_plus_us_H2_1 = data->drs_us_H2_1_plus_us_H2_1;
    double *us_H_1_plus_us_H2_1 = data->rs_us_H_1_plus_us_H2_1;
    double *rus_H_1_plus_us_H2_1 = data->drs_us_H_1_plus_us_H2_1;
    double us_H_1;
    double us_H2_1;
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
	us_H_1 = input[j];
        
        mdensity += us_H_1;
	
        if (us_H_1 < 0.0) {
            fprintf(stderr, "JNegative[%d][us_H_1] = % 0.16g [%d]\n",
                    i, us_H_1, j);
            /*us_H_1 = 0.0;*/
            us_H_1 = 1e-20;
            return 1;
        }
        j++;
        
	us_H2_1 = input[j];
        
        mdensity += us_H2_1;
	
        if (us_H2_1 < 0.0) {
            fprintf(stderr, "JNegative[%d][us_H2_1] = % 0.16g [%d]\n",
                    i, us_H2_1, j);
            /*us_H2_1 = 0.0;*/
            us_H2_1 = 1e-20;
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
        // Species: us_H_1
        //
            // us_H_1 by us_H_1
            Joutput[j] = 2*us_H_1_plus_us_H2_1[i]*us_H2_1;
	    
	    
            j++;
            // us_H2_1 by us_H_1
            Joutput[j] = -us_H_1_plus_us_H2_1[i]*us_H2_1;
	    
	    
            j++;
            // ge by us_H_1
            Joutput[j] = 0;
	    
	    Joutput[j] /= mdensity;
	    
	    
            j++;
    
        // 
        // Species: us_H2_1
        //
            // us_H_1 by us_H2_1
            Joutput[j] = 4*us_H2_1_plus_us_H2_1[i]*us_H2_1 + 2*us_H_1_plus_us_H2_1[i]*us_H_1;
	    
	    
            j++;
            // us_H2_1 by us_H2_1
            Joutput[j] = -2*us_H2_1_plus_us_H2_1[i]*us_H2_1 - us_H_1_plus_us_H2_1[i]*us_H_1;
	    
	    
            j++;
            // ge by us_H2_1
            Joutput[j] = 0;
	    
	    Joutput[j] /= mdensity;
	    
	    
            j++;
    
        // 
        // Species: ge
        //
            // us_H_1 by ge
            Joutput[j] = 0;
	    
	    
            Joutput[j] *= Tge[i];
            
            j++;
            // us_H2_1 by ge
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
    double us_H_1;
    double us_H2_1;
    double ge;
    double total_e = 0.0;
    int e_indx;

    for (i = 0; i<nstrip; i++) {
        j = i * nchem;
        us_H_1 = input[j];
        
        total_e += us_H_1 * 0;
        
        j++;
    
        us_H2_1 = input[j];
        
        total_e += us_H2_1 * 0;
        
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
    
    double us_H_1;
    double us_H2_1;
    double ge;

    for (i = 0; i<nstrip; i++) {
        j = i * nchem;
        us_H_1 = input[j];
        
        us_H_1 /= 1 * mh;
        
        /*fprintf(stderr, "us_H_1[%d] = % 0.16g\n",
                i, us_H_1);*/
        j++;
    
        us_H2_1 = input[j];
        
        us_H2_1 /= 2 * mh;
        
        /*fprintf(stderr, "us_H2_1[%d] = % 0.16g\n",
                i, us_H2_1);*/
        j++;
    
        ge = input[j];
        
        /*fprintf(stderr, "ge[%d] = % 0.16g\n",
                i, ge);*/
        j++;
    
        density = 2*us_H2_1 + us_H_1;
        strip_temperature[i] = density*ge*mh/(kb*(us_H2_1/(gamma - 1.0) + us_H_1/(gamma - 1.0)));
        if (strip_temperature[i] < 1.0)
            strip_temperature[i] = 1.0;
    }
         
}
 