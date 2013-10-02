
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

umist_data *umist_setup_data(void) {

    umist_data *data = (umist_data *) malloc(sizeof(umist_data));

    data->bounds[0] = 10.0;
    data->bounds[1] = 1000.0;
    data->nbins = 1024;
    data->dbin = (log(data->bounds[1]) - log(data->bounds[0])) / data->nbins;
    data->idbin = 1.0L / data->dbin;
    
    umist_read_rate_tables(data);
    fprintf(stderr, "Successfully read in rate tables.\n");

    umist_read_cooling_tables(data);
    fprintf(stderr, "Successfully read in cooling rate tables.\n");

    return data;

}


int umist_main(int argc, char** argv)
{
    umist_data *data = umist_setup_data();

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

    int N = 17;

    double *atol, *rtol;
    atol = (double *) alloca(N * dims * sizeof(double));
    rtol = (double *) alloca(N * dims * sizeof(double));

    double *tics = (double *) alloca(dims * sizeof(double));
    double *ics = (double *) alloca(dims * N * sizeof(double));
    double *input = (double *) alloca(dims * N * sizeof(double));
    
    unsigned int i = 0, j;
    
    fprintf(stderr, "Reading I.C. for /ge\n");
    H5LTread_dataset_double(file_id, "/ge", tics);
    for (j = 0; j < dims; j++) {
        ics[j * N + i] = tics[j] / 1.0; /* Convert to number density */
        atol[j * N + i] = tics[j] * 1e-09;
        rtol[j * N + i] = 1e-09;
        if(j==0) {
            fprintf(stderr, "ge[0] = %0.3g, atol => % 0.16g\n",
                    tics[j], atol[j]);
        }
    }
    i++;
    
    fprintf(stderr, "Reading I.C. for /us_Cm\n");
    H5LTread_dataset_double(file_id, "/us_Cm", tics);
    for (j = 0; j < dims; j++) {
        ics[j * N + i] = tics[j] / -1; /* Convert to number density */
        atol[j * N + i] = tics[j] * 1e-09;
        rtol[j * N + i] = 1e-09;
        if(j==0) {
            fprintf(stderr, "us_Cm[0] = %0.3g, atol => % 0.16g\n",
                    tics[j], atol[j]);
        }
    }
    i++;
    
    fprintf(stderr, "Reading I.C. for /us_CO\n");
    H5LTread_dataset_double(file_id, "/us_CO", tics);
    for (j = 0; j < dims; j++) {
        ics[j * N + i] = tics[j] / -1; /* Convert to number density */
        atol[j * N + i] = tics[j] * 1e-09;
        rtol[j * N + i] = 1e-09;
        if(j==0) {
            fprintf(stderr, "us_CO[0] = %0.3g, atol => % 0.16g\n",
                    tics[j], atol[j]);
        }
    }
    i++;
    
    fprintf(stderr, "Reading I.C. for /us_em\n");
    H5LTread_dataset_double(file_id, "/us_em", tics);
    for (j = 0; j < dims; j++) {
        ics[j * N + i] = tics[j] / -1; /* Convert to number density */
        atol[j * N + i] = tics[j] * 1e-09;
        rtol[j * N + i] = 1e-09;
        if(j==0) {
            fprintf(stderr, "us_em[0] = %0.3g, atol => % 0.16g\n",
                    tics[j], atol[j]);
        }
    }
    i++;
    
    fprintf(stderr, "Reading I.C. for /us_O\n");
    H5LTread_dataset_double(file_id, "/us_O", tics);
    for (j = 0; j < dims; j++) {
        ics[j * N + i] = tics[j] / -1; /* Convert to number density */
        atol[j * N + i] = tics[j] * 1e-09;
        rtol[j * N + i] = 1e-09;
        if(j==0) {
            fprintf(stderr, "us_O[0] = %0.3g, atol => % 0.16g\n",
                    tics[j], atol[j]);
        }
    }
    i++;
    
    fprintf(stderr, "Reading I.C. for /us_C\n");
    H5LTread_dataset_double(file_id, "/us_C", tics);
    for (j = 0; j < dims; j++) {
        ics[j * N + i] = tics[j] / -1; /* Convert to number density */
        atol[j * N + i] = tics[j] * 1e-09;
        rtol[j * N + i] = 1e-09;
        if(j==0) {
            fprintf(stderr, "us_C[0] = %0.3g, atol => % 0.16g\n",
                    tics[j], atol[j]);
        }
    }
    i++;
    
    fprintf(stderr, "Reading I.C. for /us_Om\n");
    H5LTread_dataset_double(file_id, "/us_Om", tics);
    for (j = 0; j < dims; j++) {
        ics[j * N + i] = tics[j] / -1; /* Convert to number density */
        atol[j * N + i] = tics[j] * 1e-09;
        rtol[j * N + i] = 1e-09;
        if(j==0) {
            fprintf(stderr, "us_Om[0] = %0.3g, atol => % 0.16g\n",
                    tics[j], atol[j]);
        }
    }
    i++;
    
    fprintf(stderr, "Reading I.C. for /us_OHm\n");
    H5LTread_dataset_double(file_id, "/us_OHm", tics);
    for (j = 0; j < dims; j++) {
        ics[j * N + i] = tics[j] / -1; /* Convert to number density */
        atol[j * N + i] = tics[j] * 1e-09;
        rtol[j * N + i] = 1e-09;
        if(j==0) {
            fprintf(stderr, "us_OHm[0] = %0.3g, atol => % 0.16g\n",
                    tics[j], atol[j]);
        }
    }
    i++;
    
    fprintf(stderr, "Reading I.C. for /us_Hm\n");
    H5LTread_dataset_double(file_id, "/us_Hm", tics);
    for (j = 0; j < dims; j++) {
        ics[j * N + i] = tics[j] / -1; /* Convert to number density */
        atol[j * N + i] = tics[j] * 1e-09;
        rtol[j * N + i] = 1e-09;
        if(j==0) {
            fprintf(stderr, "us_Hm[0] = %0.3g, atol => % 0.16g\n",
                    tics[j], atol[j]);
        }
    }
    i++;
    
    fprintf(stderr, "Reading I.C. for /us_Cp\n");
    H5LTread_dataset_double(file_id, "/us_Cp", tics);
    for (j = 0; j < dims; j++) {
        ics[j * N + i] = tics[j] / -1; /* Convert to number density */
        atol[j * N + i] = tics[j] * 1e-09;
        rtol[j * N + i] = 1e-09;
        if(j==0) {
            fprintf(stderr, "us_Cp[0] = %0.3g, atol => % 0.16g\n",
                    tics[j], atol[j]);
        }
    }
    i++;
    
    fprintf(stderr, "Reading I.C. for /us_H2\n");
    H5LTread_dataset_double(file_id, "/us_H2", tics);
    for (j = 0; j < dims; j++) {
        ics[j * N + i] = tics[j] / -1; /* Convert to number density */
        atol[j * N + i] = tics[j] * 1e-09;
        rtol[j * N + i] = 1e-09;
        if(j==0) {
            fprintf(stderr, "us_H2[0] = %0.3g, atol => % 0.16g\n",
                    tics[j], atol[j]);
        }
    }
    i++;
    
    fprintf(stderr, "Reading I.C. for /us_H2p\n");
    H5LTread_dataset_double(file_id, "/us_H2p", tics);
    for (j = 0; j < dims; j++) {
        ics[j * N + i] = tics[j] / -1; /* Convert to number density */
        atol[j * N + i] = tics[j] * 1e-09;
        rtol[j * N + i] = 1e-09;
        if(j==0) {
            fprintf(stderr, "us_H2p[0] = %0.3g, atol => % 0.16g\n",
                    tics[j], atol[j]);
        }
    }
    i++;
    
    fprintf(stderr, "Reading I.C. for /us_H\n");
    H5LTread_dataset_double(file_id, "/us_H", tics);
    for (j = 0; j < dims; j++) {
        ics[j * N + i] = tics[j] / -1; /* Convert to number density */
        atol[j * N + i] = tics[j] * 1e-09;
        rtol[j * N + i] = 1e-09;
        if(j==0) {
            fprintf(stderr, "us_H[0] = %0.3g, atol => % 0.16g\n",
                    tics[j], atol[j]);
        }
    }
    i++;
    
    fprintf(stderr, "Reading I.C. for /us_Hp\n");
    H5LTread_dataset_double(file_id, "/us_Hp", tics);
    for (j = 0; j < dims; j++) {
        ics[j * N + i] = tics[j] / -1; /* Convert to number density */
        atol[j * N + i] = tics[j] * 1e-09;
        rtol[j * N + i] = 1e-09;
        if(j==0) {
            fprintf(stderr, "us_Hp[0] = %0.3g, atol => % 0.16g\n",
                    tics[j], atol[j]);
        }
    }
    i++;
    
    fprintf(stderr, "Reading I.C. for /us_Op\n");
    H5LTread_dataset_double(file_id, "/us_Op", tics);
    for (j = 0; j < dims; j++) {
        ics[j * N + i] = tics[j] / -1; /* Convert to number density */
        atol[j * N + i] = tics[j] * 1e-09;
        rtol[j * N + i] = 1e-09;
        if(j==0) {
            fprintf(stderr, "us_Op[0] = %0.3g, atol => % 0.16g\n",
                    tics[j], atol[j]);
        }
    }
    i++;
    
    fprintf(stderr, "Reading I.C. for /us_OHp\n");
    H5LTread_dataset_double(file_id, "/us_OHp", tics);
    for (j = 0; j < dims; j++) {
        ics[j * N + i] = tics[j] / -1; /* Convert to number density */
        atol[j * N + i] = tics[j] * 1e-09;
        rtol[j * N + i] = 1e-09;
        if(j==0) {
            fprintf(stderr, "us_OHp[0] = %0.3g, atol => % 0.16g\n",
                    tics[j], atol[j]);
        }
    }
    i++;
    
    fprintf(stderr, "Reading I.C. for /us_OH\n");
    H5LTread_dataset_double(file_id, "/us_OH", tics);
    for (j = 0; j < dims; j++) {
        ics[j * N + i] = tics[j] / -1; /* Convert to number density */
        atol[j * N + i] = tics[j] * 1e-09;
        rtol[j * N + i] = 1e-09;
        if(j==0) {
            fprintf(stderr, "us_OH[0] = %0.3g, atol => % 0.16g\n",
                    tics[j], atol[j]);
        }
    }
    i++;
    

    H5Fclose(file_id);

    double dtf = 3.1557e13;
    double dt = -1.0;
    for (i = 0; i < dims * N; i++) input[i] = ics[i];
    double ttot;
    ttot = dengo_evolve_umist(dtf, dt, input, rtol, atol, dims, data);

    /* Write results to HDF5 file */
    file_id = H5Fcreate("umist_solution.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    hsize_t dimsarr[1];
    dimsarr[0] = dims;
    i = 0;
    
    double ge[dims];
    for (j = 0; j < dims; j++) {
        ge[j] = input[j * N + i] * 1.0;
    }
    fprintf(stderr, "Writing solution for /ge\n");
    H5LTmake_dataset_double(file_id, "/ge", 1, dimsarr, ge);
    i++;
    
    double us_Cm[dims];
    for (j = 0; j < dims; j++) {
        us_Cm[j] = input[j * N + i] * -1;
    }
    fprintf(stderr, "Writing solution for /us_Cm\n");
    H5LTmake_dataset_double(file_id, "/us_Cm", 1, dimsarr, us_Cm);
    i++;
    
    double us_CO[dims];
    for (j = 0; j < dims; j++) {
        us_CO[j] = input[j * N + i] * -1;
    }
    fprintf(stderr, "Writing solution for /us_CO\n");
    H5LTmake_dataset_double(file_id, "/us_CO", 1, dimsarr, us_CO);
    i++;
    
    double us_em[dims];
    for (j = 0; j < dims; j++) {
        us_em[j] = input[j * N + i] * -1;
    }
    fprintf(stderr, "Writing solution for /us_em\n");
    H5LTmake_dataset_double(file_id, "/us_em", 1, dimsarr, us_em);
    i++;
    
    double us_O[dims];
    for (j = 0; j < dims; j++) {
        us_O[j] = input[j * N + i] * -1;
    }
    fprintf(stderr, "Writing solution for /us_O\n");
    H5LTmake_dataset_double(file_id, "/us_O", 1, dimsarr, us_O);
    i++;
    
    double us_C[dims];
    for (j = 0; j < dims; j++) {
        us_C[j] = input[j * N + i] * -1;
    }
    fprintf(stderr, "Writing solution for /us_C\n");
    H5LTmake_dataset_double(file_id, "/us_C", 1, dimsarr, us_C);
    i++;
    
    double us_Om[dims];
    for (j = 0; j < dims; j++) {
        us_Om[j] = input[j * N + i] * -1;
    }
    fprintf(stderr, "Writing solution for /us_Om\n");
    H5LTmake_dataset_double(file_id, "/us_Om", 1, dimsarr, us_Om);
    i++;
    
    double us_OHm[dims];
    for (j = 0; j < dims; j++) {
        us_OHm[j] = input[j * N + i] * -1;
    }
    fprintf(stderr, "Writing solution for /us_OHm\n");
    H5LTmake_dataset_double(file_id, "/us_OHm", 1, dimsarr, us_OHm);
    i++;
    
    double us_Hm[dims];
    for (j = 0; j < dims; j++) {
        us_Hm[j] = input[j * N + i] * -1;
    }
    fprintf(stderr, "Writing solution for /us_Hm\n");
    H5LTmake_dataset_double(file_id, "/us_Hm", 1, dimsarr, us_Hm);
    i++;
    
    double us_Cp[dims];
    for (j = 0; j < dims; j++) {
        us_Cp[j] = input[j * N + i] * -1;
    }
    fprintf(stderr, "Writing solution for /us_Cp\n");
    H5LTmake_dataset_double(file_id, "/us_Cp", 1, dimsarr, us_Cp);
    i++;
    
    double us_H2[dims];
    for (j = 0; j < dims; j++) {
        us_H2[j] = input[j * N + i] * -1;
    }
    fprintf(stderr, "Writing solution for /us_H2\n");
    H5LTmake_dataset_double(file_id, "/us_H2", 1, dimsarr, us_H2);
    i++;
    
    double us_H2p[dims];
    for (j = 0; j < dims; j++) {
        us_H2p[j] = input[j * N + i] * -1;
    }
    fprintf(stderr, "Writing solution for /us_H2p\n");
    H5LTmake_dataset_double(file_id, "/us_H2p", 1, dimsarr, us_H2p);
    i++;
    
    double us_H[dims];
    for (j = 0; j < dims; j++) {
        us_H[j] = input[j * N + i] * -1;
    }
    fprintf(stderr, "Writing solution for /us_H\n");
    H5LTmake_dataset_double(file_id, "/us_H", 1, dimsarr, us_H);
    i++;
    
    double us_Hp[dims];
    for (j = 0; j < dims; j++) {
        us_Hp[j] = input[j * N + i] * -1;
    }
    fprintf(stderr, "Writing solution for /us_Hp\n");
    H5LTmake_dataset_double(file_id, "/us_Hp", 1, dimsarr, us_Hp);
    i++;
    
    double us_Op[dims];
    for (j = 0; j < dims; j++) {
        us_Op[j] = input[j * N + i] * -1;
    }
    fprintf(stderr, "Writing solution for /us_Op\n");
    H5LTmake_dataset_double(file_id, "/us_Op", 1, dimsarr, us_Op);
    i++;
    
    double us_OHp[dims];
    for (j = 0; j < dims; j++) {
        us_OHp[j] = input[j * N + i] * -1;
    }
    fprintf(stderr, "Writing solution for /us_OHp\n");
    H5LTmake_dataset_double(file_id, "/us_OHp", 1, dimsarr, us_OHp);
    i++;
    
    double us_OH[dims];
    for (j = 0; j < dims; j++) {
        us_OH[j] = input[j * N + i] * -1;
    }
    fprintf(stderr, "Writing solution for /us_OH\n");
    H5LTmake_dataset_double(file_id, "/us_OH", 1, dimsarr, us_OH);
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
 



double dengo_evolve_umist (double dtf, double &dt, double *input,
            double *rtol, double *atol, int dims, umist_data *data) {
    int i, j;
    hid_t file_id;
    fprintf(stderr, "  ncells = % 3i\n", (int) dims);

    int N = 17;
    rhs_f f = calculate_rhs_umist;
    jac_f jf = calculate_jacobian_umist;
    if (dt < 0) dt = dtf / 1e3;
    int niter = 0;
    int siter = 0;
    double ttot = 0;
    double *scale = (double *) alloca(dims * N * sizeof(double));
    double *prev = (double *) alloca(dims * N * sizeof(double));
    for (i = 0; i < dims * N; i++) scale[i] = 1.0;
    for (i = 0; i < dims * N; i++) prev[i] = input[i];
    while (ttot < dtf) {
        int rv = BE_chem_solve(f, jf, input, dt, rtol, atol, dims, N, scale, (void *) data);
        /*
        fprintf(stderr, "Return value [%d]: %i.  %0.5g / %0.5g = %0.5g (%0.5g)\n",
                niter, rv, ttot, dtf, ttot/dtf, dt);
        fprintf(stderr, "Value[80] = %0.5g %0.5g %0.5g\n",
                input[80], prev[80], ics[80]);
        */
        for (i = 0; i < dims * N; i++) {
            if (input[i] < 0) {
                rv = 1;
                break;
            }
        }
        if (rv == 0) {
	    if (siter == 49999) break;
	    siter++;
        if (siter % 1000 == 0) {
            fprintf(stderr, "Successful Iteration[%d]: (%0.4g) %0.16g / %0.16g\n",
                     siter, dt, ttot, dtf);
        }
        ttot += dt;
	    dt = DMIN(dt * 1.1, dtf - ttot);
	    
	    for (i = 0; i < dims * N; i++) prev[i] = input[i];
        } else {
            dt /= 2.0;
            for (i = 0; i < dims * N; i++) input[i] = prev[i];
            if (dt == 0.0)  {
                fprintf(stderr, "Dying!\n");
                break;
            }
        }
        niter++;
    }
    fprintf(stderr, "End: %0.5g / %0.5g (%0.5g)\n",
        ttot, dtf, dtf-ttot);

    return ttot;
}
 


void umist_read_rate_tables(umist_data *data)
{
    hid_t file_id = H5Fopen("umist_tables.h5", H5F_ACC_RDONLY, H5P_DEFAULT);
    /* Allocate the correct number of rate tables */
    H5LTread_dataset_double(file_id, "/us_C+us_Om", data->r_us_C+us_Om);
    H5LTread_dataset_double(file_id, "/us_Cm+us_Cp", data->r_us_Cm+us_Cp);
    H5LTread_dataset_double(file_id, "/us_Cm+us_Hp", data->r_us_Cm+us_Hp);
    H5LTread_dataset_double(file_id, "/us_Cm+us_O", data->r_us_Cm+us_O);
    H5LTread_dataset_double(file_id, "/us_Cm+us_Op", data->r_us_Cm+us_Op);
    H5LTread_dataset_double(file_id, "/us_H+us_H2p", data->r_us_H+us_H2p);
    H5LTread_dataset_double(file_id, "/us_H+us_Om", data->r_us_H+us_Om);
    H5LTread_dataset_double(file_id, "/us_H+us_Op", data->r_us_H+us_Op);
    H5LTread_dataset_double(file_id, "/us_H2+us_Op", data->r_us_H2+us_Op);
    H5LTread_dataset_double(file_id, "/us_H2p+us_O", data->r_us_H2p+us_O);
    H5LTread_dataset_double(file_id, "/us_H2p+us_OH", data->r_us_H2p+us_OH);
    H5LTread_dataset_double(file_id, "/us_Hm+us_Cp", data->r_us_Hm+us_Cp);
    H5LTread_dataset_double(file_id, "/us_Hm+us_Hp", data->r_us_Hm+us_Hp);
    H5LTread_dataset_double(file_id, "/us_Hm+us_O", data->r_us_Hm+us_O);
    H5LTread_dataset_double(file_id, "/us_Hm+us_Op", data->r_us_Hm+us_Op);
    H5LTread_dataset_double(file_id, "/us_Hp+us_O", data->r_us_Hp+us_O);
    H5LTread_dataset_double(file_id, "/us_Hp+us_OH", data->r_us_Hp+us_OH);
    H5LTread_dataset_double(file_id, "/us_OHm+us_Cp", data->r_us_OHm+us_Cp);
    H5LTread_dataset_double(file_id, "/us_OHm+us_Hp", data->r_us_OHm+us_Hp);
    H5LTread_dataset_double(file_id, "/us_OHm+us_Op", data->r_us_OHm+us_Op);
    H5LTread_dataset_double(file_id, "/us_Om+us_Cp", data->r_us_Om+us_Cp);
    H5LTread_dataset_double(file_id, "/us_Om+us_Hp", data->r_us_Om+us_Hp);
    H5LTread_dataset_double(file_id, "/us_Om+us_Op", data->r_us_Om+us_Op);

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
    double ge;
    double us_Cm;
    double us_CO;
    double us_em;
    double us_O;
    double us_C;
    double us_Om;
    double us_OHm;
    double us_Hm;
    double us_Cp;
    double us_H2;
    double us_H2p;
    double us_H;
    double us_Hp;
    double us_Op;
    double us_OHp;
    double us_OH;

    for (i = 0; i<nstrip; i++) {
        j = i * nchem;
        ge = input[j];
        /*fprintf(stderr, "ge[%d] = % 0.16g\n",
                i, ge);*/
        j++;
    
        us_Cm = input[j];
        /*fprintf(stderr, "us_Cm[%d] = % 0.16g\n",
                i, us_Cm);*/
        j++;
    
        us_CO = input[j];
        /*fprintf(stderr, "us_CO[%d] = % 0.16g\n",
                i, us_CO);*/
        j++;
    
        us_em = input[j];
        /*fprintf(stderr, "us_em[%d] = % 0.16g\n",
                i, us_em);*/
        j++;
    
        us_O = input[j];
        /*fprintf(stderr, "us_O[%d] = % 0.16g\n",
                i, us_O);*/
        j++;
    
        us_C = input[j];
        /*fprintf(stderr, "us_C[%d] = % 0.16g\n",
                i, us_C);*/
        j++;
    
        us_Om = input[j];
        /*fprintf(stderr, "us_Om[%d] = % 0.16g\n",
                i, us_Om);*/
        j++;
    
        us_OHm = input[j];
        /*fprintf(stderr, "us_OHm[%d] = % 0.16g\n",
                i, us_OHm);*/
        j++;
    
        us_Hm = input[j];
        /*fprintf(stderr, "us_Hm[%d] = % 0.16g\n",
                i, us_Hm);*/
        j++;
    
        us_Cp = input[j];
        /*fprintf(stderr, "us_Cp[%d] = % 0.16g\n",
                i, us_Cp);*/
        j++;
    
        us_H2 = input[j];
        /*fprintf(stderr, "us_H2[%d] = % 0.16g\n",
                i, us_H2);*/
        j++;
    
        us_H2p = input[j];
        /*fprintf(stderr, "us_H2p[%d] = % 0.16g\n",
                i, us_H2p);*/
        j++;
    
        us_H = input[j];
        /*fprintf(stderr, "us_H[%d] = % 0.16g\n",
                i, us_H);*/
        j++;
    
        us_Hp = input[j];
        /*fprintf(stderr, "us_Hp[%d] = % 0.16g\n",
                i, us_Hp);*/
        j++;
    
        us_Op = input[j];
        /*fprintf(stderr, "us_Op[%d] = % 0.16g\n",
                i, us_Op);*/
        j++;
    
        us_OHp = input[j];
        /*fprintf(stderr, "us_OHp[%d] = % 0.16g\n",
                i, us_OHp);*/
        j++;
    
        us_OH = input[j];
        /*fprintf(stderr, "us_OH[%d] = % 0.16g\n",
                i, us_OH);*/
        j++;
    
        density = -us_C - us_CO - us_Cm - us_Cp - us_H - us_H2 - us_H2p - us_Hm - us_Hp - us_O - us_OH - us_OHm - us_OHp - us_Om - us_Op - us_em;
        data->Ts[i] = density*ge*mh/(kb*(us_C/(gamma - 1.0) + us_CO/(gamma - 1.0) + us_Cm/(gamma - 1.0) + us_Cp/(gamma - 1.0) + us_H/(gamma - 1.0) + us_H2/(gamma - 1.0) + us_H2p/(gamma - 1.0) + us_Hm/(gamma - 1.0) + us_Hp/(gamma - 1.0) + us_O/(gamma - 1.0) + us_OH/(gamma - 1.0) + us_OHm/(gamma - 1.0) + us_OHp/(gamma - 1.0) + us_Om/(gamma - 1.0) + us_Op/(gamma - 1.0) + us_em/(gamma - 1.0)));
        if (data->Ts[i] < data->bounds[0]) {
            data->Ts[i] = data->bounds[0];
        } else if (data->Ts[i] > data->bounds[1]) {
            data->Ts[i] = data->bounds[1];
        }
        data->logTs[i] = log(data->Ts[i]);
	data->dTs_ge[i] = 
        density*mh/(kb*(us_C/(gamma - 1.0) + us_CO/(gamma - 1.0) + us_Cm/(gamma - 1.0) + us_Cp/(gamma - 1.0) + us_H/(gamma - 1.0) + us_H2/(gamma - 1.0) + us_H2p/(gamma - 1.0) + us_Hm/(gamma - 1.0) + us_Hp/(gamma - 1.0) + us_O/(gamma - 1.0) + us_OH/(gamma - 1.0) + us_OHm/(gamma - 1.0) + us_OHp/(gamma - 1.0) + us_Om/(gamma - 1.0) + us_Op/(gamma - 1.0) + us_em/(gamma - 1.0)));
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
    int i, bin_id;
    double lb, t1, t2;
    lb = log(data->bounds[0]);
    /*fprintf(stderr, "lb = % 0.16g, ub = % 0.16g\n", lb, ub);*/
    for (i = 0; i < nstrip; i++) {
        data->bin_id[i] = bin_id = (int) (data->idbin * (data->logTs[i] - lb));
        if (data->bin_id[i] <= 0) {
            data->bin_id[i] = 0;
        } else if (data->bin_id[i] >= data->nbins) {
            data->bin_id[i] = data->nbins - 2;
        }
        t1 = (lb + (bin_id    ) * data->dbin);
        t2 = (lb + (bin_id + 1) * data->dbin);
        data->Tdef[i] = (data->logTs[i] - t1)/(t2 - t1);
        data->dT[i] = (t2 - t1);
        /*fprintf(stderr, "INTERP: %d, bin_id = %d, dT = % 0.16g, T = % 0.16g, logT = % 0.16g\n",
                i, data->bin_id[i], data->dT[i], data->Ts[i],
                data->logTs[i]);*/
    }
    for (i = 0; i < nstrip; i++) {
        bin_id = data->bin_id[i];
        data->rs_us_C+us_Om[i] = data->r_us_C+us_Om[bin_id] +
            data->Tdef[i] * (data->r_us_C+us_Om[bin_id+1] - data->r_us_C+us_Om[bin_id]);
        data->drs_us_C+us_Om[i] = (data->r_us_C+us_Om[bin_id+1] - data->r_us_C+us_Om[bin_id]);
        data->drs_us_C+us_Om[i] /= data->dT[i];
	data->drs_us_C+us_Om[i] /= data->Ts[i];
    }
    
    for (i = 0; i < nstrip; i++) {
        bin_id = data->bin_id[i];
        data->rs_us_Cm+us_Cp[i] = data->r_us_Cm+us_Cp[bin_id] +
            data->Tdef[i] * (data->r_us_Cm+us_Cp[bin_id+1] - data->r_us_Cm+us_Cp[bin_id]);
        data->drs_us_Cm+us_Cp[i] = (data->r_us_Cm+us_Cp[bin_id+1] - data->r_us_Cm+us_Cp[bin_id]);
        data->drs_us_Cm+us_Cp[i] /= data->dT[i];
	data->drs_us_Cm+us_Cp[i] /= data->Ts[i];
    }
    
    for (i = 0; i < nstrip; i++) {
        bin_id = data->bin_id[i];
        data->rs_us_Cm+us_Hp[i] = data->r_us_Cm+us_Hp[bin_id] +
            data->Tdef[i] * (data->r_us_Cm+us_Hp[bin_id+1] - data->r_us_Cm+us_Hp[bin_id]);
        data->drs_us_Cm+us_Hp[i] = (data->r_us_Cm+us_Hp[bin_id+1] - data->r_us_Cm+us_Hp[bin_id]);
        data->drs_us_Cm+us_Hp[i] /= data->dT[i];
	data->drs_us_Cm+us_Hp[i] /= data->Ts[i];
    }
    
    for (i = 0; i < nstrip; i++) {
        bin_id = data->bin_id[i];
        data->rs_us_Cm+us_O[i] = data->r_us_Cm+us_O[bin_id] +
            data->Tdef[i] * (data->r_us_Cm+us_O[bin_id+1] - data->r_us_Cm+us_O[bin_id]);
        data->drs_us_Cm+us_O[i] = (data->r_us_Cm+us_O[bin_id+1] - data->r_us_Cm+us_O[bin_id]);
        data->drs_us_Cm+us_O[i] /= data->dT[i];
	data->drs_us_Cm+us_O[i] /= data->Ts[i];
    }
    
    for (i = 0; i < nstrip; i++) {
        bin_id = data->bin_id[i];
        data->rs_us_Cm+us_Op[i] = data->r_us_Cm+us_Op[bin_id] +
            data->Tdef[i] * (data->r_us_Cm+us_Op[bin_id+1] - data->r_us_Cm+us_Op[bin_id]);
        data->drs_us_Cm+us_Op[i] = (data->r_us_Cm+us_Op[bin_id+1] - data->r_us_Cm+us_Op[bin_id]);
        data->drs_us_Cm+us_Op[i] /= data->dT[i];
	data->drs_us_Cm+us_Op[i] /= data->Ts[i];
    }
    
    for (i = 0; i < nstrip; i++) {
        bin_id = data->bin_id[i];
        data->rs_us_H+us_H2p[i] = data->r_us_H+us_H2p[bin_id] +
            data->Tdef[i] * (data->r_us_H+us_H2p[bin_id+1] - data->r_us_H+us_H2p[bin_id]);
        data->drs_us_H+us_H2p[i] = (data->r_us_H+us_H2p[bin_id+1] - data->r_us_H+us_H2p[bin_id]);
        data->drs_us_H+us_H2p[i] /= data->dT[i];
	data->drs_us_H+us_H2p[i] /= data->Ts[i];
    }
    
    for (i = 0; i < nstrip; i++) {
        bin_id = data->bin_id[i];
        data->rs_us_H+us_Om[i] = data->r_us_H+us_Om[bin_id] +
            data->Tdef[i] * (data->r_us_H+us_Om[bin_id+1] - data->r_us_H+us_Om[bin_id]);
        data->drs_us_H+us_Om[i] = (data->r_us_H+us_Om[bin_id+1] - data->r_us_H+us_Om[bin_id]);
        data->drs_us_H+us_Om[i] /= data->dT[i];
	data->drs_us_H+us_Om[i] /= data->Ts[i];
    }
    
    for (i = 0; i < nstrip; i++) {
        bin_id = data->bin_id[i];
        data->rs_us_H+us_Op[i] = data->r_us_H+us_Op[bin_id] +
            data->Tdef[i] * (data->r_us_H+us_Op[bin_id+1] - data->r_us_H+us_Op[bin_id]);
        data->drs_us_H+us_Op[i] = (data->r_us_H+us_Op[bin_id+1] - data->r_us_H+us_Op[bin_id]);
        data->drs_us_H+us_Op[i] /= data->dT[i];
	data->drs_us_H+us_Op[i] /= data->Ts[i];
    }
    
    for (i = 0; i < nstrip; i++) {
        bin_id = data->bin_id[i];
        data->rs_us_H2+us_Op[i] = data->r_us_H2+us_Op[bin_id] +
            data->Tdef[i] * (data->r_us_H2+us_Op[bin_id+1] - data->r_us_H2+us_Op[bin_id]);
        data->drs_us_H2+us_Op[i] = (data->r_us_H2+us_Op[bin_id+1] - data->r_us_H2+us_Op[bin_id]);
        data->drs_us_H2+us_Op[i] /= data->dT[i];
	data->drs_us_H2+us_Op[i] /= data->Ts[i];
    }
    
    for (i = 0; i < nstrip; i++) {
        bin_id = data->bin_id[i];
        data->rs_us_H2p+us_O[i] = data->r_us_H2p+us_O[bin_id] +
            data->Tdef[i] * (data->r_us_H2p+us_O[bin_id+1] - data->r_us_H2p+us_O[bin_id]);
        data->drs_us_H2p+us_O[i] = (data->r_us_H2p+us_O[bin_id+1] - data->r_us_H2p+us_O[bin_id]);
        data->drs_us_H2p+us_O[i] /= data->dT[i];
	data->drs_us_H2p+us_O[i] /= data->Ts[i];
    }
    
    for (i = 0; i < nstrip; i++) {
        bin_id = data->bin_id[i];
        data->rs_us_H2p+us_OH[i] = data->r_us_H2p+us_OH[bin_id] +
            data->Tdef[i] * (data->r_us_H2p+us_OH[bin_id+1] - data->r_us_H2p+us_OH[bin_id]);
        data->drs_us_H2p+us_OH[i] = (data->r_us_H2p+us_OH[bin_id+1] - data->r_us_H2p+us_OH[bin_id]);
        data->drs_us_H2p+us_OH[i] /= data->dT[i];
	data->drs_us_H2p+us_OH[i] /= data->Ts[i];
    }
    
    for (i = 0; i < nstrip; i++) {
        bin_id = data->bin_id[i];
        data->rs_us_Hm+us_Cp[i] = data->r_us_Hm+us_Cp[bin_id] +
            data->Tdef[i] * (data->r_us_Hm+us_Cp[bin_id+1] - data->r_us_Hm+us_Cp[bin_id]);
        data->drs_us_Hm+us_Cp[i] = (data->r_us_Hm+us_Cp[bin_id+1] - data->r_us_Hm+us_Cp[bin_id]);
        data->drs_us_Hm+us_Cp[i] /= data->dT[i];
	data->drs_us_Hm+us_Cp[i] /= data->Ts[i];
    }
    
    for (i = 0; i < nstrip; i++) {
        bin_id = data->bin_id[i];
        data->rs_us_Hm+us_Hp[i] = data->r_us_Hm+us_Hp[bin_id] +
            data->Tdef[i] * (data->r_us_Hm+us_Hp[bin_id+1] - data->r_us_Hm+us_Hp[bin_id]);
        data->drs_us_Hm+us_Hp[i] = (data->r_us_Hm+us_Hp[bin_id+1] - data->r_us_Hm+us_Hp[bin_id]);
        data->drs_us_Hm+us_Hp[i] /= data->dT[i];
	data->drs_us_Hm+us_Hp[i] /= data->Ts[i];
    }
    
    for (i = 0; i < nstrip; i++) {
        bin_id = data->bin_id[i];
        data->rs_us_Hm+us_O[i] = data->r_us_Hm+us_O[bin_id] +
            data->Tdef[i] * (data->r_us_Hm+us_O[bin_id+1] - data->r_us_Hm+us_O[bin_id]);
        data->drs_us_Hm+us_O[i] = (data->r_us_Hm+us_O[bin_id+1] - data->r_us_Hm+us_O[bin_id]);
        data->drs_us_Hm+us_O[i] /= data->dT[i];
	data->drs_us_Hm+us_O[i] /= data->Ts[i];
    }
    
    for (i = 0; i < nstrip; i++) {
        bin_id = data->bin_id[i];
        data->rs_us_Hm+us_Op[i] = data->r_us_Hm+us_Op[bin_id] +
            data->Tdef[i] * (data->r_us_Hm+us_Op[bin_id+1] - data->r_us_Hm+us_Op[bin_id]);
        data->drs_us_Hm+us_Op[i] = (data->r_us_Hm+us_Op[bin_id+1] - data->r_us_Hm+us_Op[bin_id]);
        data->drs_us_Hm+us_Op[i] /= data->dT[i];
	data->drs_us_Hm+us_Op[i] /= data->Ts[i];
    }
    
    for (i = 0; i < nstrip; i++) {
        bin_id = data->bin_id[i];
        data->rs_us_Hp+us_O[i] = data->r_us_Hp+us_O[bin_id] +
            data->Tdef[i] * (data->r_us_Hp+us_O[bin_id+1] - data->r_us_Hp+us_O[bin_id]);
        data->drs_us_Hp+us_O[i] = (data->r_us_Hp+us_O[bin_id+1] - data->r_us_Hp+us_O[bin_id]);
        data->drs_us_Hp+us_O[i] /= data->dT[i];
	data->drs_us_Hp+us_O[i] /= data->Ts[i];
    }
    
    for (i = 0; i < nstrip; i++) {
        bin_id = data->bin_id[i];
        data->rs_us_Hp+us_OH[i] = data->r_us_Hp+us_OH[bin_id] +
            data->Tdef[i] * (data->r_us_Hp+us_OH[bin_id+1] - data->r_us_Hp+us_OH[bin_id]);
        data->drs_us_Hp+us_OH[i] = (data->r_us_Hp+us_OH[bin_id+1] - data->r_us_Hp+us_OH[bin_id]);
        data->drs_us_Hp+us_OH[i] /= data->dT[i];
	data->drs_us_Hp+us_OH[i] /= data->Ts[i];
    }
    
    for (i = 0; i < nstrip; i++) {
        bin_id = data->bin_id[i];
        data->rs_us_OHm+us_Cp[i] = data->r_us_OHm+us_Cp[bin_id] +
            data->Tdef[i] * (data->r_us_OHm+us_Cp[bin_id+1] - data->r_us_OHm+us_Cp[bin_id]);
        data->drs_us_OHm+us_Cp[i] = (data->r_us_OHm+us_Cp[bin_id+1] - data->r_us_OHm+us_Cp[bin_id]);
        data->drs_us_OHm+us_Cp[i] /= data->dT[i];
	data->drs_us_OHm+us_Cp[i] /= data->Ts[i];
    }
    
    for (i = 0; i < nstrip; i++) {
        bin_id = data->bin_id[i];
        data->rs_us_OHm+us_Hp[i] = data->r_us_OHm+us_Hp[bin_id] +
            data->Tdef[i] * (data->r_us_OHm+us_Hp[bin_id+1] - data->r_us_OHm+us_Hp[bin_id]);
        data->drs_us_OHm+us_Hp[i] = (data->r_us_OHm+us_Hp[bin_id+1] - data->r_us_OHm+us_Hp[bin_id]);
        data->drs_us_OHm+us_Hp[i] /= data->dT[i];
	data->drs_us_OHm+us_Hp[i] /= data->Ts[i];
    }
    
    for (i = 0; i < nstrip; i++) {
        bin_id = data->bin_id[i];
        data->rs_us_OHm+us_Op[i] = data->r_us_OHm+us_Op[bin_id] +
            data->Tdef[i] * (data->r_us_OHm+us_Op[bin_id+1] - data->r_us_OHm+us_Op[bin_id]);
        data->drs_us_OHm+us_Op[i] = (data->r_us_OHm+us_Op[bin_id+1] - data->r_us_OHm+us_Op[bin_id]);
        data->drs_us_OHm+us_Op[i] /= data->dT[i];
	data->drs_us_OHm+us_Op[i] /= data->Ts[i];
    }
    
    for (i = 0; i < nstrip; i++) {
        bin_id = data->bin_id[i];
        data->rs_us_Om+us_Cp[i] = data->r_us_Om+us_Cp[bin_id] +
            data->Tdef[i] * (data->r_us_Om+us_Cp[bin_id+1] - data->r_us_Om+us_Cp[bin_id]);
        data->drs_us_Om+us_Cp[i] = (data->r_us_Om+us_Cp[bin_id+1] - data->r_us_Om+us_Cp[bin_id]);
        data->drs_us_Om+us_Cp[i] /= data->dT[i];
	data->drs_us_Om+us_Cp[i] /= data->Ts[i];
    }
    
    for (i = 0; i < nstrip; i++) {
        bin_id = data->bin_id[i];
        data->rs_us_Om+us_Hp[i] = data->r_us_Om+us_Hp[bin_id] +
            data->Tdef[i] * (data->r_us_Om+us_Hp[bin_id+1] - data->r_us_Om+us_Hp[bin_id]);
        data->drs_us_Om+us_Hp[i] = (data->r_us_Om+us_Hp[bin_id+1] - data->r_us_Om+us_Hp[bin_id]);
        data->drs_us_Om+us_Hp[i] /= data->dT[i];
	data->drs_us_Om+us_Hp[i] /= data->Ts[i];
    }
    
    for (i = 0; i < nstrip; i++) {
        bin_id = data->bin_id[i];
        data->rs_us_Om+us_Op[i] = data->r_us_Om+us_Op[bin_id] +
            data->Tdef[i] * (data->r_us_Om+us_Op[bin_id+1] - data->r_us_Om+us_Op[bin_id]);
        data->drs_us_Om+us_Op[i] = (data->r_us_Om+us_Op[bin_id+1] - data->r_us_Om+us_Op[bin_id]);
        data->drs_us_Om+us_Op[i] /= data->dT[i];
	data->drs_us_Om+us_Op[i] /= data->Ts[i];
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
    double *us_C+us_Om = data->rs_us_C+us_Om;
    double *us_Cm+us_Cp = data->rs_us_Cm+us_Cp;
    double *us_Cm+us_Hp = data->rs_us_Cm+us_Hp;
    double *us_Cm+us_O = data->rs_us_Cm+us_O;
    double *us_Cm+us_Op = data->rs_us_Cm+us_Op;
    double *us_H+us_H2p = data->rs_us_H+us_H2p;
    double *us_H+us_Om = data->rs_us_H+us_Om;
    double *us_H+us_Op = data->rs_us_H+us_Op;
    double *us_H2+us_Op = data->rs_us_H2+us_Op;
    double *us_H2p+us_O = data->rs_us_H2p+us_O;
    double *us_H2p+us_OH = data->rs_us_H2p+us_OH;
    double *us_Hm+us_Cp = data->rs_us_Hm+us_Cp;
    double *us_Hm+us_Hp = data->rs_us_Hm+us_Hp;
    double *us_Hm+us_O = data->rs_us_Hm+us_O;
    double *us_Hm+us_Op = data->rs_us_Hm+us_Op;
    double *us_Hp+us_O = data->rs_us_Hp+us_O;
    double *us_Hp+us_OH = data->rs_us_Hp+us_OH;
    double *us_OHm+us_Cp = data->rs_us_OHm+us_Cp;
    double *us_OHm+us_Hp = data->rs_us_OHm+us_Hp;
    double *us_OHm+us_Op = data->rs_us_OHm+us_Op;
    double *us_Om+us_Cp = data->rs_us_Om+us_Cp;
    double *us_Om+us_Hp = data->rs_us_Om+us_Hp;
    double *us_Om+us_Op = data->rs_us_Om+us_Op;
    double ge;
    double us_Cm;
    double us_CO;
    double us_em;
    double us_O;
    double us_C;
    double us_Om;
    double us_OHm;
    double us_Hm;
    double us_Cp;
    double us_H2;
    double us_H2p;
    double us_H;
    double us_Hp;
    double us_Op;
    double us_OHp;
    double us_OH;

    double mh = 1.67e-24;
    double total, total_e, total_de, mdensity;
    for (i = 0; i<nstrip; i++) {
        j = i * nchem;
        total = total_e = total_de = mdensity = 0.0;
        ge = input[j];
        if (ge < 0.0) {
          fprintf(stderr, "RNegative[%d][ge] = % 0.16g [%d]\n",
            i, ge, j);
            return 1;
          ge = 1e-20;
        }
        
        j++;
    
        us_Cm = input[j];
        if (us_Cm < 0.0) {
          fprintf(stderr, "RNegative[%d][us_Cm] = % 0.16g [%d]\n",
            i, us_Cm, j);
            return 1;
          us_Cm = 1e-20;
        }
        
        
          total+=us_Cm * -1;
        
        
        j++;
    
        us_CO = input[j];
        if (us_CO < 0.0) {
          fprintf(stderr, "RNegative[%d][us_CO] = % 0.16g [%d]\n",
            i, us_CO, j);
            return 1;
          us_CO = 1e-20;
        }
        
        
          total+=us_CO * -1;
        
        
        j++;
    
        us_em = input[j];
        if (us_em < 0.0) {
          fprintf(stderr, "RNegative[%d][us_em] = % 0.16g [%d]\n",
            i, us_em, j);
            return 1;
          us_em = 1e-20;
        }
        
        
          total+=us_em * -1;
        
        
        j++;
    
        us_O = input[j];
        if (us_O < 0.0) {
          fprintf(stderr, "RNegative[%d][us_O] = % 0.16g [%d]\n",
            i, us_O, j);
            return 1;
          us_O = 1e-20;
        }
        
        
          total+=us_O * -1;
        
        
        j++;
    
        us_C = input[j];
        if (us_C < 0.0) {
          fprintf(stderr, "RNegative[%d][us_C] = % 0.16g [%d]\n",
            i, us_C, j);
            return 1;
          us_C = 1e-20;
        }
        
        
          total+=us_C * -1;
        
        
        j++;
    
        us_Om = input[j];
        if (us_Om < 0.0) {
          fprintf(stderr, "RNegative[%d][us_Om] = % 0.16g [%d]\n",
            i, us_Om, j);
            return 1;
          us_Om = 1e-20;
        }
        
        
          total+=us_Om * -1;
        
        
        j++;
    
        us_OHm = input[j];
        if (us_OHm < 0.0) {
          fprintf(stderr, "RNegative[%d][us_OHm] = % 0.16g [%d]\n",
            i, us_OHm, j);
            return 1;
          us_OHm = 1e-20;
        }
        
        
          total+=us_OHm * -1;
        
        
        j++;
    
        us_Hm = input[j];
        if (us_Hm < 0.0) {
          fprintf(stderr, "RNegative[%d][us_Hm] = % 0.16g [%d]\n",
            i, us_Hm, j);
            return 1;
          us_Hm = 1e-20;
        }
        
        
          total+=us_Hm * -1;
        
        
        j++;
    
        us_Cp = input[j];
        if (us_Cp < 0.0) {
          fprintf(stderr, "RNegative[%d][us_Cp] = % 0.16g [%d]\n",
            i, us_Cp, j);
            return 1;
          us_Cp = 1e-20;
        }
        
        
          total+=us_Cp * -1;
        
        
        j++;
    
        us_H2 = input[j];
        if (us_H2 < 0.0) {
          fprintf(stderr, "RNegative[%d][us_H2] = % 0.16g [%d]\n",
            i, us_H2, j);
            return 1;
          us_H2 = 1e-20;
        }
        
        
          total+=us_H2 * -1;
        
        
        j++;
    
        us_H2p = input[j];
        if (us_H2p < 0.0) {
          fprintf(stderr, "RNegative[%d][us_H2p] = % 0.16g [%d]\n",
            i, us_H2p, j);
            return 1;
          us_H2p = 1e-20;
        }
        
        
          total+=us_H2p * -1;
        
        
        j++;
    
        us_H = input[j];
        if (us_H < 0.0) {
          fprintf(stderr, "RNegative[%d][us_H] = % 0.16g [%d]\n",
            i, us_H, j);
            return 1;
          us_H = 1e-20;
        }
        
        
          total+=us_H * -1;
        
        
        j++;
    
        us_Hp = input[j];
        if (us_Hp < 0.0) {
          fprintf(stderr, "RNegative[%d][us_Hp] = % 0.16g [%d]\n",
            i, us_Hp, j);
            return 1;
          us_Hp = 1e-20;
        }
        
        
          total+=us_Hp * -1;
        
        
        j++;
    
        us_Op = input[j];
        if (us_Op < 0.0) {
          fprintf(stderr, "RNegative[%d][us_Op] = % 0.16g [%d]\n",
            i, us_Op, j);
            return 1;
          us_Op = 1e-20;
        }
        
        
          total+=us_Op * -1;
        
        
        j++;
    
        us_OHp = input[j];
        if (us_OHp < 0.0) {
          fprintf(stderr, "RNegative[%d][us_OHp] = % 0.16g [%d]\n",
            i, us_OHp, j);
            return 1;
          us_OHp = 1e-20;
        }
        
        
          total+=us_OHp * -1;
        
        
        j++;
    
        us_OH = input[j];
        if (us_OH < 0.0) {
          fprintf(stderr, "RNegative[%d][us_OH] = % 0.16g [%d]\n",
            i, us_OH, j);
            return 1;
          us_OH = 1e-20;
        }
        
        
          total+=us_OH * -1;
        
        
        j++;
    
        mdensity = total * mh;
        total = 0.0;
        j = i * nchem;
        // 
        // Species: ge
        // 
        rhs[j] = 0;
        
	    rhs[j] /= mdensity;
        
        
        j++;
    
        // 
        // Species: us_Cm
        // 
        rhs[j] = -us_Cm+us_Cp[i]*us_Cm*us_Cp - us_Cm+us_Hp[i]*us_Cm*us_Hp - us_Cm+us_O[i]*us_Cm*us_O - us_Cm+us_Op[i]*us_Cm*us_Op;
        
            /* Already in number density, not mass density */
            total += rhs[j] * -1;
            total_e += us_Cm * -1;
        
        
            total_de += rhs[j] * -1;
        
        j++;
    
        // 
        // Species: us_CO
        // 
        rhs[j] = us_C+us_Om[i]*us_C*us_Om + us_Cm+us_O[i]*us_Cm*us_O;
        
            /* Already in number density, not mass density */
            total += rhs[j] * -1;
            total_e += us_CO * 0;
        
        
            total_de += rhs[j] * 0;
        
        j++;
    
        // 
        // Species: us_em
        // 
        rhs[j] = us_C+us_Om[i]*us_C*us_Om + us_Cm+us_O[i]*us_Cm*us_O + us_H+us_Om[i]*us_H*us_Om + us_Hm+us_O[i]*us_Hm*us_O;
        
            /* Already in number density, not mass density */
            total += rhs[j] * -1;
            total_e += us_em * -1;
        
        
            total_de += rhs[j] * -1;
        
        j++;
    
        // 
        // Species: us_O
        // 
        rhs[j] = -us_Cm+us_O[i]*us_Cm*us_O + us_Cm+us_Op[i]*us_Cm*us_Op + us_H+us_Op[i]*us_H*us_Op - us_H2p+us_O[i]*us_H2p*us_O - us_Hm+us_O[i]*us_Hm*us_O + us_Hm+us_Op[i]*us_Hm*us_Op - us_Hp+us_O[i]*us_Hp*us_O + us_OHm+us_Op[i]*us_OHm*us_Op + us_Om+us_Cp[i]*us_Cp*us_Om + us_Om+us_Hp[i]*us_Hp*us_Om + 2*us_Om+us_Op[i]*us_Om*us_Op;
        
            /* Already in number density, not mass density */
            total += rhs[j] * -1;
            total_e += us_O * 0;
        
        
            total_de += rhs[j] * 0;
        
        j++;
    
        // 
        // Species: us_C
        // 
        rhs[j] = -us_C+us_Om[i]*us_C*us_Om + 2*us_Cm+us_Cp[i]*us_Cm*us_Cp + us_Cm+us_Hp[i]*us_Cm*us_Hp + us_Cm+us_Op[i]*us_Cm*us_Op + us_Hm+us_Cp[i]*us_Cp*us_Hm + us_OHm+us_Cp[i]*us_Cp*us_OHm + us_Om+us_Cp[i]*us_Cp*us_Om;
        
            /* Already in number density, not mass density */
            total += rhs[j] * -1;
            total_e += us_C * 0;
        
        
            total_de += rhs[j] * 0;
        
        j++;
    
        // 
        // Species: us_Om
        // 
        rhs[j] = -us_C+us_Om[i]*us_C*us_Om - us_H+us_Om[i]*us_H*us_Om - us_Om+us_Cp[i]*us_Cp*us_Om - us_Om+us_Hp[i]*us_Hp*us_Om - us_Om+us_Op[i]*us_Om*us_Op;
        
            /* Already in number density, not mass density */
            total += rhs[j] * -1;
            total_e += us_Om * -1;
        
        
            total_de += rhs[j] * -1;
        
        j++;
    
        // 
        // Species: us_OHm
        // 
        rhs[j] = -us_OHm+us_Cp[i]*us_Cp*us_OHm - us_OHm+us_Hp[i]*us_Hp*us_OHm - us_OHm+us_Op[i]*us_OHm*us_Op;
        
            /* Already in number density, not mass density */
            total += rhs[j] * -1;
            total_e += us_OHm * -1;
        
        
            total_de += rhs[j] * -1;
        
        j++;
    
        // 
        // Species: us_Hm
        // 
        rhs[j] = -us_Hm+us_Cp[i]*us_Cp*us_Hm - us_Hm+us_Hp[i]*us_Hm*us_Hp - us_Hm+us_O[i]*us_Hm*us_O - us_Hm+us_Op[i]*us_Hm*us_Op;
        
            /* Already in number density, not mass density */
            total += rhs[j] * -1;
            total_e += us_Hm * -1;
        
        
            total_de += rhs[j] * -1;
        
        j++;
    
        // 
        // Species: us_Cp
        // 
        rhs[j] = -us_Cm+us_Cp[i]*us_Cm*us_Cp - us_Hm+us_Cp[i]*us_Cp*us_Hm - us_OHm+us_Cp[i]*us_Cp*us_OHm - us_Om+us_Cp[i]*us_Cp*us_Om;
        
            /* Already in number density, not mass density */
            total += rhs[j] * -1;
            total_e += us_Cp * 1;
        
        
            total_de += rhs[j] * 1;
        
        j++;
    
        // 
        // Species: us_H2
        // 
        rhs[j] = us_H+us_H2p[i]*us_H*us_H2p - us_H2+us_Op[i]*us_H2*us_Op + us_H2p+us_OH[i]*us_H2p*us_OH;
        
            /* Already in number density, not mass density */
            total += rhs[j] * -1;
            total_e += us_H2 * 0;
        
        
            total_de += rhs[j] * 0;
        
        j++;
    
        // 
        // Species: us_H2p
        // 
        rhs[j] = -us_H+us_H2p[i]*us_H*us_H2p - us_H2p+us_OH[i]*us_H2p*us_OH - us_H2p+us_O[i]*us_H2p*us_O;
        
            /* Already in number density, not mass density */
            total += rhs[j] * -1;
            total_e += us_H2p * 1;
        
        
            total_de += rhs[j] * 1;
        
        j++;
    
        // 
        // Species: us_H
        // 
        rhs[j] = us_Cm+us_Hp[i]*us_Cm*us_Hp - us_H+us_H2p[i]*us_H*us_H2p - us_H+us_Om[i]*us_H*us_Om - us_H+us_Op[i]*us_H*us_Op + us_H2+us_Op[i]*us_H2*us_Op + us_H2p+us_O[i]*us_H2p*us_O + us_Hm+us_Cp[i]*us_Cp*us_Hm + 2*us_Hm+us_Hp[i]*us_Hm*us_Hp + us_Hm+us_Op[i]*us_Hm*us_Op + us_Hp+us_OH[i]*us_Hp*us_OH + us_Hp+us_O[i]*us_Hp*us_O + us_OHm+us_Hp[i]*us_Hp*us_OHm + us_Om+us_Hp[i]*us_Hp*us_Om;
        
            /* Already in number density, not mass density */
            total += rhs[j] * -1;
            total_e += us_H * 0;
        
        
            total_de += rhs[j] * 0;
        
        j++;
    
        // 
        // Species: us_Hp
        // 
        rhs[j] = -us_Cm+us_Hp[i]*us_Cm*us_Hp + us_H+us_H2p[i]*us_H*us_H2p + us_H+us_Op[i]*us_H*us_Op - us_Hm+us_Hp[i]*us_Hm*us_Hp - us_Hp+us_OH[i]*us_Hp*us_OH - us_Hp+us_O[i]*us_Hp*us_O - us_OHm+us_Hp[i]*us_Hp*us_OHm - us_Om+us_Hp[i]*us_Hp*us_Om;
        
            /* Already in number density, not mass density */
            total += rhs[j] * -1;
            total_e += us_Hp * 1;
        
        
            total_de += rhs[j] * 1;
        
        j++;
    
        // 
        // Species: us_Op
        // 
        rhs[j] = -us_Cm+us_Op[i]*us_Cm*us_Op - us_H+us_Op[i]*us_H*us_Op - us_H2+us_Op[i]*us_H2*us_Op - us_Hm+us_Op[i]*us_Hm*us_Op + us_Hp+us_O[i]*us_Hp*us_O - us_OHm+us_Op[i]*us_OHm*us_Op - us_Om+us_Op[i]*us_Om*us_Op;
        
            /* Already in number density, not mass density */
            total += rhs[j] * -1;
            total_e += us_Op * 1;
        
        
            total_de += rhs[j] * 1;
        
        j++;
    
        // 
        // Species: us_OHp
        // 
        rhs[j] = us_H2+us_Op[i]*us_H2*us_Op + us_H2p+us_OH[i]*us_H2p*us_OH + us_H2p+us_O[i]*us_H2p*us_O + us_Hp+us_OH[i]*us_Hp*us_OH;
        
            /* Already in number density, not mass density */
            total += rhs[j] * -1;
            total_e += us_OHp * 1;
        
        
            total_de += rhs[j] * 1;
        
        j++;
    
        // 
        // Species: us_OH
        // 
        rhs[j] = us_H+us_Om[i]*us_H*us_Om - us_H2p+us_OH[i]*us_H2p*us_OH + us_Hm+us_O[i]*us_Hm*us_O - us_Hp+us_OH[i]*us_Hp*us_OH + us_OHm+us_Cp[i]*us_Cp*us_OHm + us_OHm+us_Hp[i]*us_Hp*us_OHm + us_OHm+us_Op[i]*us_OHm*us_Op;
        
            /* Already in number density, not mass density */
            total += rhs[j] * -1;
            total_e += us_OH * 0;
        
        
            total_de += rhs[j] * 0;
        
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
    umist_calculate_temperature(data, input, nstrip, nchem);

    umist_interpolate_rates(data, nstrip);

    /* Now we set up some temporaries */
    double *Tge = data->dTs_ge;
    double *us_C+us_Om = data->rs_us_C+us_Om;
    double *rus_C+us_Om = data->drs_us_C+us_Om;
    double *us_Cm+us_Cp = data->rs_us_Cm+us_Cp;
    double *rus_Cm+us_Cp = data->drs_us_Cm+us_Cp;
    double *us_Cm+us_Hp = data->rs_us_Cm+us_Hp;
    double *rus_Cm+us_Hp = data->drs_us_Cm+us_Hp;
    double *us_Cm+us_O = data->rs_us_Cm+us_O;
    double *rus_Cm+us_O = data->drs_us_Cm+us_O;
    double *us_Cm+us_Op = data->rs_us_Cm+us_Op;
    double *rus_Cm+us_Op = data->drs_us_Cm+us_Op;
    double *us_H+us_H2p = data->rs_us_H+us_H2p;
    double *rus_H+us_H2p = data->drs_us_H+us_H2p;
    double *us_H+us_Om = data->rs_us_H+us_Om;
    double *rus_H+us_Om = data->drs_us_H+us_Om;
    double *us_H+us_Op = data->rs_us_H+us_Op;
    double *rus_H+us_Op = data->drs_us_H+us_Op;
    double *us_H2+us_Op = data->rs_us_H2+us_Op;
    double *rus_H2+us_Op = data->drs_us_H2+us_Op;
    double *us_H2p+us_O = data->rs_us_H2p+us_O;
    double *rus_H2p+us_O = data->drs_us_H2p+us_O;
    double *us_H2p+us_OH = data->rs_us_H2p+us_OH;
    double *rus_H2p+us_OH = data->drs_us_H2p+us_OH;
    double *us_Hm+us_Cp = data->rs_us_Hm+us_Cp;
    double *rus_Hm+us_Cp = data->drs_us_Hm+us_Cp;
    double *us_Hm+us_Hp = data->rs_us_Hm+us_Hp;
    double *rus_Hm+us_Hp = data->drs_us_Hm+us_Hp;
    double *us_Hm+us_O = data->rs_us_Hm+us_O;
    double *rus_Hm+us_O = data->drs_us_Hm+us_O;
    double *us_Hm+us_Op = data->rs_us_Hm+us_Op;
    double *rus_Hm+us_Op = data->drs_us_Hm+us_Op;
    double *us_Hp+us_O = data->rs_us_Hp+us_O;
    double *rus_Hp+us_O = data->drs_us_Hp+us_O;
    double *us_Hp+us_OH = data->rs_us_Hp+us_OH;
    double *rus_Hp+us_OH = data->drs_us_Hp+us_OH;
    double *us_OHm+us_Cp = data->rs_us_OHm+us_Cp;
    double *rus_OHm+us_Cp = data->drs_us_OHm+us_Cp;
    double *us_OHm+us_Hp = data->rs_us_OHm+us_Hp;
    double *rus_OHm+us_Hp = data->drs_us_OHm+us_Hp;
    double *us_OHm+us_Op = data->rs_us_OHm+us_Op;
    double *rus_OHm+us_Op = data->drs_us_OHm+us_Op;
    double *us_Om+us_Cp = data->rs_us_Om+us_Cp;
    double *rus_Om+us_Cp = data->drs_us_Om+us_Cp;
    double *us_Om+us_Hp = data->rs_us_Om+us_Hp;
    double *rus_Om+us_Hp = data->drs_us_Om+us_Hp;
    double *us_Om+us_Op = data->rs_us_Om+us_Op;
    double *rus_Om+us_Op = data->drs_us_Om+us_Op;
    double ge;
    double us_Cm;
    double us_CO;
    double us_em;
    double us_O;
    double us_C;
    double us_Om;
    double us_OHm;
    double us_Hm;
    double us_Cp;
    double us_H2;
    double us_H2p;
    double us_H;
    double us_Hp;
    double us_Op;
    double us_OHp;
    double us_OH;

    double mh = 1.67e-24;
    double total, mdensity;
    for (i = 0; i<nstrip; i++) {
        j = i * nchem;
	total = mdensity = 0.0;
	    ge = input[j];
        if (ge < 0.0) {
          fprintf(stderr, "JNegative[%d][ge] = % 0.16g [%d]\n",
            i, ge, j);
          /*ge = 0.0;*/
          ge = 1e-20;
          return 1;
        }
	
        j++;
    
	    us_Cm = input[j];
        if (us_Cm < 0.0) {
          fprintf(stderr, "JNegative[%d][us_Cm] = % 0.16g [%d]\n",
            i, us_Cm, j);
          /*us_Cm = 0.0;*/
          us_Cm = 1e-20;
          return 1;
        }
	
        
          total+=us_Cm * -1;
        
	
        j++;
    
	    us_CO = input[j];
        if (us_CO < 0.0) {
          fprintf(stderr, "JNegative[%d][us_CO] = % 0.16g [%d]\n",
            i, us_CO, j);
          /*us_CO = 0.0;*/
          us_CO = 1e-20;
          return 1;
        }
	
        
          total+=us_CO * -1;
        
	
        j++;
    
	    us_em = input[j];
        if (us_em < 0.0) {
          fprintf(stderr, "JNegative[%d][us_em] = % 0.16g [%d]\n",
            i, us_em, j);
          /*us_em = 0.0;*/
          us_em = 1e-20;
          return 1;
        }
	
        
          total+=us_em * -1;
        
	
        j++;
    
	    us_O = input[j];
        if (us_O < 0.0) {
          fprintf(stderr, "JNegative[%d][us_O] = % 0.16g [%d]\n",
            i, us_O, j);
          /*us_O = 0.0;*/
          us_O = 1e-20;
          return 1;
        }
	
        
          total+=us_O * -1;
        
	
        j++;
    
	    us_C = input[j];
        if (us_C < 0.0) {
          fprintf(stderr, "JNegative[%d][us_C] = % 0.16g [%d]\n",
            i, us_C, j);
          /*us_C = 0.0;*/
          us_C = 1e-20;
          return 1;
        }
	
        
          total+=us_C * -1;
        
	
        j++;
    
	    us_Om = input[j];
        if (us_Om < 0.0) {
          fprintf(stderr, "JNegative[%d][us_Om] = % 0.16g [%d]\n",
            i, us_Om, j);
          /*us_Om = 0.0;*/
          us_Om = 1e-20;
          return 1;
        }
	
        
          total+=us_Om * -1;
        
	
        j++;
    
	    us_OHm = input[j];
        if (us_OHm < 0.0) {
          fprintf(stderr, "JNegative[%d][us_OHm] = % 0.16g [%d]\n",
            i, us_OHm, j);
          /*us_OHm = 0.0;*/
          us_OHm = 1e-20;
          return 1;
        }
	
        
          total+=us_OHm * -1;
        
	
        j++;
    
	    us_Hm = input[j];
        if (us_Hm < 0.0) {
          fprintf(stderr, "JNegative[%d][us_Hm] = % 0.16g [%d]\n",
            i, us_Hm, j);
          /*us_Hm = 0.0;*/
          us_Hm = 1e-20;
          return 1;
        }
	
        
          total+=us_Hm * -1;
        
	
        j++;
    
	    us_Cp = input[j];
        if (us_Cp < 0.0) {
          fprintf(stderr, "JNegative[%d][us_Cp] = % 0.16g [%d]\n",
            i, us_Cp, j);
          /*us_Cp = 0.0;*/
          us_Cp = 1e-20;
          return 1;
        }
	
        
          total+=us_Cp * -1;
        
	
        j++;
    
	    us_H2 = input[j];
        if (us_H2 < 0.0) {
          fprintf(stderr, "JNegative[%d][us_H2] = % 0.16g [%d]\n",
            i, us_H2, j);
          /*us_H2 = 0.0;*/
          us_H2 = 1e-20;
          return 1;
        }
	
        
          total+=us_H2 * -1;
        
	
        j++;
    
	    us_H2p = input[j];
        if (us_H2p < 0.0) {
          fprintf(stderr, "JNegative[%d][us_H2p] = % 0.16g [%d]\n",
            i, us_H2p, j);
          /*us_H2p = 0.0;*/
          us_H2p = 1e-20;
          return 1;
        }
	
        
          total+=us_H2p * -1;
        
	
        j++;
    
	    us_H = input[j];
        if (us_H < 0.0) {
          fprintf(stderr, "JNegative[%d][us_H] = % 0.16g [%d]\n",
            i, us_H, j);
          /*us_H = 0.0;*/
          us_H = 1e-20;
          return 1;
        }
	
        
          total+=us_H * -1;
        
	
        j++;
    
	    us_Hp = input[j];
        if (us_Hp < 0.0) {
          fprintf(stderr, "JNegative[%d][us_Hp] = % 0.16g [%d]\n",
            i, us_Hp, j);
          /*us_Hp = 0.0;*/
          us_Hp = 1e-20;
          return 1;
        }
	
        
          total+=us_Hp * -1;
        
	
        j++;
    
	    us_Op = input[j];
        if (us_Op < 0.0) {
          fprintf(stderr, "JNegative[%d][us_Op] = % 0.16g [%d]\n",
            i, us_Op, j);
          /*us_Op = 0.0;*/
          us_Op = 1e-20;
          return 1;
        }
	
        
          total+=us_Op * -1;
        
	
        j++;
    
	    us_OHp = input[j];
        if (us_OHp < 0.0) {
          fprintf(stderr, "JNegative[%d][us_OHp] = % 0.16g [%d]\n",
            i, us_OHp, j);
          /*us_OHp = 0.0;*/
          us_OHp = 1e-20;
          return 1;
        }
	
        
          total+=us_OHp * -1;
        
	
        j++;
    
	    us_OH = input[j];
        if (us_OH < 0.0) {
          fprintf(stderr, "JNegative[%d][us_OH] = % 0.16g [%d]\n",
            i, us_OH, j);
          /*us_OH = 0.0;*/
          us_OH = 1e-20;
          return 1;
        }
	
        
          total+=us_OH * -1;
        
	
        j++;
    
        mdensity = total * mh;
        
        j = i * nchem * nchem;
        // 
        // Species: ge
        //
            // ge by ge
            Joutput[j] = 0;
	    
	    Joutput[j] /= mdensity;
	    
	    
            Joutput[j] *= Tge[i];
            
            j++;
            // us_Cm by ge
            Joutput[j] = -rus_Cm+us_Cp[i]*us_Cm*us_Cp - rus_Cm+us_Hp[i]*us_Cm*us_Hp - rus_Cm+us_O[i]*us_Cm*us_O - rus_Cm+us_Op[i]*us_Cm*us_Op;
	    
	    
            Joutput[j] *= Tge[i];
            
            j++;
            // us_CO by ge
            Joutput[j] = rus_C+us_Om[i]*us_C*us_Om + rus_Cm+us_O[i]*us_Cm*us_O;
	    
	    
            Joutput[j] *= Tge[i];
            
            j++;
            // us_em by ge
            Joutput[j] = rus_C+us_Om[i]*us_C*us_Om + rus_Cm+us_O[i]*us_Cm*us_O + rus_H+us_Om[i]*us_H*us_Om + rus_Hm+us_O[i]*us_Hm*us_O;
	    
	    
            Joutput[j] *= Tge[i];
            
            j++;
            // us_O by ge
            Joutput[j] = -rus_Cm+us_O[i]*us_Cm*us_O + rus_Cm+us_Op[i]*us_Cm*us_Op + rus_H+us_Op[i]*us_H*us_Op - rus_H2p+us_O[i]*us_H2p*us_O - rus_Hm+us_O[i]*us_Hm*us_O + rus_Hm+us_Op[i]*us_Hm*us_Op - rus_Hp+us_O[i]*us_Hp*us_O + rus_OHm+us_Op[i]*us_OHm*us_Op + rus_Om+us_Cp[i]*us_Cp*us_Om + rus_Om+us_Hp[i]*us_Hp*us_Om + 2*rus_Om+us_Op[i]*us_Om*us_Op;
	    
	    
            Joutput[j] *= Tge[i];
            
            j++;
            // us_C by ge
            Joutput[j] = -rus_C+us_Om[i]*us_C*us_Om + 2*rus_Cm+us_Cp[i]*us_Cm*us_Cp + rus_Cm+us_Hp[i]*us_Cm*us_Hp + rus_Cm+us_Op[i]*us_Cm*us_Op + rus_Hm+us_Cp[i]*us_Cp*us_Hm + rus_OHm+us_Cp[i]*us_Cp*us_OHm + rus_Om+us_Cp[i]*us_Cp*us_Om;
	    
	    
            Joutput[j] *= Tge[i];
            
            j++;
            // us_Om by ge
            Joutput[j] = -rus_C+us_Om[i]*us_C*us_Om - rus_H+us_Om[i]*us_H*us_Om - rus_Om+us_Cp[i]*us_Cp*us_Om - rus_Om+us_Hp[i]*us_Hp*us_Om - rus_Om+us_Op[i]*us_Om*us_Op;
	    
	    
            Joutput[j] *= Tge[i];
            
            j++;
            // us_OHm by ge
            Joutput[j] = -rus_OHm+us_Cp[i]*us_Cp*us_OHm - rus_OHm+us_Hp[i]*us_Hp*us_OHm - rus_OHm+us_Op[i]*us_OHm*us_Op;
	    
	    
            Joutput[j] *= Tge[i];
            
            j++;
            // us_Hm by ge
            Joutput[j] = -rus_Hm+us_Cp[i]*us_Cp*us_Hm - rus_Hm+us_Hp[i]*us_Hm*us_Hp - rus_Hm+us_O[i]*us_Hm*us_O - rus_Hm+us_Op[i]*us_Hm*us_Op;
	    
	    
            Joutput[j] *= Tge[i];
            
            j++;
            // us_Cp by ge
            Joutput[j] = -rus_Cm+us_Cp[i]*us_Cm*us_Cp - rus_Hm+us_Cp[i]*us_Cp*us_Hm - rus_OHm+us_Cp[i]*us_Cp*us_OHm - rus_Om+us_Cp[i]*us_Cp*us_Om;
	    
	    
            Joutput[j] *= Tge[i];
            
            j++;
            // us_H2 by ge
            Joutput[j] = rus_H+us_H2p[i]*us_H*us_H2p - rus_H2+us_Op[i]*us_H2*us_Op + rus_H2p+us_OH[i]*us_H2p*us_OH;
	    
	    
            Joutput[j] *= Tge[i];
            
            j++;
            // us_H2p by ge
            Joutput[j] = -rus_H+us_H2p[i]*us_H*us_H2p - rus_H2p+us_OH[i]*us_H2p*us_OH - rus_H2p+us_O[i]*us_H2p*us_O;
	    
	    
            Joutput[j] *= Tge[i];
            
            j++;
            // us_H by ge
            Joutput[j] = rus_Cm+us_Hp[i]*us_Cm*us_Hp - rus_H+us_H2p[i]*us_H*us_H2p - rus_H+us_Om[i]*us_H*us_Om - rus_H+us_Op[i]*us_H*us_Op + rus_H2+us_Op[i]*us_H2*us_Op + rus_H2p+us_O[i]*us_H2p*us_O + rus_Hm+us_Cp[i]*us_Cp*us_Hm + 2*rus_Hm+us_Hp[i]*us_Hm*us_Hp + rus_Hm+us_Op[i]*us_Hm*us_Op + rus_Hp+us_OH[i]*us_Hp*us_OH + rus_Hp+us_O[i]*us_Hp*us_O + rus_OHm+us_Hp[i]*us_Hp*us_OHm + rus_Om+us_Hp[i]*us_Hp*us_Om;
	    
	    
            Joutput[j] *= Tge[i];
            
            j++;
            // us_Hp by ge
            Joutput[j] = -rus_Cm+us_Hp[i]*us_Cm*us_Hp + rus_H+us_H2p[i]*us_H*us_H2p + rus_H+us_Op[i]*us_H*us_Op - rus_Hm+us_Hp[i]*us_Hm*us_Hp - rus_Hp+us_OH[i]*us_Hp*us_OH - rus_Hp+us_O[i]*us_Hp*us_O - rus_OHm+us_Hp[i]*us_Hp*us_OHm - rus_Om+us_Hp[i]*us_Hp*us_Om;
	    
	    
            Joutput[j] *= Tge[i];
            
            j++;
            // us_Op by ge
            Joutput[j] = -rus_Cm+us_Op[i]*us_Cm*us_Op - rus_H+us_Op[i]*us_H*us_Op - rus_H2+us_Op[i]*us_H2*us_Op - rus_Hm+us_Op[i]*us_Hm*us_Op + rus_Hp+us_O[i]*us_Hp*us_O - rus_OHm+us_Op[i]*us_OHm*us_Op - rus_Om+us_Op[i]*us_Om*us_Op;
	    
	    
            Joutput[j] *= Tge[i];
            
            j++;
            // us_OHp by ge
            Joutput[j] = rus_H2+us_Op[i]*us_H2*us_Op + rus_H2p+us_OH[i]*us_H2p*us_OH + rus_H2p+us_O[i]*us_H2p*us_O + rus_Hp+us_OH[i]*us_Hp*us_OH;
	    
	    
            Joutput[j] *= Tge[i];
            
            j++;
            // us_OH by ge
            Joutput[j] = rus_H+us_Om[i]*us_H*us_Om - rus_H2p+us_OH[i]*us_H2p*us_OH + rus_Hm+us_O[i]*us_Hm*us_O - rus_Hp+us_OH[i]*us_Hp*us_OH + rus_OHm+us_Cp[i]*us_Cp*us_OHm + rus_OHm+us_Hp[i]*us_Hp*us_OHm + rus_OHm+us_Op[i]*us_OHm*us_Op;
	    
	    
            Joutput[j] *= Tge[i];
            
            j++;
    
        // 
        // Species: us_Cm
        //
            // ge by us_Cm
            Joutput[j] = 0;
	    
	    Joutput[j] /= mdensity;
	    
	    
            j++;
            // us_Cm by us_Cm
            Joutput[j] = -us_Cm+us_Cp[i]*us_Cp - us_Cm+us_Hp[i]*us_Hp - us_Cm+us_O[i]*us_O - us_Cm+us_Op[i]*us_Op;
	    
	    
            j++;
            // us_CO by us_Cm
            Joutput[j] = us_Cm+us_O[i]*us_O;
	    
	    
            j++;
            // us_em by us_Cm
            Joutput[j] = us_Cm+us_O[i]*us_O;
	    
	    
            j++;
            // us_O by us_Cm
            Joutput[j] = -us_Cm+us_O[i]*us_O + us_Cm+us_Op[i]*us_Op;
	    
	    
            j++;
            // us_C by us_Cm
            Joutput[j] = 2*us_Cm+us_Cp[i]*us_Cp + us_Cm+us_Hp[i]*us_Hp + us_Cm+us_Op[i]*us_Op;
	    
	    
            j++;
            // us_Om by us_Cm
            Joutput[j] = 0;
	    
	    
            j++;
            // us_OHm by us_Cm
            Joutput[j] = 0;
	    
	    
            j++;
            // us_Hm by us_Cm
            Joutput[j] = 0;
	    
	    
            j++;
            // us_Cp by us_Cm
            Joutput[j] = -us_Cm+us_Cp[i]*us_Cp;
	    
	    
            j++;
            // us_H2 by us_Cm
            Joutput[j] = 0;
	    
	    
            j++;
            // us_H2p by us_Cm
            Joutput[j] = 0;
	    
	    
            j++;
            // us_H by us_Cm
            Joutput[j] = us_Cm+us_Hp[i]*us_Hp;
	    
	    
            j++;
            // us_Hp by us_Cm
            Joutput[j] = -us_Cm+us_Hp[i]*us_Hp;
	    
	    
            j++;
            // us_Op by us_Cm
            Joutput[j] = -us_Cm+us_Op[i]*us_Op;
	    
	    
            j++;
            // us_OHp by us_Cm
            Joutput[j] = 0;
	    
	    
            j++;
            // us_OH by us_Cm
            Joutput[j] = 0;
	    
	    
            j++;
    
        // 
        // Species: us_CO
        //
            // ge by us_CO
            Joutput[j] = 0;
	    
	    Joutput[j] /= mdensity;
	    
	    
            j++;
            // us_Cm by us_CO
            Joutput[j] = 0;
	    
	    
            j++;
            // us_CO by us_CO
            Joutput[j] = 0;
	    
	    
            j++;
            // us_em by us_CO
            Joutput[j] = 0;
	    
	    
            j++;
            // us_O by us_CO
            Joutput[j] = 0;
	    
	    
            j++;
            // us_C by us_CO
            Joutput[j] = 0;
	    
	    
            j++;
            // us_Om by us_CO
            Joutput[j] = 0;
	    
	    
            j++;
            // us_OHm by us_CO
            Joutput[j] = 0;
	    
	    
            j++;
            // us_Hm by us_CO
            Joutput[j] = 0;
	    
	    
            j++;
            // us_Cp by us_CO
            Joutput[j] = 0;
	    
	    
            j++;
            // us_H2 by us_CO
            Joutput[j] = 0;
	    
	    
            j++;
            // us_H2p by us_CO
            Joutput[j] = 0;
	    
	    
            j++;
            // us_H by us_CO
            Joutput[j] = 0;
	    
	    
            j++;
            // us_Hp by us_CO
            Joutput[j] = 0;
	    
	    
            j++;
            // us_Op by us_CO
            Joutput[j] = 0;
	    
	    
            j++;
            // us_OHp by us_CO
            Joutput[j] = 0;
	    
	    
            j++;
            // us_OH by us_CO
            Joutput[j] = 0;
	    
	    
            j++;
    
        // 
        // Species: us_em
        //
            // ge by us_em
            Joutput[j] = 0;
	    
	    Joutput[j] /= mdensity;
	    
	    
            j++;
            // us_Cm by us_em
            Joutput[j] = 0;
	    
	    
            j++;
            // us_CO by us_em
            Joutput[j] = 0;
	    
	    
            j++;
            // us_em by us_em
            Joutput[j] = 0;
	    
	    
            j++;
            // us_O by us_em
            Joutput[j] = 0;
	    
	    
            j++;
            // us_C by us_em
            Joutput[j] = 0;
	    
	    
            j++;
            // us_Om by us_em
            Joutput[j] = 0;
	    
	    
            j++;
            // us_OHm by us_em
            Joutput[j] = 0;
	    
	    
            j++;
            // us_Hm by us_em
            Joutput[j] = 0;
	    
	    
            j++;
            // us_Cp by us_em
            Joutput[j] = 0;
	    
	    
            j++;
            // us_H2 by us_em
            Joutput[j] = 0;
	    
	    
            j++;
            // us_H2p by us_em
            Joutput[j] = 0;
	    
	    
            j++;
            // us_H by us_em
            Joutput[j] = 0;
	    
	    
            j++;
            // us_Hp by us_em
            Joutput[j] = 0;
	    
	    
            j++;
            // us_Op by us_em
            Joutput[j] = 0;
	    
	    
            j++;
            // us_OHp by us_em
            Joutput[j] = 0;
	    
	    
            j++;
            // us_OH by us_em
            Joutput[j] = 0;
	    
	    
            j++;
    
        // 
        // Species: us_O
        //
            // ge by us_O
            Joutput[j] = 0;
	    
	    Joutput[j] /= mdensity;
	    
	    
            j++;
            // us_Cm by us_O
            Joutput[j] = -us_Cm+us_O[i]*us_Cm;
	    
	    
            j++;
            // us_CO by us_O
            Joutput[j] = us_Cm+us_O[i]*us_Cm;
	    
	    
            j++;
            // us_em by us_O
            Joutput[j] = us_Cm+us_O[i]*us_Cm + us_Hm+us_O[i]*us_Hm;
	    
	    
            j++;
            // us_O by us_O
            Joutput[j] = -us_Cm+us_O[i]*us_Cm - us_H2p+us_O[i]*us_H2p - us_Hm+us_O[i]*us_Hm - us_Hp+us_O[i]*us_Hp;
	    
	    
            j++;
            // us_C by us_O
            Joutput[j] = 0;
	    
	    
            j++;
            // us_Om by us_O
            Joutput[j] = 0;
	    
	    
            j++;
            // us_OHm by us_O
            Joutput[j] = 0;
	    
	    
            j++;
            // us_Hm by us_O
            Joutput[j] = -us_Hm+us_O[i]*us_Hm;
	    
	    
            j++;
            // us_Cp by us_O
            Joutput[j] = 0;
	    
	    
            j++;
            // us_H2 by us_O
            Joutput[j] = 0;
	    
	    
            j++;
            // us_H2p by us_O
            Joutput[j] = -us_H2p+us_O[i]*us_H2p;
	    
	    
            j++;
            // us_H by us_O
            Joutput[j] = us_H2p+us_O[i]*us_H2p + us_Hp+us_O[i]*us_Hp;
	    
	    
            j++;
            // us_Hp by us_O
            Joutput[j] = -us_Hp+us_O[i]*us_Hp;
	    
	    
            j++;
            // us_Op by us_O
            Joutput[j] = us_Hp+us_O[i]*us_Hp;
	    
	    
            j++;
            // us_OHp by us_O
            Joutput[j] = us_H2p+us_O[i]*us_H2p;
	    
	    
            j++;
            // us_OH by us_O
            Joutput[j] = us_Hm+us_O[i]*us_Hm;
	    
	    
            j++;
    
        // 
        // Species: us_C
        //
            // ge by us_C
            Joutput[j] = 0;
	    
	    Joutput[j] /= mdensity;
	    
	    
            j++;
            // us_Cm by us_C
            Joutput[j] = 0;
	    
	    
            j++;
            // us_CO by us_C
            Joutput[j] = us_C+us_Om[i]*us_Om;
	    
	    
            j++;
            // us_em by us_C
            Joutput[j] = us_C+us_Om[i]*us_Om;
	    
	    
            j++;
            // us_O by us_C
            Joutput[j] = 0;
	    
	    
            j++;
            // us_C by us_C
            Joutput[j] = -us_C+us_Om[i]*us_Om;
	    
	    
            j++;
            // us_Om by us_C
            Joutput[j] = -us_C+us_Om[i]*us_Om;
	    
	    
            j++;
            // us_OHm by us_C
            Joutput[j] = 0;
	    
	    
            j++;
            // us_Hm by us_C
            Joutput[j] = 0;
	    
	    
            j++;
            // us_Cp by us_C
            Joutput[j] = 0;
	    
	    
            j++;
            // us_H2 by us_C
            Joutput[j] = 0;
	    
	    
            j++;
            // us_H2p by us_C
            Joutput[j] = 0;
	    
	    
            j++;
            // us_H by us_C
            Joutput[j] = 0;
	    
	    
            j++;
            // us_Hp by us_C
            Joutput[j] = 0;
	    
	    
            j++;
            // us_Op by us_C
            Joutput[j] = 0;
	    
	    
            j++;
            // us_OHp by us_C
            Joutput[j] = 0;
	    
	    
            j++;
            // us_OH by us_C
            Joutput[j] = 0;
	    
	    
            j++;
    
        // 
        // Species: us_Om
        //
            // ge by us_Om
            Joutput[j] = 0;
	    
	    Joutput[j] /= mdensity;
	    
	    
            j++;
            // us_Cm by us_Om
            Joutput[j] = 0;
	    
	    
            j++;
            // us_CO by us_Om
            Joutput[j] = us_C+us_Om[i]*us_C;
	    
	    
            j++;
            // us_em by us_Om
            Joutput[j] = us_C+us_Om[i]*us_C + us_H+us_Om[i]*us_H;
	    
	    
            j++;
            // us_O by us_Om
            Joutput[j] = us_Om+us_Cp[i]*us_Cp + us_Om+us_Hp[i]*us_Hp + 2*us_Om+us_Op[i]*us_Op;
	    
	    
            j++;
            // us_C by us_Om
            Joutput[j] = -us_C+us_Om[i]*us_C + us_Om+us_Cp[i]*us_Cp;
	    
	    
            j++;
            // us_Om by us_Om
            Joutput[j] = -us_C+us_Om[i]*us_C - us_H+us_Om[i]*us_H - us_Om+us_Cp[i]*us_Cp - us_Om+us_Hp[i]*us_Hp - us_Om+us_Op[i]*us_Op;
	    
	    
            j++;
            // us_OHm by us_Om
            Joutput[j] = 0;
	    
	    
            j++;
            // us_Hm by us_Om
            Joutput[j] = 0;
	    
	    
            j++;
            // us_Cp by us_Om
            Joutput[j] = -us_Om+us_Cp[i]*us_Cp;
	    
	    
            j++;
            // us_H2 by us_Om
            Joutput[j] = 0;
	    
	    
            j++;
            // us_H2p by us_Om
            Joutput[j] = 0;
	    
	    
            j++;
            // us_H by us_Om
            Joutput[j] = -us_H+us_Om[i]*us_H + us_Om+us_Hp[i]*us_Hp;
	    
	    
            j++;
            // us_Hp by us_Om
            Joutput[j] = -us_Om+us_Hp[i]*us_Hp;
	    
	    
            j++;
            // us_Op by us_Om
            Joutput[j] = -us_Om+us_Op[i]*us_Op;
	    
	    
            j++;
            // us_OHp by us_Om
            Joutput[j] = 0;
	    
	    
            j++;
            // us_OH by us_Om
            Joutput[j] = us_H+us_Om[i]*us_H;
	    
	    
            j++;
    
        // 
        // Species: us_OHm
        //
            // ge by us_OHm
            Joutput[j] = 0;
	    
	    Joutput[j] /= mdensity;
	    
	    
            j++;
            // us_Cm by us_OHm
            Joutput[j] = 0;
	    
	    
            j++;
            // us_CO by us_OHm
            Joutput[j] = 0;
	    
	    
            j++;
            // us_em by us_OHm
            Joutput[j] = 0;
	    
	    
            j++;
            // us_O by us_OHm
            Joutput[j] = us_OHm+us_Op[i]*us_Op;
	    
	    
            j++;
            // us_C by us_OHm
            Joutput[j] = us_OHm+us_Cp[i]*us_Cp;
	    
	    
            j++;
            // us_Om by us_OHm
            Joutput[j] = 0;
	    
	    
            j++;
            // us_OHm by us_OHm
            Joutput[j] = -us_OHm+us_Cp[i]*us_Cp - us_OHm+us_Hp[i]*us_Hp - us_OHm+us_Op[i]*us_Op;
	    
	    
            j++;
            // us_Hm by us_OHm
            Joutput[j] = 0;
	    
	    
            j++;
            // us_Cp by us_OHm
            Joutput[j] = -us_OHm+us_Cp[i]*us_Cp;
	    
	    
            j++;
            // us_H2 by us_OHm
            Joutput[j] = 0;
	    
	    
            j++;
            // us_H2p by us_OHm
            Joutput[j] = 0;
	    
	    
            j++;
            // us_H by us_OHm
            Joutput[j] = us_OHm+us_Hp[i]*us_Hp;
	    
	    
            j++;
            // us_Hp by us_OHm
            Joutput[j] = -us_OHm+us_Hp[i]*us_Hp;
	    
	    
            j++;
            // us_Op by us_OHm
            Joutput[j] = -us_OHm+us_Op[i]*us_Op;
	    
	    
            j++;
            // us_OHp by us_OHm
            Joutput[j] = 0;
	    
	    
            j++;
            // us_OH by us_OHm
            Joutput[j] = us_OHm+us_Cp[i]*us_Cp + us_OHm+us_Hp[i]*us_Hp + us_OHm+us_Op[i]*us_Op;
	    
	    
            j++;
    
        // 
        // Species: us_Hm
        //
            // ge by us_Hm
            Joutput[j] = 0;
	    
	    Joutput[j] /= mdensity;
	    
	    
            j++;
            // us_Cm by us_Hm
            Joutput[j] = 0;
	    
	    
            j++;
            // us_CO by us_Hm
            Joutput[j] = 0;
	    
	    
            j++;
            // us_em by us_Hm
            Joutput[j] = us_Hm+us_O[i]*us_O;
	    
	    
            j++;
            // us_O by us_Hm
            Joutput[j] = -us_Hm+us_O[i]*us_O + us_Hm+us_Op[i]*us_Op;
	    
	    
            j++;
            // us_C by us_Hm
            Joutput[j] = us_Hm+us_Cp[i]*us_Cp;
	    
	    
            j++;
            // us_Om by us_Hm
            Joutput[j] = 0;
	    
	    
            j++;
            // us_OHm by us_Hm
            Joutput[j] = 0;
	    
	    
            j++;
            // us_Hm by us_Hm
            Joutput[j] = -us_Hm+us_Cp[i]*us_Cp - us_Hm+us_Hp[i]*us_Hp - us_Hm+us_O[i]*us_O - us_Hm+us_Op[i]*us_Op;
	    
	    
            j++;
            // us_Cp by us_Hm
            Joutput[j] = -us_Hm+us_Cp[i]*us_Cp;
	    
	    
            j++;
            // us_H2 by us_Hm
            Joutput[j] = 0;
	    
	    
            j++;
            // us_H2p by us_Hm
            Joutput[j] = 0;
	    
	    
            j++;
            // us_H by us_Hm
            Joutput[j] = us_Hm+us_Cp[i]*us_Cp + 2*us_Hm+us_Hp[i]*us_Hp + us_Hm+us_Op[i]*us_Op;
	    
	    
            j++;
            // us_Hp by us_Hm
            Joutput[j] = -us_Hm+us_Hp[i]*us_Hp;
	    
	    
            j++;
            // us_Op by us_Hm
            Joutput[j] = -us_Hm+us_Op[i]*us_Op;
	    
	    
            j++;
            // us_OHp by us_Hm
            Joutput[j] = 0;
	    
	    
            j++;
            // us_OH by us_Hm
            Joutput[j] = us_Hm+us_O[i]*us_O;
	    
	    
            j++;
    
        // 
        // Species: us_Cp
        //
            // ge by us_Cp
            Joutput[j] = 0;
	    
	    Joutput[j] /= mdensity;
	    
	    
            j++;
            // us_Cm by us_Cp
            Joutput[j] = -us_Cm+us_Cp[i]*us_Cm;
	    
	    
            j++;
            // us_CO by us_Cp
            Joutput[j] = 0;
	    
	    
            j++;
            // us_em by us_Cp
            Joutput[j] = 0;
	    
	    
            j++;
            // us_O by us_Cp
            Joutput[j] = us_Om+us_Cp[i]*us_Om;
	    
	    
            j++;
            // us_C by us_Cp
            Joutput[j] = 2*us_Cm+us_Cp[i]*us_Cm + us_Hm+us_Cp[i]*us_Hm + us_OHm+us_Cp[i]*us_OHm + us_Om+us_Cp[i]*us_Om;
	    
	    
            j++;
            // us_Om by us_Cp
            Joutput[j] = -us_Om+us_Cp[i]*us_Om;
	    
	    
            j++;
            // us_OHm by us_Cp
            Joutput[j] = -us_OHm+us_Cp[i]*us_OHm;
	    
	    
            j++;
            // us_Hm by us_Cp
            Joutput[j] = -us_Hm+us_Cp[i]*us_Hm;
	    
	    
            j++;
            // us_Cp by us_Cp
            Joutput[j] = -us_Cm+us_Cp[i]*us_Cm - us_Hm+us_Cp[i]*us_Hm - us_OHm+us_Cp[i]*us_OHm - us_Om+us_Cp[i]*us_Om;
	    
	    
            j++;
            // us_H2 by us_Cp
            Joutput[j] = 0;
	    
	    
            j++;
            // us_H2p by us_Cp
            Joutput[j] = 0;
	    
	    
            j++;
            // us_H by us_Cp
            Joutput[j] = us_Hm+us_Cp[i]*us_Hm;
	    
	    
            j++;
            // us_Hp by us_Cp
            Joutput[j] = 0;
	    
	    
            j++;
            // us_Op by us_Cp
            Joutput[j] = 0;
	    
	    
            j++;
            // us_OHp by us_Cp
            Joutput[j] = 0;
	    
	    
            j++;
            // us_OH by us_Cp
            Joutput[j] = us_OHm+us_Cp[i]*us_OHm;
	    
	    
            j++;
    
        // 
        // Species: us_H2
        //
            // ge by us_H2
            Joutput[j] = 0;
	    
	    Joutput[j] /= mdensity;
	    
	    
            j++;
            // us_Cm by us_H2
            Joutput[j] = 0;
	    
	    
            j++;
            // us_CO by us_H2
            Joutput[j] = 0;
	    
	    
            j++;
            // us_em by us_H2
            Joutput[j] = 0;
	    
	    
            j++;
            // us_O by us_H2
            Joutput[j] = 0;
	    
	    
            j++;
            // us_C by us_H2
            Joutput[j] = 0;
	    
	    
            j++;
            // us_Om by us_H2
            Joutput[j] = 0;
	    
	    
            j++;
            // us_OHm by us_H2
            Joutput[j] = 0;
	    
	    
            j++;
            // us_Hm by us_H2
            Joutput[j] = 0;
	    
	    
            j++;
            // us_Cp by us_H2
            Joutput[j] = 0;
	    
	    
            j++;
            // us_H2 by us_H2
            Joutput[j] = -us_H2+us_Op[i]*us_Op;
	    
	    
            j++;
            // us_H2p by us_H2
            Joutput[j] = 0;
	    
	    
            j++;
            // us_H by us_H2
            Joutput[j] = us_H2+us_Op[i]*us_Op;
	    
	    
            j++;
            // us_Hp by us_H2
            Joutput[j] = 0;
	    
	    
            j++;
            // us_Op by us_H2
            Joutput[j] = -us_H2+us_Op[i]*us_Op;
	    
	    
            j++;
            // us_OHp by us_H2
            Joutput[j] = us_H2+us_Op[i]*us_Op;
	    
	    
            j++;
            // us_OH by us_H2
            Joutput[j] = 0;
	    
	    
            j++;
    
        // 
        // Species: us_H2p
        //
            // ge by us_H2p
            Joutput[j] = 0;
	    
	    Joutput[j] /= mdensity;
	    
	    
            j++;
            // us_Cm by us_H2p
            Joutput[j] = 0;
	    
	    
            j++;
            // us_CO by us_H2p
            Joutput[j] = 0;
	    
	    
            j++;
            // us_em by us_H2p
            Joutput[j] = 0;
	    
	    
            j++;
            // us_O by us_H2p
            Joutput[j] = -us_H2p+us_O[i]*us_O;
	    
	    
            j++;
            // us_C by us_H2p
            Joutput[j] = 0;
	    
	    
            j++;
            // us_Om by us_H2p
            Joutput[j] = 0;
	    
	    
            j++;
            // us_OHm by us_H2p
            Joutput[j] = 0;
	    
	    
            j++;
            // us_Hm by us_H2p
            Joutput[j] = 0;
	    
	    
            j++;
            // us_Cp by us_H2p
            Joutput[j] = 0;
	    
	    
            j++;
            // us_H2 by us_H2p
            Joutput[j] = us_H+us_H2p[i]*us_H + us_H2p+us_OH[i]*us_OH;
	    
	    
            j++;
            // us_H2p by us_H2p
            Joutput[j] = -us_H+us_H2p[i]*us_H - us_H2p+us_OH[i]*us_OH - us_H2p+us_O[i]*us_O;
	    
	    
            j++;
            // us_H by us_H2p
            Joutput[j] = -us_H+us_H2p[i]*us_H + us_H2p+us_O[i]*us_O;
	    
	    
            j++;
            // us_Hp by us_H2p
            Joutput[j] = us_H+us_H2p[i]*us_H;
	    
	    
            j++;
            // us_Op by us_H2p
            Joutput[j] = 0;
	    
	    
            j++;
            // us_OHp by us_H2p
            Joutput[j] = us_H2p+us_OH[i]*us_OH + us_H2p+us_O[i]*us_O;
	    
	    
            j++;
            // us_OH by us_H2p
            Joutput[j] = -us_H2p+us_OH[i]*us_OH;
	    
	    
            j++;
    
        // 
        // Species: us_H
        //
            // ge by us_H
            Joutput[j] = 0;
	    
	    Joutput[j] /= mdensity;
	    
	    
            j++;
            // us_Cm by us_H
            Joutput[j] = 0;
	    
	    
            j++;
            // us_CO by us_H
            Joutput[j] = 0;
	    
	    
            j++;
            // us_em by us_H
            Joutput[j] = us_H+us_Om[i]*us_Om;
	    
	    
            j++;
            // us_O by us_H
            Joutput[j] = us_H+us_Op[i]*us_Op;
	    
	    
            j++;
            // us_C by us_H
            Joutput[j] = 0;
	    
	    
            j++;
            // us_Om by us_H
            Joutput[j] = -us_H+us_Om[i]*us_Om;
	    
	    
            j++;
            // us_OHm by us_H
            Joutput[j] = 0;
	    
	    
            j++;
            // us_Hm by us_H
            Joutput[j] = 0;
	    
	    
            j++;
            // us_Cp by us_H
            Joutput[j] = 0;
	    
	    
            j++;
            // us_H2 by us_H
            Joutput[j] = us_H+us_H2p[i]*us_H2p;
	    
	    
            j++;
            // us_H2p by us_H
            Joutput[j] = -us_H+us_H2p[i]*us_H2p;
	    
	    
            j++;
            // us_H by us_H
            Joutput[j] = -us_H+us_H2p[i]*us_H2p - us_H+us_Om[i]*us_Om - us_H+us_Op[i]*us_Op;
	    
	    
            j++;
            // us_Hp by us_H
            Joutput[j] = us_H+us_H2p[i]*us_H2p + us_H+us_Op[i]*us_Op;
	    
	    
            j++;
            // us_Op by us_H
            Joutput[j] = -us_H+us_Op[i]*us_Op;
	    
	    
            j++;
            // us_OHp by us_H
            Joutput[j] = 0;
	    
	    
            j++;
            // us_OH by us_H
            Joutput[j] = us_H+us_Om[i]*us_Om;
	    
	    
            j++;
    
        // 
        // Species: us_Hp
        //
            // ge by us_Hp
            Joutput[j] = 0;
	    
	    Joutput[j] /= mdensity;
	    
	    
            j++;
            // us_Cm by us_Hp
            Joutput[j] = -us_Cm+us_Hp[i]*us_Cm;
	    
	    
            j++;
            // us_CO by us_Hp
            Joutput[j] = 0;
	    
	    
            j++;
            // us_em by us_Hp
            Joutput[j] = 0;
	    
	    
            j++;
            // us_O by us_Hp
            Joutput[j] = -us_Hp+us_O[i]*us_O + us_Om+us_Hp[i]*us_Om;
	    
	    
            j++;
            // us_C by us_Hp
            Joutput[j] = us_Cm+us_Hp[i]*us_Cm;
	    
	    
            j++;
            // us_Om by us_Hp
            Joutput[j] = -us_Om+us_Hp[i]*us_Om;
	    
	    
            j++;
            // us_OHm by us_Hp
            Joutput[j] = -us_OHm+us_Hp[i]*us_OHm;
	    
	    
            j++;
            // us_Hm by us_Hp
            Joutput[j] = -us_Hm+us_Hp[i]*us_Hm;
	    
	    
            j++;
            // us_Cp by us_Hp
            Joutput[j] = 0;
	    
	    
            j++;
            // us_H2 by us_Hp
            Joutput[j] = 0;
	    
	    
            j++;
            // us_H2p by us_Hp
            Joutput[j] = 0;
	    
	    
            j++;
            // us_H by us_Hp
            Joutput[j] = us_Cm+us_Hp[i]*us_Cm + 2*us_Hm+us_Hp[i]*us_Hm + us_Hp+us_OH[i]*us_OH + us_Hp+us_O[i]*us_O + us_OHm+us_Hp[i]*us_OHm + us_Om+us_Hp[i]*us_Om;
	    
	    
            j++;
            // us_Hp by us_Hp
            Joutput[j] = -us_Cm+us_Hp[i]*us_Cm - us_Hm+us_Hp[i]*us_Hm - us_Hp+us_OH[i]*us_OH - us_Hp+us_O[i]*us_O - us_OHm+us_Hp[i]*us_OHm - us_Om+us_Hp[i]*us_Om;
	    
	    
            j++;
            // us_Op by us_Hp
            Joutput[j] = us_Hp+us_O[i]*us_O;
	    
	    
            j++;
            // us_OHp by us_Hp
            Joutput[j] = us_Hp+us_OH[i]*us_OH;
	    
	    
            j++;
            // us_OH by us_Hp
            Joutput[j] = -us_Hp+us_OH[i]*us_OH + us_OHm+us_Hp[i]*us_OHm;
	    
	    
            j++;
    
        // 
        // Species: us_Op
        //
            // ge by us_Op
            Joutput[j] = 0;
	    
	    Joutput[j] /= mdensity;
	    
	    
            j++;
            // us_Cm by us_Op
            Joutput[j] = -us_Cm+us_Op[i]*us_Cm;
	    
	    
            j++;
            // us_CO by us_Op
            Joutput[j] = 0;
	    
	    
            j++;
            // us_em by us_Op
            Joutput[j] = 0;
	    
	    
            j++;
            // us_O by us_Op
            Joutput[j] = us_Cm+us_Op[i]*us_Cm + us_H+us_Op[i]*us_H + us_Hm+us_Op[i]*us_Hm + us_OHm+us_Op[i]*us_OHm + 2*us_Om+us_Op[i]*us_Om;
	    
	    
            j++;
            // us_C by us_Op
            Joutput[j] = us_Cm+us_Op[i]*us_Cm;
	    
	    
            j++;
            // us_Om by us_Op
            Joutput[j] = -us_Om+us_Op[i]*us_Om;
	    
	    
            j++;
            // us_OHm by us_Op
            Joutput[j] = -us_OHm+us_Op[i]*us_OHm;
	    
	    
            j++;
            // us_Hm by us_Op
            Joutput[j] = -us_Hm+us_Op[i]*us_Hm;
	    
	    
            j++;
            // us_Cp by us_Op
            Joutput[j] = 0;
	    
	    
            j++;
            // us_H2 by us_Op
            Joutput[j] = -us_H2+us_Op[i]*us_H2;
	    
	    
            j++;
            // us_H2p by us_Op
            Joutput[j] = 0;
	    
	    
            j++;
            // us_H by us_Op
            Joutput[j] = -us_H+us_Op[i]*us_H + us_H2+us_Op[i]*us_H2 + us_Hm+us_Op[i]*us_Hm;
	    
	    
            j++;
            // us_Hp by us_Op
            Joutput[j] = us_H+us_Op[i]*us_H;
	    
	    
            j++;
            // us_Op by us_Op
            Joutput[j] = -us_Cm+us_Op[i]*us_Cm - us_H+us_Op[i]*us_H - us_H2+us_Op[i]*us_H2 - us_Hm+us_Op[i]*us_Hm - us_OHm+us_Op[i]*us_OHm - us_Om+us_Op[i]*us_Om;
	    
	    
            j++;
            // us_OHp by us_Op
            Joutput[j] = us_H2+us_Op[i]*us_H2;
	    
	    
            j++;
            // us_OH by us_Op
            Joutput[j] = us_OHm+us_Op[i]*us_OHm;
	    
	    
            j++;
    
        // 
        // Species: us_OHp
        //
            // ge by us_OHp
            Joutput[j] = 0;
	    
	    Joutput[j] /= mdensity;
	    
	    
            j++;
            // us_Cm by us_OHp
            Joutput[j] = 0;
	    
	    
            j++;
            // us_CO by us_OHp
            Joutput[j] = 0;
	    
	    
            j++;
            // us_em by us_OHp
            Joutput[j] = 0;
	    
	    
            j++;
            // us_O by us_OHp
            Joutput[j] = 0;
	    
	    
            j++;
            // us_C by us_OHp
            Joutput[j] = 0;
	    
	    
            j++;
            // us_Om by us_OHp
            Joutput[j] = 0;
	    
	    
            j++;
            // us_OHm by us_OHp
            Joutput[j] = 0;
	    
	    
            j++;
            // us_Hm by us_OHp
            Joutput[j] = 0;
	    
	    
            j++;
            // us_Cp by us_OHp
            Joutput[j] = 0;
	    
	    
            j++;
            // us_H2 by us_OHp
            Joutput[j] = 0;
	    
	    
            j++;
            // us_H2p by us_OHp
            Joutput[j] = 0;
	    
	    
            j++;
            // us_H by us_OHp
            Joutput[j] = 0;
	    
	    
            j++;
            // us_Hp by us_OHp
            Joutput[j] = 0;
	    
	    
            j++;
            // us_Op by us_OHp
            Joutput[j] = 0;
	    
	    
            j++;
            // us_OHp by us_OHp
            Joutput[j] = 0;
	    
	    
            j++;
            // us_OH by us_OHp
            Joutput[j] = 0;
	    
	    
            j++;
    
        // 
        // Species: us_OH
        //
            // ge by us_OH
            Joutput[j] = 0;
	    
	    Joutput[j] /= mdensity;
	    
	    
            j++;
            // us_Cm by us_OH
            Joutput[j] = 0;
	    
	    
            j++;
            // us_CO by us_OH
            Joutput[j] = 0;
	    
	    
            j++;
            // us_em by us_OH
            Joutput[j] = 0;
	    
	    
            j++;
            // us_O by us_OH
            Joutput[j] = 0;
	    
	    
            j++;
            // us_C by us_OH
            Joutput[j] = 0;
	    
	    
            j++;
            // us_Om by us_OH
            Joutput[j] = 0;
	    
	    
            j++;
            // us_OHm by us_OH
            Joutput[j] = 0;
	    
	    
            j++;
            // us_Hm by us_OH
            Joutput[j] = 0;
	    
	    
            j++;
            // us_Cp by us_OH
            Joutput[j] = 0;
	    
	    
            j++;
            // us_H2 by us_OH
            Joutput[j] = us_H2p+us_OH[i]*us_H2p;
	    
	    
            j++;
            // us_H2p by us_OH
            Joutput[j] = -us_H2p+us_OH[i]*us_H2p;
	    
	    
            j++;
            // us_H by us_OH
            Joutput[j] = us_Hp+us_OH[i]*us_Hp;
	    
	    
            j++;
            // us_Hp by us_OH
            Joutput[j] = -us_Hp+us_OH[i]*us_Hp;
	    
	    
            j++;
            // us_Op by us_OH
            Joutput[j] = 0;
	    
	    
            j++;
            // us_OHp by us_OH
            Joutput[j] = us_H2p+us_OH[i]*us_H2p + us_Hp+us_OH[i]*us_Hp;
	    
	    
            j++;
            // us_OH by us_OH
            Joutput[j] = -us_H2p+us_OH[i]*us_H2p - us_Hp+us_OH[i]*us_Hp;
	    
	    
            j++;
    
    }

    return 0;
    
}
