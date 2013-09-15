
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


#include "oxygen_solver.h"

oxygen_data *oxygen_setup_data(void) {

    oxygen_data *data = (oxygen_data *) malloc(sizeof(oxygen_data));

    data->bounds[0] = 1.0;
    data->bounds[1] = 100000000.0;
    data->nbins = 1024;
    data->dbin = (log(data->bounds[1]) - log(data->bounds[0])) / data->nbins;
    data->idbin = 1.0L / data->dbin;
    
    oxygen_read_rate_tables(data);
    fprintf(stderr, "Successfully read in rate tables.\n");

    oxygen_read_cooling_tables(data);
    fprintf(stderr, "Successfully read in cooling rate tables.\n");

    return data;

}


int oxygen_main(int argc, char** argv)
{
    oxygen_data *data = oxygen_setup_data();

    /* Initial conditions */

    hid_t file_id = H5Fopen("oxygen_initial_conditions.h5", H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file_id < 0) {fprintf(stderr, "Failed to open "
        "oxygen_initial_conditions.h5 so dying.\n");
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

    int N = 11;

    double *atol, *rtol;
    atol = (double *) alloca(N * dims * sizeof(double));
    rtol = (double *) alloca(N * dims * sizeof(double));

    double *tics = (double *) alloca(dims * sizeof(double));
    double *ics = (double *) alloca(dims * N * sizeof(double));
    double *input = (double *) alloca(dims * N * sizeof(double));
    
    unsigned int i = 0, j;
    
    fprintf(stderr, "Reading I.C. for /de\n");
    H5LTread_dataset_double(file_id, "/de", tics);
    for (j = 0; j < dims; j++) {
        ics[j * N + i] = tics[j] / 1.0; /* Convert to number density */
        atol[j * N + i] = tics[j] * 1e-11;
        rtol[j * N + i] = 1e-11;
        if(j==0) {
            fprintf(stderr, "de[0] = %0.3g, atol => % 0.16g\n",
                    tics[j], atol[j]);
        }
    }
    i++;
    
    fprintf(stderr, "Reading I.C. for /ge\n");
    H5LTread_dataset_double(file_id, "/ge", tics);
    for (j = 0; j < dims; j++) {
        ics[j * N + i] = tics[j] / 1.0; /* Convert to number density */
        atol[j * N + i] = tics[j] * 1e-11;
        rtol[j * N + i] = 1e-11;
        if(j==0) {
            fprintf(stderr, "ge[0] = %0.3g, atol => % 0.16g\n",
                    tics[j], atol[j]);
        }
    }
    i++;
    
    fprintf(stderr, "Reading I.C. for /OI\n");
    H5LTread_dataset_double(file_id, "/OI", tics);
    for (j = 0; j < dims; j++) {
        ics[j * N + i] = tics[j] / 16; /* Convert to number density */
        atol[j * N + i] = tics[j] * 1e-11;
        rtol[j * N + i] = 1e-11;
        if(j==0) {
            fprintf(stderr, "OI[0] = %0.3g, atol => % 0.16g\n",
                    tics[j], atol[j]);
        }
    }
    i++;
    
    fprintf(stderr, "Reading I.C. for /OII\n");
    H5LTread_dataset_double(file_id, "/OII", tics);
    for (j = 0; j < dims; j++) {
        ics[j * N + i] = tics[j] / 16; /* Convert to number density */
        atol[j * N + i] = tics[j] * 1e-11;
        rtol[j * N + i] = 1e-11;
        if(j==0) {
            fprintf(stderr, "OII[0] = %0.3g, atol => % 0.16g\n",
                    tics[j], atol[j]);
        }
    }
    i++;
    
    fprintf(stderr, "Reading I.C. for /OIII\n");
    H5LTread_dataset_double(file_id, "/OIII", tics);
    for (j = 0; j < dims; j++) {
        ics[j * N + i] = tics[j] / 16; /* Convert to number density */
        atol[j * N + i] = tics[j] * 1e-11;
        rtol[j * N + i] = 1e-11;
        if(j==0) {
            fprintf(stderr, "OIII[0] = %0.3g, atol => % 0.16g\n",
                    tics[j], atol[j]);
        }
    }
    i++;
    
    fprintf(stderr, "Reading I.C. for /OIV\n");
    H5LTread_dataset_double(file_id, "/OIV", tics);
    for (j = 0; j < dims; j++) {
        ics[j * N + i] = tics[j] / 16; /* Convert to number density */
        atol[j * N + i] = tics[j] * 1e-11;
        rtol[j * N + i] = 1e-11;
        if(j==0) {
            fprintf(stderr, "OIV[0] = %0.3g, atol => % 0.16g\n",
                    tics[j], atol[j]);
        }
    }
    i++;
    
    fprintf(stderr, "Reading I.C. for /OV\n");
    H5LTread_dataset_double(file_id, "/OV", tics);
    for (j = 0; j < dims; j++) {
        ics[j * N + i] = tics[j] / 16; /* Convert to number density */
        atol[j * N + i] = tics[j] * 1e-11;
        rtol[j * N + i] = 1e-11;
        if(j==0) {
            fprintf(stderr, "OV[0] = %0.3g, atol => % 0.16g\n",
                    tics[j], atol[j]);
        }
    }
    i++;
    
    fprintf(stderr, "Reading I.C. for /OVI\n");
    H5LTread_dataset_double(file_id, "/OVI", tics);
    for (j = 0; j < dims; j++) {
        ics[j * N + i] = tics[j] / 16; /* Convert to number density */
        atol[j * N + i] = tics[j] * 1e-11;
        rtol[j * N + i] = 1e-11;
        if(j==0) {
            fprintf(stderr, "OVI[0] = %0.3g, atol => % 0.16g\n",
                    tics[j], atol[j]);
        }
    }
    i++;
    
    fprintf(stderr, "Reading I.C. for /OVII\n");
    H5LTread_dataset_double(file_id, "/OVII", tics);
    for (j = 0; j < dims; j++) {
        ics[j * N + i] = tics[j] / 16; /* Convert to number density */
        atol[j * N + i] = tics[j] * 1e-11;
        rtol[j * N + i] = 1e-11;
        if(j==0) {
            fprintf(stderr, "OVII[0] = %0.3g, atol => % 0.16g\n",
                    tics[j], atol[j]);
        }
    }
    i++;
    
    fprintf(stderr, "Reading I.C. for /OVIII\n");
    H5LTread_dataset_double(file_id, "/OVIII", tics);
    for (j = 0; j < dims; j++) {
        ics[j * N + i] = tics[j] / 16; /* Convert to number density */
        atol[j * N + i] = tics[j] * 1e-11;
        rtol[j * N + i] = 1e-11;
        if(j==0) {
            fprintf(stderr, "OVIII[0] = %0.3g, atol => % 0.16g\n",
                    tics[j], atol[j]);
        }
    }
    i++;
    
    fprintf(stderr, "Reading I.C. for /OIX\n");
    H5LTread_dataset_double(file_id, "/OIX", tics);
    for (j = 0; j < dims; j++) {
        ics[j * N + i] = tics[j] / 16; /* Convert to number density */
        atol[j * N + i] = tics[j] * 1e-11;
        rtol[j * N + i] = 1e-11;
        if(j==0) {
            fprintf(stderr, "OIX[0] = %0.3g, atol => % 0.16g\n",
                    tics[j], atol[j]);
        }
    }
    i++;
    

    H5Fclose(file_id);

    double dtf = 3.1557e16;
    double dt = -1.0;
    for (i = 0; i < dims * N; i++) input[i] = ics[i];
    double ttot;
    ttot = dengo_evolve_oxygen(dtf, dt, input, rtol, atol, dims, data);

    /* Write results to HDF5 file */
    file_id = H5Fcreate("oxygen_solution.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    hsize_t dimsarr[1];
    dimsarr[0] = dims;
    i = 0;
    
    double de[dims];
    for (j = 0; j < dims; j++) {
        de[j] = input[j * N + i] * 1.0;
    }
    fprintf(stderr, "Writing solution for /de\n");
    H5LTmake_dataset_double(file_id, "/de", 1, dimsarr, de);
    i++;
    
    double ge[dims];
    for (j = 0; j < dims; j++) {
        ge[j] = input[j * N + i] * 1.0;
    }
    fprintf(stderr, "Writing solution for /ge\n");
    H5LTmake_dataset_double(file_id, "/ge", 1, dimsarr, ge);
    i++;
    
    double OI[dims];
    for (j = 0; j < dims; j++) {
        OI[j] = input[j * N + i] * 16;
    }
    fprintf(stderr, "Writing solution for /OI\n");
    H5LTmake_dataset_double(file_id, "/OI", 1, dimsarr, OI);
    i++;
    
    double OII[dims];
    for (j = 0; j < dims; j++) {
        OII[j] = input[j * N + i] * 16;
    }
    fprintf(stderr, "Writing solution for /OII\n");
    H5LTmake_dataset_double(file_id, "/OII", 1, dimsarr, OII);
    i++;
    
    double OIII[dims];
    for (j = 0; j < dims; j++) {
        OIII[j] = input[j * N + i] * 16;
    }
    fprintf(stderr, "Writing solution for /OIII\n");
    H5LTmake_dataset_double(file_id, "/OIII", 1, dimsarr, OIII);
    i++;
    
    double OIV[dims];
    for (j = 0; j < dims; j++) {
        OIV[j] = input[j * N + i] * 16;
    }
    fprintf(stderr, "Writing solution for /OIV\n");
    H5LTmake_dataset_double(file_id, "/OIV", 1, dimsarr, OIV);
    i++;
    
    double OV[dims];
    for (j = 0; j < dims; j++) {
        OV[j] = input[j * N + i] * 16;
    }
    fprintf(stderr, "Writing solution for /OV\n");
    H5LTmake_dataset_double(file_id, "/OV", 1, dimsarr, OV);
    i++;
    
    double OVI[dims];
    for (j = 0; j < dims; j++) {
        OVI[j] = input[j * N + i] * 16;
    }
    fprintf(stderr, "Writing solution for /OVI\n");
    H5LTmake_dataset_double(file_id, "/OVI", 1, dimsarr, OVI);
    i++;
    
    double OVII[dims];
    for (j = 0; j < dims; j++) {
        OVII[j] = input[j * N + i] * 16;
    }
    fprintf(stderr, "Writing solution for /OVII\n");
    H5LTmake_dataset_double(file_id, "/OVII", 1, dimsarr, OVII);
    i++;
    
    double OVIII[dims];
    for (j = 0; j < dims; j++) {
        OVIII[j] = input[j * N + i] * 16;
    }
    fprintf(stderr, "Writing solution for /OVIII\n");
    H5LTmake_dataset_double(file_id, "/OVIII", 1, dimsarr, OVIII);
    i++;
    
    double OIX[dims];
    for (j = 0; j < dims; j++) {
        OIX[j] = input[j * N + i] * 16;
    }
    fprintf(stderr, "Writing solution for /OIX\n");
    H5LTmake_dataset_double(file_id, "/OIX", 1, dimsarr, OIX);
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
 



double dengo_evolve_oxygen (double dtf, double &dt, double *input,
            double *rtol, double *atol, int dims, oxygen_data *data) {
    int i, j;
    hid_t file_id;
    fprintf(stderr, "  ncells = % 3i\n", (int) dims);

    int N = 11;
    rhs_f f = calculate_rhs_oxygen;
    jac_f jf = calculate_jacobian_oxygen;
    if (dt < 0) dt = dtf / 10000.0;
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
	    if (siter == 9999) break;
	    siter++;
        if (siter % 1000 == 0) {
            fprintf(stderr, "Successful Iteration[%d]: (%0.4g) %0.16g / %0.16g\n",
                     siter, dt, ttot, dtf);
        }
        ttot += dt;
	    dt = DMIN(dt * 1.1, dtf - ttot);
	    
	    /* Write intermediate  results to HDF5 file */
	    char imfilename[255];
	    snprintf(imfilename, 255, "oxygen_intermediate_%06d.h5", siter);
    	    file_id = H5Fcreate(imfilename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    	    hsize_t dimsarr[1];
    	    dimsarr[0] = dims;
    	    i = 0;
    	    
    	    double de[dims];
    	    for (j = 0; j < dims; j++) {
            	de[j] = prev[j * N + i];
    	    }
    	    //fprintf(stderr, "Writing solution for /de\n");
    	    H5LTmake_dataset_double(file_id, "/de", 1, dimsarr, de);
    	    i++;
    	    
    	    double ge[dims];
    	    for (j = 0; j < dims; j++) {
            	ge[j] = prev[j * N + i];
    	    }
    	    //fprintf(stderr, "Writing solution for /ge\n");
    	    H5LTmake_dataset_double(file_id, "/ge", 1, dimsarr, ge);
    	    i++;
    	    
    	    double OI[dims];
    	    for (j = 0; j < dims; j++) {
            	OI[j] = prev[j * N + i];
    	    }
    	    //fprintf(stderr, "Writing solution for /OI\n");
    	    H5LTmake_dataset_double(file_id, "/OI", 1, dimsarr, OI);
    	    i++;
    	    
    	    double OII[dims];
    	    for (j = 0; j < dims; j++) {
            	OII[j] = prev[j * N + i];
    	    }
    	    //fprintf(stderr, "Writing solution for /OII\n");
    	    H5LTmake_dataset_double(file_id, "/OII", 1, dimsarr, OII);
    	    i++;
    	    
    	    double OIII[dims];
    	    for (j = 0; j < dims; j++) {
            	OIII[j] = prev[j * N + i];
    	    }
    	    //fprintf(stderr, "Writing solution for /OIII\n");
    	    H5LTmake_dataset_double(file_id, "/OIII", 1, dimsarr, OIII);
    	    i++;
    	    
    	    double OIV[dims];
    	    for (j = 0; j < dims; j++) {
            	OIV[j] = prev[j * N + i];
    	    }
    	    //fprintf(stderr, "Writing solution for /OIV\n");
    	    H5LTmake_dataset_double(file_id, "/OIV", 1, dimsarr, OIV);
    	    i++;
    	    
    	    double OV[dims];
    	    for (j = 0; j < dims; j++) {
            	OV[j] = prev[j * N + i];
    	    }
    	    //fprintf(stderr, "Writing solution for /OV\n");
    	    H5LTmake_dataset_double(file_id, "/OV", 1, dimsarr, OV);
    	    i++;
    	    
    	    double OVI[dims];
    	    for (j = 0; j < dims; j++) {
            	OVI[j] = prev[j * N + i];
    	    }
    	    //fprintf(stderr, "Writing solution for /OVI\n");
    	    H5LTmake_dataset_double(file_id, "/OVI", 1, dimsarr, OVI);
    	    i++;
    	    
    	    double OVII[dims];
    	    for (j = 0; j < dims; j++) {
            	OVII[j] = prev[j * N + i];
    	    }
    	    //fprintf(stderr, "Writing solution for /OVII\n");
    	    H5LTmake_dataset_double(file_id, "/OVII", 1, dimsarr, OVII);
    	    i++;
    	    
    	    double OVIII[dims];
    	    for (j = 0; j < dims; j++) {
            	OVIII[j] = prev[j * N + i];
    	    }
    	    //fprintf(stderr, "Writing solution for /OVIII\n");
    	    H5LTmake_dataset_double(file_id, "/OVIII", 1, dimsarr, OVIII);
    	    i++;
    	    
    	    double OIX[dims];
    	    for (j = 0; j < dims; j++) {
            	OIX[j] = prev[j * N + i];
    	    }
    	    //fprintf(stderr, "Writing solution for /OIX\n");
    	    H5LTmake_dataset_double(file_id, "/OIX", 1, dimsarr, OIX);
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
 


void oxygen_read_rate_tables(oxygen_data *data)
{
    hid_t file_id = H5Fopen("oxygen_tables.h5", H5F_ACC_RDONLY, H5P_DEFAULT);
    /* Allocate the correct number of rate tables */
    H5LTread_dataset_double(file_id, "/OI_i", data->r_OI_i);
    H5LTread_dataset_double(file_id, "/OII_i", data->r_OII_i);
    H5LTread_dataset_double(file_id, "/OII_r", data->r_OII_r);
    H5LTread_dataset_double(file_id, "/OIII_i", data->r_OIII_i);
    H5LTread_dataset_double(file_id, "/OIII_r", data->r_OIII_r);
    H5LTread_dataset_double(file_id, "/OIV_i", data->r_OIV_i);
    H5LTread_dataset_double(file_id, "/OIV_r", data->r_OIV_r);
    H5LTread_dataset_double(file_id, "/OIX_r", data->r_OIX_r);
    H5LTread_dataset_double(file_id, "/OV_i", data->r_OV_i);
    H5LTread_dataset_double(file_id, "/OV_r", data->r_OV_r);
    H5LTread_dataset_double(file_id, "/OVI_i", data->r_OVI_i);
    H5LTread_dataset_double(file_id, "/OVI_r", data->r_OVI_r);
    H5LTread_dataset_double(file_id, "/OVII_i", data->r_OVII_i);
    H5LTread_dataset_double(file_id, "/OVII_r", data->r_OVII_r);
    H5LTread_dataset_double(file_id, "/OVIII_i", data->r_OVIII_i);
    H5LTread_dataset_double(file_id, "/OVIII_r", data->r_OVIII_r);

    H5Fclose(file_id);
}

void oxygen_read_cooling_tables(oxygen_data *data)
{

    hid_t file_id = H5Fopen("oxygen_tables.h5", H5F_ACC_RDONLY, H5P_DEFAULT);
    /* Allocate the correct number of rate tables */
    H5LTread_dataset_double(file_id, "/OI_c_OI_c",
                            data->c_OI_c_OI_c);
    H5LTread_dataset_double(file_id, "/OII_c_OII_c",
                            data->c_OII_c_OII_c);
    H5LTread_dataset_double(file_id, "/OIII_c_OIII_c",
                            data->c_OIII_c_OIII_c);
    H5LTread_dataset_double(file_id, "/OIV_c_OIV_c",
                            data->c_OIV_c_OIV_c);
    H5LTread_dataset_double(file_id, "/OIX_c_OIX_c",
                            data->c_OIX_c_OIX_c);
    H5LTread_dataset_double(file_id, "/OV_c_OV_c",
                            data->c_OV_c_OV_c);
    H5LTread_dataset_double(file_id, "/OVI_c_OVI_c",
                            data->c_OVI_c_OVI_c);
    H5LTread_dataset_double(file_id, "/OVII_c_OVII_c",
                            data->c_OVII_c_OVII_c);
    H5LTread_dataset_double(file_id, "/OVIII_c_OVIII_c",
                            data->c_OVIII_c_OVIII_c);

    H5Fclose(file_id);
}

 


void oxygen_calculate_temperature(oxygen_data *data,
                        double *input, int nstrip, int nchem)
{
    int i, j;
    double density;
    double kb = 1.3806504e-16; // Boltzmann constant [erg/K]
    double mh = 1.67e-24;
    double gamma = 5.e0/3.e0;
    
    /* Calculate total density */
    double de;
    double ge;
    double OI;
    double OII;
    double OIII;
    double OIV;
    double OV;
    double OVI;
    double OVII;
    double OVIII;
    double OIX;

    for (i = 0; i<nstrip; i++) {
        j = i * nchem;
        de = input[j];
        /*fprintf(stderr, "de[%d] = % 0.16g\n",
                i, de);*/
        j++;
    
        ge = input[j];
        /*fprintf(stderr, "ge[%d] = % 0.16g\n",
                i, ge);*/
        j++;
    
        OI = input[j];
        /*fprintf(stderr, "OI[%d] = % 0.16g\n",
                i, OI);*/
        j++;
    
        OII = input[j];
        /*fprintf(stderr, "OII[%d] = % 0.16g\n",
                i, OII);*/
        j++;
    
        OIII = input[j];
        /*fprintf(stderr, "OIII[%d] = % 0.16g\n",
                i, OIII);*/
        j++;
    
        OIV = input[j];
        /*fprintf(stderr, "OIV[%d] = % 0.16g\n",
                i, OIV);*/
        j++;
    
        OV = input[j];
        /*fprintf(stderr, "OV[%d] = % 0.16g\n",
                i, OV);*/
        j++;
    
        OVI = input[j];
        /*fprintf(stderr, "OVI[%d] = % 0.16g\n",
                i, OVI);*/
        j++;
    
        OVII = input[j];
        /*fprintf(stderr, "OVII[%d] = % 0.16g\n",
                i, OVII);*/
        j++;
    
        OVIII = input[j];
        /*fprintf(stderr, "OVIII[%d] = % 0.16g\n",
                i, OVIII);*/
        j++;
    
        OIX = input[j];
        /*fprintf(stderr, "OIX[%d] = % 0.16g\n",
                i, OIX);*/
        j++;
    
        density = 16*OI + 16*OII + 16*OIII + 16*OIV + 16*OIX + 16*OV + 16*OVI + 16*OVII + 16*OVIII;
        data->Ts[i] = density*ge*mh/(kb*(OI/(gamma - 1.0) + OII/(gamma - 1.0) + OIII/(gamma - 1.0) + OIV/(gamma - 1.0) + OIX/(gamma - 1.0) + OV/(gamma - 1.0) + OVI/(gamma - 1.0) + OVII/(gamma - 1.0) + OVIII/(gamma - 1.0) + de/(gamma - 1.0)));
        if (data->Ts[i] < data->bounds[0]) {
            data->Ts[i] = data->bounds[0];
        } else if (data->Ts[i] > data->bounds[1]) {
            data->Ts[i] = data->bounds[1];
        }
        data->logTs[i] = log(data->Ts[i]);
	data->dTs_ge[i] = 
        density*mh/(kb*(OI/(gamma - 1.0) + OII/(gamma - 1.0) + OIII/(gamma - 1.0) + OIV/(gamma - 1.0) + OIX/(gamma - 1.0) + OV/(gamma - 1.0) + OVI/(gamma - 1.0) + OVII/(gamma - 1.0) + OVIII/(gamma - 1.0) + de/(gamma - 1.0)));
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
void oxygen_interpolate_rates(oxygen_data *data,
                    int nstrip)
{
    int i, bin_id;
    double lb, t1, t2;
    lb = log(data->bounds[0]);
    /*fprintf(stderr, "lb = % 0.16g, ub = % 0.16g\n", lb, ub);*/
    for (i = 0; i < nstrip; i++) {
        data->bin_id[i] = bin_id = (int) (data->idbin * (data->logTs[i] - lb));
        if (data->bin_id[i] <= 0) {
            data->bin_id[i] = 1;
        } else if (data->bin_id[i] >= data->nbins) {
            data->bin_id[i] = data->nbins - 1;
        }
        t1 = (lb + (bin_id - 1) * data->dbin);
        t2 = (lb + (bin_id    ) * data->dbin);
        data->Tdef[i] = (data->logTs[i] - t1)/(t2 - t1);
        data->dT[i] = (t2 - t1);
        /*fprintf(stderr, "INTERP: %d, bin_id = %d, dT = % 0.16g, T = % 0.16g, logT = % 0.16g\n",
                i, data->bin_id[i], data->dT[i], data->Ts[i],
                data->logTs[i]);*/
    }
    for (i = 0; i < nstrip; i++) {
        bin_id = data->bin_id[i];
        data->rs_OI_i[i] = data->r_OI_i[bin_id] +
            data->Tdef[i] * (data->r_OI_i[bin_id+1] - data->r_OI_i[bin_id]);
        data->drs_OI_i[i] = (data->r_OI_i[bin_id+1] - data->r_OI_i[bin_id]);
        data->drs_OI_i[i] /= data->dT[i];
	data->drs_OI_i[i] /= data->Ts[i];
    }
    
    for (i = 0; i < nstrip; i++) {
        bin_id = data->bin_id[i];
        data->rs_OII_i[i] = data->r_OII_i[bin_id] +
            data->Tdef[i] * (data->r_OII_i[bin_id+1] - data->r_OII_i[bin_id]);
        data->drs_OII_i[i] = (data->r_OII_i[bin_id+1] - data->r_OII_i[bin_id]);
        data->drs_OII_i[i] /= data->dT[i];
	data->drs_OII_i[i] /= data->Ts[i];
    }
    
    for (i = 0; i < nstrip; i++) {
        bin_id = data->bin_id[i];
        data->rs_OII_r[i] = data->r_OII_r[bin_id] +
            data->Tdef[i] * (data->r_OII_r[bin_id+1] - data->r_OII_r[bin_id]);
        data->drs_OII_r[i] = (data->r_OII_r[bin_id+1] - data->r_OII_r[bin_id]);
        data->drs_OII_r[i] /= data->dT[i];
	data->drs_OII_r[i] /= data->Ts[i];
    }
    
    for (i = 0; i < nstrip; i++) {
        bin_id = data->bin_id[i];
        data->rs_OIII_i[i] = data->r_OIII_i[bin_id] +
            data->Tdef[i] * (data->r_OIII_i[bin_id+1] - data->r_OIII_i[bin_id]);
        data->drs_OIII_i[i] = (data->r_OIII_i[bin_id+1] - data->r_OIII_i[bin_id]);
        data->drs_OIII_i[i] /= data->dT[i];
	data->drs_OIII_i[i] /= data->Ts[i];
    }
    
    for (i = 0; i < nstrip; i++) {
        bin_id = data->bin_id[i];
        data->rs_OIII_r[i] = data->r_OIII_r[bin_id] +
            data->Tdef[i] * (data->r_OIII_r[bin_id+1] - data->r_OIII_r[bin_id]);
        data->drs_OIII_r[i] = (data->r_OIII_r[bin_id+1] - data->r_OIII_r[bin_id]);
        data->drs_OIII_r[i] /= data->dT[i];
	data->drs_OIII_r[i] /= data->Ts[i];
    }
    
    for (i = 0; i < nstrip; i++) {
        bin_id = data->bin_id[i];
        data->rs_OIV_i[i] = data->r_OIV_i[bin_id] +
            data->Tdef[i] * (data->r_OIV_i[bin_id+1] - data->r_OIV_i[bin_id]);
        data->drs_OIV_i[i] = (data->r_OIV_i[bin_id+1] - data->r_OIV_i[bin_id]);
        data->drs_OIV_i[i] /= data->dT[i];
	data->drs_OIV_i[i] /= data->Ts[i];
    }
    
    for (i = 0; i < nstrip; i++) {
        bin_id = data->bin_id[i];
        data->rs_OIV_r[i] = data->r_OIV_r[bin_id] +
            data->Tdef[i] * (data->r_OIV_r[bin_id+1] - data->r_OIV_r[bin_id]);
        data->drs_OIV_r[i] = (data->r_OIV_r[bin_id+1] - data->r_OIV_r[bin_id]);
        data->drs_OIV_r[i] /= data->dT[i];
	data->drs_OIV_r[i] /= data->Ts[i];
    }
    
    for (i = 0; i < nstrip; i++) {
        bin_id = data->bin_id[i];
        data->rs_OIX_r[i] = data->r_OIX_r[bin_id] +
            data->Tdef[i] * (data->r_OIX_r[bin_id+1] - data->r_OIX_r[bin_id]);
        data->drs_OIX_r[i] = (data->r_OIX_r[bin_id+1] - data->r_OIX_r[bin_id]);
        data->drs_OIX_r[i] /= data->dT[i];
	data->drs_OIX_r[i] /= data->Ts[i];
    }
    
    for (i = 0; i < nstrip; i++) {
        bin_id = data->bin_id[i];
        data->rs_OV_i[i] = data->r_OV_i[bin_id] +
            data->Tdef[i] * (data->r_OV_i[bin_id+1] - data->r_OV_i[bin_id]);
        data->drs_OV_i[i] = (data->r_OV_i[bin_id+1] - data->r_OV_i[bin_id]);
        data->drs_OV_i[i] /= data->dT[i];
	data->drs_OV_i[i] /= data->Ts[i];
    }
    
    for (i = 0; i < nstrip; i++) {
        bin_id = data->bin_id[i];
        data->rs_OV_r[i] = data->r_OV_r[bin_id] +
            data->Tdef[i] * (data->r_OV_r[bin_id+1] - data->r_OV_r[bin_id]);
        data->drs_OV_r[i] = (data->r_OV_r[bin_id+1] - data->r_OV_r[bin_id]);
        data->drs_OV_r[i] /= data->dT[i];
	data->drs_OV_r[i] /= data->Ts[i];
    }
    
    for (i = 0; i < nstrip; i++) {
        bin_id = data->bin_id[i];
        data->rs_OVI_i[i] = data->r_OVI_i[bin_id] +
            data->Tdef[i] * (data->r_OVI_i[bin_id+1] - data->r_OVI_i[bin_id]);
        data->drs_OVI_i[i] = (data->r_OVI_i[bin_id+1] - data->r_OVI_i[bin_id]);
        data->drs_OVI_i[i] /= data->dT[i];
	data->drs_OVI_i[i] /= data->Ts[i];
    }
    
    for (i = 0; i < nstrip; i++) {
        bin_id = data->bin_id[i];
        data->rs_OVI_r[i] = data->r_OVI_r[bin_id] +
            data->Tdef[i] * (data->r_OVI_r[bin_id+1] - data->r_OVI_r[bin_id]);
        data->drs_OVI_r[i] = (data->r_OVI_r[bin_id+1] - data->r_OVI_r[bin_id]);
        data->drs_OVI_r[i] /= data->dT[i];
	data->drs_OVI_r[i] /= data->Ts[i];
    }
    
    for (i = 0; i < nstrip; i++) {
        bin_id = data->bin_id[i];
        data->rs_OVII_i[i] = data->r_OVII_i[bin_id] +
            data->Tdef[i] * (data->r_OVII_i[bin_id+1] - data->r_OVII_i[bin_id]);
        data->drs_OVII_i[i] = (data->r_OVII_i[bin_id+1] - data->r_OVII_i[bin_id]);
        data->drs_OVII_i[i] /= data->dT[i];
	data->drs_OVII_i[i] /= data->Ts[i];
    }
    
    for (i = 0; i < nstrip; i++) {
        bin_id = data->bin_id[i];
        data->rs_OVII_r[i] = data->r_OVII_r[bin_id] +
            data->Tdef[i] * (data->r_OVII_r[bin_id+1] - data->r_OVII_r[bin_id]);
        data->drs_OVII_r[i] = (data->r_OVII_r[bin_id+1] - data->r_OVII_r[bin_id]);
        data->drs_OVII_r[i] /= data->dT[i];
	data->drs_OVII_r[i] /= data->Ts[i];
    }
    
    for (i = 0; i < nstrip; i++) {
        bin_id = data->bin_id[i];
        data->rs_OVIII_i[i] = data->r_OVIII_i[bin_id] +
            data->Tdef[i] * (data->r_OVIII_i[bin_id+1] - data->r_OVIII_i[bin_id]);
        data->drs_OVIII_i[i] = (data->r_OVIII_i[bin_id+1] - data->r_OVIII_i[bin_id]);
        data->drs_OVIII_i[i] /= data->dT[i];
	data->drs_OVIII_i[i] /= data->Ts[i];
    }
    
    for (i = 0; i < nstrip; i++) {
        bin_id = data->bin_id[i];
        data->rs_OVIII_r[i] = data->r_OVIII_r[bin_id] +
            data->Tdef[i] * (data->r_OVIII_r[bin_id+1] - data->r_OVIII_r[bin_id]);
        data->drs_OVIII_r[i] = (data->r_OVIII_r[bin_id+1] - data->r_OVIII_r[bin_id]);
        data->drs_OVIII_r[i] /= data->dT[i];
	data->drs_OVIII_r[i] /= data->Ts[i];
    }
    
    for (i = 0; i < nstrip; i++) {
        bin_id = data->bin_id[i];
        data->cs_OI_c_OI_c[i] = data->c_OI_c_OI_c[bin_id] +
            data->Tdef[i] * (data->c_OI_c_OI_c[bin_id+1] - data->c_OI_c_OI_c[bin_id]);
        data->dcs_OI_c_OI_c[i] = (data->c_OI_c_OI_c[bin_id+1] - data->c_OI_c_OI_c[bin_id]);;
        data->dcs_OI_c_OI_c[i] /= data->dT[i];
	data->dcs_OI_c_OI_c[i] /= data->Ts[i];
    }
    
    
    for (i = 0; i < nstrip; i++) {
        bin_id = data->bin_id[i];
        data->cs_OII_c_OII_c[i] = data->c_OII_c_OII_c[bin_id] +
            data->Tdef[i] * (data->c_OII_c_OII_c[bin_id+1] - data->c_OII_c_OII_c[bin_id]);
        data->dcs_OII_c_OII_c[i] = (data->c_OII_c_OII_c[bin_id+1] - data->c_OII_c_OII_c[bin_id]);;
        data->dcs_OII_c_OII_c[i] /= data->dT[i];
	data->dcs_OII_c_OII_c[i] /= data->Ts[i];
    }
    
    
    for (i = 0; i < nstrip; i++) {
        bin_id = data->bin_id[i];
        data->cs_OIII_c_OIII_c[i] = data->c_OIII_c_OIII_c[bin_id] +
            data->Tdef[i] * (data->c_OIII_c_OIII_c[bin_id+1] - data->c_OIII_c_OIII_c[bin_id]);
        data->dcs_OIII_c_OIII_c[i] = (data->c_OIII_c_OIII_c[bin_id+1] - data->c_OIII_c_OIII_c[bin_id]);;
        data->dcs_OIII_c_OIII_c[i] /= data->dT[i];
	data->dcs_OIII_c_OIII_c[i] /= data->Ts[i];
    }
    
    
    for (i = 0; i < nstrip; i++) {
        bin_id = data->bin_id[i];
        data->cs_OIV_c_OIV_c[i] = data->c_OIV_c_OIV_c[bin_id] +
            data->Tdef[i] * (data->c_OIV_c_OIV_c[bin_id+1] - data->c_OIV_c_OIV_c[bin_id]);
        data->dcs_OIV_c_OIV_c[i] = (data->c_OIV_c_OIV_c[bin_id+1] - data->c_OIV_c_OIV_c[bin_id]);;
        data->dcs_OIV_c_OIV_c[i] /= data->dT[i];
	data->dcs_OIV_c_OIV_c[i] /= data->Ts[i];
    }
    
    
    for (i = 0; i < nstrip; i++) {
        bin_id = data->bin_id[i];
        data->cs_OIX_c_OIX_c[i] = data->c_OIX_c_OIX_c[bin_id] +
            data->Tdef[i] * (data->c_OIX_c_OIX_c[bin_id+1] - data->c_OIX_c_OIX_c[bin_id]);
        data->dcs_OIX_c_OIX_c[i] = (data->c_OIX_c_OIX_c[bin_id+1] - data->c_OIX_c_OIX_c[bin_id]);;
        data->dcs_OIX_c_OIX_c[i] /= data->dT[i];
	data->dcs_OIX_c_OIX_c[i] /= data->Ts[i];
    }
    
    
    for (i = 0; i < nstrip; i++) {
        bin_id = data->bin_id[i];
        data->cs_OV_c_OV_c[i] = data->c_OV_c_OV_c[bin_id] +
            data->Tdef[i] * (data->c_OV_c_OV_c[bin_id+1] - data->c_OV_c_OV_c[bin_id]);
        data->dcs_OV_c_OV_c[i] = (data->c_OV_c_OV_c[bin_id+1] - data->c_OV_c_OV_c[bin_id]);;
        data->dcs_OV_c_OV_c[i] /= data->dT[i];
	data->dcs_OV_c_OV_c[i] /= data->Ts[i];
    }
    
    
    for (i = 0; i < nstrip; i++) {
        bin_id = data->bin_id[i];
        data->cs_OVI_c_OVI_c[i] = data->c_OVI_c_OVI_c[bin_id] +
            data->Tdef[i] * (data->c_OVI_c_OVI_c[bin_id+1] - data->c_OVI_c_OVI_c[bin_id]);
        data->dcs_OVI_c_OVI_c[i] = (data->c_OVI_c_OVI_c[bin_id+1] - data->c_OVI_c_OVI_c[bin_id]);;
        data->dcs_OVI_c_OVI_c[i] /= data->dT[i];
	data->dcs_OVI_c_OVI_c[i] /= data->Ts[i];
    }
    
    
    for (i = 0; i < nstrip; i++) {
        bin_id = data->bin_id[i];
        data->cs_OVII_c_OVII_c[i] = data->c_OVII_c_OVII_c[bin_id] +
            data->Tdef[i] * (data->c_OVII_c_OVII_c[bin_id+1] - data->c_OVII_c_OVII_c[bin_id]);
        data->dcs_OVII_c_OVII_c[i] = (data->c_OVII_c_OVII_c[bin_id+1] - data->c_OVII_c_OVII_c[bin_id]);;
        data->dcs_OVII_c_OVII_c[i] /= data->dT[i];
	data->dcs_OVII_c_OVII_c[i] /= data->Ts[i];
    }
    
    
    for (i = 0; i < nstrip; i++) {
        bin_id = data->bin_id[i];
        data->cs_OVIII_c_OVIII_c[i] = data->c_OVIII_c_OVIII_c[bin_id] +
            data->Tdef[i] * (data->c_OVIII_c_OVIII_c[bin_id+1] - data->c_OVIII_c_OVIII_c[bin_id]);
        data->dcs_OVIII_c_OVIII_c[i] = (data->c_OVIII_c_OVIII_c[bin_id+1] - data->c_OVIII_c_OVIII_c[bin_id]);;
        data->dcs_OVIII_c_OVIII_c[i] /= data->dT[i];
	data->dcs_OVIII_c_OVIII_c[i] /= data->Ts[i];
    }
    
    

}
 



int calculate_rhs_oxygen(double *input, double *rhs, int nstrip,
                  int nchem, void *sdata)
{
    /* We iterate over all of the rates */
    /* Calculate temperature first */
    oxygen_data *data = (oxygen_data*)sdata;
    int i, j;
    oxygen_calculate_temperature(data, input, nstrip, nchem);

    oxygen_interpolate_rates(data, nstrip);

    /* Now we set up some temporaries */
    double *OI_i = data->rs_OI_i;
    double *OII_i = data->rs_OII_i;
    double *OII_r = data->rs_OII_r;
    double *OIII_i = data->rs_OIII_i;
    double *OIII_r = data->rs_OIII_r;
    double *OIV_i = data->rs_OIV_i;
    double *OIV_r = data->rs_OIV_r;
    double *OIX_r = data->rs_OIX_r;
    double *OV_i = data->rs_OV_i;
    double *OV_r = data->rs_OV_r;
    double *OVI_i = data->rs_OVI_i;
    double *OVI_r = data->rs_OVI_r;
    double *OVII_i = data->rs_OVII_i;
    double *OVII_r = data->rs_OVII_r;
    double *OVIII_i = data->rs_OVIII_i;
    double *OVIII_r = data->rs_OVIII_r;
    double *OI_c_OI_c = data->cs_OI_c_OI_c;
    double *OII_c_OII_c = data->cs_OII_c_OII_c;
    double *OIII_c_OIII_c = data->cs_OIII_c_OIII_c;
    double *OIV_c_OIV_c = data->cs_OIV_c_OIV_c;
    double *OIX_c_OIX_c = data->cs_OIX_c_OIX_c;
    double *OV_c_OV_c = data->cs_OV_c_OV_c;
    double *OVI_c_OVI_c = data->cs_OVI_c_OVI_c;
    double *OVII_c_OVII_c = data->cs_OVII_c_OVII_c;
    double *OVIII_c_OVIII_c = data->cs_OVIII_c_OVIII_c;
    double de;
    double ge;
    double OI;
    double OII;
    double OIII;
    double OIV;
    double OV;
    double OVI;
    double OVII;
    double OVIII;
    double OIX;

    double mh = 1.67e-24;
    double total, total_e, total_de, mdensity;
    for (i = 0; i<nstrip; i++) {
        j = i * nchem;
        total = total_e = total_de = mdensity = 0.0;
        de = input[j];
        if (de < 0.0) {
          fprintf(stderr, "RNegative[%d][de] = % 0.16g [%d]\n",
            i, de, j);
            return 1;
          de = 1e-20;
        }
        
        
        
        j++;
    
        ge = input[j];
        if (ge < 0.0) {
          fprintf(stderr, "RNegative[%d][ge] = % 0.16g [%d]\n",
            i, ge, j);
            return 1;
          ge = 1e-20;
        }
        
        j++;
    
        OI = input[j];
        if (OI < 0.0) {
          fprintf(stderr, "RNegative[%d][OI] = % 0.16g [%d]\n",
            i, OI, j);
            return 1;
          OI = 1e-20;
        }
        
        
          total+=OI * 16;
        
        
        j++;
    
        OII = input[j];
        if (OII < 0.0) {
          fprintf(stderr, "RNegative[%d][OII] = % 0.16g [%d]\n",
            i, OII, j);
            return 1;
          OII = 1e-20;
        }
        
        
          total+=OII * 16;
        
        
        j++;
    
        OIII = input[j];
        if (OIII < 0.0) {
          fprintf(stderr, "RNegative[%d][OIII] = % 0.16g [%d]\n",
            i, OIII, j);
            return 1;
          OIII = 1e-20;
        }
        
        
          total+=OIII * 16;
        
        
        j++;
    
        OIV = input[j];
        if (OIV < 0.0) {
          fprintf(stderr, "RNegative[%d][OIV] = % 0.16g [%d]\n",
            i, OIV, j);
            return 1;
          OIV = 1e-20;
        }
        
        
          total+=OIV * 16;
        
        
        j++;
    
        OV = input[j];
        if (OV < 0.0) {
          fprintf(stderr, "RNegative[%d][OV] = % 0.16g [%d]\n",
            i, OV, j);
            return 1;
          OV = 1e-20;
        }
        
        
          total+=OV * 16;
        
        
        j++;
    
        OVI = input[j];
        if (OVI < 0.0) {
          fprintf(stderr, "RNegative[%d][OVI] = % 0.16g [%d]\n",
            i, OVI, j);
            return 1;
          OVI = 1e-20;
        }
        
        
          total+=OVI * 16;
        
        
        j++;
    
        OVII = input[j];
        if (OVII < 0.0) {
          fprintf(stderr, "RNegative[%d][OVII] = % 0.16g [%d]\n",
            i, OVII, j);
            return 1;
          OVII = 1e-20;
        }
        
        
          total+=OVII * 16;
        
        
        j++;
    
        OVIII = input[j];
        if (OVIII < 0.0) {
          fprintf(stderr, "RNegative[%d][OVIII] = % 0.16g [%d]\n",
            i, OVIII, j);
            return 1;
          OVIII = 1e-20;
        }
        
        
          total+=OVIII * 16;
        
        
        j++;
    
        OIX = input[j];
        if (OIX < 0.0) {
          fprintf(stderr, "RNegative[%d][OIX] = % 0.16g [%d]\n",
            i, OIX, j);
            return 1;
          OIX = 1e-20;
        }
        
        
          total+=OIX * 16;
        
        
        j++;
    
        mdensity = total * mh;
        total = 0.0;
        j = i * nchem;
        // 
        // Species: de
        // 
        rhs[j] = OIII_i[i]*OIII*de - OIII_r[i]*OIII*de + OII_i[i]*OII*de - OII_r[i]*OII*de + OIV_i[i]*OIV*de - OIV_r[i]*OIV*de - OIX_r[i]*OIX*de + OI_i[i]*OI*de + OVIII_i[i]*OVIII*de - OVIII_r[i]*OVIII*de + OVII_i[i]*OVII*de - OVII_r[i]*OVII*de + OVI_i[i]*OVI*de - OVI_r[i]*OVI*de + OV_i[i]*OV*de - OV_r[i]*OV*de;
        
        
            total_de += -rhs[j];
        
        j++;
    
        // 
        // Species: ge
        // 
        rhs[j] = -OI*OI_c_OI_c[i]*de - OII*OII_c_OII_c[i]*de - OIII*OIII_c_OIII_c[i]*de - OIV*OIV_c_OIV_c[i]*de - OIX*OIX_c_OIX_c[i]*de - OV*OV_c_OV_c[i]*de - OVI*OVI_c_OVI_c[i]*de - OVII*OVII_c_OVII_c[i]*de - OVIII*OVIII_c_OVIII_c[i]*de;
        
	    rhs[j] /= mdensity;
        
        
        j++;
    
        // 
        // Species: OI
        // 
        rhs[j] = OII_r[i]*OII*de - OI_i[i]*OI*de;
        
            /* Already in number density, not mass density */
            total += rhs[j] * 16;
            total_e += OI * 0;
        
        
            total_de += rhs[j] * 0;
        
        j++;
    
        // 
        // Species: OII
        // 
        rhs[j] = OIII_r[i]*OIII*de - OII_i[i]*OII*de - OII_r[i]*OII*de + OI_i[i]*OI*de;
        
            /* Already in number density, not mass density */
            total += rhs[j] * 16;
            total_e += OII * 1;
        
        
            total_de += rhs[j] * 1;
        
        j++;
    
        // 
        // Species: OIII
        // 
        rhs[j] = -OIII_i[i]*OIII*de - OIII_r[i]*OIII*de + OII_i[i]*OII*de + OIV_r[i]*OIV*de;
        
            /* Already in number density, not mass density */
            total += rhs[j] * 16;
            total_e += OIII * 2;
        
        
            total_de += rhs[j] * 2;
        
        j++;
    
        // 
        // Species: OIV
        // 
        rhs[j] = OIII_i[i]*OIII*de - OIV_i[i]*OIV*de - OIV_r[i]*OIV*de + OV_r[i]*OV*de;
        
            /* Already in number density, not mass density */
            total += rhs[j] * 16;
            total_e += OIV * 3;
        
        
            total_de += rhs[j] * 3;
        
        j++;
    
        // 
        // Species: OV
        // 
        rhs[j] = OIV_i[i]*OIV*de + OVI_r[i]*OVI*de - OV_i[i]*OV*de - OV_r[i]*OV*de;
        
            /* Already in number density, not mass density */
            total += rhs[j] * 16;
            total_e += OV * 4;
        
        
            total_de += rhs[j] * 4;
        
        j++;
    
        // 
        // Species: OVI
        // 
        rhs[j] = OVII_r[i]*OVII*de - OVI_i[i]*OVI*de - OVI_r[i]*OVI*de + OV_i[i]*OV*de;
        
            /* Already in number density, not mass density */
            total += rhs[j] * 16;
            total_e += OVI * 5;
        
        
            total_de += rhs[j] * 5;
        
        j++;
    
        // 
        // Species: OVII
        // 
        rhs[j] = OVIII_r[i]*OVIII*de - OVII_i[i]*OVII*de - OVII_r[i]*OVII*de + OVI_i[i]*OVI*de;
        
            /* Already in number density, not mass density */
            total += rhs[j] * 16;
            total_e += OVII * 6;
        
        
            total_de += rhs[j] * 6;
        
        j++;
    
        // 
        // Species: OVIII
        // 
        rhs[j] = OIX_r[i]*OIX*de - OVIII_i[i]*OVIII*de - OVIII_r[i]*OVIII*de + OVII_i[i]*OVII*de;
        
            /* Already in number density, not mass density */
            total += rhs[j] * 16;
            total_e += OVIII * 7;
        
        
            total_de += rhs[j] * 7;
        
        j++;
    
        // 
        // Species: OIX
        // 
        rhs[j] = -OIX_r[i]*OIX*de + OVIII_i[i]*OVIII*de;
        
            /* Already in number density, not mass density */
            total += rhs[j] * 16;
            total_e += OIX * 8;
        
        
            total_de += rhs[j] * 8;
        
        j++;
    
    }  
    return 0;
}




int calculate_jacobian_oxygen(double *input, double *Joutput,
        int nstrip, int nchem, void *sdata)
{
    /* We iterate over all of the rates */
    /* Calculate temperature first */
    oxygen_data *data = (oxygen_data*)sdata;

    int i, j;
    oxygen_calculate_temperature(data, input, nstrip, nchem);

    oxygen_interpolate_rates(data, nstrip);

    /* Now we set up some temporaries */
    double *Tge = data->dTs_ge;
    double *OI_i = data->rs_OI_i;
    double *rOI_i = data->drs_OI_i;
    double *OII_i = data->rs_OII_i;
    double *rOII_i = data->drs_OII_i;
    double *OII_r = data->rs_OII_r;
    double *rOII_r = data->drs_OII_r;
    double *OIII_i = data->rs_OIII_i;
    double *rOIII_i = data->drs_OIII_i;
    double *OIII_r = data->rs_OIII_r;
    double *rOIII_r = data->drs_OIII_r;
    double *OIV_i = data->rs_OIV_i;
    double *rOIV_i = data->drs_OIV_i;
    double *OIV_r = data->rs_OIV_r;
    double *rOIV_r = data->drs_OIV_r;
    double *OIX_r = data->rs_OIX_r;
    double *rOIX_r = data->drs_OIX_r;
    double *OV_i = data->rs_OV_i;
    double *rOV_i = data->drs_OV_i;
    double *OV_r = data->rs_OV_r;
    double *rOV_r = data->drs_OV_r;
    double *OVI_i = data->rs_OVI_i;
    double *rOVI_i = data->drs_OVI_i;
    double *OVI_r = data->rs_OVI_r;
    double *rOVI_r = data->drs_OVI_r;
    double *OVII_i = data->rs_OVII_i;
    double *rOVII_i = data->drs_OVII_i;
    double *OVII_r = data->rs_OVII_r;
    double *rOVII_r = data->drs_OVII_r;
    double *OVIII_i = data->rs_OVIII_i;
    double *rOVIII_i = data->drs_OVIII_i;
    double *OVIII_r = data->rs_OVIII_r;
    double *rOVIII_r = data->drs_OVIII_r;
    double *OI_c_OI_c = data->cs_OI_c_OI_c;
    double *rOI_c_OI_c = data->dcs_OI_c_OI_c;
    double *OII_c_OII_c = data->cs_OII_c_OII_c;
    double *rOII_c_OII_c = data->dcs_OII_c_OII_c;
    double *OIII_c_OIII_c = data->cs_OIII_c_OIII_c;
    double *rOIII_c_OIII_c = data->dcs_OIII_c_OIII_c;
    double *OIV_c_OIV_c = data->cs_OIV_c_OIV_c;
    double *rOIV_c_OIV_c = data->dcs_OIV_c_OIV_c;
    double *OIX_c_OIX_c = data->cs_OIX_c_OIX_c;
    double *rOIX_c_OIX_c = data->dcs_OIX_c_OIX_c;
    double *OV_c_OV_c = data->cs_OV_c_OV_c;
    double *rOV_c_OV_c = data->dcs_OV_c_OV_c;
    double *OVI_c_OVI_c = data->cs_OVI_c_OVI_c;
    double *rOVI_c_OVI_c = data->dcs_OVI_c_OVI_c;
    double *OVII_c_OVII_c = data->cs_OVII_c_OVII_c;
    double *rOVII_c_OVII_c = data->dcs_OVII_c_OVII_c;
    double *OVIII_c_OVIII_c = data->cs_OVIII_c_OVIII_c;
    double *rOVIII_c_OVIII_c = data->dcs_OVIII_c_OVIII_c;
    double de;
    double ge;
    double OI;
    double OII;
    double OIII;
    double OIV;
    double OV;
    double OVI;
    double OVII;
    double OVIII;
    double OIX;

    double mh = 1.67e-24;
    double total, mdensity;
    for (i = 0; i<nstrip; i++) {
        j = i * nchem;
	total = mdensity = 0.0;
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
    
	    OI = input[j];
        if (OI < 0.0) {
          fprintf(stderr, "JNegative[%d][OI] = % 0.16g [%d]\n",
            i, OI, j);
          /*OI = 0.0;*/
          OI = 1e-20;
          return 1;
        }
	
        
          total+=OI * 16;
        
	
        j++;
    
	    OII = input[j];
        if (OII < 0.0) {
          fprintf(stderr, "JNegative[%d][OII] = % 0.16g [%d]\n",
            i, OII, j);
          /*OII = 0.0;*/
          OII = 1e-20;
          return 1;
        }
	
        
          total+=OII * 16;
        
	
        j++;
    
	    OIII = input[j];
        if (OIII < 0.0) {
          fprintf(stderr, "JNegative[%d][OIII] = % 0.16g [%d]\n",
            i, OIII, j);
          /*OIII = 0.0;*/
          OIII = 1e-20;
          return 1;
        }
	
        
          total+=OIII * 16;
        
	
        j++;
    
	    OIV = input[j];
        if (OIV < 0.0) {
          fprintf(stderr, "JNegative[%d][OIV] = % 0.16g [%d]\n",
            i, OIV, j);
          /*OIV = 0.0;*/
          OIV = 1e-20;
          return 1;
        }
	
        
          total+=OIV * 16;
        
	
        j++;
    
	    OV = input[j];
        if (OV < 0.0) {
          fprintf(stderr, "JNegative[%d][OV] = % 0.16g [%d]\n",
            i, OV, j);
          /*OV = 0.0;*/
          OV = 1e-20;
          return 1;
        }
	
        
          total+=OV * 16;
        
	
        j++;
    
	    OVI = input[j];
        if (OVI < 0.0) {
          fprintf(stderr, "JNegative[%d][OVI] = % 0.16g [%d]\n",
            i, OVI, j);
          /*OVI = 0.0;*/
          OVI = 1e-20;
          return 1;
        }
	
        
          total+=OVI * 16;
        
	
        j++;
    
	    OVII = input[j];
        if (OVII < 0.0) {
          fprintf(stderr, "JNegative[%d][OVII] = % 0.16g [%d]\n",
            i, OVII, j);
          /*OVII = 0.0;*/
          OVII = 1e-20;
          return 1;
        }
	
        
          total+=OVII * 16;
        
	
        j++;
    
	    OVIII = input[j];
        if (OVIII < 0.0) {
          fprintf(stderr, "JNegative[%d][OVIII] = % 0.16g [%d]\n",
            i, OVIII, j);
          /*OVIII = 0.0;*/
          OVIII = 1e-20;
          return 1;
        }
	
        
          total+=OVIII * 16;
        
	
        j++;
    
	    OIX = input[j];
        if (OIX < 0.0) {
          fprintf(stderr, "JNegative[%d][OIX] = % 0.16g [%d]\n",
            i, OIX, j);
          /*OIX = 0.0;*/
          OIX = 1e-20;
          return 1;
        }
	
        
          total+=OIX * 16;
        
	
        j++;
    
        mdensity = total * mh;
        
        j = i * nchem * nchem;
        // 
        // Species: de
        //
            // de by de
            Joutput[j] = OIII_i[i]*OIII - OIII_r[i]*OIII + OII_i[i]*OII - OII_r[i]*OII + OIV_i[i]*OIV - OIV_r[i]*OIV - OIX_r[i]*OIX + OI_i[i]*OI + OVIII_i[i]*OVIII - OVIII_r[i]*OVIII + OVII_i[i]*OVII - OVII_r[i]*OVII + OVI_i[i]*OVI - OVI_r[i]*OVI + OV_i[i]*OV - OV_r[i]*OV;
	    
	    
            j++;
            // ge by de
            Joutput[j] = -OI*OI_c_OI_c[i] - OII*OII_c_OII_c[i] - OIII*OIII_c_OIII_c[i] - OIV*OIV_c_OIV_c[i] - OIX*OIX_c_OIX_c[i] - OV*OV_c_OV_c[i] - OVI*OVI_c_OVI_c[i] - OVII*OVII_c_OVII_c[i] - OVIII*OVIII_c_OVIII_c[i];
	    
	    Joutput[j] /= mdensity;
	    
	    
            j++;
            // OI by de
            Joutput[j] = OII_r[i]*OII - OI_i[i]*OI;
	    
	    
            j++;
            // OII by de
            Joutput[j] = OIII_r[i]*OIII - OII_i[i]*OII - OII_r[i]*OII + OI_i[i]*OI;
	    
	    
            j++;
            // OIII by de
            Joutput[j] = -OIII_i[i]*OIII - OIII_r[i]*OIII + OII_i[i]*OII + OIV_r[i]*OIV;
	    
	    
            j++;
            // OIV by de
            Joutput[j] = OIII_i[i]*OIII - OIV_i[i]*OIV - OIV_r[i]*OIV + OV_r[i]*OV;
	    
	    
            j++;
            // OV by de
            Joutput[j] = OIV_i[i]*OIV + OVI_r[i]*OVI - OV_i[i]*OV - OV_r[i]*OV;
	    
	    
            j++;
            // OVI by de
            Joutput[j] = OVII_r[i]*OVII - OVI_i[i]*OVI - OVI_r[i]*OVI + OV_i[i]*OV;
	    
	    
            j++;
            // OVII by de
            Joutput[j] = OVIII_r[i]*OVIII - OVII_i[i]*OVII - OVII_r[i]*OVII + OVI_i[i]*OVI;
	    
	    
            j++;
            // OVIII by de
            Joutput[j] = OIX_r[i]*OIX - OVIII_i[i]*OVIII - OVIII_r[i]*OVIII + OVII_i[i]*OVII;
	    
	    
            j++;
            // OIX by de
            Joutput[j] = -OIX_r[i]*OIX + OVIII_i[i]*OVIII;
	    
	    
            j++;
    
        // 
        // Species: ge
        //
            // de by ge
            Joutput[j] = OI*de*rOI_i[i] + OII*de*rOII_i[i] - OII*de*rOII_r[i] + OIII*de*rOIII_i[i] - OIII*de*rOIII_r[i] + OIV*de*rOIV_i[i] - OIV*de*rOIV_r[i] - OIX*de*rOIX_r[i] + OV*de*rOV_i[i] - OV*de*rOV_r[i] + OVI*de*rOVI_i[i] - OVI*de*rOVI_r[i] + OVII*de*rOVII_i[i] - OVII*de*rOVII_r[i] + OVIII*de*rOVIII_i[i] - OVIII*de*rOVIII_r[i];
	    
	    
            Joutput[j] *= Tge[i];
            
            j++;
            // ge by ge
            Joutput[j] = 0;
	    
	    Joutput[j] /= mdensity;
	    
	    
            Joutput[j] *= Tge[i];
            
            j++;
            // OI by ge
            Joutput[j] = -OI*de*rOI_i[i] + OII*de*rOII_r[i];
	    
	    
            Joutput[j] *= Tge[i];
            
            j++;
            // OII by ge
            Joutput[j] = OI*de*rOI_i[i] - OII*de*rOII_i[i] - OII*de*rOII_r[i] + OIII*de*rOIII_r[i];
	    
	    
            Joutput[j] *= Tge[i];
            
            j++;
            // OIII by ge
            Joutput[j] = OII*de*rOII_i[i] - OIII*de*rOIII_i[i] - OIII*de*rOIII_r[i] + OIV*de*rOIV_r[i];
	    
	    
            Joutput[j] *= Tge[i];
            
            j++;
            // OIV by ge
            Joutput[j] = OIII*de*rOIII_i[i] - OIV*de*rOIV_i[i] - OIV*de*rOIV_r[i] + OV*de*rOV_r[i];
	    
	    
            Joutput[j] *= Tge[i];
            
            j++;
            // OV by ge
            Joutput[j] = OIV*de*rOIV_i[i] - OV*de*rOV_i[i] - OV*de*rOV_r[i] + OVI*de*rOVI_r[i];
	    
	    
            Joutput[j] *= Tge[i];
            
            j++;
            // OVI by ge
            Joutput[j] = OV*de*rOV_i[i] - OVI*de*rOVI_i[i] - OVI*de*rOVI_r[i] + OVII*de*rOVII_r[i];
	    
	    
            Joutput[j] *= Tge[i];
            
            j++;
            // OVII by ge
            Joutput[j] = OVI*de*rOVI_i[i] - OVII*de*rOVII_i[i] - OVII*de*rOVII_r[i] + OVIII*de*rOVIII_r[i];
	    
	    
            Joutput[j] *= Tge[i];
            
            j++;
            // OVIII by ge
            Joutput[j] = OIX*de*rOIX_r[i] + OVII*de*rOVII_i[i] - OVIII*de*rOVIII_i[i] - OVIII*de*rOVIII_r[i];
	    
	    
            Joutput[j] *= Tge[i];
            
            j++;
            // OIX by ge
            Joutput[j] = -OIX*de*rOIX_r[i] + OVIII*de*rOVIII_i[i];
	    
	    
            Joutput[j] *= Tge[i];
            
            j++;
    
        // 
        // Species: OI
        //
            // de by OI
            Joutput[j] = OI_i[i]*de;
	    
	    
            j++;
            // ge by OI
            Joutput[j] = -OI_c_OI_c[i]*de;
	    
	    Joutput[j] /= mdensity;
	    
	    
            j++;
            // OI by OI
            Joutput[j] = -OI_i[i]*de;
	    
	    
            j++;
            // OII by OI
            Joutput[j] = OI_i[i]*de;
	    
	    
            j++;
            // OIII by OI
            Joutput[j] = 0;
	    
	    
            j++;
            // OIV by OI
            Joutput[j] = 0;
	    
	    
            j++;
            // OV by OI
            Joutput[j] = 0;
	    
	    
            j++;
            // OVI by OI
            Joutput[j] = 0;
	    
	    
            j++;
            // OVII by OI
            Joutput[j] = 0;
	    
	    
            j++;
            // OVIII by OI
            Joutput[j] = 0;
	    
	    
            j++;
            // OIX by OI
            Joutput[j] = 0;
	    
	    
            j++;
    
        // 
        // Species: OII
        //
            // de by OII
            Joutput[j] = OII_i[i]*de - OII_r[i]*de;
	    
	    
            j++;
            // ge by OII
            Joutput[j] = -OII_c_OII_c[i]*de;
	    
	    Joutput[j] /= mdensity;
	    
	    
            j++;
            // OI by OII
            Joutput[j] = OII_r[i]*de;
	    
	    
            j++;
            // OII by OII
            Joutput[j] = -OII_i[i]*de - OII_r[i]*de;
	    
	    
            j++;
            // OIII by OII
            Joutput[j] = OII_i[i]*de;
	    
	    
            j++;
            // OIV by OII
            Joutput[j] = 0;
	    
	    
            j++;
            // OV by OII
            Joutput[j] = 0;
	    
	    
            j++;
            // OVI by OII
            Joutput[j] = 0;
	    
	    
            j++;
            // OVII by OII
            Joutput[j] = 0;
	    
	    
            j++;
            // OVIII by OII
            Joutput[j] = 0;
	    
	    
            j++;
            // OIX by OII
            Joutput[j] = 0;
	    
	    
            j++;
    
        // 
        // Species: OIII
        //
            // de by OIII
            Joutput[j] = OIII_i[i]*de - OIII_r[i]*de;
	    
	    
            j++;
            // ge by OIII
            Joutput[j] = -OIII_c_OIII_c[i]*de;
	    
	    Joutput[j] /= mdensity;
	    
	    
            j++;
            // OI by OIII
            Joutput[j] = 0;
	    
	    
            j++;
            // OII by OIII
            Joutput[j] = OIII_r[i]*de;
	    
	    
            j++;
            // OIII by OIII
            Joutput[j] = -OIII_i[i]*de - OIII_r[i]*de;
	    
	    
            j++;
            // OIV by OIII
            Joutput[j] = OIII_i[i]*de;
	    
	    
            j++;
            // OV by OIII
            Joutput[j] = 0;
	    
	    
            j++;
            // OVI by OIII
            Joutput[j] = 0;
	    
	    
            j++;
            // OVII by OIII
            Joutput[j] = 0;
	    
	    
            j++;
            // OVIII by OIII
            Joutput[j] = 0;
	    
	    
            j++;
            // OIX by OIII
            Joutput[j] = 0;
	    
	    
            j++;
    
        // 
        // Species: OIV
        //
            // de by OIV
            Joutput[j] = OIV_i[i]*de - OIV_r[i]*de;
	    
	    
            j++;
            // ge by OIV
            Joutput[j] = -OIV_c_OIV_c[i]*de;
	    
	    Joutput[j] /= mdensity;
	    
	    
            j++;
            // OI by OIV
            Joutput[j] = 0;
	    
	    
            j++;
            // OII by OIV
            Joutput[j] = 0;
	    
	    
            j++;
            // OIII by OIV
            Joutput[j] = OIV_r[i]*de;
	    
	    
            j++;
            // OIV by OIV
            Joutput[j] = -OIV_i[i]*de - OIV_r[i]*de;
	    
	    
            j++;
            // OV by OIV
            Joutput[j] = OIV_i[i]*de;
	    
	    
            j++;
            // OVI by OIV
            Joutput[j] = 0;
	    
	    
            j++;
            // OVII by OIV
            Joutput[j] = 0;
	    
	    
            j++;
            // OVIII by OIV
            Joutput[j] = 0;
	    
	    
            j++;
            // OIX by OIV
            Joutput[j] = 0;
	    
	    
            j++;
    
        // 
        // Species: OV
        //
            // de by OV
            Joutput[j] = OV_i[i]*de - OV_r[i]*de;
	    
	    
            j++;
            // ge by OV
            Joutput[j] = -OV_c_OV_c[i]*de;
	    
	    Joutput[j] /= mdensity;
	    
	    
            j++;
            // OI by OV
            Joutput[j] = 0;
	    
	    
            j++;
            // OII by OV
            Joutput[j] = 0;
	    
	    
            j++;
            // OIII by OV
            Joutput[j] = 0;
	    
	    
            j++;
            // OIV by OV
            Joutput[j] = OV_r[i]*de;
	    
	    
            j++;
            // OV by OV
            Joutput[j] = -OV_i[i]*de - OV_r[i]*de;
	    
	    
            j++;
            // OVI by OV
            Joutput[j] = OV_i[i]*de;
	    
	    
            j++;
            // OVII by OV
            Joutput[j] = 0;
	    
	    
            j++;
            // OVIII by OV
            Joutput[j] = 0;
	    
	    
            j++;
            // OIX by OV
            Joutput[j] = 0;
	    
	    
            j++;
    
        // 
        // Species: OVI
        //
            // de by OVI
            Joutput[j] = OVI_i[i]*de - OVI_r[i]*de;
	    
	    
            j++;
            // ge by OVI
            Joutput[j] = -OVI_c_OVI_c[i]*de;
	    
	    Joutput[j] /= mdensity;
	    
	    
            j++;
            // OI by OVI
            Joutput[j] = 0;
	    
	    
            j++;
            // OII by OVI
            Joutput[j] = 0;
	    
	    
            j++;
            // OIII by OVI
            Joutput[j] = 0;
	    
	    
            j++;
            // OIV by OVI
            Joutput[j] = 0;
	    
	    
            j++;
            // OV by OVI
            Joutput[j] = OVI_r[i]*de;
	    
	    
            j++;
            // OVI by OVI
            Joutput[j] = -OVI_i[i]*de - OVI_r[i]*de;
	    
	    
            j++;
            // OVII by OVI
            Joutput[j] = OVI_i[i]*de;
	    
	    
            j++;
            // OVIII by OVI
            Joutput[j] = 0;
	    
	    
            j++;
            // OIX by OVI
            Joutput[j] = 0;
	    
	    
            j++;
    
        // 
        // Species: OVII
        //
            // de by OVII
            Joutput[j] = OVII_i[i]*de - OVII_r[i]*de;
	    
	    
            j++;
            // ge by OVII
            Joutput[j] = -OVII_c_OVII_c[i]*de;
	    
	    Joutput[j] /= mdensity;
	    
	    
            j++;
            // OI by OVII
            Joutput[j] = 0;
	    
	    
            j++;
            // OII by OVII
            Joutput[j] = 0;
	    
	    
            j++;
            // OIII by OVII
            Joutput[j] = 0;
	    
	    
            j++;
            // OIV by OVII
            Joutput[j] = 0;
	    
	    
            j++;
            // OV by OVII
            Joutput[j] = 0;
	    
	    
            j++;
            // OVI by OVII
            Joutput[j] = OVII_r[i]*de;
	    
	    
            j++;
            // OVII by OVII
            Joutput[j] = -OVII_i[i]*de - OVII_r[i]*de;
	    
	    
            j++;
            // OVIII by OVII
            Joutput[j] = OVII_i[i]*de;
	    
	    
            j++;
            // OIX by OVII
            Joutput[j] = 0;
	    
	    
            j++;
    
        // 
        // Species: OVIII
        //
            // de by OVIII
            Joutput[j] = OVIII_i[i]*de - OVIII_r[i]*de;
	    
	    
            j++;
            // ge by OVIII
            Joutput[j] = -OVIII_c_OVIII_c[i]*de;
	    
	    Joutput[j] /= mdensity;
	    
	    
            j++;
            // OI by OVIII
            Joutput[j] = 0;
	    
	    
            j++;
            // OII by OVIII
            Joutput[j] = 0;
	    
	    
            j++;
            // OIII by OVIII
            Joutput[j] = 0;
	    
	    
            j++;
            // OIV by OVIII
            Joutput[j] = 0;
	    
	    
            j++;
            // OV by OVIII
            Joutput[j] = 0;
	    
	    
            j++;
            // OVI by OVIII
            Joutput[j] = 0;
	    
	    
            j++;
            // OVII by OVIII
            Joutput[j] = OVIII_r[i]*de;
	    
	    
            j++;
            // OVIII by OVIII
            Joutput[j] = -OVIII_i[i]*de - OVIII_r[i]*de;
	    
	    
            j++;
            // OIX by OVIII
            Joutput[j] = OVIII_i[i]*de;
	    
	    
            j++;
    
        // 
        // Species: OIX
        //
            // de by OIX
            Joutput[j] = -OIX_r[i]*de;
	    
	    
            j++;
            // ge by OIX
            Joutput[j] = -OIX_c_OIX_c[i]*de;
	    
	    Joutput[j] /= mdensity;
	    
	    
            j++;
            // OI by OIX
            Joutput[j] = 0;
	    
	    
            j++;
            // OII by OIX
            Joutput[j] = 0;
	    
	    
            j++;
            // OIII by OIX
            Joutput[j] = 0;
	    
	    
            j++;
            // OIV by OIX
            Joutput[j] = 0;
	    
	    
            j++;
            // OV by OIX
            Joutput[j] = 0;
	    
	    
            j++;
            // OVI by OIX
            Joutput[j] = 0;
	    
	    
            j++;
            // OVII by OIX
            Joutput[j] = 0;
	    
	    
            j++;
            // OVIII by OIX
            Joutput[j] = OIX_r[i]*de;
	    
	    
            j++;
            // OIX by OIX
            Joutput[j] = -OIX_r[i]*de;
	    
	    
            j++;
    
    }

    return 0;
    
}
