#include "primordial_solver.h"
#include "hdf5.h"
#include "hdf5_hl.h"

#define kb      RCONST(1.3806504e-16)
#define mh      RCONST(1.67e-24)
#define gamma   RCONST(5.0/3.0)
#define nchem   RCONST(10)

void initialize_long_ydata(double *ydata, unsigned long NSYSTEM, double density, double temperature, double h2fraction, double efraction)
{
    double tiny_fraction = 1.0e-20;
    for (unsigned long i = 0; i < NSYSTEM; i++)
    {
        // H2I
        ydata[i*nchem]   = density*h2fraction*0.76;
        // H2II
        ydata[i*nchem+1] = density*tiny_fraction;
        // HI
        ydata[i*nchem+2] = density*0.76*(1-h2fraction);
        // HII
        ydata[i*nchem+3] = density*efraction;
        // H-
        ydata[i*nchem+4] = density*tiny_fraction;
        // HeI
        ydata[i*nchem+5] = density*0.24;
        // HeII
        ydata[i*nchem+6] = density*tiny_fraction;
        // HeIII
        ydata[i*nchem+7] = density*tiny_fraction;
        // de
        ydata[i*nchem+8] = density*efraction;
        // ge
        ydata[i*nchem+9] = 3./2.*kb* temperature / mh; 

    }
}


int run_dengo_struct(double density, double temperature, double h2fraction, double efraction, unsigned long dims, double *output)
{
    primordial_data *data = primordial_setup_data("primordial_tables.h5", NULL, NULL);

    double *input, *rtol, *atol, *temp; 
    input = new double[dims*nchem];
    rtol  = new double[1];
    atol  = new double[dims*nchem];
    temp  = new double[dims];

    initialize_long_ydata(input, dims, density, temperature, h2fraction, efraction);
    rtol[0] = 1.0e-5;

    for (int i = 0; i<nchem; i++) printf("input[%d] = %0.5g\n", i, input[i]);

    // main evolution;
    double dtf = pow(6.67e-8*mh*density, -0.5);
    double z   = 0.0;
    double dt   = 0.0;

    int flag = dengo_evolve_primordial(dtf, dt, z, input, rtol, atol, dims, data, temp);


    for (int i = 0; i<nchem; i++)
    {
        printf("output[%d] = %0.5g\n", i, input[i]);
        output[i] = input[i];
    }
    // supposingly this would return the output in "input" array

    /*
       double diff;
       for (unsigned long i = 0; i<dims; i++)
       {
       diff = (input[i] - input[i%nchem]) / input[i%nchem];
       if (fabs(diff) > rtol[0]){
       printf("output %lu diff[%d] = %0.5g; ref = %0.5g; out = %0.5g\n", i,  i%nchem, diff, input[i], input[i%nchem]);
       }
       }
       */

}

int test_scaling_dims()
{
    unsigned long dims = 4096;
    double density, temperature, h2fraction, efraction;
    double start, end;
    double cpu_time_used;

    float milliseconds = 0;

    int ntimes = 8;
    double *timelapsed = new double[ntimes];
    long *ncells     = new long[ntimes];
    double *output   = new double[ntimes*nchem];

    density = 1.0;
    temperature = 1.0;
    h2fraction = 1e-4;
    efraction = 1e-4;

    for (int i = 0; i < ntimes; i++)
    {
        start = omp_get_wtime();

        run_dengo_struct(density, temperature, h2fraction, efraction, dims, &output[i*nchem]);

        dims *= 2;
        end = omp_get_wtime();
        cpu_time_used = end - start;

        printf("took %f seconds to execute %lu \n", cpu_time_used, dims); 
        // measured the time lapsed
        // and the cells needed
        ncells[i] = dims;
        timelapsed[i] = cpu_time_used;
    }
    // create a file
    hid_t file_id, dataset;
    hid_t datatype, dataspace;
    hsize_t dimsf[1];
    herr_t status;
    dimsf[0] = ntimes;

    file_id = H5Fcreate("scaling.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    H5LTmake_dataset_double(file_id, "/time", 1, dimsf, timelapsed);
    H5LTmake_dataset_long(file_id, "/ncells", 1, dimsf, ncells);
    H5Fclose(file_id);


}

int grid_performance()
{
    unsigned long dims = 16384;
    double density, temperature, h2fraction, efraction;
    double d, T;
    double start, end;
    double cpu_time_used;


    density = 1.0;
    temperature = 10.0;
    h2fraction = 1e-4;
    efraction = 1e-4;

    const int nd = 9;
    const int nT = 9;
    double *timelapsed = new double[nd*nT];
    double *darray      = new double[nd*nT];
    double *Tarray      = new double[nd*nT];

    double *output = new double[nchem*nd*nT];

    for (int i = 0; i < nd; i++)
    {
        for (int j = 0; j < nT; j++)
        {
            d = density* pow(10., i);
            T = temperature* pow(2.2, j);

            start = omp_get_wtime();
            run_dengo_struct(d, T, h2fraction, efraction, dims, &output[(i*nT+j)*nchem]);
            end = omp_get_wtime();
            cpu_time_used = end - start;

            printf("d=%0.5g; T=%0.5g;  %f seconds %lu cells; %0.5g percell \n", d, T, cpu_time_used, dims, cpu_time_used/ dims); 

            //printf("took %f milliseconds to execute %lu; %0.5g percell \n", milliseconds, dims, milliseconds/ dims); 

            timelapsed[i*nT+j] = cpu_time_used;
            darray    [i*nT+j] = d;
            Tarray    [i*nT+j] = T;

        }
    }

    // create a file
    hid_t file_id, dataset;
    hid_t datatype, dataspace;
    hsize_t dimsf[1], dimsOut[1];
    herr_t status;
    dimsf[0] = nd*nT;
    dimsOut[0] = nd*nT*nchem;

    file_id = H5Fcreate("performance.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    H5LTmake_dataset_double(file_id, "/time", 1, dimsf, timelapsed);
    H5LTmake_dataset_double(file_id, "/density", 1, dimsf, darray);
    H5LTmake_dataset_double(file_id, "/temperature", 1, dimsf, Tarray);
    H5LTmake_dataset_double(file_id, "/output", 1, dimsOut, output);

    H5Fclose(file_id);

}

int writefile()
{
    double *Tarray = new double[10];
    double *darray = new double[10];
    
    for (int i = 0; i < 10; i++){
        Tarray[i] = i;
        darray[i] = i;
    }

    // create a file
    hid_t file_id, dataset;
    hid_t datatype, dataspace;
    hsize_t dimsf[1];
    herr_t status;
    dimsf[0] = 10;

    file_id = H5Fcreate("performance.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    //H5LTmake_dataset_double(file_id, "/time", 1, dimsf, timelapsed);
    H5LTmake_dataset_double(file_id, "/density", 1, dimsf, darray);
    H5LTmake_dataset_double(file_id, "/temperature", 1, dimsf, Tarray);

    H5Fclose(file_id);
}

int main()
{
    //run_dengo_struct(1, 1000, 1e-4, 1e-4,16384);
    //
    //writefile();
    test_scaling_dims();
    grid_performance();
}
