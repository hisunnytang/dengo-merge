#include "primordial_cuda_solver.h"
#include <time.h>

int PrintDeviceInfo() {
    int nDevices;

    cudaGetDeviceCount(&nDevices);
    for (int i = 0; i < nDevices; i++) {
        cudaDeviceProp prop;
        cudaGetDeviceProperties(&prop, i);
        printf("Device Number: %d\n", i);
        printf("  Device name: %s\n", prop.name);
        printf("  Memory Clock Rate (KHz): %d\n",
                prop.memoryClockRate);
        printf("  Memory Bus Width (bits): %d\n",
                prop.memoryBusWidth);
        printf("  Peak Memory Bandwidth (GB/s): %f\n\n",
                2.0*prop.memoryClockRate*(prop.memoryBusWidth/8)/1.0e6);

    }

}

int main(int argc, char* argv[])
{
    // read the rate data
    PrintDeviceInfo();

    cudaDeviceSynchronize();

    // test how the performance scales with
    // number of cells
    test_scaling_dims();
    // run the solver on grids of density and temperature
    grid_performance();

    primordial_cuda_data data = primordial_cuda_setup_data(NULL, NULL);
    primordial_cuda_read_cooling_tables( &data);
    primordial_cuda_read_rate_tables( &data);


    // experimental
    // FInding the optimal block config
    printf("\nlaunch temperature\n");
    launchTemperatureKernel(&data);
    printf("\nlaunch interpolation\n");
    launchInterpolationKernel(&data);
    printf("\nlaunch rhs\n");
    launchRhsKernel(&data);;
    printf("\nlaunch jacobian\n");
    launchJacobianKernel(&data);

    // printf("rk01 = %0.5g\n", data.r_k22[213]);
    // printf("h2mheat = %0.5g\n", data.c_h2formation_h2mheat[1020]);

    // test interpolation first
    // test_interpolation_kernel(data);

    // test temperature kerenel
    //    test_temperature_kernel(data);

    // test the rhs function
    // test_rhs_function(data);

    // initialize initial conditions
    // create a y_vec that holds NSYSTEM  * nchem elements
    // test_jacobian_function(data);

    // initialize yvec to see if we can have it print out accurate ydot

    // run_solver(NULL, NULL);


    // launch the test with DengoSolver
    double density, temperature, efraction, h2fraction;
    unsigned long dims;
    if (argc < 5){
        printf("error you have to specify input\n");
        return 1;
    }else{
        density = atof(argv[1]);
        temperature = atof(argv[2]);
        efraction   = atof(argv[3]);
        h2fraction  = atof(argv[4]);
        dims = (unsigned long)atol(argv[5]);
    }
    run_dengo_solver(density, temperature, efraction, h2fraction, dims);

    cudaEvent_t startCuda, stopCuda;
    cudaEventCreate(&startCuda);
    cudaEventCreate(&stopCuda);
    float milliseconds = 0;

    cudaEventRecord(startCuda);
    //run_dengo_struct(density, temperature, efraction, h2fraction, dims);

    cudaEventRecord(stopCuda);
    cudaEventSynchronize(stopCuda);
    cudaEventElapsedTime(&milliseconds, startCuda, stopCuda);
    printf("took %f milliseconds\n", milliseconds);
    cudaDeviceSynchronize();

}
