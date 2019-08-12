/**
 * \file solver_main.cu
 * \brief the generic main file for all GPU solvers
 *
 * \author Nicholas Curtis
 * \date 03/09/2015
 *
 * Contains main function, setup, initialization, logging, timing and driver functions
 */


/* Include common code. */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>
#include <complex.h>
#include <sys/types.h>
#include <sys/stat.h>
#ifdef DEBUG
#include <fenv.h>
#endif

/* Include CUDA libraries. */
#include <cuda.h>
#include <cuda_runtime.h>
#include <helper_cuda.h>
#include <cuComplex.h>

//our code
#include "header.cuh"
#include "timer.h"
//get our solver stuff
#include "solver.cuh"
#include "gpu_memory.cuh"
#include "read_initial_conditions.cuh"
#include "launch_bounds.cuh"
#include "mechanism.cuh"

#ifdef DIVERGENCE_TEST
    #include <assert.h>
    //! If #DIVERGENCE_TEST is defined, this creates a device array for tracking
    //  internal integrator steps per thread
    __device__ int integrator_steps[DIVERGENCE_TEST] = {0};
#endif



/**
 * \brief Writes state vectors to file
 * \param[in]           NUM                 the number of state vectors to write
 * \param[in]           t                   the current system time
 * \param[in]           y_host              the current state vectors
 * \param[in]           pFile               the opened binary file object
 *
 * The resulting file is updated as:
 * system time\n
 * temperature, mass fractions (State #1)\n
 * temperature, mass fractions (State #2)...
 */

/*
void write_log(int NUM, double t, const double* y_host, FILE* pFile)
{
    fwrite(&t, sizeof(double), 1, pFile);
    double buffer[NN];
    for (int j = 0; j < NUM; j++)
    {
        double Y_N = 1.0;
        buffer[0] = y_host[j];
        for (int i = 1; i < NSP; ++i)
        {
            buffer[i] = y_host[NUM * i + j];
            Y_N -= buffer[i];
        }
        #if NN == NSP + 1 //pyjac
        buffer[NSP] = Y_N;
        #endif
        apply_reverse_mask(&buffer[1]);
        fwrite(buffer, sizeof(double), NN, pFile);
    }
}
*/

/**
 * \brief A convienience method to copy memory between host pointers of different pitches, widths and heights.
 *        Enables easier use of CUDA's cudaMemcpy2D functions.
 *
 * \param[out]              dst             The destination array
 * \param[in]               pitch_dst       The width (in bytes) of the destination array.
                                            This corresponds to the padded number of IVPs to be solved.
 * \param[in]               src             The source pointer
 * \param[in]               pitch_src       The width (in bytes)  of the source array.
                                            This corresponds to the (non-padded) number of IVPs read by read_initial_conditions
 * \param[in]               offset          The offset within the source array (IVP index) to copy from.
                                            This is useful in the case (for large models) where the solver and state vector memory will not fit in device memory
                                            and the integration must be split into multiple kernel calls.
 * \param[in]               width           The size (in bytes) of memory to copy for each entry in the state vector
 * \param[in]               height          The number of entries in the state vector
 */
inline void memcpy2D_in(double* dst, const int pitch_dst, double const * src, const int pitch_src,
                                     const int offset, const size_t width, const int height) {
    for (int i = 0; i < height; ++i)
    {
        memcpy(dst, &src[offset], width);
        dst += pitch_dst;
        src += pitch_src;
    }
}

/**
 * \brief A convienience method to copy memory between host pointers of different pitches, widths and heights.
 *        Enables easier use of CUDA's cudaMemcpy2D functions.
 *
 * \param[out]              dst             The destination array
 * \param[in]               pitch_dst       The width (in bytes)  of the source array.
                                            This corresponds to the (non-padded) number of IVPs read by read_initial_conditions
 * \param[in]               src             The source pointer
 * \param[in]               pitch_src       The width (in bytes) of the destination array.
                                            This corresponds to the padded number of IVPs to be solved.
 * \param[in]               offset          The offset within the destination array (IVP index) to copy to.
                                            This is useful in the case (for large models) where the solver and state vector memory will not fit in device memory
                                            and the integration must be split into multiple kernel calls.
 * \param[in]               width           The size (in bytes) of memory to copy for each entry in the state vector
 * \param[in]               height          The number of entries in the state vector
 */
inline void memcpy2D_out(double* dst, const int pitch_dst, double const * src, const int pitch_src,
                                      const int offset, const size_t width, const int height) {
    for (int i = 0; i < height; ++i)
    {
        memcpy(&dst[offset], src, width);
        dst += pitch_dst;
        src += pitch_src;
    }
}

//////////////////////////////////////////////////////////////////////////////

/** Main function
 *
 * \param[in]       argc    command line argument count
 * \param[in]       argv    command line argument vector
 *
 * This allows running the integrators from the command line.  The syntax is as follows:\n
 * `./solver-name [num_threads] [num_IVPs]`\n
 * *  num_threads  [Optional, Default:1]
 *      *  The number OpenMP threads to utilize
 *      *  The number of threads cannot be greater than recognized by OpenMP via `omp_get_max_threads()`
 * *  num_IVPs     [Optional, Default:1]
 *      *  The number of initial value problems to solve.
 *      *  This must be less than the number of conditions in the data file if #SAME_IC is not defined.
 *      *  If #SAME_IC is defined, then the initial conditions in the mechanism files will be used.
 *
 */
int main (int argc, char *argv[])
{

//enable signaling NAN and other bad numerics tracking for easier debugging
#ifdef DEBUG
    feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#endif
    /** Number of independent systems */
    int NUM = 8192;
    double density = 1.0e10; // cm**-3
    double temperature = 1000.0; // K
    double fH2  = 1.0e-6;
    
    // check for problem size given as command line option
    if (argc > 2)
    {
        int problemsize = NUM;
        if (sscanf(argv[1], "%i", &problemsize) != 1 || (problemsize <= 0))
        {
            printf("Error: Problem size not in correct range\n");
            printf("Provide number greater than 0\n");
            exit(1);
        }
        NUM = problemsize;

        double rho = density;
        rho = atof(argv[2]);
        if ( rho <= 0)
        {
            printf("Error: density not in correct range\n");
            printf("Provide number greater than 0\n");
            exit(1);
        }
        density = rho;

        double T   = temperature;
        T  = atof(argv[3]);
        if (  T <= 0 )
        {
            printf("Error: temperature not in correct range\n");
            printf("Provide number greater than 0\n");
            exit(1);
        }
        temperature = T;

        double H2frac   = fH2;
        H2frac  = atof(argv[4]);
        if (  H2frac <= 0 || H2frac > 1.0 )
        {
            printf("Error: fH2 not in correct range\n");
            printf("Provide number between 1 and 0\n");
            exit(1);
        }
        fH2 = H2frac;

        printf("using user_specified test parameter\n");
        printf("NUM = %d \n",NUM);
        printf("density = %0.5g \n", density);
        printf("temperature = %0.5g\n", temperature);
        printf("fH2         = %0.5g\n", fH2);

    } else{ 
        printf("no input specified...\n");
        printf("using default test parameter\n");
        printf("NUM = %d \n",NUM);
        printf("density = %0.5g \n", density);
        printf("temperature = %0.5g\n", temperature);
        printf("fH2         = %0.5g\n", fH2);

    }

    // set & initialize device using command line argument (if any)
    cudaDeviceProp devProp;
    cudaErrorCheck (cudaSetDevice (0) );
    cudaErrorCheck (cudaGetDeviceProperties(&devProp, 0));
    

    cudaErrorCheck( cudaDeviceReset() );
    cudaErrorCheck( cudaPeekAtLastError() );
    cudaErrorCheck( cudaDeviceSynchronize() );

    #ifdef DIVERGENCE_TEST
        NUM = DIVERGENCE_TEST;
        assert(NUM % 32 == 0);
    #endif
    //bump up shared mem bank size
    cudaErrorCheck(cudaDeviceSetSharedMemConfig(cudaSharedMemBankSizeEightByte));
    //and L1 size
    cudaErrorCheck(cudaDeviceSetCacheConfig(cudaFuncCachePreferL1));

    size_t size_per_thread = required_mechanism_size() + required_solver_size();
    size_t free_mem = 0;
    size_t total_mem = 0;
    cudaErrorCheck( cudaMemGetInfo (&free_mem, &total_mem) );

    //conservatively estimate the maximum allowable threads
    int max_threads = int(floor(0.8 * ((double)free_mem -sizeof(cvklu_data) ) / ((double)size_per_thread)));
    int padded = min(NUM, max_threads);
    //padded is next factor of block size up
    padded = int(ceil(padded / float(TARGET_BLOCK_SIZE)) * TARGET_BLOCK_SIZE);
    if (padded == 0)
    {
        printf("Mechanism is too large to fit into global CUDA memory... exiting.");
        exit(-1);
    }

    solver_memory* host_solver, *device_solver;
    mechanism_memory* host_mech, *device_mech;
    host_solver = (solver_memory*)malloc(sizeof(solver_memory));
    host_mech = (mechanism_memory*)malloc(sizeof(mechanism_memory));

    const char * FileLocation = "cvklu_tables.h5";
    cvklu_data *host_rateData = cvklu_setup_data( FileLocation, NULL, NULL);


    initialize_gpu_memory(padded, &host_mech, &device_mech, &host_rateData);
    initialize_solver(padded, &host_solver, &device_solver);

    double* y_host;
    double* var_host;
    dengo_set_initial_conditions( density, temperature, fH2, NUM, &y_host, &var_host);
    double t0 = 2.992e15;
    double end_time = t0 / sqrt( density );
    double t_step = end_time;
    // printf("end_time: %0.5g \n", end_time );


    double *temperature_array;
    double *density_array;
    double *h2_optical_depth_array;
    dengo_set_additional_constant(density, temperature, NUM, &y_host, 
                                 &temperature_array, &density_array, &h2_optical_depth_array ); 
    double* temp_temp = 0;
    temp_temp = (double*)malloc(padded * sizeof(double));
    double* density_temp = 0;
    density_temp = (double*)malloc(padded * sizeof(double));
    double* h2_optical_depth_temp = 0;
    h2_optical_depth_temp = (double*)malloc(padded * sizeof(double));



    #ifdef SCALE_INPUT
    double* _y_host;
    _y_host = (double*) malloc(NUM*NSP*sizeof(double));
    
    for (int i = 0; i < NSP*NUM; i++){
       _y_host[i] = 1.0 / y_host[i];
    }
    #endif


    for (int i =0; i< NSP; i++){
        printf("y_init[%d] = %0.5g \n", i, y_host[i*NUM]);
    }


    //grid sizes
    dim3 dimBlock(TARGET_BLOCK_SIZE, 1);
    dim3 dimGrid (padded / TARGET_BLOCK_SIZE, 1 );
    int* result_flag = (int*)malloc(padded * sizeof(int));

    // print number of threads and block size
    printf ("grid size: %d \t block size: %d\n", padded / TARGET_BLOCK_SIZE, TARGET_BLOCK_SIZE);
    printf ("max_threads: %d \n", max_threads);
    printf ("padded: %d\n ", padded);
    
    double* y_temp = 0;
    y_temp = (double*)malloc(padded * NSP * sizeof(double));
    #ifdef SCALE_INPUT
    double* _y_temp = 0;
    double* ones = 0;
    _y_temp = (double*)malloc(padded * NSP * sizeof(double));
    ones    = (double*)malloc(padded * NSP * sizeof(double));

    for (int i = 0; i < NSP*padded; i++){
      ones[i] = 1.0;
    }
    #endif

    //////////////////////////////
    // start timer
    StartTimer();
    //////////////////////////////

    // set initial time
    double t = 0;
    double t_next = fmin(end_time, t_step);
    int numSteps = 0;

    int *rhs_call=0;
    int *jac_call=0;

    // time integration loop
    while (t + EPS < end_time)
    {
        numSteps++;

        int num_solved = 0;
        while (num_solved < NUM)
        {
            int num_cond = min(NUM - num_solved, padded);

            cudaErrorCheck( cudaMemcpy (host_mech->var, &var_host[num_solved],
                                        num_cond * sizeof(double), cudaMemcpyHostToDevice));

             // copy our memory into y_temp
            memcpy2D_in(y_temp, padded, y_host, NUM,
                            num_solved, num_cond * sizeof(double), NSP);

            memcpy2D_in(temp_temp, padded, temperature_array, NUM,
                            num_solved, num_cond * sizeof(double), 1);
            memcpy2D_in(density_temp, padded, density_array, NUM,
                            num_solved, num_cond * sizeof(double), 1);
            memcpy2D_in(h2_optical_depth_temp, padded, h2_optical_depth_array, NUM,
                            num_solved, num_cond * sizeof(double), 1);


	    // 

            #ifdef SCALE_INPUT
            // copy additional inv_scale into _y_temp
            memcpy2D_in( _y_temp, padded, _y_host, NUM,
                            num_solved, num_cond * sizeof(double), NSP);

           // initialize y as ones
            cudaErrorCheck( cudaMemcpy2D ( host_mech->y, padded * sizeof(double),
                                            ones, padded * sizeof(double),
                                            num_cond * sizeof(double), NSP,
                                            cudaMemcpyHostToDevice) );           
            cudaErrorCheck( cudaMemcpy2D ( host_mech->scale, padded * sizeof(double),
                                            y_temp, padded * sizeof(double),
                                            num_cond * sizeof(double), NSP,
                                            cudaMemcpyHostToDevice) );
            cudaErrorCheck( cudaMemcpy2D ( host_mech->inv_scale, padded * sizeof(double),
                                            _y_temp, padded * sizeof(double),
                                            num_cond * sizeof(double), NSP,
                                            cudaMemcpyHostToDevice) );

            cudaErrorCheck( cudaMemcpy2D (host_mech->temperature, padded * sizeof(double),
                                            temp_temp, padded * sizeof(double),
                                            num_cond * sizeof(double), 1,
                                            cudaMemcpyHostToDevice) );
            cudaErrorCheck( cudaMemcpy2D (host_mech->density, padded * sizeof(double),
                                            density_temp, padded * sizeof(double),
                                            num_cond * sizeof(double), 1,
                                            cudaMemcpyHostToDevice) );
            cudaErrorCheck( cudaMemcpy2D (host_mech->h2_optical_depth_approx, padded * sizeof(double),
                                            h2_optical_depth_temp, padded * sizeof(double),
                                            num_cond * sizeof(double), 1,
                                            cudaMemcpyHostToDevice) );


         
            cudaErrorCheck( cudaMemcpy( device_mech, host_mech, sizeof(mechanism_memory), cudaMemcpyHostToDevice) );


/*
            cudaErrorCheck( cudaMemcpy(device_mech->scale, y_temp, num_cond* NSP * sizeof(double), cudaMemcpyHostToDevice) );
            cudaErrorCheck( cudaMemcpy(device_mech->inv_scale, _y_temp, num_cond* NSP * sizeof(double), cudaMemcpyHostToDevice) );
*/            
            #else
            // transfer memory to GPU
            cudaErrorCheck( cudaMemcpy2D (host_mech->y, padded * sizeof(double),
                                            y_temp, padded * sizeof(double),
                                            num_cond * sizeof(double), NSP,
                                            cudaMemcpyHostToDevice) );

            cudaErrorCheck( cudaMemcpy2D (host_mech->temperature, padded * sizeof(double),
                                            temp_temp, padded * sizeof(double),
                                            num_cond * sizeof(double), 1,
                                            cudaMemcpyHostToDevice) );
            cudaErrorCheck( cudaMemcpy2D (host_mech->density, padded * sizeof(double),
                                            density_temp, padded * sizeof(double),
                                            num_cond * sizeof(double), 1,
                                            cudaMemcpyHostToDevice) );
            cudaErrorCheck( cudaMemcpy2D (host_mech->h2_optical_depth_approx, padded * sizeof(double),
                                            h2_optical_depth_temp, padded * sizeof(double),
                                            num_cond * sizeof(double), 1,
                                            cudaMemcpyHostToDevice) );

            cudaErrorCheck( cudaMemcpy( device_mech, host_mech, sizeof(mechanism_memory), cudaMemcpyHostToDevice) );

            #endif

            intDriver <<< dimGrid, dimBlock, SHARED_SIZE >>> (num_cond, t, t_next, host_mech->var, host_mech->y, device_mech, device_solver);
            cudaErrorCheck( cudaPeekAtLastError() );
            cudaErrorCheck( cudaDeviceSynchronize() );

            // copy the result flag back
            cudaErrorCheck( cudaMemcpy(result_flag, host_solver->result, num_cond * sizeof(int), cudaMemcpyDeviceToHost) );
            check_error(num_cond, result_flag);


            // transfer memory back to CPU
            #ifdef SCALE_INPUT
            cudaErrorCheck( cudaMemcpy2D (_y_temp, padded * sizeof(double),
                                            host_mech->y, padded * sizeof(double),
                                            num_cond * sizeof(double), NSP,
                                            cudaMemcpyDeviceToHost) );
            // now _y_temp holds the scaled output
            // and it is multiplied by the orginal scale y_temp
            // and scaled back to original values
            for (int i = 0; i < num_cond*NSP; i++ ){
                y_temp[i] *= _y_temp[i];
            }
            memcpy2D_out(y_host, NUM, y_temp, padded,
                            num_solved, num_cond * sizeof(double), NSP);

            #else
            cudaErrorCheck( cudaMemcpy2D (y_temp, padded * sizeof(double),
                                            host_mech->y, padded * sizeof(double),
                                            num_cond * sizeof(double), NSP,
                                            cudaMemcpyDeviceToHost) );

            memcpy2D_out(y_host, NUM, y_temp, padded,
                            num_solved, num_cond * sizeof(double), NSP);
            #endif
            num_solved += num_cond;


        }

        t = t_next;
       // t_next = fmin(end_time, pow(1.1, (numSteps + 1 )) * t_step);
        t_next = fmin(end_time, (numSteps + 1) * t_step);
        
        for (int i = 0; i < 10; i++){
           printf("y_final[%d] = %0.5g\n", i, y_host[i*NUM]);
        }
	printf("----------\n");

/*
        double perr;
        for (int sp = 0; sp < 10; sp++){
          for (int jj = 0; jj < NUM; jj++){
            // printf("%.5le\t%.15le\n", t, y_host[sp*NUM+jj]);
            perr =  fabs(y_host[sp*NUM+jj] - y_host[sp*NUM] ) / y_host[sp*NUM];
            if ( perr > RTOL){
              printf(" thread[%d] sp[%d] has error of %0.5g\n ", jj, sp, perr);
            }
          }
          printf("\n");
        }
        
        // check if within bounds
        if (y_host[0] != y_host[0] ){
          printf("Error, out of bounds.\n");
          printf("Time: %e, ind %d val %e\n", t, 0, y_host[0]);
          return 1;
        }

*/

    }
    // end timer
    double runtime = GetTimer();
    /////////////////////////////////



    runtime /= 1000.0;
    printf ("Time: %.15e sec\n", runtime);
    runtime = runtime / ((double)(numSteps));
    printf ("Time per step: %e (s)\t%.15e (s/thread)\n", runtime, runtime / NUM);

    free_gpu_memory(&host_mech, &device_mech);
    cleanup_solver(&host_solver, &device_solver);
    free(y_temp);
    free(host_mech);
    free(host_solver);
    free(result_flag);
    
    cudaErrorCheck( cudaDeviceReset() );
    cudaErrorCheck( cudaPeekAtLastError() );
 
    return 0;
}

