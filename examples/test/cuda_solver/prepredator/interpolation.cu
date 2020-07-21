// includes, system
#include <cstdlib>
#include <math.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <stdio.h>
// includes, cuda
#include <cuda.h>
#include <cuda_runtime.h>

// dengo specific header
#include "predator_prey_solver.h"


using namespace std;

#define BLOCK_SIZE 256

texture<float, 1, cudaReadModeElementType> data_d_texture_filtering;
texture<float, 1> data_d_texture;

/******************/
/* ERROR CHECKING */
/******************/
#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__);  }
inline void gpuAssert(cudaError_t code, char *file, int line, bool abort=true)
{
    if (code != cudaSuccess)
    {
        fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
        if (abort) { getchar(); exit(code);  }

    }

}

/************/
/* linspace */
/************/
// --- generates n equally spaced, increasing points between a and b and stores them in x
void linspace(float* x, float a, float b, int n) {
    float delta_x=(b-a)/(float)n;
    x[0]=a;
    for(int k=1;k<n;k++) x[k]=x[k-1]+delta_x;

}

/*************/
/* randspace */
/*************/
// --- generates n randomly spaced, increasing points between a and b and stores them in x
void randspace(float* x, float a, float b, int n) {
    float delta_x=(b-a)/(float)n;
    x[0]=a;
    for(int k=1;k<n;k++) x[k]=x[k-1]+delta_x+(((float)rand()/(float)RAND_MAX-0.5)*(1./(float)n));

}

/******************/
/* DATA GENERATOR */
/******************/
// --- Generates N complex random data points, with real and imaginary parts ranging in (0.f,1.f)
void Data_Generator(float* data, int N) {
    for(int k=0;k<N;k++) {
        data[k]=(float)rand()/(float)RAND_MAX;

    }

}

/* LINEAR INTERPOLATION FUNCTION - CPU */
/***************************************/
void linear_interpolation_function_CPU(float* result_GPU, float* data, float* x_in, float* x_out, int M, int N){

    float a;
    for(int j=0; j<N; j++){
        int k = floor(x_out[j]+M/2);
        a = x_out[j]+M/2-floor(x_out[j]+M/2);
        result_GPU[j] = a * data[k+1] + (-data[k] * a + data[k]);

    }


}

/***************************************************************/
/* LINEAR INTERPOLATION KERNEL FUNCTION - GPU - TEXTURE MEMORY */
/***************************************************************/
__global__ void linear_interpolation_kernel_function_GPU_texture(float* __restrict__ result_d, const float* __restrict__ x_out_d, const int M, const int N)
{
    int j = threadIdx.x + blockDim.x * blockIdx.x;

    if(j<N)
    {
        float reg_x_out = x_out_d[j]+M/2;
        int k = __float2int_rz(reg_x_out);
        float a = reg_x_out - truncf(reg_x_out);
        float dk = tex1Dfetch(data_d_texture,k);
        float dkp1 = tex1Dfetch(data_d_texture,k+1);
        result_d[j] = a * dkp1 + (-dk * a + dk);

    }

}

/********************************************************/
/* LINEAR INTERPOLATION FUNCTION - GPU - TEXTURE MEMORY */
/********************************************************/
void linear_interpolation_function_GPU_texture(float* result_d, float* data_d, float* x_in_d, float* x_out_d, int M, int N){

    cudaBindTexture(NULL, data_d_texture, data_d, M*sizeof(float));

    dim3 dimBlock(BLOCK_SIZE,1); dim3 dimGrid(N/BLOCK_SIZE + (N%BLOCK_SIZE == 0 ? 0:1),1);
    linear_interpolation_kernel_function_GPU_texture<<<dimGrid,dimBlock>>>(result_d, x_out_d, M, N);
    gpuErrchk(cudaPeekAtLastError());
    gpuErrchk(cudaDeviceSynchronize());

}

__global__
void interpolation_kernel(double *rates_out, double *rate_in, int *bin_id, double *Tdef, int N)
{

    int j = threadIdx.x + blockDim.x*blockIdx.x;
    int k;

    // https://forums.developer.nvidia.com/t/accuracy-of-1d-linear-interpolation-by-cuda-texture-interpolation/28013/9
    // so the interpolation can be mapped to exactly two FMAs?
    if (j < N)
    {
        k = bin_id[j];
        rates_out[j] = Tdef[k]*rate_in[k+1] + (-rate_in[k]*Tdef[k] + rate_in[k]);
    }
}

__global__ 
void find_binID_Tdef(int *bin_id, double *Tdef, double *logTs, double dbin, double idbin, double lb)
{
    int j = threadIdx.x + blockDim.x*blockIdx.x;
    double t1;

    if (j < N)
    {
        int bId = __float2int_rz(idbin *logTs[j] - lb);
        t1 = lb + bId*dbin;
        bin_id[j] = bId;
        Tdef  [j] = (logTs[j] - t1) * idbin;
    }

}

void linear_interpolation_function_GPU_texture(* result_d, predator_prey_data* data_d, float* x_in_d, float* x_out_d, int M, int N){

    cudaBindTexture(NULL, data_d_texture, data_d, M*sizeof(float));

    dim3 dimBlock(BLOCK_SIZE,1); dim3 dimGrid(N/BLOCK_SIZE + (N%BLOCK_SIZE == 0 ? 0:1),1);
    linear_interpolation_kernel_function_GPU_texture<<<dimGrid,dimBlock>>>(result_d, x_out_d, M, N);
    gpuErrchk(cudaPeekAtLastError());
    gpuErrchk(cudaDeviceSynchronize());

}


// now something simple to read in the data
// and do interpolation to create a basline

int main()
{

    // read the data in from the hdf5;
    predator_prey_data *data = predator_prey_setup_data(NULL, NULL, NULL);

    /*
    // something temporary for now to keep track of the temperature
    int nstrip = 32;

    for (int i = 0; i < nstrip; i++)
    {
        data->Ts[0][i] = 100.0;
        data->logTs[0][i] = log(data->Ts[0][i]);
        data->invTs[0][i] = 1.0 / data->Ts[0][i];
    }

    // interpolate rates typical CPU code
    predator_prey_interpolate_rates(data, nstrip);

    for (int i = 0; i < nstrip; i++)
    {
        fprintf(stderr, "data->rs_exp_growth_prey = %0.5g\n", data->rs_exp_growth_prey[0][i]);
    }
     */

    // now use interpolation from CUDA;

    // first downcast the double to float
    int NRATE = 1024;
    float* temp =  malloc(NRATE*sizeof(float));
    for (int i = 0; i < NRATE; i++)
    {
        temp[i] = (float) data->r_exp_growth_prey[i];
    }
    cudaBindTexture(NULL, data_d_texture, temp, NRATE*sizeof(float));


    dim3 dimBlock(BLOCK_SIZE,1); dim3 dimGrid(N/BLOCK_SIZE + (N%BLOCK_SIZE == 0 ? 0:1),1);

    // create a working memory in bin_id and Tdef, logTs
    find_binID_Tdef<<<dimGrid,dimBlock>>>(int *bin_id, double *Tdef, double *logTs, double dbin, double idbin, double lb)

}




/********/
/* MAIN */
/********/


/*
int main()
{

    int M=1024;                // --- Number of input points

    int N=1024;                // --- Number of output points

    int Nit = 1000;            // --- Number of computations for time measurement

    // --- Input sampling
    float* x_in; gpuErrchk(cudaMallocManaged(&x_in,sizeof(float)*M));

    // --- Input data
    float *data;        gpuErrchk(cudaMallocManaged(&data,(M+1)*sizeof(float))); Data_Generator(data,M); data[M]=0.;

    // --- Output sampling
    float* x_out;        gpuErrchk(cudaMallocManaged((void**)&x_out,sizeof(float)*N)); randspace(x_out,-M/2.,M/2.,N);

    // --- Result allocation
    float *result_CPU;                            result_CPU=(float*)malloc(N*sizeof(float));
    float *result_d;                            gpuErrchk(cudaMallocManaged(&result_d,sizeof(float)*N));
    float *result_d_texture;                    gpuErrchk(cudaMallocManaged(&result_d_texture,sizeof(float)*N));
    float *result_d_texture_filtering;            gpuErrchk(cudaMallocManaged(&result_d_texture_filtering,sizeof(float)*N));

    // --- Reference interpolation result as evaluated on the CPU
    linear_interpolation_function_CPU(result_CPU, data, x_in, x_out, M, N);

    float time;
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start, 0);
    for (int k=0; k<Nit; k++) linear_interpolation_function_GPU_texture(result_d_texture, data, x_in, x_out, M, N);
    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&time, start, stop);
    cout << "GPU Texture [ms]: " << setprecision (10) << time/Nit << endl;

    float diff_norm_texture=0.f, norm=0.f;
    for(int j=0; j<N; j++) {
        norm      = norm      + result_CPU[j]*result_CPU[j];
        diff_norm_texture = diff_norm_texture + (result_CPU[j]-result_d_texture[j])*(result_CPU[j]-result_d_texture[j]);

    }

    cout << "Error texture [percentage] = " << setprecision(10) << 100.*sqrt(diff_norm_texture/norm) << endl;

    cudaDeviceReset();

    return 0;

}
*/

