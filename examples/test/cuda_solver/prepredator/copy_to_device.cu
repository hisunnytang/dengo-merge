
#include "hdf5.h"
#include "hdf5_hl.h"
// copy the data thing to the device


typedef struct {
    double* r_exp_growth_prey;
    double* rs_exp_growth_prey;
    float* float_rate;

    // temperature-related piece
    double nbins;
    double dbin;
    double idbin;
    double lb;

    double *logTs;
    double *Tdef;
} rateData;

typedef struct 
{
    double nbins;
    double ib;
    double dbin;
    double idbin;

    double *logTs;
    double *Tdef;
    double *dTs_ge;

    // cooling & chemical tables
    double *r_exp_growth_prey;
    double *rs_exp_growth_prey;
    double *drs_exp_growth_prey;
    
    double *r_natural_death_predator;
    double *rs_natural_death_predator;
    double *drs_natural_death_predator;
    
    double *r_predation;
    double *rs_predation;
    double *drs_predation;
} predator_prey_data;

rateData use_unified_memory(void)
{
    
    // read data from hdf5 file
    rateData UnifiedRateData;

    // initialize memory space for rates;
    cudaMallocManaged(&UnifiedRateData.r_exp_growth_prey, 
            sizeof(double)*1024);

    // temperature related pieces
    cudaMallocManaged(&UnifiedRateData.logTs, 
            sizeof(double)*256);
    cudaMallocManaged(&UnifiedRateData.Tdef, 
            sizeof(double)*256);

    // to store the rates
    cudaMallocManaged(&UnifiedRateData.rs_exp_growth_prey, 
            sizeof(double)*256);

    for (int i = 0; i< 256; i++)
    {
        UnifiedRateData.logTs[i] = (double )i/256*1.001;
    }

    UnifiedRateData.nbins = 1023;
    UnifiedRateData.dbin = (log(100000000.0)-log(1.0)) / 1023;
    UnifiedRateData.idbin = 1.0 / UnifiedRateData.dbin;
    UnifiedRateData.lb   = log(1.0);

    // read data from h5
    // can copy it to our cpu array
    const char *filedir;
    filedir = "primordial_tables.h5";
    hid_t file_id = H5Fopen(filedir, H5F_ACC_RDONLY, H5P_DEFAULT);
    H5LTread_dataset_double(file_id, "/k22", UnifiedRateData.r_exp_growth_prey);
    H5Fclose(file_id);

    cudaMallocManaged(&UnifiedRateData.float_rate, 
            sizeof(float)*1024);
    convert_double2float( UnifiedRateData.float_rate, UnifiedRateData.r_exp_growth_prey, 1024);

    return UnifiedRateData;
}

predator_prey_data init_rate_data()
{
    predator_prey_data ratedata;

    ratedata.nbins = 1023;
    ratedata.dbin = (log(100000000.0)-log(1.0)) / 1023;
    ratedata.idbin = 1.0 / ratedata.dbin;
    ratedata.lb   = log(1.0);

    // initialize memory space for rates
    cudaMallocManaged(&ratedata.r_exp_growth_prey,        sizeof(double)*1024);
    cudaMallocManaged(&ratedata.r_predation,              sizeof(double)*1024);
    cudaMallocManaged(&ratedata.r_natural_death_predator, sizeof(double)*1024);

    cudaMallocManaged(&ratedata.rs_exp_growth_prey,       sizeof(double)*256);
    cudaMallocManaged(&ratedata.drs_exp_growth_prey,      sizeof(double)*256);
    cudaMallocManaged(&ratedata.rs_predation,             sizeof(double)*256);
    cudaMallocManaged(&ratedata.drs_predation,            sizeof(double)*256);
    cudaMallocManaged(&ratedata.rs_natural_death_predator,sizeof(double)*256);
    cudaMallocManaged(&ratedata.drs_natural_death_predator,sizeof(double)*256);

    // read rate data from h5
    // can copy it to our cpu array
    const char *filedir = "predator_prey_tables.h5";
    hid_t file_id = H5Fopen(filedir, H5F_ACC_RDONLY, H5P_DEFAULT);
    H5LTread_dataset_double(file_id, "/exp_growth_prey",        ratedata.r_exp_growth_prey);
    H5LTread_dataset_double(file_id, "/predation",              ratedata.r_predation);
    H5LTread_dataset_double(file_id, "/natural_depth_predator", ratedata.r_natural_depth_predator);
    H5Fclose(file_id);


    // temperature related pieces
    cudaMallocManaged(&UnifiedRateData.logTs, sizeof(double)*256);
    cudaMallocManaged(&UnifiedRateData.Tdef,  sizeof(double)*256);

    return ratedata;
}



// use global memory
__global__
void interpolate_kernel(rateData data)
{
    int j = threadIdx.x + blockDim.x* blockIdx.x;

    int k;
    double Tdef, t1;

    double *rates_out = data.rs_exp_growth_prey;
    double *rates_in  = data.r_exp_growth_prey;
    if(j< 256)
    {
        k = __float2int_rz(data.idbin*data.logTs[j] - data.lb);
        t1 = data.lb + k*data.dbin;
        Tdef = (data.logTs[j] - t1) * data.idbin;

        printf("Tdef[%d] = %0.5g; k = %d\n", j, data.logTs[j], k);
        
        rates_out[j] = Tdef*rates_in[k+1] + (-rates_in[k]*Tdef + rates_in[k]);
        printf("rates_out[%d] = %0.5g; @ Tdef = %0.5g\n", j, rates_out[j], Tdef);
    }
}

__global__
void checkIndex(void) {
     printf("threadIdx:(%d, %d, %d) blockIdx:(%d, %d, %d) blockDim:(%d, %d, %d) "
              "gridDim:(%d, %d, %d)\n", threadIdx.x, threadIdx.y, threadIdx.z,
               blockIdx.x, blockIdx.y, blockIdx.z, blockDim.x, blockDim.y, blockDim.z,
                gridDim.x,gridDim.y,gridDim.z);
     
}


texture<float, 1> rate_tex;


__global__ void linear_interpolation_texture(rateData data)
{
    int j = threadIdx.x + blockDim.x * blockIdx.x;
    double *rates_out = data.rs_exp_growth_prey;

    if (j < 256)
    {
        
        int k = __float2int_rz(data.idbin*data.logTs[j] - data.lb);
        float t1 = data.lb + k*data.dbin;
        float Tdef = (data.logTs[j] - t1) * data.idbin;

        float rk   = tex1Dfetch(rate_tex, k);
        float rkp1 = tex1Dfetch(rate_tex, k+1);
        //printf("rk[%d] = %0.5g \n", k, rk);
        rates_out[j] = Tdef*rkp1 + (-rk*Tdef + rk);
        printf("rates_out[%d] = %0.5g; @Tdef = %0.5g \n", j, rates_out[j], Tdef);
    }
}

void copy2texture(const rateData data)
{
    ////////////////////////////////////
    // Copy data to device
    // load data to device memory
    
    int N = 256;
    int BLOCK_SIZE = 8;
    dim3 block(BLOCK_SIZE);
    dim3 grid( (N+block.x -1) / BLOCK_SIZE);

    printf("just before kernel\n");
    // check grid and block dimension from host side
     printf("grid.x %d grid.y %d grid.z %d\n",grid.x, grid.y, grid.z);
    //  printf("block.x %d block.y %d block.z %d\n",block.x, block.y, block.z);
    ///////////////////////////////////
    // load data to device memory
    
    printf("begore float");

    ///////////////////////////////////
    // bind double array to texture
    //////////////////////////////////
    //cudaMalloc()
    //cudaMemcpy()


    
    cudaBindTexture(NULL, rate_tex, data.float_rate, sizeof(float)*1024);
    
    printf("done bining texture\n");
    linear_interpolation_texture<<<grid, block>>>(data);

    cudaDeviceReset();
}


void copy2device(const rateData data)
{
    // load data to device memory
    rateData dA;
    dA.nbins = 1023;
    dA.dbin = (log(100000000.0)-log(1.0)) / 1023;
    dA.idbin = 1.0 / dA.dbin;
    dA.lb   = log(1.0);
    
    // allocate memory in CUDA
    cudaMalloc(&dA.r_exp_growth_prey, 
            sizeof(double)*1024);

    int N = 256;
    int BLOCK_SIZE = 8;

    // to hold interpolated rates
    cudaMalloc(&dA.rs_exp_growth_prey, 
            sizeof(double)*N);
    cudaMalloc(&dA.logTs, 
            sizeof(double)*N);
    cudaMalloc(&dA.Tdef, 
            sizeof(double)*N);


    // copy it to device such that they can be accessed
    // from the kernel
    cudaMemcpy(dA.r_exp_growth_prey, 
            data.r_exp_growth_prey,
            sizeof(double) *N,
            cudaMemcpyHostToDevice);

    // copy the temperature to device
    cudaMemcpy(dA.logTs, 
            data.logTs, 
            sizeof(double) *N,
            cudaMemcpyHostToDevice);
    cudaMemcpy(dA.Tdef, 
            data.Tdef, 
            sizeof(double) *N,
            cudaMemcpyHostToDevice);

    dim3 block(BLOCK_SIZE);
    dim3 grid( (N+block.x -1) / BLOCK_SIZE);

    printf("just before kernel\n");
    // check grid and block dimension from host side
     printf("grid.x %d grid.y %d grid.z %d\n",grid.x, grid.y, grid.z);
    //  printf("block.x %d block.y %d block.z %d\n",block.x, block.y, block.z);

    interpolate_kernel<<<grid, block>>>(dA);

//checkIndex <<<grid, block>>> ();

    cudaDeviceReset();

    // should now copy the device array back to host and print it;


}

void convert_double2float(float *outrate, double *inrate, int N)
{
    for (int i = 0; i< N; i++)
    {
        outrate[i] = (float) inrate[i];    
    }
}

rateData use_unified_memory(void)
{
    
    // read data from hdf5 file
    rateData UnifiedRateData;

    // initialize memory space for rates;
    cudaMallocManaged(&UnifiedRateData.r_exp_growth_prey, 
            sizeof(double)*1024);

    // temperature related pieces
    cudaMallocManaged(&UnifiedRateData.logTs, 
            sizeof(double)*256);
    cudaMallocManaged(&UnifiedRateData.Tdef, 
            sizeof(double)*256);

    // to store the rates
    cudaMallocManaged(&UnifiedRateData.rs_exp_growth_prey, 
            sizeof(double)*256);

    for (int i = 0; i< 256; i++)
    {
        UnifiedRateData.logTs[i] = (double )i/256*1.001;
    }

    UnifiedRateData.nbins = 1023;
    UnifiedRateData.dbin = (log(100000000.0)-log(1.0)) / 1023;
    UnifiedRateData.idbin = 1.0 / UnifiedRateData.dbin;
    UnifiedRateData.lb   = log(1.0);

    // read data from h5
    // can copy it to our cpu array
    const char *filedir;
    filedir = "primordial_tables.h5";
    hid_t file_id = H5Fopen(filedir, H5F_ACC_RDONLY, H5P_DEFAULT);
    H5LTread_dataset_double(file_id, "/k22", UnifiedRateData.r_exp_growth_prey);
    H5Fclose(file_id);

    cudaMallocManaged(&UnifiedRateData.float_rate, 
            sizeof(float)*1024);
    convert_double2float( UnifiedRateData.float_rate, UnifiedRateData.r_exp_growth_prey, 1024);

    return UnifiedRateData;
}


static int rhs(double t, double *y, double *ydot, void* rateData)
{

}

typedef


int main()
{

    // read data from hdf5 file
    rateData hostRateData;

    // initialize memory space for rates;
    hostRateData.r_exp_growth_prey = (double *) malloc(sizeof(double)* 1024);
    hostRateData.logTs = (double *) malloc(sizeof(double)*256);
    for (int i = 0; i< 256; i++)
    {
        hostRateData.logTs[i] = (double )i/256*1.001;
    }

    // read data from h5
    // can copy it to our cpu array
    const char *filedir;
    filedir = "primordial_tables.h5";
    hid_t file_id = H5Fopen(filedir, H5F_ACC_RDONLY, H5P_DEFAULT);
    H5LTread_dataset_double(file_id, "/k22", hostRateData.r_exp_growth_prey);
    H5Fclose(file_id);

    // to verify it's actually there
    for (int i = 0; i< 100; i++)
    {
        fprintf(stderr, "rate[%d] = %0.5g\n", i*10, hostRateData.r_exp_growth_prey[i*10]);
    }

    copy2device(hostRateData);

    rateData uRateData = use_unified_memory();
    copy2texture(uRateData);

}

