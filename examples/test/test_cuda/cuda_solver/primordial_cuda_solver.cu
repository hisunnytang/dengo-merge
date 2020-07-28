#include "primordial_cuda_solver.h"
#include <time.h>

class DengoSolver{

    private:
        primordial_cuda_data rate_data;
        void*    cvode_mem;
        cusparseHandle_t   cusp_handle;
        cusolverSpHandle_t cusol_handle;
        realtype *ydata, *abstol_data;
        N_Vector y, abstol;
        SUNMatrix A;
        SUNLinearSolver LS;
        unsigned long d;
        unsigned long neq;
        int initialize_cvode_solver();

    public:
        DengoSolver()
        {
            neq = BATCHSIZE*nchem;
            int flag = initialize_cvode_solver();
        }

        ~DengoSolver()
        {
            /* Print some final statistics */
            PrintFinalStats(cvode_mem, LS);
            /* Free y and abstol vectors */
            N_VDestroy(y);
            N_VDestroy(abstol);
            /* Free integrator memory */
            CVodeFree(&cvode_mem);
            /* Free the linear solver memory */
            SUNLinSolFree(LS);
            /* Free the matrix memory */
            SUNMatDestroy(A);
            /* Destroy the cuSOLVER and cuSPARSE handles */
            cusparseDestroy(cusp_handle);
            cusolverSpDestroy(cusol_handle);
        };

        int EvolveChemistry(double dtf, double z, double *input, double *rtol, double *atol, unsigned long dims);

};

int DengoSolver::initialize_cvode_solver(){
    rate_data = primordial_cuda_setup_data(NULL, NULL);
    primordial_cuda_read_cooling_tables( &rate_data);
    primordial_cuda_read_rate_tables( &rate_data);

    int retval, iout;
    // create CUDA vector
    y = abstol = NULL;
    A = NULL;
    LS = NULL;
    cvode_mem = NULL;

    neq = BATCHSIZE*nchem;

    // Initialize cuSolver and suSPARSE
    cusparseCreate  (&cusp_handle);
    cusolverSpCreate(&cusol_handle);

    /////////////////////////////////////////////////////////////////
    // Initialize memory spaces for cvode solver //
    /////////////////////////////////////////////////////////////////
    // Create CUDA vector of length neq for I.C. and abstol
    y = N_VNew_Cuda(neq);
    if (check_retval ((void *)y, "N_VNew_Cuda", 0)) return(1);
    abstol = N_VNew_Cuda(neq);
    if (check_retval ((void *)abstol, "N_VNew_Cuda", 0)) return(1);

    ydata       = N_VGetHostArrayPointer_Cuda(y);
    abstol_data = N_VGetHostArrayPointer_Cuda(abstol);
    // Initialize cuda vector
    for (unsigned long i = 0; i < neq; i++){
        ydata[i] = 1.0;
        abstol_data[i] = 1.0e-6;
    }
    N_VCopyToDevice_Cuda(y);
    N_VCopyToDevice_Cuda(abstol);

    /* Call CVodeCreate to create the solver memory and specify the
     * Backward Differentiation Formula */
    cvode_mem = CVodeCreate(CV_BDF);
    if (check_retval((void *)cvode_mem, "CVodeCreate", 0)) return(1);

    /* Call CVodeInit to initialize the integrator memory and specify the
     * user's right hand side function in y'=f(t,y), the inital time T0, and
     * the initial dependent variable vector y. */
    retval = CVodeInit(cvode_mem, calculate_rhs_primordial_cuda, 0.0, y);
    if (check_retval(&retval, "CVodeInit", 1)) return(1);

    /* Call CVodeSetUserData to attach the user data structure */
    retval = CVodeSetUserData(cvode_mem, &rate_data);
    if (check_retval(&retval, "CVodeSetUserData", 1)) return(1);

    /* Call CVodeSVtolerances to specify the scalar relative tolerance
     * and vector absolute tolerances */
    retval = CVodeSVtolerances(cvode_mem, 1.0e-5, abstol);
    if (check_retval(&retval, "CVodeSVtolerances", 1)) return(1);

#ifndef USESPARSE
    /* Create sparse SUNMatrix for use in linear solves */
    A = SUNMatrix_cuSparse_NewBlockCSR(BATCHSIZE, nchem, nchem, nchem*nchem, cusp_handle);
    if(check_retval((void *)A, "SUNMatrix_cuSparse_NewBlockCSR", 0)) return(1);

    /* Set the sparsity pattern to be fixed so that the row pointers
     * and column indicies are not zeroed out by SUNMatZero */
    retval = SUNMatrix_cuSparse_SetFixedPattern(A, 1);

    /* Initialiize the Jacobian with its fixed sparsity pattern */
    blockJacInit(A);
#else
    int nnz = 64;
    /* Create sparse SUNMatrix for use in linear solves */
    A = SUNMatrix_cuSparse_NewBlockCSR(BATCHSIZE, nchem, nchem, nnz, cusp_handle);
    if(check_retval((void *)A, "SUNMatrix_cuSparse_NewBlockCSR", 0)) return(1);
    /* Set the sparsity pattern to be fixed so that the row pointers
     * and column indicies are not zeroed out by SUNMatZero */
    retval = SUNMatrix_cuSparse_SetFixedPattern(A, 1);
    blockSparseJacInit(A);
#endif

    /* Create the SUNLinearSolver object for use by CVode */
    LS = SUNLinSol_cuSolverSp_batchQR(y, A, cusol_handle);
    if(check_retval((void *)LS, "SUNLinSol_cuSolverSp_batchQR", 0)) return(1);

    /* Call CVodeSetLinearSolver to attach the matrix and linear solver to CVode */
    retval = CVodeSetLinearSolver(cvode_mem, LS, A);
    if(check_retval(&retval, "CVodeSetLinearSolver", 1)) return(1);

    /* Set the user-supplied Jacobian routine Jac */

#ifndef USESPARSE
    retval = CVodeSetJacFn(cvode_mem, calculate_jacobian_primordial_cuda);
    if(check_retval(&retval, "CVodeSetJacFn", 1)) return(1);
#else
    retval = CVodeSetJacFn(cvode_mem, calculate_sparse_jacobian_primordial_cuda);
    if(check_retval(&retval, "CVodeSetJacFn", 1)) return(1);
#endif

    return 0;
}


int DengoSolver::EvolveChemistry(double dtf, double z, double *input, double *rtol, double *atol, unsigned long dims)
{
    unsigned long ntimes = dims/BATCHSIZE;
    unsigned long eqlast = dims % BATCHSIZE;
    unsigned long count  = 0;
    int retval;
    realtype t;
    realtype reltol = rtol[0];

    // split input by batchsize
    for (count = 0; count < ntimes; count++)
    {
        for (int i = 0; i < neq; i++)
        {
            ydata[i]       = input[count*neq+i];
            abstol_data[i] = reltol*reltol*input[count*neq+i];
        }
        //for (int i = 0; i < nchem; i++) printf("ydata[%d] = %0.5g\n", i, ydata[i]);

        N_VCopyToDevice_Cuda(y);
        N_VCopyToDevice_Cuda(abstol);

        retval = CVodeReInit(cvode_mem, 0.0, y);
        retval = CVode(cvode_mem, dtf, y, &t, CV_NORMAL);
        N_VCopyFromDevice_Cuda(y);

        // copy cvode output from cuda vector
        for (int i = 0; i < neq; i++) input[count*neq+i] = ydata[i];

        //for (int i = 0; i < nchem; i++) printf("outydata[%d] = %0.5g\n", i, ydata[i]);
    }
    return 0;

}


// Initialize a data object that stores the reaction/ cooling rate data
primordial_cuda_data primordial_cuda_setup_data(int *NumberOfFields, char ***FieldNames)
{

    //-----------------------------------------------------
    // Function : primordial_cuda_setup_data
    // Description: Initialize a data object that stores the reaction/ cooling rate data 
    //-----------------------------------------------------


    // let's not re-scale the data yet...
    primordial_cuda_data ratedata;

    ratedata.nbins = 1024;
    ratedata.dbin = (log( 100000000.0)-log(1.0)) / 1023;
    ratedata.idbin = 1.0 / ratedata.dbin;
    ratedata.lb   = log(1.0);
    ratedata.ub   = log(100000000.0);

    /* Redshift-related pieces */
    /*
       data->z_bounds[0] = 0.0;
       data->z_bounds[1] = 10.0;
       data->n_zbins = 0 - 1;
       data->d_zbin = (log(data->z_bounds[1] + 1.0) - log(data->z_bounds[0] + 1.0)) / data->n_zbins;
       data->id_zbin = 1.0L / data->d_zbin;
     */


    // initialize memory space for reaction rates and cooling rates
    // we use managed data, so the pointer can simultaneously be accessed from device and the host
    cudaMallocManaged(&ratedata.r_k01, sizeof(double)*1024);
    cudaMallocManaged(&ratedata.rs_k01, sizeof(double)*1024);
    cudaMallocManaged(&ratedata.drs_k01, sizeof(double)*1024);
    cudaMallocManaged(&ratedata.r_k02, sizeof(double)*1024);
    cudaMallocManaged(&ratedata.rs_k02, sizeof(double)*1024);
    cudaMallocManaged(&ratedata.drs_k02, sizeof(double)*1024);
    cudaMallocManaged(&ratedata.r_k03, sizeof(double)*1024);
    cudaMallocManaged(&ratedata.rs_k03, sizeof(double)*1024);
    cudaMallocManaged(&ratedata.drs_k03, sizeof(double)*1024);
    cudaMallocManaged(&ratedata.r_k04, sizeof(double)*1024);
    cudaMallocManaged(&ratedata.rs_k04, sizeof(double)*1024);
    cudaMallocManaged(&ratedata.drs_k04, sizeof(double)*1024);
    cudaMallocManaged(&ratedata.r_k05, sizeof(double)*1024);
    cudaMallocManaged(&ratedata.rs_k05, sizeof(double)*1024);
    cudaMallocManaged(&ratedata.drs_k05, sizeof(double)*1024);
    cudaMallocManaged(&ratedata.r_k06, sizeof(double)*1024);
    cudaMallocManaged(&ratedata.rs_k06, sizeof(double)*1024);
    cudaMallocManaged(&ratedata.drs_k06, sizeof(double)*1024);
    cudaMallocManaged(&ratedata.r_k07, sizeof(double)*1024);
    cudaMallocManaged(&ratedata.rs_k07, sizeof(double)*1024);
    cudaMallocManaged(&ratedata.drs_k07, sizeof(double)*1024);
    cudaMallocManaged(&ratedata.r_k08, sizeof(double)*1024);
    cudaMallocManaged(&ratedata.rs_k08, sizeof(double)*1024);
    cudaMallocManaged(&ratedata.drs_k08, sizeof(double)*1024);
    cudaMallocManaged(&ratedata.r_k09, sizeof(double)*1024);
    cudaMallocManaged(&ratedata.rs_k09, sizeof(double)*1024);
    cudaMallocManaged(&ratedata.drs_k09, sizeof(double)*1024);
    cudaMallocManaged(&ratedata.r_k10, sizeof(double)*1024);
    cudaMallocManaged(&ratedata.rs_k10, sizeof(double)*1024);
    cudaMallocManaged(&ratedata.drs_k10, sizeof(double)*1024);
    cudaMallocManaged(&ratedata.r_k11, sizeof(double)*1024);
    cudaMallocManaged(&ratedata.rs_k11, sizeof(double)*1024);
    cudaMallocManaged(&ratedata.drs_k11, sizeof(double)*1024);
    cudaMallocManaged(&ratedata.r_k12, sizeof(double)*1024);
    cudaMallocManaged(&ratedata.rs_k12, sizeof(double)*1024);
    cudaMallocManaged(&ratedata.drs_k12, sizeof(double)*1024);
    cudaMallocManaged(&ratedata.r_k13, sizeof(double)*1024);
    cudaMallocManaged(&ratedata.rs_k13, sizeof(double)*1024);
    cudaMallocManaged(&ratedata.drs_k13, sizeof(double)*1024);
    cudaMallocManaged(&ratedata.r_k14, sizeof(double)*1024);
    cudaMallocManaged(&ratedata.rs_k14, sizeof(double)*1024);
    cudaMallocManaged(&ratedata.drs_k14, sizeof(double)*1024);
    cudaMallocManaged(&ratedata.r_k15, sizeof(double)*1024);
    cudaMallocManaged(&ratedata.rs_k15, sizeof(double)*1024);
    cudaMallocManaged(&ratedata.drs_k15, sizeof(double)*1024);
    cudaMallocManaged(&ratedata.r_k16, sizeof(double)*1024);
    cudaMallocManaged(&ratedata.rs_k16, sizeof(double)*1024);
    cudaMallocManaged(&ratedata.drs_k16, sizeof(double)*1024);
    cudaMallocManaged(&ratedata.r_k17, sizeof(double)*1024);
    cudaMallocManaged(&ratedata.rs_k17, sizeof(double)*1024);
    cudaMallocManaged(&ratedata.drs_k17, sizeof(double)*1024);
    cudaMallocManaged(&ratedata.r_k18, sizeof(double)*1024);
    cudaMallocManaged(&ratedata.rs_k18, sizeof(double)*1024);
    cudaMallocManaged(&ratedata.drs_k18, sizeof(double)*1024);
    cudaMallocManaged(&ratedata.r_k19, sizeof(double)*1024);
    cudaMallocManaged(&ratedata.rs_k19, sizeof(double)*1024);
    cudaMallocManaged(&ratedata.drs_k19, sizeof(double)*1024);
    cudaMallocManaged(&ratedata.r_k21, sizeof(double)*1024);
    cudaMallocManaged(&ratedata.rs_k21, sizeof(double)*1024);
    cudaMallocManaged(&ratedata.drs_k21, sizeof(double)*1024);
    cudaMallocManaged(&ratedata.r_k22, sizeof(double)*1024);
    cudaMallocManaged(&ratedata.rs_k22, sizeof(double)*1024);
    cudaMallocManaged(&ratedata.drs_k22, sizeof(double)*1024);
    cudaMallocManaged(&ratedata.r_k23, sizeof(double)*1024);
    cudaMallocManaged(&ratedata.rs_k23, sizeof(double)*1024);
    cudaMallocManaged(&ratedata.drs_k23, sizeof(double)*1024);

    // Cooling Rates
    cudaMallocManaged(&ratedata.c_brem_brem, sizeof(double)* 1024 );
    cudaMallocManaged(&ratedata.cs_brem_brem, sizeof(double)* 1024 );
    cudaMallocManaged(&ratedata.dcs_brem_brem, sizeof(double)* 1024 );
    cudaMallocManaged(&ratedata.c_ceHeI_ceHeI, sizeof(double)* 1024 );
    cudaMallocManaged(&ratedata.cs_ceHeI_ceHeI, sizeof(double)* 1024 );
    cudaMallocManaged(&ratedata.dcs_ceHeI_ceHeI, sizeof(double)* 1024 );
    cudaMallocManaged(&ratedata.c_ceHeII_ceHeII, sizeof(double)* 1024 );
    cudaMallocManaged(&ratedata.cs_ceHeII_ceHeII, sizeof(double)* 1024 );
    cudaMallocManaged(&ratedata.dcs_ceHeII_ceHeII, sizeof(double)* 1024 );
    cudaMallocManaged(&ratedata.c_ceHI_ceHI, sizeof(double)* 1024 );
    cudaMallocManaged(&ratedata.cs_ceHI_ceHI, sizeof(double)* 1024 );
    cudaMallocManaged(&ratedata.dcs_ceHI_ceHI, sizeof(double)* 1024 );
    cudaMallocManaged(&ratedata.c_cie_cooling_cieco, sizeof(double)* 1024 );
    cudaMallocManaged(&ratedata.cs_cie_cooling_cieco, sizeof(double)* 1024 );
    cudaMallocManaged(&ratedata.dcs_cie_cooling_cieco, sizeof(double)* 1024 );
    cudaMallocManaged(&ratedata.c_ciHeI_ciHeI, sizeof(double)* 1024 );
    cudaMallocManaged(&ratedata.cs_ciHeI_ciHeI, sizeof(double)* 1024 );
    cudaMallocManaged(&ratedata.dcs_ciHeI_ciHeI, sizeof(double)* 1024 );
    cudaMallocManaged(&ratedata.c_ciHeII_ciHeII, sizeof(double)* 1024 );
    cudaMallocManaged(&ratedata.cs_ciHeII_ciHeII, sizeof(double)* 1024 );
    cudaMallocManaged(&ratedata.dcs_ciHeII_ciHeII, sizeof(double)* 1024 );
    cudaMallocManaged(&ratedata.c_ciHeIS_ciHeIS, sizeof(double)* 1024 );
    cudaMallocManaged(&ratedata.cs_ciHeIS_ciHeIS, sizeof(double)* 1024 );
    cudaMallocManaged(&ratedata.dcs_ciHeIS_ciHeIS, sizeof(double)* 1024 );
    cudaMallocManaged(&ratedata.c_ciHI_ciHI, sizeof(double)* 1024 );
    cudaMallocManaged(&ratedata.cs_ciHI_ciHI, sizeof(double)* 1024 );
    cudaMallocManaged(&ratedata.dcs_ciHI_ciHI, sizeof(double)* 1024 );
    cudaMallocManaged(&ratedata.c_compton_comp_, sizeof(double)* 1024 );
    cudaMallocManaged(&ratedata.cs_compton_comp_, sizeof(double)* 1024 );
    cudaMallocManaged(&ratedata.dcs_compton_comp_, sizeof(double)* 1024 );
    cudaMallocManaged(&ratedata.c_gammah_gammah, sizeof(double)* 1024 );
    cudaMallocManaged(&ratedata.cs_gammah_gammah, sizeof(double)* 1024 );
    cudaMallocManaged(&ratedata.dcs_gammah_gammah, sizeof(double)* 1024 );
    cudaMallocManaged(&ratedata.c_gloverabel08_gael, sizeof(double)* 1024 );
    cudaMallocManaged(&ratedata.cs_gloverabel08_gael, sizeof(double)* 1024 );
    cudaMallocManaged(&ratedata.dcs_gloverabel08_gael, sizeof(double)* 1024 );
    cudaMallocManaged(&ratedata.c_gloverabel08_gaH2, sizeof(double)* 1024 );
    cudaMallocManaged(&ratedata.cs_gloverabel08_gaH2, sizeof(double)* 1024 );
    cudaMallocManaged(&ratedata.dcs_gloverabel08_gaH2, sizeof(double)* 1024 );
    cudaMallocManaged(&ratedata.c_gloverabel08_gaHe, sizeof(double)* 1024 );
    cudaMallocManaged(&ratedata.cs_gloverabel08_gaHe, sizeof(double)* 1024 );
    cudaMallocManaged(&ratedata.dcs_gloverabel08_gaHe, sizeof(double)* 1024 );
    cudaMallocManaged(&ratedata.c_gloverabel08_gaHI, sizeof(double)* 1024 );
    cudaMallocManaged(&ratedata.cs_gloverabel08_gaHI, sizeof(double)* 1024 );
    cudaMallocManaged(&ratedata.dcs_gloverabel08_gaHI, sizeof(double)* 1024 );
    cudaMallocManaged(&ratedata.c_gloverabel08_gaHp, sizeof(double)* 1024 );
    cudaMallocManaged(&ratedata.cs_gloverabel08_gaHp, sizeof(double)* 1024 );
    cudaMallocManaged(&ratedata.dcs_gloverabel08_gaHp, sizeof(double)* 1024 );
    cudaMallocManaged(&ratedata.c_gloverabel08_gphdl, sizeof(double)* 1024 );
    cudaMallocManaged(&ratedata.cs_gloverabel08_gphdl, sizeof(double)* 1024 );
    cudaMallocManaged(&ratedata.dcs_gloverabel08_gphdl, sizeof(double)* 1024 );
    cudaMallocManaged(&ratedata.c_gloverabel08_gpldl, sizeof(double)* 1024 );
    cudaMallocManaged(&ratedata.cs_gloverabel08_gpldl, sizeof(double)* 1024 );
    cudaMallocManaged(&ratedata.dcs_gloverabel08_gpldl, sizeof(double)* 1024 );
    cudaMallocManaged(&ratedata.c_gloverabel08_h2lte, sizeof(double)* 1024 );
    cudaMallocManaged(&ratedata.cs_gloverabel08_h2lte, sizeof(double)* 1024 );
    cudaMallocManaged(&ratedata.dcs_gloverabel08_h2lte, sizeof(double)* 1024 );
    cudaMallocManaged(&ratedata.c_h2formation_h2mcool, sizeof(double)* 1024 );
    cudaMallocManaged(&ratedata.cs_h2formation_h2mcool, sizeof(double)* 1024 );
    cudaMallocManaged(&ratedata.dcs_h2formation_h2mcool, sizeof(double)* 1024 );
    cudaMallocManaged(&ratedata.c_h2formation_h2mheat, sizeof(double)* 1024 );
    cudaMallocManaged(&ratedata.cs_h2formation_h2mheat, sizeof(double)* 1024 );
    cudaMallocManaged(&ratedata.dcs_h2formation_h2mheat, sizeof(double)* 1024 );
    cudaMallocManaged(&ratedata.c_h2formation_ncrd1, sizeof(double)* 1024 );
    cudaMallocManaged(&ratedata.cs_h2formation_ncrd1, sizeof(double)* 1024 );
    cudaMallocManaged(&ratedata.dcs_h2formation_ncrd1, sizeof(double)* 1024 );
    cudaMallocManaged(&ratedata.c_h2formation_ncrd2, sizeof(double)* 1024 );
    cudaMallocManaged(&ratedata.cs_h2formation_ncrd2, sizeof(double)* 1024 );
    cudaMallocManaged(&ratedata.dcs_h2formation_ncrd2, sizeof(double)* 1024 );
    cudaMallocManaged(&ratedata.c_h2formation_ncrn, sizeof(double)* 1024 );
    cudaMallocManaged(&ratedata.cs_h2formation_ncrn, sizeof(double)* 1024 );
    cudaMallocManaged(&ratedata.dcs_h2formation_ncrn, sizeof(double)* 1024 );
    cudaMallocManaged(&ratedata.c_reHeII1_reHeII1, sizeof(double)* 1024 );
    cudaMallocManaged(&ratedata.cs_reHeII1_reHeII1, sizeof(double)* 1024 );
    cudaMallocManaged(&ratedata.dcs_reHeII1_reHeII1, sizeof(double)* 1024 );
    cudaMallocManaged(&ratedata.c_reHeII2_reHeII2, sizeof(double)* 1024 );
    cudaMallocManaged(&ratedata.cs_reHeII2_reHeII2, sizeof(double)* 1024 );
    cudaMallocManaged(&ratedata.dcs_reHeII2_reHeII2, sizeof(double)* 1024 );
    cudaMallocManaged(&ratedata.c_reHeIII_reHeIII, sizeof(double)* 1024 );
    cudaMallocManaged(&ratedata.cs_reHeIII_reHeIII, sizeof(double)* 1024 );
    cudaMallocManaged(&ratedata.dcs_reHeIII_reHeIII, sizeof(double)* 1024 );
    cudaMallocManaged(&ratedata.c_reHII_reHII, sizeof(double)* 1024 );
    cudaMallocManaged(&ratedata.cs_reHII_reHII, sizeof(double)* 1024 );
    cudaMallocManaged(&ratedata.dcs_reHII_reHII, sizeof(double)* 1024 );

    // initialize memory space for the temperature-related pieces
    cudaMallocManaged(&ratedata.Ts, sizeof(double)* BATCHSIZE);
    cudaMallocManaged(&ratedata.logTs, sizeof(double)* BATCHSIZE);
    cudaMallocManaged(&ratedata.Tdef,  sizeof(double)* BATCHSIZE);
    cudaMallocManaged(&ratedata.dTs_ge,  sizeof(double)* BATCHSIZE);
    cudaMallocManaged(&ratedata.Tge,  sizeof(double)* BATCHSIZE);

    // gamma as a function of temperature
    /*
       cudaMallocManaged(&ratedata.g_gammaH2_1, sizeof(double)* BATCHSIZE);
       cudaMallocManaged(&ratedata.g_dgammaH2_1_dT,  sizeof(double)* BATCHSIZE);
       cudaMallocManaged(&ratedata.gammaH2_1,  sizeof(double)* BATCHSIZE);
       cudaMallocManaged(&ratedata.dgamma_dTH2_1,  sizeof(double)* BATCHSIZE);
       cudaMallocManaged(&ratedata._gammaH2_1_dT, sizeof(double)*BATCHSIZE);
       cudaMallocManaged(&ratedata.g_gammaH2_2, sizeof(double)* BATCHSIZE);
       cudaMallocManaged(&ratedata.g_dgammaH2_2_dT,  sizeof(double)* BATCHSIZE);
       cudaMallocManaged(&ratedata.gammaH2_2,  sizeof(double)* BATCHSIZE);
       cudaMallocManaged(&ratedata.dgamma_dTH2_2,  sizeof(double)* BATCHSIZE);
       cudaMallocManaged(&ratedata._gammaH2_2_dT, sizeof(double)*BATCHSIZE);

    // maybe we can calculate the density on the fly
    // space to store the mass density
    cudaMallocManaged(&ratedata.mdensity, sizeof(double)* BATCHSIZE);
    cudaMallocManaged(&ratedata.inv_mdensity, sizeof(double)* BATCHSIZE);
     */
    // extra stuff like the density-dependent cooling rate
    cudaMallocManaged(&ratedata.cie_optical_depth_approx, sizeof(double)* BATCHSIZE);
    cudaMallocManaged(&ratedata.h2_optical_depth_approx, sizeof(double)* BATCHSIZE);

    return ratedata;
}


void primordial_cuda_read_rate_tables(primordial_cuda_data *data)
{
    hid_t file_id = H5Fopen("primordial_cuda_tables.h5", H5F_ACC_RDONLY, H5P_DEFAULT);
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


void primordial_cuda_read_cooling_tables(primordial_cuda_data *data)
{

    hid_t file_id = H5Fopen("primordial_cuda_tables.h5", H5F_ACC_RDONLY, H5P_DEFAULT);
    /* Allocate the correct number of rate tables */
    H5LTread_dataset_double(file_id, "/brem_brem",
            data->c_brem_brem);
    H5LTread_dataset_double(file_id, "/ceHeI_ceHeI",
            data->c_ceHeI_ceHeI);
    H5LTread_dataset_double(file_id, "/ceHeII_ceHeII",
            data->c_ceHeII_ceHeII);
    H5LTread_dataset_double(file_id, "/ceHI_ceHI",
            data->c_ceHI_ceHI);
    H5LTread_dataset_double(file_id, "/cie_cooling_cieco",
            data->c_cie_cooling_cieco);
    H5LTread_dataset_double(file_id, "/ciHeI_ciHeI",
            data->c_ciHeI_ciHeI);
    H5LTread_dataset_double(file_id, "/ciHeII_ciHeII",
            data->c_ciHeII_ciHeII);
    H5LTread_dataset_double(file_id, "/ciHeIS_ciHeIS",
            data->c_ciHeIS_ciHeIS);
    H5LTread_dataset_double(file_id, "/ciHI_ciHI",
            data->c_ciHI_ciHI);
    H5LTread_dataset_double(file_id, "/compton_comp_",
            data->c_compton_comp_);
    H5LTread_dataset_double(file_id, "/gammah_gammah",
            data->c_gammah_gammah);
    H5LTread_dataset_double(file_id, "/gloverabel08_gael",
            data->c_gloverabel08_gael);
    H5LTread_dataset_double(file_id, "/gloverabel08_gaH2",
            data->c_gloverabel08_gaH2);
    H5LTread_dataset_double(file_id, "/gloverabel08_gaHe",
            data->c_gloverabel08_gaHe);
    H5LTread_dataset_double(file_id, "/gloverabel08_gaHI",
            data->c_gloverabel08_gaHI);
    H5LTread_dataset_double(file_id, "/gloverabel08_gaHp",
            data->c_gloverabel08_gaHp);
    H5LTread_dataset_double(file_id, "/gloverabel08_gphdl",
            data->c_gloverabel08_gphdl);
    H5LTread_dataset_double(file_id, "/gloverabel08_gpldl",
            data->c_gloverabel08_gpldl);
    H5LTread_dataset_double(file_id, "/gloverabel08_h2lte",
            data->c_gloverabel08_h2lte);
    H5LTread_dataset_double(file_id, "/h2formation_h2mcool",
            data->c_h2formation_h2mcool);
    H5LTread_dataset_double(file_id, "/h2formation_h2mheat",
            data->c_h2formation_h2mheat);
    H5LTread_dataset_double(file_id, "/h2formation_ncrd1",
            data->c_h2formation_ncrd1);
    H5LTread_dataset_double(file_id, "/h2formation_ncrd2",
            data->c_h2formation_ncrd2);
    H5LTread_dataset_double(file_id, "/h2formation_ncrn",
            data->c_h2formation_ncrn);
    H5LTread_dataset_double(file_id, "/reHeII1_reHeII1",
            data->c_reHeII1_reHeII1);
    H5LTread_dataset_double(file_id, "/reHeII2_reHeII2",
            data->c_reHeII2_reHeII2);
    H5LTread_dataset_double(file_id, "/reHeIII_reHeIII",
            data->c_reHeIII_reHeIII);
    H5LTread_dataset_double(file_id, "/reHII_reHII",
            data->c_reHII_reHII);

    H5Fclose(file_id);
}

/*
   void primordial_cuda_read_gamma(primordial_cuda_data *data)
   {

   hid_t file_id = H5Fopen("primordial_cuda_tables.h5", H5F_ACC_RDONLY, H5P_DEFAULT);
// Allocate the correct number of rate tables
H5LTread_dataset_double(file_id, "/gammaH2_1",
data->g_gammaH2_1 );
H5LTread_dataset_double(file_id, "/dgammaH2_1_dT",
data->g_dgammaH2_1_dT );   

H5LTread_dataset_double(file_id, "/gammaH2_2",
data->g_gammaH2_2 );
H5LTread_dataset_double(file_id, "/dgammaH2_2_dT",
data->g_dgammaH2_2_dT );   


H5Fclose(file_id);

}
 
 */


// interpolation kernel
// ideally, we should use texture to do interpolation,
// let's ignore it for now, cos i guess most time is spent in doing the matrix thingy

    __global__
void linear_interpolation_kernel(primordial_cuda_data data)
{
    int j = threadIdx.x + blockDim.x* blockIdx.x;

    int k;
    double Tdef, t1;
    double *k01 = data.r_k01;
    double *rs_k01  = data.rs_k01;
    double *drs_k01 = data.drs_k01;
    double *k02 = data.r_k02;
    double *rs_k02  = data.rs_k02;
    double *drs_k02 = data.drs_k02;
    double *k03 = data.r_k03;
    double *rs_k03  = data.rs_k03;
    double *drs_k03 = data.drs_k03;
    double *k04 = data.r_k04;
    double *rs_k04  = data.rs_k04;
    double *drs_k04 = data.drs_k04;
    double *k05 = data.r_k05;
    double *rs_k05  = data.rs_k05;
    double *drs_k05 = data.drs_k05;
    double *k06 = data.r_k06;
    double *rs_k06  = data.rs_k06;
    double *drs_k06 = data.drs_k06;
    double *k07 = data.r_k07;
    double *rs_k07  = data.rs_k07;
    double *drs_k07 = data.drs_k07;
    double *k08 = data.r_k08;
    double *rs_k08  = data.rs_k08;
    double *drs_k08 = data.drs_k08;
    double *k09 = data.r_k09;
    double *rs_k09  = data.rs_k09;
    double *drs_k09 = data.drs_k09;
    double *k10 = data.r_k10;
    double *rs_k10  = data.rs_k10;
    double *drs_k10 = data.drs_k10;
    double *k11 = data.r_k11;
    double *rs_k11  = data.rs_k11;
    double *drs_k11 = data.drs_k11;
    double *k12 = data.r_k12;
    double *rs_k12  = data.rs_k12;
    double *drs_k12 = data.drs_k12;
    double *k13 = data.r_k13;
    double *rs_k13  = data.rs_k13;
    double *drs_k13 = data.drs_k13;
    double *k14 = data.r_k14;
    double *rs_k14  = data.rs_k14;
    double *drs_k14 = data.drs_k14;
    double *k15 = data.r_k15;
    double *rs_k15  = data.rs_k15;
    double *drs_k15 = data.drs_k15;
    double *k16 = data.r_k16;
    double *rs_k16  = data.rs_k16;
    double *drs_k16 = data.drs_k16;
    double *k17 = data.r_k17;
    double *rs_k17  = data.rs_k17;
    double *drs_k17 = data.drs_k17;
    double *k18 = data.r_k18;
    double *rs_k18  = data.rs_k18;
    double *drs_k18 = data.drs_k18;
    double *k19 = data.r_k19;
    double *rs_k19  = data.rs_k19;
    double *drs_k19 = data.drs_k19;
    double *k21 = data.r_k21;
    double *rs_k21  = data.rs_k21;
    double *drs_k21 = data.drs_k21;
    double *k22 = data.r_k22;
    double *rs_k22  = data.rs_k22;
    double *drs_k22 = data.drs_k22;
    double *k23 = data.r_k23;
    double *rs_k23  = data.rs_k23;
    double *drs_k23 = data.drs_k23;
    double *brem_brem = data.c_brem_brem;
    double *cs_brem_brem  = data.cs_brem_brem;
    double *dcs_brem_brem = data.dcs_brem_brem;
    double *ceHeI_ceHeI = data.c_ceHeI_ceHeI;
    double *cs_ceHeI_ceHeI  = data.cs_ceHeI_ceHeI;
    double *dcs_ceHeI_ceHeI = data.dcs_ceHeI_ceHeI;
    double *ceHeII_ceHeII = data.c_ceHeII_ceHeII;
    double *cs_ceHeII_ceHeII  = data.cs_ceHeII_ceHeII;
    double *dcs_ceHeII_ceHeII = data.dcs_ceHeII_ceHeII;
    double *ceHI_ceHI = data.c_ceHI_ceHI;
    double *cs_ceHI_ceHI  = data.cs_ceHI_ceHI;
    double *dcs_ceHI_ceHI = data.dcs_ceHI_ceHI;
    double *cie_cooling_cieco = data.c_cie_cooling_cieco;
    double *cs_cie_cooling_cieco  = data.cs_cie_cooling_cieco;
    double *dcs_cie_cooling_cieco = data.dcs_cie_cooling_cieco;
    double *ciHeI_ciHeI = data.c_ciHeI_ciHeI;
    double *cs_ciHeI_ciHeI  = data.cs_ciHeI_ciHeI;
    double *dcs_ciHeI_ciHeI = data.dcs_ciHeI_ciHeI;
    double *ciHeII_ciHeII = data.c_ciHeII_ciHeII;
    double *cs_ciHeII_ciHeII  = data.cs_ciHeII_ciHeII;
    double *dcs_ciHeII_ciHeII = data.dcs_ciHeII_ciHeII;
    double *ciHeIS_ciHeIS = data.c_ciHeIS_ciHeIS;
    double *cs_ciHeIS_ciHeIS  = data.cs_ciHeIS_ciHeIS;
    double *dcs_ciHeIS_ciHeIS = data.dcs_ciHeIS_ciHeIS;
    double *ciHI_ciHI = data.c_ciHI_ciHI;
    double *cs_ciHI_ciHI  = data.cs_ciHI_ciHI;
    double *dcs_ciHI_ciHI = data.dcs_ciHI_ciHI;
    double *compton_comp_ = data.c_compton_comp_;
    double *cs_compton_comp_  = data.cs_compton_comp_;
    double *dcs_compton_comp_ = data.dcs_compton_comp_;
    double *gammah_gammah = data.c_gammah_gammah;
    double *cs_gammah_gammah  = data.cs_gammah_gammah;
    double *dcs_gammah_gammah = data.dcs_gammah_gammah;
    double *gloverabel08_gael = data.c_gloverabel08_gael;
    double *cs_gloverabel08_gael  = data.cs_gloverabel08_gael;
    double *dcs_gloverabel08_gael = data.dcs_gloverabel08_gael;
    double *gloverabel08_gaH2 = data.c_gloverabel08_gaH2;
    double *cs_gloverabel08_gaH2  = data.cs_gloverabel08_gaH2;
    double *dcs_gloverabel08_gaH2 = data.dcs_gloverabel08_gaH2;
    double *gloverabel08_gaHe = data.c_gloverabel08_gaHe;
    double *cs_gloverabel08_gaHe  = data.cs_gloverabel08_gaHe;
    double *dcs_gloverabel08_gaHe = data.dcs_gloverabel08_gaHe;
    double *gloverabel08_gaHI = data.c_gloverabel08_gaHI;
    double *cs_gloverabel08_gaHI  = data.cs_gloverabel08_gaHI;
    double *dcs_gloverabel08_gaHI = data.dcs_gloverabel08_gaHI;
    double *gloverabel08_gaHp = data.c_gloverabel08_gaHp;
    double *cs_gloverabel08_gaHp  = data.cs_gloverabel08_gaHp;
    double *dcs_gloverabel08_gaHp = data.dcs_gloverabel08_gaHp;
    double *gloverabel08_gphdl = data.c_gloverabel08_gphdl;
    double *cs_gloverabel08_gphdl  = data.cs_gloverabel08_gphdl;
    double *dcs_gloverabel08_gphdl = data.dcs_gloverabel08_gphdl;
    double *gloverabel08_gpldl = data.c_gloverabel08_gpldl;
    double *cs_gloverabel08_gpldl  = data.cs_gloverabel08_gpldl;
    double *dcs_gloverabel08_gpldl = data.dcs_gloverabel08_gpldl;
    double *gloverabel08_h2lte = data.c_gloverabel08_h2lte;
    double *cs_gloverabel08_h2lte  = data.cs_gloverabel08_h2lte;
    double *dcs_gloverabel08_h2lte = data.dcs_gloverabel08_h2lte;
    double *h2formation_h2mcool = data.c_h2formation_h2mcool;
    double *cs_h2formation_h2mcool  = data.cs_h2formation_h2mcool;
    double *dcs_h2formation_h2mcool = data.dcs_h2formation_h2mcool;
    double *h2formation_h2mheat = data.c_h2formation_h2mheat;
    double *cs_h2formation_h2mheat  = data.cs_h2formation_h2mheat;
    double *dcs_h2formation_h2mheat = data.dcs_h2formation_h2mheat;
    double *h2formation_ncrd1 = data.c_h2formation_ncrd1;
    double *cs_h2formation_ncrd1  = data.cs_h2formation_ncrd1;
    double *dcs_h2formation_ncrd1 = data.dcs_h2formation_ncrd1;
    double *h2formation_ncrd2 = data.c_h2formation_ncrd2;
    double *cs_h2formation_ncrd2  = data.cs_h2formation_ncrd2;
    double *dcs_h2formation_ncrd2 = data.dcs_h2formation_ncrd2;
    double *h2formation_ncrn = data.c_h2formation_ncrn;
    double *cs_h2formation_ncrn  = data.cs_h2formation_ncrn;
    double *dcs_h2formation_ncrn = data.dcs_h2formation_ncrn;
    double *reHeII1_reHeII1 = data.c_reHeII1_reHeII1;
    double *cs_reHeII1_reHeII1  = data.cs_reHeII1_reHeII1;
    double *dcs_reHeII1_reHeII1 = data.dcs_reHeII1_reHeII1;
    double *reHeII2_reHeII2 = data.c_reHeII2_reHeII2;
    double *cs_reHeII2_reHeII2  = data.cs_reHeII2_reHeII2;
    double *dcs_reHeII2_reHeII2 = data.dcs_reHeII2_reHeII2;
    double *reHeIII_reHeIII = data.c_reHeIII_reHeIII;
    double *cs_reHeIII_reHeIII  = data.cs_reHeIII_reHeIII;
    double *dcs_reHeIII_reHeIII = data.dcs_reHeIII_reHeIII;
    double *reHII_reHII = data.c_reHII_reHII;
    double *cs_reHII_reHII  = data.cs_reHII_reHII;
    double *dcs_reHII_reHII = data.dcs_reHII_reHII;

    if (j < BATCHSIZE)
    {
        k    = __float2int_rz(data.idbin*data.logTs[j] - data.lb);
        t1   = data.lb + k*data.dbin;
        Tdef = (data.logTs[j] - t1) * data.idbin;
        rs_k01[j] = Tdef*k01[k+1] + (-k01[k]*Tdef + k01[k]);
        rs_k02[j] = Tdef*k02[k+1] + (-k02[k]*Tdef + k02[k]);
        rs_k03[j] = Tdef*k03[k+1] + (-k03[k]*Tdef + k03[k]);
        rs_k04[j] = Tdef*k04[k+1] + (-k04[k]*Tdef + k04[k]);
        rs_k05[j] = Tdef*k05[k+1] + (-k05[k]*Tdef + k05[k]);
        rs_k06[j] = Tdef*k06[k+1] + (-k06[k]*Tdef + k06[k]);
        rs_k07[j] = Tdef*k07[k+1] + (-k07[k]*Tdef + k07[k]);
        rs_k08[j] = Tdef*k08[k+1] + (-k08[k]*Tdef + k08[k]);
        rs_k09[j] = Tdef*k09[k+1] + (-k09[k]*Tdef + k09[k]);
        rs_k10[j] = Tdef*k10[k+1] + (-k10[k]*Tdef + k10[k]);
        rs_k11[j] = Tdef*k11[k+1] + (-k11[k]*Tdef + k11[k]);
        rs_k12[j] = Tdef*k12[k+1] + (-k12[k]*Tdef + k12[k]);
        rs_k13[j] = Tdef*k13[k+1] + (-k13[k]*Tdef + k13[k]);
        rs_k14[j] = Tdef*k14[k+1] + (-k14[k]*Tdef + k14[k]);
        rs_k15[j] = Tdef*k15[k+1] + (-k15[k]*Tdef + k15[k]);
        rs_k16[j] = Tdef*k16[k+1] + (-k16[k]*Tdef + k16[k]);
        rs_k17[j] = Tdef*k17[k+1] + (-k17[k]*Tdef + k17[k]);
        rs_k18[j] = Tdef*k18[k+1] + (-k18[k]*Tdef + k18[k]);
        rs_k19[j] = Tdef*k19[k+1] + (-k19[k]*Tdef + k19[k]);
        rs_k21[j] = Tdef*k21[k+1] + (-k21[k]*Tdef + k21[k]);
        rs_k22[j] = Tdef*k22[k+1] + (-k22[k]*Tdef + k22[k]);
        rs_k23[j] = Tdef*k23[k+1] + (-k23[k]*Tdef + k23[k]);
        cs_brem_brem[j] = Tdef*brem_brem[k+1] + (-brem_brem[k]*Tdef + brem_brem[k]);
        cs_ceHeI_ceHeI[j] = Tdef*ceHeI_ceHeI[k+1] + (-ceHeI_ceHeI[k]*Tdef + ceHeI_ceHeI[k]);
        cs_ceHeII_ceHeII[j] = Tdef*ceHeII_ceHeII[k+1] + (-ceHeII_ceHeII[k]*Tdef + ceHeII_ceHeII[k]);
        cs_ceHI_ceHI[j] = Tdef*ceHI_ceHI[k+1] + (-ceHI_ceHI[k]*Tdef + ceHI_ceHI[k]);
        cs_cie_cooling_cieco[j] = Tdef*cie_cooling_cieco[k+1] + (-cie_cooling_cieco[k]*Tdef + cie_cooling_cieco[k]);
        cs_ciHeI_ciHeI[j] = Tdef*ciHeI_ciHeI[k+1] + (-ciHeI_ciHeI[k]*Tdef + ciHeI_ciHeI[k]);
        cs_ciHeII_ciHeII[j] = Tdef*ciHeII_ciHeII[k+1] + (-ciHeII_ciHeII[k]*Tdef + ciHeII_ciHeII[k]);
        cs_ciHeIS_ciHeIS[j] = Tdef*ciHeIS_ciHeIS[k+1] + (-ciHeIS_ciHeIS[k]*Tdef + ciHeIS_ciHeIS[k]);
        cs_ciHI_ciHI[j] = Tdef*ciHI_ciHI[k+1] + (-ciHI_ciHI[k]*Tdef + ciHI_ciHI[k]);
        cs_compton_comp_[j] = Tdef*compton_comp_[k+1] + (-compton_comp_[k]*Tdef + compton_comp_[k]);
        cs_gammah_gammah[j] = Tdef*gammah_gammah[k+1] + (-gammah_gammah[k]*Tdef + gammah_gammah[k]);
        cs_gloverabel08_gael[j] = Tdef*gloverabel08_gael[k+1] + (-gloverabel08_gael[k]*Tdef + gloverabel08_gael[k]);
        cs_gloverabel08_gaH2[j] = Tdef*gloverabel08_gaH2[k+1] + (-gloverabel08_gaH2[k]*Tdef + gloverabel08_gaH2[k]);
        cs_gloverabel08_gaHe[j] = Tdef*gloverabel08_gaHe[k+1] + (-gloverabel08_gaHe[k]*Tdef + gloverabel08_gaHe[k]);
        cs_gloverabel08_gaHI[j] = Tdef*gloverabel08_gaHI[k+1] + (-gloverabel08_gaHI[k]*Tdef + gloverabel08_gaHI[k]);
        cs_gloverabel08_gaHp[j] = Tdef*gloverabel08_gaHp[k+1] + (-gloverabel08_gaHp[k]*Tdef + gloverabel08_gaHp[k]);
        cs_gloverabel08_gphdl[j] = Tdef*gloverabel08_gphdl[k+1] + (-gloverabel08_gphdl[k]*Tdef + gloverabel08_gphdl[k]);
        cs_gloverabel08_gpldl[j] = Tdef*gloverabel08_gpldl[k+1] + (-gloverabel08_gpldl[k]*Tdef + gloverabel08_gpldl[k]);
        cs_gloverabel08_h2lte[j] = Tdef*gloverabel08_h2lte[k+1] + (-gloverabel08_h2lte[k]*Tdef + gloverabel08_h2lte[k]);
        cs_h2formation_h2mcool[j] = Tdef*h2formation_h2mcool[k+1] + (-h2formation_h2mcool[k]*Tdef + h2formation_h2mcool[k]);
        cs_h2formation_h2mheat[j] = Tdef*h2formation_h2mheat[k+1] + (-h2formation_h2mheat[k]*Tdef + h2formation_h2mheat[k]);
        cs_h2formation_ncrd1[j] = Tdef*h2formation_ncrd1[k+1] + (-h2formation_ncrd1[k]*Tdef + h2formation_ncrd1[k]);
        cs_h2formation_ncrd2[j] = Tdef*h2formation_ncrd2[k+1] + (-h2formation_ncrd2[k]*Tdef + h2formation_ncrd2[k]);
        cs_h2formation_ncrn[j] = Tdef*h2formation_ncrn[k+1] + (-h2formation_ncrn[k]*Tdef + h2formation_ncrn[k]);
        cs_reHeII1_reHeII1[j] = Tdef*reHeII1_reHeII1[k+1] + (-reHeII1_reHeII1[k]*Tdef + reHeII1_reHeII1[k]);
        cs_reHeII2_reHeII2[j] = Tdef*reHeII2_reHeII2[k+1] + (-reHeII2_reHeII2[k]*Tdef + reHeII2_reHeII2[k]);
        cs_reHeIII_reHeIII[j] = Tdef*reHeIII_reHeIII[k+1] + (-reHeIII_reHeIII[k]*Tdef + reHeIII_reHeIII[k]);
        cs_reHII_reHII[j] = Tdef*reHII_reHII[k+1] + (-reHII_reHII[k]*Tdef + reHII_reHII[k]);

    }
}


    __global__
static void rhs_kernel(double y, double *ydata, double *ydotdata, primordial_cuda_data data)
{
    int i = blockIdx.x* blockDim.x + threadIdx.x;

    int groupi = i * nchem; 

    // get rate pointer
    double *k01 = data.rs_k01;
    double *k02 = data.rs_k02;
    double *k03 = data.rs_k03;
    double *k04 = data.rs_k04;
    double *k05 = data.rs_k05;
    double *k06 = data.rs_k06;
    double *k07 = data.rs_k07;
    double *k08 = data.rs_k08;
    double *k09 = data.rs_k09;
    double *k10 = data.rs_k10;
    double *k11 = data.rs_k11;
    double *k12 = data.rs_k12;
    double *k13 = data.rs_k13;
    double *k14 = data.rs_k14;
    double *k15 = data.rs_k15;
    double *k16 = data.rs_k16;
    double *k17 = data.rs_k17;
    double *k18 = data.rs_k18;
    double *k19 = data.rs_k19;
    double *k21 = data.rs_k21;
    double *k22 = data.rs_k22;
    double *k23 = data.rs_k23;
    double *brem_brem = data.cs_brem_brem;
    double *ceHeI_ceHeI = data.cs_ceHeI_ceHeI;
    double *ceHeII_ceHeII = data.cs_ceHeII_ceHeII;
    double *ceHI_ceHI = data.cs_ceHI_ceHI;
    double *cie_cooling_cieco = data.cs_cie_cooling_cieco;
    double *ciHeI_ciHeI = data.cs_ciHeI_ciHeI;
    double *ciHeII_ciHeII = data.cs_ciHeII_ciHeII;
    double *ciHeIS_ciHeIS = data.cs_ciHeIS_ciHeIS;
    double *ciHI_ciHI = data.cs_ciHI_ciHI;
    double *compton_comp_ = data.cs_compton_comp_;
    double *gammah_gammah = data.cs_gammah_gammah;
    double *gloverabel08_gael = data.cs_gloverabel08_gael;
    double *gloverabel08_gaH2 = data.cs_gloverabel08_gaH2;
    double *gloverabel08_gaHe = data.cs_gloverabel08_gaHe;
    double *gloverabel08_gaHI = data.cs_gloverabel08_gaHI;
    double *gloverabel08_gaHp = data.cs_gloverabel08_gaHp;
    double *gloverabel08_gphdl = data.cs_gloverabel08_gphdl;
    double *gloverabel08_gpldl = data.cs_gloverabel08_gpldl;
    double *gloverabel08_h2lte = data.cs_gloverabel08_h2lte;
    double *h2formation_h2mcool = data.cs_h2formation_h2mcool;
    double *h2formation_h2mheat = data.cs_h2formation_h2mheat;
    double *h2formation_ncrd1 = data.cs_h2formation_ncrd1;
    double *h2formation_ncrd2 = data.cs_h2formation_ncrd2;
    double *h2formation_ncrn = data.cs_h2formation_ncrn;
    double *reHeII1_reHeII1 = data.cs_reHeII1_reHeII1;
    double *reHeII2_reHeII2 = data.cs_reHeII2_reHeII2;
    double *reHeIII_reHeIII = data.cs_reHeIII_reHeIII;
    double *reHII_reHII = data.cs_reHII_reHII;


    
    double h2_optical_depth_approx;    
    
    
    double cie_optical_depth_approx;
    

    int j;
    double z, T, mdensity, inv_mdensity;

    if (i < BATCHSIZE)
    {
        T = data.Ts[i];
        z = data.current_z;

        
        h2_optical_depth_approx = 1.0; // data.h2_optical_depth_approx[i];
        
        
        cie_optical_depth_approx = 1.0; // data.cie_optical_depth_approx[i];
        
        double H2_1 = ydata[groupi+0];
        double H2_2 = ydata[groupi+1];
        double H_1 = ydata[groupi+2];
        double H_2 = ydata[groupi+3];
        double H_m0 = ydata[groupi+4];
        double He_1 = ydata[groupi+5];
        double He_2 = ydata[groupi+6];
        double He_3 = ydata[groupi+7];
        double de = ydata[groupi+8];
        double ge = ydata[groupi+9];

        mdensity     = mh*(2.0*H2_1 + 2.0*H2_2 + 1.0079400000000001*H_1 + 1.0079400000000001*H_2 + 1.0079400000000001*H_m0 + 4.0026020000000004*He_1 + 4.0026020000000004*He_2 + 4.0026020000000004*He_3);
        inv_mdensity = 1.0/mdensity;
        //
        // Species: H2_1
        //
        j = 0;
        ydotdata[groupi+j] = k08[i]*H_1*H_m0 + k10[i]*H2_2*H_1 - k11[i]*H2_1*H_2 - k12[i]*H2_1*de - k13[i]*H2_1*H_1 + k19[i]*H2_2*H_m0 + k21[i]*H2_1*H_1*H_1 + k22[i]*H_1*H_1*H_1 - k23[i]*H2_1*H2_1;
        
        j++;
        //
        // Species: H2_2
        //
        j = 1;
        ydotdata[groupi+j] = k09[i]*H_1*H_2 - k10[i]*H2_2*H_1 + k11[i]*H2_1*H_2 + k17[i]*H_2*H_m0 - k18[i]*H2_2*de - k19[i]*H2_2*H_m0;
        
        j++;
        //
        // Species: H_1
        //
        j = 2;
        ydotdata[groupi+j] = -k01[i]*H_1*de + k02[i]*H_2*de - k07[i]*H_1*de - k08[i]*H_1*H_m0 - k09[i]*H_1*H_2 - k10[i]*H2_2*H_1 + k11[i]*H2_1*H_2 + 2*k12[i]*H2_1*de + 2*k13[i]*H2_1*H_1 + k14[i]*H_m0*de + k15[i]*H_1*H_m0 + 2*k16[i]*H_2*H_m0 + 2*k18[i]*H2_2*de + k19[i]*H2_2*H_m0 - 2*k21[i]*H2_1*H_1*H_1 - 2*k22[i]*H_1*H_1*H_1 + 2*k23[i]*H2_1*H2_1;
        
        j++;
        //
        // Species: H_2
        //
        j = 3;
        ydotdata[groupi+j] = k01[i]*H_1*de - k02[i]*H_2*de - k09[i]*H_1*H_2 + k10[i]*H2_2*H_1 - k11[i]*H2_1*H_2 - k16[i]*H_2*H_m0 - k17[i]*H_2*H_m0;
        
        j++;
        //
        // Species: H_m0
        //
        j = 4;
        ydotdata[groupi+j] = k07[i]*H_1*de - k08[i]*H_1*H_m0 - k14[i]*H_m0*de - k15[i]*H_1*H_m0 - k16[i]*H_2*H_m0 - k17[i]*H_2*H_m0 - k19[i]*H2_2*H_m0;
        
        j++;
        //
        // Species: He_1
        //
        j = 5;
        ydotdata[groupi+j] = -k03[i]*He_1*de + k04[i]*He_2*de;
        
        j++;
        //
        // Species: He_2
        //
        j = 6;
        ydotdata[groupi+j] = k03[i]*He_1*de - k04[i]*He_2*de - k05[i]*He_2*de + k06[i]*He_3*de;
        
        j++;
        //
        // Species: He_3
        //
        j = 7;
        ydotdata[groupi+j] = k05[i]*He_2*de - k06[i]*He_3*de;
        
        j++;
        //
        // Species: de
        //
        j = 8;
        ydotdata[groupi+j] = k01[i]*H_1*de - k02[i]*H_2*de + k03[i]*He_1*de - k04[i]*He_2*de + k05[i]*He_2*de - k06[i]*He_3*de - k07[i]*H_1*de + k08[i]*H_1*H_m0 + k14[i]*H_m0*de + k15[i]*H_1*H_m0 + k17[i]*H_2*H_m0 - k18[i]*H2_2*de;
        
        j++;
        //
        // Species: ge
        //
        j = 9;
        ydotdata[groupi+j] = -brem_brem[i]*cie_optical_depth_approx*de*(H_2 + He_2 + 4.0*He_3) - ceHI_ceHI[i]*H_1*cie_optical_depth_approx*de - ceHeII_ceHeII[i]*He_2*cie_optical_depth_approx*de - ceHeI_ceHeI[i]*He_2*cie_optical_depth_approx*pow(de, 2) - ciHI_ciHI[i]*H_1*cie_optical_depth_approx*de - ciHeII_ciHeII[i]*He_2*cie_optical_depth_approx*de - ciHeIS_ciHeIS[i]*He_2*cie_optical_depth_approx*pow(de, 2) - ciHeI_ciHeI[i]*He_1*cie_optical_depth_approx*de - 2.0158800000000001*cie_cooling_cieco[i]*H2_1*cie_optical_depth_approx*mdensity - compton_comp_[i]*cie_optical_depth_approx*de*pow(z + 1.0, 4)*(T - 2.73*z - 2.73) - gloverabel08_h2lte[i]*H2_1*cie_optical_depth_approx*h2_optical_depth_approx/(gloverabel08_h2lte[i]/(gloverabel08_gaH2[i]*H2_1 + gloverabel08_gaHI[i]*H_1 + gloverabel08_gaHe[i]*He_1 + gloverabel08_gaHp[i]*H_2 + gloverabel08_gael[i]*de) + 1.0) - reHII_reHII[i]*H_2*cie_optical_depth_approx*de - reHeII1_reHeII1[i]*He_2*cie_optical_depth_approx*de - reHeII2_reHeII2[i]*He_2*cie_optical_depth_approx*de - reHeIII_reHeIII[i]*He_3*cie_optical_depth_approx*de + 0.5*1.0/(h2formation_ncrn[i]/(h2formation_ncrd1[i]*H_1 + h2formation_ncrd2[i]*H2_1) + 1.0)*(-h2formation_h2mcool[i]*H2_1*H_1 + h2formation_h2mheat[i]*pow(H_1, 3));
        
        ydotdata[groupi+j] *= inv_mdensity;
        
        j++;

    }

    /*
       for (int ii = 0; ii < 10; ii++)
       {
       printf("from %d: ydot[%d] = %0.5g\n", i, ii, ydotdata[groupi+ii]);
       }
     */

}


    __global__
void temperature_kernel(double* ydata, primordial_cuda_data data)
{
    int i = blockIdx.x* blockDim.x + threadIdx.x;
    int groupi = i * nchem; 

    double *temperature = data.Ts;
    double *logTs      = data.logTs;
    double *Tge        = data.Tge;

    double gammaH2_1 = 7./5.;
    double gammaH2_2 = 7./5.;
    // as of now just do not use any "temperature-dependent" gamma
    // which simplifies my life, and not having the need to iterate it to convergence

    if (i < BATCHSIZE)
    {
        double H2_1 = ydata[groupi+0];
        double H2_2 = ydata[groupi+1];
        double H_1 = ydata[groupi+2];
        double H_2 = ydata[groupi+3];
        double H_m0 = ydata[groupi+4];
        double He_1 = ydata[groupi+5];
        double He_2 = ydata[groupi+6];
        double He_3 = ydata[groupi+7];
        double de = ydata[groupi+8];
        double ge = ydata[groupi+9];
        double density = 2.0*H2_1 + 2.0*H2_2 + 1.0079400000000001*H_1 + 1.0079400000000001*H_2 + 1.0079400000000001*H_m0 + 4.0026020000000004*He_1 + 4.0026020000000004*He_2 + 4.0026020000000004*He_3;
        temperature[i] = density*ge*mh/(kb*(H2_1/(gammaH2_1 - 1.0) + H2_2/(gammaH2_2 - 1.0) + H_1*_gamma_m1 + H_2*_gamma_m1 + H_m0*_gamma_m1 + He_1*_gamma_m1 + He_2*_gamma_m1 + He_3*_gamma_m1 + _gamma_m1*de));
        logTs      [i] = log(temperature[i]);
        Tge        [i] = 0.0; //TODO: update this to dT_dge;
    }
}

// Function Called by the solver
static int calculate_rhs_primordial_cuda(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
    primordial_cuda_data *udata = (primordial_cuda_data *) user_data;
    double *ydata    = N_VGetDeviceArrayPointer_Cuda(y);
    double *ydotdata = N_VGetDeviceArrayPointer_Cuda(ydot);

    // calculate temperature kernel
    temperature_kernel<<<GRIDSIZE, BLOCKSIZE>>> (ydata, *udata);
    // interpolate the rates with updated temperature
    linear_interpolation_kernel<<<GRIDSIZE, BLOCKSIZE>>>(*udata);
    // update ydot with the kernel function
    rhs_kernel<<<GRIDSIZE, BLOCKSIZE>>>(t, ydata, ydotdata, *udata);

    cudaDeviceSynchronize();
    cudaError_t cuerr = cudaGetLastError();
    if (cuerr != cudaSuccess) {
        fprintf(stderr,
                ">>> ERROR in f: cudaGetLastError returned %s\n",
                cudaGetErrorName(cuerr));
        return(-1);
    }


    return 0;

}


// write jacobian



/*
 * taken from cvRoberts_block_cusolversp_batchqr.cu
 *
 * Jacobian initialization routine. This sets the sparisty pattern of
 * the blocks of the Jacobian J(t,y) = df/dy. This is performed on the CPU,
 * and only occurs at the beginning of the simulation.
 */
static int blockSparseJacInit(SUNMatrix J)
{
    int nnz = 64;

    int rowptrs[nchem+1];
    int colvals[nnz];

    SUNMatZero(J);
    rowptrs[0] = 0;
    rowptrs[1] = 7;
    rowptrs[2] = 14;
    rowptrs[3] = 21;
    rowptrs[4] = 28;
    rowptrs[5] = 34;
    rowptrs[6] = 38;
    rowptrs[7] = 43;
    rowptrs[8] = 47;
    rowptrs[9] = 56;
    rowptrs[nchem] = nnz;
        // H2_1 by H2_1
        colvals[0] = 0;
        // H2_1 by H2_2
        colvals[1] = 1;
        // H2_1 by H_1
        colvals[2] = 2;
        // H2_1 by H_2
        colvals[3] = 3;
        // H2_1 by H_m0
        colvals[4] = 4;
        // H2_1 by de
        colvals[5] = 8;
        // H2_1 by ge
        colvals[6] = 9;
        // H2_2 by H2_1
        colvals[7] = 0;
        // H2_2 by H2_2
        colvals[8] = 1;
        // H2_2 by H_1
        colvals[9] = 2;
        // H2_2 by H_2
        colvals[10] = 3;
        // H2_2 by H_m0
        colvals[11] = 4;
        // H2_2 by de
        colvals[12] = 8;
        // H2_2 by ge
        colvals[13] = 9;
        // H_1 by H2_1
        colvals[14] = 0;
        // H_1 by H2_2
        colvals[15] = 1;
        // H_1 by H_1
        colvals[16] = 2;
        // H_1 by H_2
        colvals[17] = 3;
        // H_1 by H_m0
        colvals[18] = 4;
        // H_1 by de
        colvals[19] = 8;
        // H_1 by ge
        colvals[20] = 9;
        // H_2 by H2_1
        colvals[21] = 0;
        // H_2 by H2_2
        colvals[22] = 1;
        // H_2 by H_1
        colvals[23] = 2;
        // H_2 by H_2
        colvals[24] = 3;
        // H_2 by H_m0
        colvals[25] = 4;
        // H_2 by de
        colvals[26] = 8;
        // H_2 by ge
        colvals[27] = 9;
        // H_m0 by H2_2
        colvals[28] = 1;
        // H_m0 by H_1
        colvals[29] = 2;
        // H_m0 by H_2
        colvals[30] = 3;
        // H_m0 by H_m0
        colvals[31] = 4;
        // H_m0 by de
        colvals[32] = 8;
        // H_m0 by ge
        colvals[33] = 9;
        // He_1 by He_1
        colvals[34] = 5;
        // He_1 by He_2
        colvals[35] = 6;
        // He_1 by de
        colvals[36] = 8;
        // He_1 by ge
        colvals[37] = 9;
        // He_2 by He_1
        colvals[38] = 5;
        // He_2 by He_2
        colvals[39] = 6;
        // He_2 by He_3
        colvals[40] = 7;
        // He_2 by de
        colvals[41] = 8;
        // He_2 by ge
        colvals[42] = 9;
        // He_3 by He_2
        colvals[43] = 6;
        // He_3 by He_3
        colvals[44] = 7;
        // He_3 by de
        colvals[45] = 8;
        // He_3 by ge
        colvals[46] = 9;
        // de by H2_2
        colvals[47] = 1;
        // de by H_1
        colvals[48] = 2;
        // de by H_2
        colvals[49] = 3;
        // de by H_m0
        colvals[50] = 4;
        // de by He_1
        colvals[51] = 5;
        // de by He_2
        colvals[52] = 6;
        // de by He_3
        colvals[53] = 7;
        // de by de
        colvals[54] = 8;
        // de by ge
        colvals[55] = 9;
        // ge by H2_1
        colvals[56] = 0;
        // ge by H_1
        colvals[57] = 2;
        // ge by H_2
        colvals[58] = 3;
        // ge by He_1
        colvals[59] = 5;
        // ge by He_2
        colvals[60] = 6;
        // ge by He_3
        colvals[61] = 7;
        // ge by de
        colvals[62] = 8;
        // ge by ge
        colvals[63] = 9;

    // copy rowptrs, colvals to the device
    SUNMatrix_cuSparse_CopyToDevice(J, NULL, rowptrs, colvals);
    cudaDeviceSynchronize();
    return (0);
}

static int blockJacInit(SUNMatrix J)
{
    int nnz = nchem*nchem;

    int rowptrs[nchem+1];
    int colvals[nnz];

    SUNMatZero(J);
    for (int r = 0; r < nchem+1; r++)
    {
        rowptrs[r] = r*nchem;
        // printf("rowptrs[%d] = %d\n", r, rowptrs[r]);
    }

    int bIdx;
    for (int c = 0; c < nnz; c++)
    {
        bIdx = c /nnz; 
        colvals[c] = bIdx*nchem + c%nchem;
        // printf("colvals[%d] = %d\n", c, colvals[c]);
    }
    // copy rowptrs, colvals to the device
    SUNMatrix_cuSparse_CopyToDevice(J, NULL, rowptrs, colvals);
    cudaDeviceSynchronize();
    return (0);
}


static int JacInit(SUNMatrix J)
{
    int rowptrs[10+1];
    int colvals[10 * 10];

    /* Zero out the Jacobian */
    SUNMatZero(J);

    /* there are 10 entries per row*/
    rowptrs[0] = 0 * 10;
    rowptrs[1] = 1 * 10;
    rowptrs[2] = 2 * 10;
    rowptrs[3] = 3 * 10;
    rowptrs[4] = 4 * 10;
    rowptrs[5] = 5 * 10;
    rowptrs[6] = 6 * 10;
    rowptrs[7] = 7 * 10;
    rowptrs[8] = 8 * 10;
    rowptrs[9] = 9 * 10;
    
    // 0 row of block
    colvals[0] = 0;
    colvals[1] = 1;
    colvals[2] = 2;
    colvals[3] = 3;
    colvals[4] = 4;
    colvals[5] = 5;
    colvals[6] = 6;
    colvals[7] = 7;
    colvals[8] = 8;
    colvals[9] = 9;
    
    // 1 row of block
    colvals[10] = 0;
    colvals[11] = 1;
    colvals[12] = 2;
    colvals[13] = 3;
    colvals[14] = 4;
    colvals[15] = 5;
    colvals[16] = 6;
    colvals[17] = 7;
    colvals[18] = 8;
    colvals[19] = 9;
    
    // 2 row of block
    colvals[20] = 0;
    colvals[21] = 1;
    colvals[22] = 2;
    colvals[23] = 3;
    colvals[24] = 4;
    colvals[25] = 5;
    colvals[26] = 6;
    colvals[27] = 7;
    colvals[28] = 8;
    colvals[29] = 9;
    
    // 3 row of block
    colvals[30] = 0;
    colvals[31] = 1;
    colvals[32] = 2;
    colvals[33] = 3;
    colvals[34] = 4;
    colvals[35] = 5;
    colvals[36] = 6;
    colvals[37] = 7;
    colvals[38] = 8;
    colvals[39] = 9;
    
    // 4 row of block
    colvals[40] = 0;
    colvals[41] = 1;
    colvals[42] = 2;
    colvals[43] = 3;
    colvals[44] = 4;
    colvals[45] = 5;
    colvals[46] = 6;
    colvals[47] = 7;
    colvals[48] = 8;
    colvals[49] = 9;
    
    // 5 row of block
    colvals[50] = 0;
    colvals[51] = 1;
    colvals[52] = 2;
    colvals[53] = 3;
    colvals[54] = 4;
    colvals[55] = 5;
    colvals[56] = 6;
    colvals[57] = 7;
    colvals[58] = 8;
    colvals[59] = 9;
    
    // 6 row of block
    colvals[60] = 0;
    colvals[61] = 1;
    colvals[62] = 2;
    colvals[63] = 3;
    colvals[64] = 4;
    colvals[65] = 5;
    colvals[66] = 6;
    colvals[67] = 7;
    colvals[68] = 8;
    colvals[69] = 9;
    
    // 7 row of block
    colvals[70] = 0;
    colvals[71] = 1;
    colvals[72] = 2;
    colvals[73] = 3;
    colvals[74] = 4;
    colvals[75] = 5;
    colvals[76] = 6;
    colvals[77] = 7;
    colvals[78] = 8;
    colvals[79] = 9;
    
    // 8 row of block
    colvals[80] = 0;
    colvals[81] = 1;
    colvals[82] = 2;
    colvals[83] = 3;
    colvals[84] = 4;
    colvals[85] = 5;
    colvals[86] = 6;
    colvals[87] = 7;
    colvals[88] = 8;
    colvals[89] = 9;
    
    // 9 row of block
    colvals[90] = 0;
    colvals[91] = 1;
    colvals[92] = 2;
    colvals[93] = 3;
    colvals[94] = 4;
    colvals[95] = 5;
    colvals[96] = 6;
    colvals[97] = 7;
    colvals[98] = 8;
    colvals[99] = 9;

    // copy rowptrs, colvals to the device
    SUNMatrix_cuSparse_CopyToDevice(J, NULL, rowptrs, colvals);
    cudaDeviceSynchronize();
    return (0);
}


/* Jacobian evaluation with GPU */
    __global__
static void jacobian_kernel(realtype *ydata, realtype *Jdata, primordial_cuda_data data)
{
    int groupj;

    // temporary:
    int nnzper = nchem*nchem;
    int i;
    double *Tge = data.Tge;
    double z, T;

    
    double *h2_optical_depth_approx_arr = data.h2_optical_depth_approx;
    
    
    double *cie_optical_depth_approx_arr = data.cie_optical_depth_approx;
    
    groupj = blockIdx.x*blockDim.x + threadIdx.x; 

    T = 1000.0;
    z = 0.0;

    if (groupj < BATCHSIZE)
    {
        i = groupj;
        
        double h2_optical_depth_approx = 1.0; //h2_optical_depth_approx_arr[i];
        
        
        double cie_optical_depth_approx = 1.0; //cie_optical_depth_approx_arr[i];
        
        // pulled the species data
        double H2_1 = ydata[nchem*groupj+0];
        double H2_2 = ydata[nchem*groupj+1];
        double H_1 = ydata[nchem*groupj+2];
        double H_2 = ydata[nchem*groupj+3];
        double H_m0 = ydata[nchem*groupj+4];
        double He_1 = ydata[nchem*groupj+5];
        double He_2 = ydata[nchem*groupj+6];
        double He_3 = ydata[nchem*groupj+7];
        double de = ydata[nchem*groupj+8];
        double ge = ydata[nchem*groupj+9];
        double mdensity = mh * (2.0*H2_1 + 2.0*H2_2 + 1.0079400000000001*H_1 + 1.0079400000000001*H_2 + 1.0079400000000001*H_m0 + 4.0026020000000004*He_1 + 4.0026020000000004*He_2 + 4.0026020000000004*He_3);
        double inv_mdensity = 1.0/ mdensity;
        double *k01 = data.rs_k01;
        double *rk01= data.drs_k01;
        double *k02 = data.rs_k02;
        double *rk02= data.drs_k02;
        double *k03 = data.rs_k03;
        double *rk03= data.drs_k03;
        double *k04 = data.rs_k04;
        double *rk04= data.drs_k04;
        double *k05 = data.rs_k05;
        double *rk05= data.drs_k05;
        double *k06 = data.rs_k06;
        double *rk06= data.drs_k06;
        double *k07 = data.rs_k07;
        double *rk07= data.drs_k07;
        double *k08 = data.rs_k08;
        double *rk08= data.drs_k08;
        double *k09 = data.rs_k09;
        double *rk09= data.drs_k09;
        double *k10 = data.rs_k10;
        double *rk10= data.drs_k10;
        double *k11 = data.rs_k11;
        double *rk11= data.drs_k11;
        double *k12 = data.rs_k12;
        double *rk12= data.drs_k12;
        double *k13 = data.rs_k13;
        double *rk13= data.drs_k13;
        double *k14 = data.rs_k14;
        double *rk14= data.drs_k14;
        double *k15 = data.rs_k15;
        double *rk15= data.drs_k15;
        double *k16 = data.rs_k16;
        double *rk16= data.drs_k16;
        double *k17 = data.rs_k17;
        double *rk17= data.drs_k17;
        double *k18 = data.rs_k18;
        double *rk18= data.drs_k18;
        double *k19 = data.rs_k19;
        double *rk19= data.drs_k19;
        double *k21 = data.rs_k21;
        double *rk21= data.drs_k21;
        double *k22 = data.rs_k22;
        double *rk22= data.drs_k22;
        double *k23 = data.rs_k23;
        double *rk23= data.drs_k23;
        double *brem_brem = data.cs_brem_brem;
        double *rbrem_brem = data.dcs_brem_brem;
        double *ceHeI_ceHeI = data.cs_ceHeI_ceHeI;
        double *rceHeI_ceHeI = data.dcs_ceHeI_ceHeI;
        double *ceHeII_ceHeII = data.cs_ceHeII_ceHeII;
        double *rceHeII_ceHeII = data.dcs_ceHeII_ceHeII;
        double *ceHI_ceHI = data.cs_ceHI_ceHI;
        double *rceHI_ceHI = data.dcs_ceHI_ceHI;
        double *cie_cooling_cieco = data.cs_cie_cooling_cieco;
        double *rcie_cooling_cieco = data.dcs_cie_cooling_cieco;
        double *ciHeI_ciHeI = data.cs_ciHeI_ciHeI;
        double *rciHeI_ciHeI = data.dcs_ciHeI_ciHeI;
        double *ciHeII_ciHeII = data.cs_ciHeII_ciHeII;
        double *rciHeII_ciHeII = data.dcs_ciHeII_ciHeII;
        double *ciHeIS_ciHeIS = data.cs_ciHeIS_ciHeIS;
        double *rciHeIS_ciHeIS = data.dcs_ciHeIS_ciHeIS;
        double *ciHI_ciHI = data.cs_ciHI_ciHI;
        double *rciHI_ciHI = data.dcs_ciHI_ciHI;
        double *compton_comp_ = data.cs_compton_comp_;
        double *rcompton_comp_ = data.dcs_compton_comp_;
        double *gammah_gammah = data.cs_gammah_gammah;
        double *rgammah_gammah = data.dcs_gammah_gammah;
        double *gloverabel08_gael = data.cs_gloverabel08_gael;
        double *rgloverabel08_gael = data.dcs_gloverabel08_gael;
        double *gloverabel08_gaH2 = data.cs_gloverabel08_gaH2;
        double *rgloverabel08_gaH2 = data.dcs_gloverabel08_gaH2;
        double *gloverabel08_gaHe = data.cs_gloverabel08_gaHe;
        double *rgloverabel08_gaHe = data.dcs_gloverabel08_gaHe;
        double *gloverabel08_gaHI = data.cs_gloverabel08_gaHI;
        double *rgloverabel08_gaHI = data.dcs_gloverabel08_gaHI;
        double *gloverabel08_gaHp = data.cs_gloverabel08_gaHp;
        double *rgloverabel08_gaHp = data.dcs_gloverabel08_gaHp;
        double *gloverabel08_gphdl = data.cs_gloverabel08_gphdl;
        double *rgloverabel08_gphdl = data.dcs_gloverabel08_gphdl;
        double *gloverabel08_gpldl = data.cs_gloverabel08_gpldl;
        double *rgloverabel08_gpldl = data.dcs_gloverabel08_gpldl;
        double *gloverabel08_h2lte = data.cs_gloverabel08_h2lte;
        double *rgloverabel08_h2lte = data.dcs_gloverabel08_h2lte;
        double *h2formation_h2mcool = data.cs_h2formation_h2mcool;
        double *rh2formation_h2mcool = data.dcs_h2formation_h2mcool;
        double *h2formation_h2mheat = data.cs_h2formation_h2mheat;
        double *rh2formation_h2mheat = data.dcs_h2formation_h2mheat;
        double *h2formation_ncrd1 = data.cs_h2formation_ncrd1;
        double *rh2formation_ncrd1 = data.dcs_h2formation_ncrd1;
        double *h2formation_ncrd2 = data.cs_h2formation_ncrd2;
        double *rh2formation_ncrd2 = data.dcs_h2formation_ncrd2;
        double *h2formation_ncrn = data.cs_h2formation_ncrn;
        double *rh2formation_ncrn = data.dcs_h2formation_ncrn;
        double *reHeII1_reHeII1 = data.cs_reHeII1_reHeII1;
        double *rreHeII1_reHeII1 = data.dcs_reHeII1_reHeII1;
        double *reHeII2_reHeII2 = data.cs_reHeII2_reHeII2;
        double *rreHeII2_reHeII2 = data.dcs_reHeII2_reHeII2;
        double *reHeIII_reHeIII = data.cs_reHeIII_reHeIII;
        double *rreHeIII_reHeIII = data.dcs_reHeIII_reHeIII;
        double *reHII_reHII = data.cs_reHII_reHII;
        double *rreHII_reHII = data.dcs_reHII_reHII;
        //
        // Species: H2_1
        //
        
        
        // H2_1 by H2_1

        
        Jdata[nnzper*groupj + 0*nchem+ 0] = -k11[i]*H_2 - k12[i]*de - k13[i]*H_1 + k21[i]*pow(H_1, 2) - 2*k23[i]*H2_1;
        

        

        
        
        // H2_1 by H2_2

        
        Jdata[nnzper*groupj + 0*nchem+ 1] = k10[i]*H_1 + k19[i]*H_m0;
        

        

        
        
        // H2_1 by H_1

        
        Jdata[nnzper*groupj + 0*nchem+ 2] = k08[i]*H_m0 + k10[i]*H2_2 - k13[i]*H2_1 + 2*k21[i]*H2_1*H_1 + 3*k22[i]*pow(H_1, 2);
        

        

        
        
        // H2_1 by H_2

        
        Jdata[nnzper*groupj + 0*nchem+ 3] = -k11[i]*H2_1;
        

        

        
        
        // H2_1 by H_m0

        
        Jdata[nnzper*groupj + 0*nchem+ 4] = k08[i]*H_1 + k19[i]*H2_2;
        

        

        
        
        // H2_1 by He_1

        
        // because the Jacobian is initialized to zeros by default
        Jdata[nnzper*groupj+ 0*nchem+ 5] = ZERO;
        

        

        
        
        // H2_1 by He_2

        
        // because the Jacobian is initialized to zeros by default
        Jdata[nnzper*groupj+ 0*nchem+ 6] = ZERO;
        

        

        
        
        // H2_1 by He_3

        
        // because the Jacobian is initialized to zeros by default
        Jdata[nnzper*groupj+ 0*nchem+ 7] = ZERO;
        

        

        
        
        // H2_1 by de

        
        Jdata[nnzper*groupj + 0*nchem+ 8] = -k12[i]*H2_1;
        

        

        
        
        // H2_1 by ge

        
        Jdata[nnzper*groupj + 0*nchem+ 9] = -pow(H2_1, 2)*rk23[i] + H2_1*pow(H_1, 2)*rk21[i] - H2_1*H_1*rk13[i] - H2_1*H_2*rk11[i] - H2_1*de*rk12[i] + H2_2*H_1*rk10[i] + H2_2*H_m0*rk19[i] + pow(H_1, 3)*rk22[i] + H_1*H_m0*rk08[i];
        

        

        
        Jdata[nnzper*groupj+ 0*nchem+ 9] *= Tge[i];
        
        //
        // Species: H2_2
        //
        
        
        // H2_2 by H2_1

        
        Jdata[nnzper*groupj + 1*nchem+ 0] = k11[i]*H_2;
        

        

        
        
        // H2_2 by H2_2

        
        Jdata[nnzper*groupj + 1*nchem+ 1] = -k10[i]*H_1 - k18[i]*de - k19[i]*H_m0;
        

        

        
        
        // H2_2 by H_1

        
        Jdata[nnzper*groupj + 1*nchem+ 2] = k09[i]*H_2 - k10[i]*H2_2;
        

        

        
        
        // H2_2 by H_2

        
        Jdata[nnzper*groupj + 1*nchem+ 3] = k09[i]*H_1 + k11[i]*H2_1 + k17[i]*H_m0;
        

        

        
        
        // H2_2 by H_m0

        
        Jdata[nnzper*groupj + 1*nchem+ 4] = k17[i]*H_2 - k19[i]*H2_2;
        

        

        
        
        // H2_2 by He_1

        
        // because the Jacobian is initialized to zeros by default
        Jdata[nnzper*groupj+ 1*nchem+ 5] = ZERO;
        

        

        
        
        // H2_2 by He_2

        
        // because the Jacobian is initialized to zeros by default
        Jdata[nnzper*groupj+ 1*nchem+ 6] = ZERO;
        

        

        
        
        // H2_2 by He_3

        
        // because the Jacobian is initialized to zeros by default
        Jdata[nnzper*groupj+ 1*nchem+ 7] = ZERO;
        

        

        
        
        // H2_2 by de

        
        Jdata[nnzper*groupj + 1*nchem+ 8] = -k18[i]*H2_2;
        

        

        
        
        // H2_2 by ge

        
        Jdata[nnzper*groupj + 1*nchem+ 9] = H2_1*H_2*rk11[i] - H2_2*H_1*rk10[i] - H2_2*H_m0*rk19[i] - H2_2*de*rk18[i] + H_1*H_2*rk09[i] + H_2*H_m0*rk17[i];
        

        

        
        Jdata[nnzper*groupj+ 1*nchem+ 9] *= Tge[i];
        
        //
        // Species: H_1
        //
        
        
        // H_1 by H2_1

        
        Jdata[nnzper*groupj + 2*nchem+ 0] = k11[i]*H_2 + 2*k12[i]*de + 2*k13[i]*H_1 - 2*k21[i]*pow(H_1, 2) + 4*k23[i]*H2_1;
        

        

        
        
        // H_1 by H2_2

        
        Jdata[nnzper*groupj + 2*nchem+ 1] = -k10[i]*H_1 + 2*k18[i]*de + k19[i]*H_m0;
        

        

        
        
        // H_1 by H_1

        
        Jdata[nnzper*groupj + 2*nchem+ 2] = -k01[i]*de - k07[i]*de - k08[i]*H_m0 - k09[i]*H_2 - k10[i]*H2_2 + 2*k13[i]*H2_1 + k15[i]*H_m0 - 4*k21[i]*H2_1*H_1 - 6*k22[i]*pow(H_1, 2);
        

        

        
        
        // H_1 by H_2

        
        Jdata[nnzper*groupj + 2*nchem+ 3] = k02[i]*de - k09[i]*H_1 + k11[i]*H2_1 + 2*k16[i]*H_m0;
        

        

        
        
        // H_1 by H_m0

        
        Jdata[nnzper*groupj + 2*nchem+ 4] = -k08[i]*H_1 + k14[i]*de + k15[i]*H_1 + 2*k16[i]*H_2 + k19[i]*H2_2;
        

        

        
        
        // H_1 by He_1

        
        // because the Jacobian is initialized to zeros by default
        Jdata[nnzper*groupj+ 2*nchem+ 5] = ZERO;
        

        

        
        
        // H_1 by He_2

        
        // because the Jacobian is initialized to zeros by default
        Jdata[nnzper*groupj+ 2*nchem+ 6] = ZERO;
        

        

        
        
        // H_1 by He_3

        
        // because the Jacobian is initialized to zeros by default
        Jdata[nnzper*groupj+ 2*nchem+ 7] = ZERO;
        

        

        
        
        // H_1 by de

        
        Jdata[nnzper*groupj + 2*nchem+ 8] = -k01[i]*H_1 + k02[i]*H_2 - k07[i]*H_1 + 2*k12[i]*H2_1 + k14[i]*H_m0 + 2*k18[i]*H2_2;
        

        

        
        
        // H_1 by ge

        
        Jdata[nnzper*groupj + 2*nchem+ 9] = 2*pow(H2_1, 2)*rk23[i] - 2*H2_1*pow(H_1, 2)*rk21[i] + 2*H2_1*H_1*rk13[i] + H2_1*H_2*rk11[i] + 2*H2_1*de*rk12[i] - H2_2*H_1*rk10[i] + H2_2*H_m0*rk19[i] + 2*H2_2*de*rk18[i] - 2*pow(H_1, 3)*rk22[i] - H_1*H_2*rk09[i] - H_1*H_m0*rk08[i] + H_1*H_m0*rk15[i] - H_1*de*rk01[i] - H_1*de*rk07[i] + 2*H_2*H_m0*rk16[i] + H_2*de*rk02[i] + H_m0*de*rk14[i];
        

        

        
        Jdata[nnzper*groupj+ 2*nchem+ 9] *= Tge[i];
        
        //
        // Species: H_2
        //
        
        
        // H_2 by H2_1

        
        Jdata[nnzper*groupj + 3*nchem+ 0] = -k11[i]*H_2;
        

        

        
        
        // H_2 by H2_2

        
        Jdata[nnzper*groupj + 3*nchem+ 1] = k10[i]*H_1;
        

        

        
        
        // H_2 by H_1

        
        Jdata[nnzper*groupj + 3*nchem+ 2] = k01[i]*de - k09[i]*H_2 + k10[i]*H2_2;
        

        

        
        
        // H_2 by H_2

        
        Jdata[nnzper*groupj + 3*nchem+ 3] = -k02[i]*de - k09[i]*H_1 - k11[i]*H2_1 - k16[i]*H_m0 - k17[i]*H_m0;
        

        

        
        
        // H_2 by H_m0

        
        Jdata[nnzper*groupj + 3*nchem+ 4] = -k16[i]*H_2 - k17[i]*H_2;
        

        

        
        
        // H_2 by He_1

        
        // because the Jacobian is initialized to zeros by default
        Jdata[nnzper*groupj+ 3*nchem+ 5] = ZERO;
        

        

        
        
        // H_2 by He_2

        
        // because the Jacobian is initialized to zeros by default
        Jdata[nnzper*groupj+ 3*nchem+ 6] = ZERO;
        

        

        
        
        // H_2 by He_3

        
        // because the Jacobian is initialized to zeros by default
        Jdata[nnzper*groupj+ 3*nchem+ 7] = ZERO;
        

        

        
        
        // H_2 by de

        
        Jdata[nnzper*groupj + 3*nchem+ 8] = k01[i]*H_1 - k02[i]*H_2;
        

        

        
        
        // H_2 by ge

        
        Jdata[nnzper*groupj + 3*nchem+ 9] = -H2_1*H_2*rk11[i] + H2_2*H_1*rk10[i] - H_1*H_2*rk09[i] + H_1*de*rk01[i] - H_2*H_m0*rk16[i] - H_2*H_m0*rk17[i] - H_2*de*rk02[i];
        

        

        
        Jdata[nnzper*groupj+ 3*nchem+ 9] *= Tge[i];
        
        //
        // Species: H_m0
        //
        
        
        // H_m0 by H2_1

        
        // because the Jacobian is initialized to zeros by default
        Jdata[nnzper*groupj+ 4*nchem+ 0] = ZERO;
        

        

        
        
        // H_m0 by H2_2

        
        Jdata[nnzper*groupj + 4*nchem+ 1] = -k19[i]*H_m0;
        

        

        
        
        // H_m0 by H_1

        
        Jdata[nnzper*groupj + 4*nchem+ 2] = k07[i]*de - k08[i]*H_m0 - k15[i]*H_m0;
        

        

        
        
        // H_m0 by H_2

        
        Jdata[nnzper*groupj + 4*nchem+ 3] = -k16[i]*H_m0 - k17[i]*H_m0;
        

        

        
        
        // H_m0 by H_m0

        
        Jdata[nnzper*groupj + 4*nchem+ 4] = -k08[i]*H_1 - k14[i]*de - k15[i]*H_1 - k16[i]*H_2 - k17[i]*H_2 - k19[i]*H2_2;
        

        

        
        
        // H_m0 by He_1

        
        // because the Jacobian is initialized to zeros by default
        Jdata[nnzper*groupj+ 4*nchem+ 5] = ZERO;
        

        

        
        
        // H_m0 by He_2

        
        // because the Jacobian is initialized to zeros by default
        Jdata[nnzper*groupj+ 4*nchem+ 6] = ZERO;
        

        

        
        
        // H_m0 by He_3

        
        // because the Jacobian is initialized to zeros by default
        Jdata[nnzper*groupj+ 4*nchem+ 7] = ZERO;
        

        

        
        
        // H_m0 by de

        
        Jdata[nnzper*groupj + 4*nchem+ 8] = k07[i]*H_1 - k14[i]*H_m0;
        

        

        
        
        // H_m0 by ge

        
        Jdata[nnzper*groupj + 4*nchem+ 9] = -H2_2*H_m0*rk19[i] - H_1*H_m0*rk08[i] - H_1*H_m0*rk15[i] + H_1*de*rk07[i] - H_2*H_m0*rk16[i] - H_2*H_m0*rk17[i] - H_m0*de*rk14[i];
        

        

        
        Jdata[nnzper*groupj+ 4*nchem+ 9] *= Tge[i];
        
        //
        // Species: He_1
        //
        
        
        // He_1 by H2_1

        
        // because the Jacobian is initialized to zeros by default
        Jdata[nnzper*groupj+ 5*nchem+ 0] = ZERO;
        

        

        
        
        // He_1 by H2_2

        
        // because the Jacobian is initialized to zeros by default
        Jdata[nnzper*groupj+ 5*nchem+ 1] = ZERO;
        

        

        
        
        // He_1 by H_1

        
        // because the Jacobian is initialized to zeros by default
        Jdata[nnzper*groupj+ 5*nchem+ 2] = ZERO;
        

        

        
        
        // He_1 by H_2

        
        // because the Jacobian is initialized to zeros by default
        Jdata[nnzper*groupj+ 5*nchem+ 3] = ZERO;
        

        

        
        
        // He_1 by H_m0

        
        // because the Jacobian is initialized to zeros by default
        Jdata[nnzper*groupj+ 5*nchem+ 4] = ZERO;
        

        

        
        
        // He_1 by He_1

        
        Jdata[nnzper*groupj + 5*nchem+ 5] = -k03[i]*de;
        

        

        
        
        // He_1 by He_2

        
        Jdata[nnzper*groupj + 5*nchem+ 6] = k04[i]*de;
        

        

        
        
        // He_1 by He_3

        
        // because the Jacobian is initialized to zeros by default
        Jdata[nnzper*groupj+ 5*nchem+ 7] = ZERO;
        

        

        
        
        // He_1 by de

        
        Jdata[nnzper*groupj + 5*nchem+ 8] = -k03[i]*He_1 + k04[i]*He_2;
        

        

        
        
        // He_1 by ge

        
        Jdata[nnzper*groupj + 5*nchem+ 9] = -He_1*de*rk03[i] + He_2*de*rk04[i];
        

        

        
        Jdata[nnzper*groupj+ 5*nchem+ 9] *= Tge[i];
        
        //
        // Species: He_2
        //
        
        
        // He_2 by H2_1

        
        // because the Jacobian is initialized to zeros by default
        Jdata[nnzper*groupj+ 6*nchem+ 0] = ZERO;
        

        

        
        
        // He_2 by H2_2

        
        // because the Jacobian is initialized to zeros by default
        Jdata[nnzper*groupj+ 6*nchem+ 1] = ZERO;
        

        

        
        
        // He_2 by H_1

        
        // because the Jacobian is initialized to zeros by default
        Jdata[nnzper*groupj+ 6*nchem+ 2] = ZERO;
        

        

        
        
        // He_2 by H_2

        
        // because the Jacobian is initialized to zeros by default
        Jdata[nnzper*groupj+ 6*nchem+ 3] = ZERO;
        

        

        
        
        // He_2 by H_m0

        
        // because the Jacobian is initialized to zeros by default
        Jdata[nnzper*groupj+ 6*nchem+ 4] = ZERO;
        

        

        
        
        // He_2 by He_1

        
        Jdata[nnzper*groupj + 6*nchem+ 5] = k03[i]*de;
        

        

        
        
        // He_2 by He_2

        
        Jdata[nnzper*groupj + 6*nchem+ 6] = -k04[i]*de - k05[i]*de;
        

        

        
        
        // He_2 by He_3

        
        Jdata[nnzper*groupj + 6*nchem+ 7] = k06[i]*de;
        

        

        
        
        // He_2 by de

        
        Jdata[nnzper*groupj + 6*nchem+ 8] = k03[i]*He_1 - k04[i]*He_2 - k05[i]*He_2 + k06[i]*He_3;
        

        

        
        
        // He_2 by ge

        
        Jdata[nnzper*groupj + 6*nchem+ 9] = He_1*de*rk03[i] - He_2*de*rk04[i] - He_2*de*rk05[i] + He_3*de*rk06[i];
        

        

        
        Jdata[nnzper*groupj+ 6*nchem+ 9] *= Tge[i];
        
        //
        // Species: He_3
        //
        
        
        // He_3 by H2_1

        
        // because the Jacobian is initialized to zeros by default
        Jdata[nnzper*groupj+ 7*nchem+ 0] = ZERO;
        

        

        
        
        // He_3 by H2_2

        
        // because the Jacobian is initialized to zeros by default
        Jdata[nnzper*groupj+ 7*nchem+ 1] = ZERO;
        

        

        
        
        // He_3 by H_1

        
        // because the Jacobian is initialized to zeros by default
        Jdata[nnzper*groupj+ 7*nchem+ 2] = ZERO;
        

        

        
        
        // He_3 by H_2

        
        // because the Jacobian is initialized to zeros by default
        Jdata[nnzper*groupj+ 7*nchem+ 3] = ZERO;
        

        

        
        
        // He_3 by H_m0

        
        // because the Jacobian is initialized to zeros by default
        Jdata[nnzper*groupj+ 7*nchem+ 4] = ZERO;
        

        

        
        
        // He_3 by He_1

        
        // because the Jacobian is initialized to zeros by default
        Jdata[nnzper*groupj+ 7*nchem+ 5] = ZERO;
        

        

        
        
        // He_3 by He_2

        
        Jdata[nnzper*groupj + 7*nchem+ 6] = k05[i]*de;
        

        

        
        
        // He_3 by He_3

        
        Jdata[nnzper*groupj + 7*nchem+ 7] = -k06[i]*de;
        

        

        
        
        // He_3 by de

        
        Jdata[nnzper*groupj + 7*nchem+ 8] = k05[i]*He_2 - k06[i]*He_3;
        

        

        
        
        // He_3 by ge

        
        Jdata[nnzper*groupj + 7*nchem+ 9] = He_2*de*rk05[i] - He_3*de*rk06[i];
        

        

        
        Jdata[nnzper*groupj+ 7*nchem+ 9] *= Tge[i];
        
        //
        // Species: de
        //
        
        
        // de by H2_1

        
        // because the Jacobian is initialized to zeros by default
        Jdata[nnzper*groupj+ 8*nchem+ 0] = ZERO;
        

        

        
        
        // de by H2_2

        
        Jdata[nnzper*groupj + 8*nchem+ 1] = -k18[i]*de;
        

        

        
        
        // de by H_1

        
        Jdata[nnzper*groupj + 8*nchem+ 2] = k01[i]*de - k07[i]*de + k08[i]*H_m0 + k15[i]*H_m0;
        

        

        
        
        // de by H_2

        
        Jdata[nnzper*groupj + 8*nchem+ 3] = -k02[i]*de + k17[i]*H_m0;
        

        

        
        
        // de by H_m0

        
        Jdata[nnzper*groupj + 8*nchem+ 4] = k08[i]*H_1 + k14[i]*de + k15[i]*H_1 + k17[i]*H_2;
        

        

        
        
        // de by He_1

        
        Jdata[nnzper*groupj + 8*nchem+ 5] = k03[i]*de;
        

        

        
        
        // de by He_2

        
        Jdata[nnzper*groupj + 8*nchem+ 6] = -k04[i]*de + k05[i]*de;
        

        

        
        
        // de by He_3

        
        Jdata[nnzper*groupj + 8*nchem+ 7] = -k06[i]*de;
        

        

        
        
        // de by de

        
        Jdata[nnzper*groupj + 8*nchem+ 8] = k01[i]*H_1 - k02[i]*H_2 + k03[i]*He_1 - k04[i]*He_2 + k05[i]*He_2 - k06[i]*He_3 - k07[i]*H_1 + k14[i]*H_m0 - k18[i]*H2_2;
        

        

        
        
        // de by ge

        
        Jdata[nnzper*groupj + 8*nchem+ 9] = -H2_2*de*rk18[i] + H_1*H_m0*rk08[i] + H_1*H_m0*rk15[i] + H_1*de*rk01[i] - H_1*de*rk07[i] + H_2*H_m0*rk17[i] - H_2*de*rk02[i] + H_m0*de*rk14[i] + He_1*de*rk03[i] - He_2*de*rk04[i] + He_2*de*rk05[i] - He_3*de*rk06[i];
        

        

        
        Jdata[nnzper*groupj+ 8*nchem+ 9] *= Tge[i];
        
        //
        // Species: ge
        //
        
        
        // ge by H2_1

        
        Jdata[nnzper*groupj + 9*nchem+ 0] = -2.0158800000000001*cie_cooling_cieco[i]*mdensity - gloverabel08_gaH2[i]*pow(gloverabel08_h2lte[i], 2)*H2_1*h2_optical_depth_approx/(pow(gloverabel08_h2lte[i]/(gloverabel08_gaH2[i]*H2_1 + gloverabel08_gaHI[i]*H_1 + gloverabel08_gaHe[i]*He_1 + gloverabel08_gaHp[i]*H_2 + gloverabel08_gael[i]*de) + 1.0, 2)*pow(gloverabel08_gaH2[i]*H2_1 + gloverabel08_gaHI[i]*H_1 + gloverabel08_gaHe[i]*He_1 + gloverabel08_gaHp[i]*H_2 + gloverabel08_gael[i]*de, 2)) - gloverabel08_h2lte[i]*h2_optical_depth_approx/(gloverabel08_h2lte[i]/(gloverabel08_gaH2[i]*H2_1 + gloverabel08_gaHI[i]*H_1 + gloverabel08_gaHe[i]*He_1 + gloverabel08_gaHp[i]*H_2 + gloverabel08_gael[i]*de) + 1.0) - 0.5*h2formation_h2mcool[i]*H_1*1.0/(h2formation_ncrn[i]/(h2formation_ncrd1[i]*H_1 + h2formation_ncrd2[i]*H2_1) + 1.0) + 0.5*h2formation_ncrd2[i]*h2formation_ncrn[i]*pow(h2formation_ncrn[i]/(h2formation_ncrd1[i]*H_1 + h2formation_ncrd2[i]*H2_1) + 1.0, -2.0)*(-h2formation_h2mcool[i]*H2_1*H_1 + h2formation_h2mheat[i]*pow(H_1, 3))/pow(h2formation_ncrd1[i]*H_1 + h2formation_ncrd2[i]*H2_1, 2);
        

        
        Jdata[nnzper*groupj+ 9*nchem+ 0] *= inv_mdensity;
        

        
        
        // ge by H2_2

        
        // because the Jacobian is initialized to zeros by default
        Jdata[nnzper*groupj+ 9*nchem+ 1] = ZERO;
        

        
        Jdata[nnzper*groupj+ 9*nchem+ 1] *= inv_mdensity;
        

        
        
        // ge by H_1

        
        Jdata[nnzper*groupj + 9*nchem+ 2] = -ceHI_ceHI[i]*de - ciHI_ciHI[i]*de - gloverabel08_gaHI[i]*pow(gloverabel08_h2lte[i], 2)*H2_1*h2_optical_depth_approx/(pow(gloverabel08_h2lte[i]/(gloverabel08_gaH2[i]*H2_1 + gloverabel08_gaHI[i]*H_1 + gloverabel08_gaHe[i]*He_1 + gloverabel08_gaHp[i]*H_2 + gloverabel08_gael[i]*de) + 1.0, 2)*pow(gloverabel08_gaH2[i]*H2_1 + gloverabel08_gaHI[i]*H_1 + gloverabel08_gaHe[i]*He_1 + gloverabel08_gaHp[i]*H_2 + gloverabel08_gael[i]*de, 2)) + 0.5*h2formation_ncrd1[i]*h2formation_ncrn[i]*pow(h2formation_ncrn[i]/(h2formation_ncrd1[i]*H_1 + h2formation_ncrd2[i]*H2_1) + 1.0, -2.0)*(-h2formation_h2mcool[i]*H2_1*H_1 + h2formation_h2mheat[i]*pow(H_1, 3))/pow(h2formation_ncrd1[i]*H_1 + h2formation_ncrd2[i]*H2_1, 2) + 0.5*(-h2formation_h2mcool[i]*H2_1 + 3*h2formation_h2mheat[i]*pow(H_1, 2))*1.0/(h2formation_ncrn[i]/(h2formation_ncrd1[i]*H_1 + h2formation_ncrd2[i]*H2_1) + 1.0);
        

        
        Jdata[nnzper*groupj+ 9*nchem+ 2] *= inv_mdensity;
        

        
        
        // ge by H_2

        
        Jdata[nnzper*groupj + 9*nchem+ 3] = -brem_brem[i]*de - gloverabel08_gaHp[i]*pow(gloverabel08_h2lte[i], 2)*H2_1*h2_optical_depth_approx/(pow(gloverabel08_h2lte[i]/(gloverabel08_gaH2[i]*H2_1 + gloverabel08_gaHI[i]*H_1 + gloverabel08_gaHe[i]*He_1 + gloverabel08_gaHp[i]*H_2 + gloverabel08_gael[i]*de) + 1.0, 2)*pow(gloverabel08_gaH2[i]*H2_1 + gloverabel08_gaHI[i]*H_1 + gloverabel08_gaHe[i]*He_1 + gloverabel08_gaHp[i]*H_2 + gloverabel08_gael[i]*de, 2)) - reHII_reHII[i]*de;
        

        
        Jdata[nnzper*groupj+ 9*nchem+ 3] *= inv_mdensity;
        

        
        
        // ge by H_m0

        
        // because the Jacobian is initialized to zeros by default
        Jdata[nnzper*groupj+ 9*nchem+ 4] = ZERO;
        

        
        Jdata[nnzper*groupj+ 9*nchem+ 4] *= inv_mdensity;
        

        
        
        // ge by He_1

        
        Jdata[nnzper*groupj + 9*nchem+ 5] = -ciHeI_ciHeI[i]*de - gloverabel08_gaHe[i]*pow(gloverabel08_h2lte[i], 2)*H2_1*h2_optical_depth_approx/(pow(gloverabel08_h2lte[i]/(gloverabel08_gaH2[i]*H2_1 + gloverabel08_gaHI[i]*H_1 + gloverabel08_gaHe[i]*He_1 + gloverabel08_gaHp[i]*H_2 + gloverabel08_gael[i]*de) + 1.0, 2)*pow(gloverabel08_gaH2[i]*H2_1 + gloverabel08_gaHI[i]*H_1 + gloverabel08_gaHe[i]*He_1 + gloverabel08_gaHp[i]*H_2 + gloverabel08_gael[i]*de, 2));
        

        
        Jdata[nnzper*groupj+ 9*nchem+ 5] *= inv_mdensity;
        

        
        
        // ge by He_2

        
        Jdata[nnzper*groupj + 9*nchem+ 6] = -brem_brem[i]*de - ceHeII_ceHeII[i]*de - ceHeI_ceHeI[i]*pow(de, 2) - ciHeII_ciHeII[i]*de - ciHeIS_ciHeIS[i]*pow(de, 2) - reHeII1_reHeII1[i]*de - reHeII2_reHeII2[i]*de;
        

        
        Jdata[nnzper*groupj+ 9*nchem+ 6] *= inv_mdensity;
        

        
        
        // ge by He_3

        
        Jdata[nnzper*groupj + 9*nchem+ 7] = -4.0*brem_brem[i]*de - reHeIII_reHeIII[i]*de;
        

        
        Jdata[nnzper*groupj+ 9*nchem+ 7] *= inv_mdensity;
        

        
        
        // ge by de

        
        Jdata[nnzper*groupj + 9*nchem+ 8] = brem_brem[i]*(-H_2 - He_2 - 4.0*He_3) - ceHI_ceHI[i]*H_1 - ceHeII_ceHeII[i]*He_2 - 2*ceHeI_ceHeI[i]*He_2*de - ciHI_ciHI[i]*H_1 - ciHeII_ciHeII[i]*He_2 - 2*ciHeIS_ciHeIS[i]*He_2*de - ciHeI_ciHeI[i]*He_1 - compton_comp_[i]*pow(z + 1.0, 4)*(T - 2.73*z - 2.73) - gloverabel08_gael[i]*pow(gloverabel08_h2lte[i], 2)*H2_1*h2_optical_depth_approx/(pow(gloverabel08_h2lte[i]/(gloverabel08_gaH2[i]*H2_1 + gloverabel08_gaHI[i]*H_1 + gloverabel08_gaHe[i]*He_1 + gloverabel08_gaHp[i]*H_2 + gloverabel08_gael[i]*de) + 1.0, 2)*pow(gloverabel08_gaH2[i]*H2_1 + gloverabel08_gaHI[i]*H_1 + gloverabel08_gaHe[i]*He_1 + gloverabel08_gaHp[i]*H_2 + gloverabel08_gael[i]*de, 2)) - reHII_reHII[i]*H_2 - reHeII1_reHeII1[i]*He_2 - reHeII2_reHeII2[i]*He_2 - reHeIII_reHeIII[i]*He_3;
        

        
        Jdata[nnzper*groupj+ 9*nchem+ 8] *= inv_mdensity;
        

        
        
        // ge by ge

        
        Jdata[nnzper*groupj + 9*nchem+ 9] = -gloverabel08_h2lte[i]*H2_1*h2_optical_depth_approx*(-gloverabel08_h2lte[i]*(-H2_1*rgloverabel08_gaH2[i] - H_1*rgloverabel08_gaHI[i] - H_2*rgloverabel08_gaHp[i] - He_1*rgloverabel08_gaHe[i] - de*rgloverabel08_gael[i])/pow(gloverabel08_gaH2[i]*H2_1 + gloverabel08_gaHI[i]*H_1 + gloverabel08_gaHe[i]*He_1 + gloverabel08_gaHp[i]*H_2 + gloverabel08_gael[i]*de, 2) - rgloverabel08_h2lte[i]/(gloverabel08_gaH2[i]*H2_1 + gloverabel08_gaHI[i]*H_1 + gloverabel08_gaHe[i]*He_1 + gloverabel08_gaHp[i]*H_2 + gloverabel08_gael[i]*de))/pow(gloverabel08_h2lte[i]/(gloverabel08_gaH2[i]*H2_1 + gloverabel08_gaHI[i]*H_1 + gloverabel08_gaHe[i]*He_1 + gloverabel08_gaHp[i]*H_2 + gloverabel08_gael[i]*de) + 1.0, 2) - H2_1*h2_optical_depth_approx*rgloverabel08_h2lte[i]/(gloverabel08_h2lte[i]/(gloverabel08_gaH2[i]*H2_1 + gloverabel08_gaHI[i]*H_1 + gloverabel08_gaHe[i]*He_1 + gloverabel08_gaHp[i]*H_2 + gloverabel08_gael[i]*de) + 1.0) - 2.0158800000000001*H2_1*mdensity*rcie_cooling_cieco[i] - H_1*de*rceHI_ceHI[i] - H_1*de*rciHI_ciHI[i] - H_2*de*rreHII_reHII[i] - He_1*de*rciHeI_ciHeI[i] - He_2*pow(de, 2)*rceHeI_ceHeI[i] - He_2*pow(de, 2)*rciHeIS_ciHeIS[i] - He_2*de*rceHeII_ceHeII[i] - He_2*de*rciHeII_ciHeII[i] - He_2*de*rreHeII1_reHeII1[i] - He_2*de*rreHeII2_reHeII2[i] - He_3*de*rreHeIII_reHeIII[i] - de*rbrem_brem[i]*(H_2 + He_2 + 4.0*He_3) - de*rcompton_comp_[i]*pow(z + 1.0, 4)*(T - 2.73*z - 2.73) + 0.5*pow(h2formation_ncrn[i]/(h2formation_ncrd1[i]*H_1 + h2formation_ncrd2[i]*H2_1) + 1.0, -2.0)*(-h2formation_h2mcool[i]*H2_1*H_1 + h2formation_h2mheat[i]*pow(H_1, 3))*(-1.0*h2formation_ncrn[i]*(-H2_1*rh2formation_ncrd2[i] - H_1*rh2formation_ncrd1[i])/pow(h2formation_ncrd1[i]*H_1 + h2formation_ncrd2[i]*H2_1, 2) - 1.0*rh2formation_ncrn[i]/(h2formation_ncrd1[i]*H_1 + h2formation_ncrd2[i]*H2_1)) + 0.5*1.0/(h2formation_ncrn[i]/(h2formation_ncrd1[i]*H_1 + h2formation_ncrd2[i]*H2_1) + 1.0)*(-H2_1*H_1*rh2formation_h2mcool[i] + pow(H_1, 3)*rh2formation_h2mheat[i]);
        

        
        Jdata[nnzper*groupj+ 9*nchem+ 9] *= inv_mdensity;
        

        
        Jdata[nnzper*groupj+ 9*nchem+ 9] *= Tge[i];
        

    }
    /*
       if (groupj < 1){
       for (int i =0; i < 100; i++){
       printf("from %d: Jdata[%d] = %0.5g\n", groupj, i, Jdata[nnzper*groupj+i]);
       }
       printf("\n");
       }
     */
}

    __global__
static void sparse_jacobian_kernel(realtype *ydata, realtype *Jdata, primordial_cuda_data data)
{
    int groupj;

    // temporary:
    int nnzper = 64;
    int i;
    double *Tge = data.Tge;
    double z, T;
    double *h2_optical_depth_approx_arr = data.h2_optical_depth_approx;
    double *cie_optical_depth_approx_arr = data.cie_optical_depth_approx;
    groupj = blockIdx.x*blockDim.x + threadIdx.x; 

    T = 1000.0;
    z = 0.0;

    if (groupj < BATCHSIZE)
    {
        i = groupj;
        double h2_optical_depth_approx = 1.0; //h2_optical_depth_approx_arr[i];
        double cie_optical_depth_approx = 1.0; //cie_optical_depth_approx_arr[i];
        // pulled the species data
        double H2_1 = ydata[nchem*groupj+0];
        double H2_2 = ydata[nchem*groupj+1];
        double H_1 = ydata[nchem*groupj+2];
        double H_2 = ydata[nchem*groupj+3];
        double H_m0 = ydata[nchem*groupj+4];
        double He_1 = ydata[nchem*groupj+5];
        double He_2 = ydata[nchem*groupj+6];
        double He_3 = ydata[nchem*groupj+7];
        double de = ydata[nchem*groupj+8];
        double ge = ydata[nchem*groupj+9];
        double mdensity = mh * (2.0*H2_1 + 2.0*H2_2 + 1.0079400000000001*H_1 + 1.0079400000000001*H_2 + 1.0079400000000001*H_m0 + 4.0026020000000004*He_1 + 4.0026020000000004*He_2 + 4.0026020000000004*He_3);
        double inv_mdensity = 1.0/ mdensity;
        double *k01 = data.rs_k01;
        double *rk01= data.drs_k01;
        double *k02 = data.rs_k02;
        double *rk02= data.drs_k02;
        double *k03 = data.rs_k03;
        double *rk03= data.drs_k03;
        double *k04 = data.rs_k04;
        double *rk04= data.drs_k04;
        double *k05 = data.rs_k05;
        double *rk05= data.drs_k05;
        double *k06 = data.rs_k06;
        double *rk06= data.drs_k06;
        double *k07 = data.rs_k07;
        double *rk07= data.drs_k07;
        double *k08 = data.rs_k08;
        double *rk08= data.drs_k08;
        double *k09 = data.rs_k09;
        double *rk09= data.drs_k09;
        double *k10 = data.rs_k10;
        double *rk10= data.drs_k10;
        double *k11 = data.rs_k11;
        double *rk11= data.drs_k11;
        double *k12 = data.rs_k12;
        double *rk12= data.drs_k12;
        double *k13 = data.rs_k13;
        double *rk13= data.drs_k13;
        double *k14 = data.rs_k14;
        double *rk14= data.drs_k14;
        double *k15 = data.rs_k15;
        double *rk15= data.drs_k15;
        double *k16 = data.rs_k16;
        double *rk16= data.drs_k16;
        double *k17 = data.rs_k17;
        double *rk17= data.drs_k17;
        double *k18 = data.rs_k18;
        double *rk18= data.drs_k18;
        double *k19 = data.rs_k19;
        double *rk19= data.drs_k19;
        double *k21 = data.rs_k21;
        double *rk21= data.drs_k21;
        double *k22 = data.rs_k22;
        double *rk22= data.drs_k22;
        double *k23 = data.rs_k23;
        double *rk23= data.drs_k23;
        double *brem_brem = data.cs_brem_brem;
        double *rbrem_brem = data.dcs_brem_brem;
        double *ceHeI_ceHeI = data.cs_ceHeI_ceHeI;
        double *rceHeI_ceHeI = data.dcs_ceHeI_ceHeI;
        double *ceHeII_ceHeII = data.cs_ceHeII_ceHeII;
        double *rceHeII_ceHeII = data.dcs_ceHeII_ceHeII;
        double *ceHI_ceHI = data.cs_ceHI_ceHI;
        double *rceHI_ceHI = data.dcs_ceHI_ceHI;
        double *cie_cooling_cieco = data.cs_cie_cooling_cieco;
        double *rcie_cooling_cieco = data.dcs_cie_cooling_cieco;
        double *ciHeI_ciHeI = data.cs_ciHeI_ciHeI;
        double *rciHeI_ciHeI = data.dcs_ciHeI_ciHeI;
        double *ciHeII_ciHeII = data.cs_ciHeII_ciHeII;
        double *rciHeII_ciHeII = data.dcs_ciHeII_ciHeII;
        double *ciHeIS_ciHeIS = data.cs_ciHeIS_ciHeIS;
        double *rciHeIS_ciHeIS = data.dcs_ciHeIS_ciHeIS;
        double *ciHI_ciHI = data.cs_ciHI_ciHI;
        double *rciHI_ciHI = data.dcs_ciHI_ciHI;
        double *compton_comp_ = data.cs_compton_comp_;
        double *rcompton_comp_ = data.dcs_compton_comp_;
        double *gammah_gammah = data.cs_gammah_gammah;
        double *rgammah_gammah = data.dcs_gammah_gammah;
        double *gloverabel08_gael = data.cs_gloverabel08_gael;
        double *rgloverabel08_gael = data.dcs_gloverabel08_gael;
        double *gloverabel08_gaH2 = data.cs_gloverabel08_gaH2;
        double *rgloverabel08_gaH2 = data.dcs_gloverabel08_gaH2;
        double *gloverabel08_gaHe = data.cs_gloverabel08_gaHe;
        double *rgloverabel08_gaHe = data.dcs_gloverabel08_gaHe;
        double *gloverabel08_gaHI = data.cs_gloverabel08_gaHI;
        double *rgloverabel08_gaHI = data.dcs_gloverabel08_gaHI;
        double *gloverabel08_gaHp = data.cs_gloverabel08_gaHp;
        double *rgloverabel08_gaHp = data.dcs_gloverabel08_gaHp;
        double *gloverabel08_gphdl = data.cs_gloverabel08_gphdl;
        double *rgloverabel08_gphdl = data.dcs_gloverabel08_gphdl;
        double *gloverabel08_gpldl = data.cs_gloverabel08_gpldl;
        double *rgloverabel08_gpldl = data.dcs_gloverabel08_gpldl;
        double *gloverabel08_h2lte = data.cs_gloverabel08_h2lte;
        double *rgloverabel08_h2lte = data.dcs_gloverabel08_h2lte;
        double *h2formation_h2mcool = data.cs_h2formation_h2mcool;
        double *rh2formation_h2mcool = data.dcs_h2formation_h2mcool;
        double *h2formation_h2mheat = data.cs_h2formation_h2mheat;
        double *rh2formation_h2mheat = data.dcs_h2formation_h2mheat;
        double *h2formation_ncrd1 = data.cs_h2formation_ncrd1;
        double *rh2formation_ncrd1 = data.dcs_h2formation_ncrd1;
        double *h2formation_ncrd2 = data.cs_h2formation_ncrd2;
        double *rh2formation_ncrd2 = data.dcs_h2formation_ncrd2;
        double *h2formation_ncrn = data.cs_h2formation_ncrn;
        double *rh2formation_ncrn = data.dcs_h2formation_ncrn;
        double *reHeII1_reHeII1 = data.cs_reHeII1_reHeII1;
        double *rreHeII1_reHeII1 = data.dcs_reHeII1_reHeII1;
        double *reHeII2_reHeII2 = data.cs_reHeII2_reHeII2;
        double *rreHeII2_reHeII2 = data.dcs_reHeII2_reHeII2;
        double *reHeIII_reHeIII = data.cs_reHeIII_reHeIII;
        double *rreHeIII_reHeIII = data.dcs_reHeIII_reHeIII;
        double *reHII_reHII = data.cs_reHII_reHII;
        double *rreHII_reHII = data.dcs_reHII_reHII;
        // H2_1 by H2_1
        Jdata[nnzper*groupj+0] = -k11[i]*H_2 - k12[i]*de - k13[i]*H_1 + k21[i]*pow(H_1, 2) - 2*k23[i]*H2_1;

        
        // H2_1 by H2_2
        Jdata[nnzper*groupj+1] = k10[i]*H_1 + k19[i]*H_m0;

        
        // H2_1 by H_1
        Jdata[nnzper*groupj+2] = k08[i]*H_m0 + k10[i]*H2_2 - k13[i]*H2_1 + 2*k21[i]*H2_1*H_1 + 3*k22[i]*pow(H_1, 2);

        
        // H2_1 by H_2
        Jdata[nnzper*groupj+3] = -k11[i]*H2_1;

        
        // H2_1 by H_m0
        Jdata[nnzper*groupj+4] = k08[i]*H_1 + k19[i]*H2_2;

        
        // H2_1 by de
        Jdata[nnzper*groupj+5] = -k12[i]*H2_1;

        
        // H2_1 by ge
        Jdata[nnzper*groupj+6] = -pow(H2_1, 2)*rk23[i] + H2_1*pow(H_1, 2)*rk21[i] - H2_1*H_1*rk13[i] - H2_1*H_2*rk11[i] - H2_1*de*rk12[i] + H2_2*H_1*rk10[i] + H2_2*H_m0*rk19[i] + pow(H_1, 3)*rk22[i] + H_1*H_m0*rk08[i];

        
        Jdata[nnzper*groupj+6] *= Tge[i];
        // H2_2 by H2_1
        Jdata[nnzper*groupj+7] = k11[i]*H_2;

        
        // H2_2 by H2_2
        Jdata[nnzper*groupj+8] = -k10[i]*H_1 - k18[i]*de - k19[i]*H_m0;

        
        // H2_2 by H_1
        Jdata[nnzper*groupj+9] = k09[i]*H_2 - k10[i]*H2_2;

        
        // H2_2 by H_2
        Jdata[nnzper*groupj+10] = k09[i]*H_1 + k11[i]*H2_1 + k17[i]*H_m0;

        
        // H2_2 by H_m0
        Jdata[nnzper*groupj+11] = k17[i]*H_2 - k19[i]*H2_2;

        
        // H2_2 by de
        Jdata[nnzper*groupj+12] = -k18[i]*H2_2;

        
        // H2_2 by ge
        Jdata[nnzper*groupj+13] = H2_1*H_2*rk11[i] - H2_2*H_1*rk10[i] - H2_2*H_m0*rk19[i] - H2_2*de*rk18[i] + H_1*H_2*rk09[i] + H_2*H_m0*rk17[i];

        
        Jdata[nnzper*groupj+13] *= Tge[i];
        // H_1 by H2_1
        Jdata[nnzper*groupj+14] = k11[i]*H_2 + 2*k12[i]*de + 2*k13[i]*H_1 - 2*k21[i]*pow(H_1, 2) + 4*k23[i]*H2_1;

        
        // H_1 by H2_2
        Jdata[nnzper*groupj+15] = -k10[i]*H_1 + 2*k18[i]*de + k19[i]*H_m0;

        
        // H_1 by H_1
        Jdata[nnzper*groupj+16] = -k01[i]*de - k07[i]*de - k08[i]*H_m0 - k09[i]*H_2 - k10[i]*H2_2 + 2*k13[i]*H2_1 + k15[i]*H_m0 - 4*k21[i]*H2_1*H_1 - 6*k22[i]*pow(H_1, 2);

        
        // H_1 by H_2
        Jdata[nnzper*groupj+17] = k02[i]*de - k09[i]*H_1 + k11[i]*H2_1 + 2*k16[i]*H_m0;

        
        // H_1 by H_m0
        Jdata[nnzper*groupj+18] = -k08[i]*H_1 + k14[i]*de + k15[i]*H_1 + 2*k16[i]*H_2 + k19[i]*H2_2;

        
        // H_1 by de
        Jdata[nnzper*groupj+19] = -k01[i]*H_1 + k02[i]*H_2 - k07[i]*H_1 + 2*k12[i]*H2_1 + k14[i]*H_m0 + 2*k18[i]*H2_2;

        
        // H_1 by ge
        Jdata[nnzper*groupj+20] = 2*pow(H2_1, 2)*rk23[i] - 2*H2_1*pow(H_1, 2)*rk21[i] + 2*H2_1*H_1*rk13[i] + H2_1*H_2*rk11[i] + 2*H2_1*de*rk12[i] - H2_2*H_1*rk10[i] + H2_2*H_m0*rk19[i] + 2*H2_2*de*rk18[i] - 2*pow(H_1, 3)*rk22[i] - H_1*H_2*rk09[i] - H_1*H_m0*rk08[i] + H_1*H_m0*rk15[i] - H_1*de*rk01[i] - H_1*de*rk07[i] + 2*H_2*H_m0*rk16[i] + H_2*de*rk02[i] + H_m0*de*rk14[i];

        
        Jdata[nnzper*groupj+20] *= Tge[i];
        // H_2 by H2_1
        Jdata[nnzper*groupj+21] = -k11[i]*H_2;

        
        // H_2 by H2_2
        Jdata[nnzper*groupj+22] = k10[i]*H_1;

        
        // H_2 by H_1
        Jdata[nnzper*groupj+23] = k01[i]*de - k09[i]*H_2 + k10[i]*H2_2;

        
        // H_2 by H_2
        Jdata[nnzper*groupj+24] = -k02[i]*de - k09[i]*H_1 - k11[i]*H2_1 - k16[i]*H_m0 - k17[i]*H_m0;

        
        // H_2 by H_m0
        Jdata[nnzper*groupj+25] = -k16[i]*H_2 - k17[i]*H_2;

        
        // H_2 by de
        Jdata[nnzper*groupj+26] = k01[i]*H_1 - k02[i]*H_2;

        
        // H_2 by ge
        Jdata[nnzper*groupj+27] = -H2_1*H_2*rk11[i] + H2_2*H_1*rk10[i] - H_1*H_2*rk09[i] + H_1*de*rk01[i] - H_2*H_m0*rk16[i] - H_2*H_m0*rk17[i] - H_2*de*rk02[i];

        
        Jdata[nnzper*groupj+27] *= Tge[i];
        // H_m0 by H2_2
        Jdata[nnzper*groupj+28] = -k19[i]*H_m0;

        
        // H_m0 by H_1
        Jdata[nnzper*groupj+29] = k07[i]*de - k08[i]*H_m0 - k15[i]*H_m0;

        
        // H_m0 by H_2
        Jdata[nnzper*groupj+30] = -k16[i]*H_m0 - k17[i]*H_m0;

        
        // H_m0 by H_m0
        Jdata[nnzper*groupj+31] = -k08[i]*H_1 - k14[i]*de - k15[i]*H_1 - k16[i]*H_2 - k17[i]*H_2 - k19[i]*H2_2;

        
        // H_m0 by de
        Jdata[nnzper*groupj+32] = k07[i]*H_1 - k14[i]*H_m0;

        
        // H_m0 by ge
        Jdata[nnzper*groupj+33] = -H2_2*H_m0*rk19[i] - H_1*H_m0*rk08[i] - H_1*H_m0*rk15[i] + H_1*de*rk07[i] - H_2*H_m0*rk16[i] - H_2*H_m0*rk17[i] - H_m0*de*rk14[i];

        
        Jdata[nnzper*groupj+33] *= Tge[i];
        // He_1 by He_1
        Jdata[nnzper*groupj+34] = -k03[i]*de;

        
        // He_1 by He_2
        Jdata[nnzper*groupj+35] = k04[i]*de;

        
        // He_1 by de
        Jdata[nnzper*groupj+36] = -k03[i]*He_1 + k04[i]*He_2;

        
        // He_1 by ge
        Jdata[nnzper*groupj+37] = -He_1*de*rk03[i] + He_2*de*rk04[i];

        
        Jdata[nnzper*groupj+37] *= Tge[i];
        // He_2 by He_1
        Jdata[nnzper*groupj+38] = k03[i]*de;

        
        // He_2 by He_2
        Jdata[nnzper*groupj+39] = -k04[i]*de - k05[i]*de;

        
        // He_2 by He_3
        Jdata[nnzper*groupj+40] = k06[i]*de;

        
        // He_2 by de
        Jdata[nnzper*groupj+41] = k03[i]*He_1 - k04[i]*He_2 - k05[i]*He_2 + k06[i]*He_3;

        
        // He_2 by ge
        Jdata[nnzper*groupj+42] = He_1*de*rk03[i] - He_2*de*rk04[i] - He_2*de*rk05[i] + He_3*de*rk06[i];

        
        Jdata[nnzper*groupj+42] *= Tge[i];
        // He_3 by He_2
        Jdata[nnzper*groupj+43] = k05[i]*de;

        
        // He_3 by He_3
        Jdata[nnzper*groupj+44] = -k06[i]*de;

        
        // He_3 by de
        Jdata[nnzper*groupj+45] = k05[i]*He_2 - k06[i]*He_3;

        
        // He_3 by ge
        Jdata[nnzper*groupj+46] = He_2*de*rk05[i] - He_3*de*rk06[i];

        
        Jdata[nnzper*groupj+46] *= Tge[i];
        // de by H2_2
        Jdata[nnzper*groupj+47] = -k18[i]*de;

        
        // de by H_1
        Jdata[nnzper*groupj+48] = k01[i]*de - k07[i]*de + k08[i]*H_m0 + k15[i]*H_m0;

        
        // de by H_2
        Jdata[nnzper*groupj+49] = -k02[i]*de + k17[i]*H_m0;

        
        // de by H_m0
        Jdata[nnzper*groupj+50] = k08[i]*H_1 + k14[i]*de + k15[i]*H_1 + k17[i]*H_2;

        
        // de by He_1
        Jdata[nnzper*groupj+51] = k03[i]*de;

        
        // de by He_2
        Jdata[nnzper*groupj+52] = -k04[i]*de + k05[i]*de;

        
        // de by He_3
        Jdata[nnzper*groupj+53] = -k06[i]*de;

        
        // de by de
        Jdata[nnzper*groupj+54] = k01[i]*H_1 - k02[i]*H_2 + k03[i]*He_1 - k04[i]*He_2 + k05[i]*He_2 - k06[i]*He_3 - k07[i]*H_1 + k14[i]*H_m0 - k18[i]*H2_2;

        
        // de by ge
        Jdata[nnzper*groupj+55] = -H2_2*de*rk18[i] + H_1*H_m0*rk08[i] + H_1*H_m0*rk15[i] + H_1*de*rk01[i] - H_1*de*rk07[i] + H_2*H_m0*rk17[i] - H_2*de*rk02[i] + H_m0*de*rk14[i] + He_1*de*rk03[i] - He_2*de*rk04[i] + He_2*de*rk05[i] - He_3*de*rk06[i];

        
        Jdata[nnzper*groupj+55] *= Tge[i];
        // ge by H2_1
        Jdata[nnzper*groupj+56] = -2.0158800000000001*cie_cooling_cieco[i]*mdensity - gloverabel08_gaH2[i]*pow(gloverabel08_h2lte[i], 2)*H2_1*h2_optical_depth_approx/(pow(gloverabel08_h2lte[i]/(gloverabel08_gaH2[i]*H2_1 + gloverabel08_gaHI[i]*H_1 + gloverabel08_gaHe[i]*He_1 + gloverabel08_gaHp[i]*H_2 + gloverabel08_gael[i]*de) + 1.0, 2)*pow(gloverabel08_gaH2[i]*H2_1 + gloverabel08_gaHI[i]*H_1 + gloverabel08_gaHe[i]*He_1 + gloverabel08_gaHp[i]*H_2 + gloverabel08_gael[i]*de, 2)) - gloverabel08_h2lte[i]*h2_optical_depth_approx/(gloverabel08_h2lte[i]/(gloverabel08_gaH2[i]*H2_1 + gloverabel08_gaHI[i]*H_1 + gloverabel08_gaHe[i]*He_1 + gloverabel08_gaHp[i]*H_2 + gloverabel08_gael[i]*de) + 1.0) - 0.5*h2formation_h2mcool[i]*H_1*1.0/(h2formation_ncrn[i]/(h2formation_ncrd1[i]*H_1 + h2formation_ncrd2[i]*H2_1) + 1.0) + 0.5*h2formation_ncrd2[i]*h2formation_ncrn[i]*pow(h2formation_ncrn[i]/(h2formation_ncrd1[i]*H_1 + h2formation_ncrd2[i]*H2_1) + 1.0, -2.0)*(-h2formation_h2mcool[i]*H2_1*H_1 + h2formation_h2mheat[i]*pow(H_1, 3))/pow(h2formation_ncrd1[i]*H_1 + h2formation_ncrd2[i]*H2_1, 2);

        
        Jdata[nnzper*groupj+56] *= inv_mdensity;
        // ge by H_1
        Jdata[nnzper*groupj+57] = -ceHI_ceHI[i]*de - ciHI_ciHI[i]*de - gloverabel08_gaHI[i]*pow(gloverabel08_h2lte[i], 2)*H2_1*h2_optical_depth_approx/(pow(gloverabel08_h2lte[i]/(gloverabel08_gaH2[i]*H2_1 + gloverabel08_gaHI[i]*H_1 + gloverabel08_gaHe[i]*He_1 + gloverabel08_gaHp[i]*H_2 + gloverabel08_gael[i]*de) + 1.0, 2)*pow(gloverabel08_gaH2[i]*H2_1 + gloverabel08_gaHI[i]*H_1 + gloverabel08_gaHe[i]*He_1 + gloverabel08_gaHp[i]*H_2 + gloverabel08_gael[i]*de, 2)) + 0.5*h2formation_ncrd1[i]*h2formation_ncrn[i]*pow(h2formation_ncrn[i]/(h2formation_ncrd1[i]*H_1 + h2formation_ncrd2[i]*H2_1) + 1.0, -2.0)*(-h2formation_h2mcool[i]*H2_1*H_1 + h2formation_h2mheat[i]*pow(H_1, 3))/pow(h2formation_ncrd1[i]*H_1 + h2formation_ncrd2[i]*H2_1, 2) + 0.5*(-h2formation_h2mcool[i]*H2_1 + 3*h2formation_h2mheat[i]*pow(H_1, 2))*1.0/(h2formation_ncrn[i]/(h2formation_ncrd1[i]*H_1 + h2formation_ncrd2[i]*H2_1) + 1.0);

        
        Jdata[nnzper*groupj+57] *= inv_mdensity;
        // ge by H_2
        Jdata[nnzper*groupj+58] = -brem_brem[i]*de - gloverabel08_gaHp[i]*pow(gloverabel08_h2lte[i], 2)*H2_1*h2_optical_depth_approx/(pow(gloverabel08_h2lte[i]/(gloverabel08_gaH2[i]*H2_1 + gloverabel08_gaHI[i]*H_1 + gloverabel08_gaHe[i]*He_1 + gloverabel08_gaHp[i]*H_2 + gloverabel08_gael[i]*de) + 1.0, 2)*pow(gloverabel08_gaH2[i]*H2_1 + gloverabel08_gaHI[i]*H_1 + gloverabel08_gaHe[i]*He_1 + gloverabel08_gaHp[i]*H_2 + gloverabel08_gael[i]*de, 2)) - reHII_reHII[i]*de;

        
        Jdata[nnzper*groupj+58] *= inv_mdensity;
        // ge by He_1
        Jdata[nnzper*groupj+59] = -ciHeI_ciHeI[i]*de - gloverabel08_gaHe[i]*pow(gloverabel08_h2lte[i], 2)*H2_1*h2_optical_depth_approx/(pow(gloverabel08_h2lte[i]/(gloverabel08_gaH2[i]*H2_1 + gloverabel08_gaHI[i]*H_1 + gloverabel08_gaHe[i]*He_1 + gloverabel08_gaHp[i]*H_2 + gloverabel08_gael[i]*de) + 1.0, 2)*pow(gloverabel08_gaH2[i]*H2_1 + gloverabel08_gaHI[i]*H_1 + gloverabel08_gaHe[i]*He_1 + gloverabel08_gaHp[i]*H_2 + gloverabel08_gael[i]*de, 2));

        
        Jdata[nnzper*groupj+59] *= inv_mdensity;
        // ge by He_2
        Jdata[nnzper*groupj+60] = -brem_brem[i]*de - ceHeII_ceHeII[i]*de - ceHeI_ceHeI[i]*pow(de, 2) - ciHeII_ciHeII[i]*de - ciHeIS_ciHeIS[i]*pow(de, 2) - reHeII1_reHeII1[i]*de - reHeII2_reHeII2[i]*de;

        
        Jdata[nnzper*groupj+60] *= inv_mdensity;
        // ge by He_3
        Jdata[nnzper*groupj+61] = -4.0*brem_brem[i]*de - reHeIII_reHeIII[i]*de;

        
        Jdata[nnzper*groupj+61] *= inv_mdensity;
        // ge by de
        Jdata[nnzper*groupj+62] = brem_brem[i]*(-H_2 - He_2 - 4.0*He_3) - ceHI_ceHI[i]*H_1 - ceHeII_ceHeII[i]*He_2 - 2*ceHeI_ceHeI[i]*He_2*de - ciHI_ciHI[i]*H_1 - ciHeII_ciHeII[i]*He_2 - 2*ciHeIS_ciHeIS[i]*He_2*de - ciHeI_ciHeI[i]*He_1 - compton_comp_[i]*pow(z + 1.0, 4)*(T - 2.73*z - 2.73) - gloverabel08_gael[i]*pow(gloverabel08_h2lte[i], 2)*H2_1*h2_optical_depth_approx/(pow(gloverabel08_h2lte[i]/(gloverabel08_gaH2[i]*H2_1 + gloverabel08_gaHI[i]*H_1 + gloverabel08_gaHe[i]*He_1 + gloverabel08_gaHp[i]*H_2 + gloverabel08_gael[i]*de) + 1.0, 2)*pow(gloverabel08_gaH2[i]*H2_1 + gloverabel08_gaHI[i]*H_1 + gloverabel08_gaHe[i]*He_1 + gloverabel08_gaHp[i]*H_2 + gloverabel08_gael[i]*de, 2)) - reHII_reHII[i]*H_2 - reHeII1_reHeII1[i]*He_2 - reHeII2_reHeII2[i]*He_2 - reHeIII_reHeIII[i]*He_3;

        
        Jdata[nnzper*groupj+62] *= inv_mdensity;
        // ge by ge
        Jdata[nnzper*groupj+63] = -gloverabel08_h2lte[i]*H2_1*h2_optical_depth_approx*(-gloverabel08_h2lte[i]*(-H2_1*rgloverabel08_gaH2[i] - H_1*rgloverabel08_gaHI[i] - H_2*rgloverabel08_gaHp[i] - He_1*rgloverabel08_gaHe[i] - de*rgloverabel08_gael[i])/pow(gloverabel08_gaH2[i]*H2_1 + gloverabel08_gaHI[i]*H_1 + gloverabel08_gaHe[i]*He_1 + gloverabel08_gaHp[i]*H_2 + gloverabel08_gael[i]*de, 2) - rgloverabel08_h2lte[i]/(gloverabel08_gaH2[i]*H2_1 + gloverabel08_gaHI[i]*H_1 + gloverabel08_gaHe[i]*He_1 + gloverabel08_gaHp[i]*H_2 + gloverabel08_gael[i]*de))/pow(gloverabel08_h2lte[i]/(gloverabel08_gaH2[i]*H2_1 + gloverabel08_gaHI[i]*H_1 + gloverabel08_gaHe[i]*He_1 + gloverabel08_gaHp[i]*H_2 + gloverabel08_gael[i]*de) + 1.0, 2) - H2_1*h2_optical_depth_approx*rgloverabel08_h2lte[i]/(gloverabel08_h2lte[i]/(gloverabel08_gaH2[i]*H2_1 + gloverabel08_gaHI[i]*H_1 + gloverabel08_gaHe[i]*He_1 + gloverabel08_gaHp[i]*H_2 + gloverabel08_gael[i]*de) + 1.0) - 2.0158800000000001*H2_1*mdensity*rcie_cooling_cieco[i] - H_1*de*rceHI_ceHI[i] - H_1*de*rciHI_ciHI[i] - H_2*de*rreHII_reHII[i] - He_1*de*rciHeI_ciHeI[i] - He_2*pow(de, 2)*rceHeI_ceHeI[i] - He_2*pow(de, 2)*rciHeIS_ciHeIS[i] - He_2*de*rceHeII_ceHeII[i] - He_2*de*rciHeII_ciHeII[i] - He_2*de*rreHeII1_reHeII1[i] - He_2*de*rreHeII2_reHeII2[i] - He_3*de*rreHeIII_reHeIII[i] - de*rbrem_brem[i]*(H_2 + He_2 + 4.0*He_3) - de*rcompton_comp_[i]*pow(z + 1.0, 4)*(T - 2.73*z - 2.73) + 0.5*pow(h2formation_ncrn[i]/(h2formation_ncrd1[i]*H_1 + h2formation_ncrd2[i]*H2_1) + 1.0, -2.0)*(-h2formation_h2mcool[i]*H2_1*H_1 + h2formation_h2mheat[i]*pow(H_1, 3))*(-1.0*h2formation_ncrn[i]*(-H2_1*rh2formation_ncrd2[i] - H_1*rh2formation_ncrd1[i])/pow(h2formation_ncrd1[i]*H_1 + h2formation_ncrd2[i]*H2_1, 2) - 1.0*rh2formation_ncrn[i]/(h2formation_ncrd1[i]*H_1 + h2formation_ncrd2[i]*H2_1)) + 0.5*1.0/(h2formation_ncrn[i]/(h2formation_ncrd1[i]*H_1 + h2formation_ncrd2[i]*H2_1) + 1.0)*(-H2_1*H_1*rh2formation_h2mcool[i] + pow(H_1, 3)*rh2formation_h2mheat[i]);

        
        Jdata[nnzper*groupj+63] *= inv_mdensity;
        Jdata[nnzper*groupj+63] *= Tge[i];
    }
    /*
       if (groupj < 1){
       for (int i =0; i < 100; i++){
       printf("from %d: Jdata[%d] = %0.5g\n", groupj, i, Jdata[nnzper*groupj+i]);
       }
       printf("\n");
       }
     */
}

/*
 * Jacobian routine. COmpute J(t,y) = df/dy.
 * This is done on the GPU.
 */
static int calculate_jacobian_primordial_cuda(realtype t, N_Vector y, N_Vector fy, SUNMatrix J,
        void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
    primordial_cuda_data *data = (primordial_cuda_data*)user_data;

    int nnzper;
    realtype *Jdata, *ydata;
    nnzper = 10* 10;
    ydata = N_VGetDeviceArrayPointer_Cuda(y);
    Jdata = SUNMatrix_cuSparse_Data(J);

    jacobian_kernel<<<GRIDSIZE, BLOCKSIZE>>>(ydata, Jdata, *data);

    cudaDeviceSynchronize();
    cudaError_t cuerr = cudaGetLastError();
    if (cuerr != cudaSuccess) {
        fprintf(stderr, ">>> ERROR in Jac: cudaGetLastError returned %s\n",
                cudaGetErrorName(cuerr));
        return(-1);
    }

    return(0);

}

static int calculate_sparse_jacobian_primordial_cuda(realtype t, N_Vector y, N_Vector fy, SUNMatrix J,
        void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
    primordial_cuda_data *data = (primordial_cuda_data*)user_data;

    realtype *Jdata, *ydata;
    ydata = N_VGetDeviceArrayPointer_Cuda(y);
    Jdata = SUNMatrix_cuSparse_Data(J);

    sparse_jacobian_kernel<<<GRIDSIZE, BLOCKSIZE>>>(ydata, Jdata, *data);

    cudaDeviceSynchronize();
    cudaError_t cuerr = cudaGetLastError();
    if (cuerr != cudaSuccess) {
        fprintf(stderr, ">>> ERROR in Jac: cudaGetLastError returned %s\n",
                cudaGetErrorName(cuerr));
        return(-1);
    }

    return(0);

}

// now write tests kit

void test_interpolation_kernel(primordial_cuda_data data)
{
    // initialize temperature;
    for (int i = 0; i < BATCHSIZE; i++)
    {
        data.Ts[i] = (double) 3000.0 * (i+10)/BATCHSIZE;
        data.logTs[i] = log(data.Ts[i]);
    }

    float time;
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start, 0);
    for (int j = 0; j < 8192; j++){
        linear_interpolation_kernel<<<GRIDSIZE, BLOCKSIZE>>>(data);
    }
    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&time, start, stop);
    printf("Time to generate:  %3.1f ms \n", time);

    cudaDeviceSynchronize();
    cudaError_t cuerr = cudaGetLastError();
    if (cuerr != cudaSuccess) {
        fprintf(stderr,
                ">>> ERROR in interpolation_kernel: cudaGetLastError returned %s\n",
                cudaGetErrorName(cuerr));
    }
}

void initialize_ydata(double *ydata, int NSYSTEM)
{
    double h2fraction = 1.0e-5;
    double efraction  = 1.0e-5;
    double density    = 1.0e0;
    double temperature = 1000.0;
    double tiny_fraction = 1.0e-20;
    for (int i = 0; i < NSYSTEM; i++)
    {
        // H2I
        ydata[i*nchem]   = density*h2fraction*0.76 /2.;
        // H2II
        ydata[i*nchem+1] = density*tiny_fraction;
        // HI
        ydata[i*nchem+2] = density*0.76*(1-h2fraction);
        // HII
        ydata[i*nchem+3] = density*efraction;
        // H-
        ydata[i*nchem+4] = density*tiny_fraction;
        // HeI
        ydata[i*nchem+5] = density*0.24 / 4.;
        // HeII
        ydata[i*nchem+6] = density*tiny_fraction;
        // HeIII
        ydata[i*nchem+7] = density*tiny_fraction;
        // de
        ydata[i*nchem+8] = density*efraction;
        // ge
        ydata[i*nchem+9] = 3./2.*kb* temperature / mh / density ;

    }
}


void test_temperature_kernel(primordial_cuda_data data)
{
    int neq = BATCHSIZE*nchem;

    N_Vector y = N_VNew_Cuda(neq);
    double *ydata;
    ydata = N_VGetHostArrayPointer_Cuda(y);
    initialize_ydata(ydata, BATCHSIZE);
    N_VCopyToDevice_Cuda(y);


    ydata = N_VGetDeviceArrayPointer_Cuda(y);
    temperature_kernel<<<GRIDSIZE,BLOCKSIZE>>>(ydata, data);

    for (int i = 0; i<BATCHSIZE; i++){
        printf("temperature[%d] = %0.5g\n", i, data.Ts[i]);
    }

    cudaDeviceSynchronize();
    cudaError_t cuerr = cudaGetLastError();
    if (cuerr != cudaSuccess) {
        fprintf(stderr,
                ">>> ERROR in temperature kernel: cudaGetLastError returned %s\n",
                cudaGetErrorName(cuerr));
    }
}



void test_rhs_function(primordial_cuda_data data)
{
    double t = 1.0;
    int neq = BATCHSIZE*nchem;

    N_Vector y = N_VNew_Cuda(neq);
    double *ydata;
    ydata = N_VGetHostArrayPointer_Cuda(y);
    initialize_ydata(ydata, BATCHSIZE);
    N_VCopyToDevice_Cuda(y);


    ydata = N_VGetDeviceArrayPointer_Cuda(y);
    N_Vector ydot = N_VNew_Cuda(neq);

    calculate_rhs_primordial_cuda(t, y, ydot, &data);
    //f(realtype t, N_Vector y, N_Vector ydot, void *user_data)
}


void test_jacobian_function(primordial_cuda_data data)
{
    double t = 1.0;
    int neq = BATCHSIZE*nchem;

    N_Vector y = N_VNew_Cuda(neq);
    double *ydata;
    ydata = N_VGetHostArrayPointer_Cuda(y);
    initialize_ydata(ydata, BATCHSIZE);
    N_VCopyToDevice_Cuda(y);

    ydata = N_VGetDeviceArrayPointer_Cuda(y);
    N_Vector ydot = N_VNew_Cuda(neq);

    // also need to initialize jacobian data space

    /* Create sparse SUNMatrix for use in linear solves */
    SUNMatrix A;
    A = NULL;
    cusparseHandle_t cusp_handle;
    cusparseCreate(&cusp_handle);
    A = SUNMatrix_cuSparse_NewBlockCSR(BATCHSIZE, nchem, nchem, nchem*nchem, cusp_handle);
    /* Initialiize the Jacobian with its fixed sparsity pattern */
    blockJacInit(A);

    calculate_jacobian_primordial_cuda(t, y, y, A, &data, y, y, y);
    //f(realtype t, N_Vector y, N_Vector ydot, void *user_data)
}

/*
 * Private Helper Function
 * Get and print some final statistics
 */

/*
 * Check function return value...
 *   opt == 0 means SUNDIALS function allocates memory so check if
 *            returned NULL pointer
 *   opt == 1 means SUNDIALS function returns an integer value so check if
 *            retval < 0
 *   opt == 2 means function allocates memory so check if returned
 *            NULL pointer
 */

static int check_retval(void *returnvalue, const char *funcname, int opt)
{
    int *retval;

    /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
    if (opt == 0 && returnvalue == NULL) {
        fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
                funcname);
        return(1); }

    /* Check if retval < 0 */
    else if (opt == 1) {
        retval = (int *) returnvalue;
        if (*retval < 0) {
            fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with retval = %d\n\n",
                    funcname, *retval);
            return(1); }}

    /* Check if function returned NULL pointer - no memory allocated */
    else if (opt == 2 && returnvalue == NULL) {
        fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
                funcname);
        return(1); }

    return(0);
}

static void PrintFinalStats(void *cvode_mem, SUNLinearSolver LS)
{
    long int nst, nfe, nsetups, nje, nni, ncfn, netf, nge;
    size_t cuSpInternalSize, cuSpWorkSize;
    int retval;

    retval = CVodeGetNumSteps(cvode_mem, &nst);
    check_retval(&retval, "CVodeGetNumSteps", 1);
    retval = CVodeGetNumRhsEvals(cvode_mem, &nfe);
    check_retval(&retval, "CVodeGetNumRhsEvals", 1);
    retval = CVodeGetNumLinSolvSetups(cvode_mem, &nsetups);
    check_retval(&retval, "CVodeGetNumLinSolvSetups", 1);
    retval = CVodeGetNumErrTestFails(cvode_mem, &netf);
    check_retval(&retval, "CVodeGetNumErrTestFails", 1);
    retval = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
    check_retval(&retval, "CVodeGetNumNonlinSolvIters", 1);
    retval = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
    check_retval(&retval, "CVodeGetNumNonlinSolvConvFails", 1);

    retval = CVodeGetNumJacEvals(cvode_mem, &nje);
    check_retval(&retval, "CVodeGetNumJacEvals", 1);

    retval = CVodeGetNumGEvals(cvode_mem, &nge);
    check_retval(&retval, "CVodeGetNumGEvals", 1);

    SUNLinSol_cuSolverSp_batchQR_GetDeviceSpace(LS, &cuSpInternalSize, &cuSpWorkSize);

    printf("\nFinal Statistics:\n");
    printf("nst = %-6ld nfe  = %-6ld nsetups = %-6ld nje = %ld\n",
            nst, nfe, nsetups, nje);
    printf("nni = %-6ld ncfn = %-6ld netf = %-6ld    nge = %ld\n \n",
            nni, ncfn, netf, nge);
    printf("cuSolverSp numerical factorization workspace size (in bytes) = %ld\n", cuSpWorkSize);
    printf("cuSolverSp internal Q, R buffer size (in bytes) = %ld\n", cuSpInternalSize);
}


int run_solver(int argc, char *argv[])
{
    realtype reltol, t, tout;
    realtype *ydata, *abstol_data;
    N_Vector y, abstol;
    SUNMatrix A;
    SUNLinearSolver LS;
    void *cvode_mem;
    int retval, iout;
    int neq, ngroups, groupj;
    primordial_cuda_data data = primordial_cuda_setup_data(NULL, NULL);
    primordial_cuda_read_cooling_tables( &data);
    primordial_cuda_read_rate_tables( &data);

    cusparseHandle_t cusp_handle;
    cusolverSpHandle_t cusol_handle;

    y = abstol = NULL;
    A = NULL;
    LS = NULL;
    cvode_mem = NULL;

    /* Parse command line arguments */
    ngroups = BATCHSIZE;
    neq = ngroups* nchem;

    reltol = 1.0e-6;
    /* Initialize cuSOLVER and cuSPARSE handles */
    cusparseCreate(&cusp_handle);
    cusolverSpCreate(&cusol_handle);

    /* Create CUDA vector of length neq for I.C. and abstol */
    y = N_VNew_Cuda(neq);
    if (check_retval((void *)y, "N_VNew_Cuda", 0)) return(1);
    abstol = N_VNew_Cuda(neq);
    if (check_retval((void *)abstol, "N_VNew_Cuda", 0)) return(1);

    ydata = N_VGetHostArrayPointer_Cuda(y);
    abstol_data = N_VGetHostArrayPointer_Cuda(abstol);

    /* Initialize */
    initialize_ydata(ydata, BATCHSIZE);
    for (int i = 0; i < neq; i++){
        abstol_data[i] = 1.0e-25;
    }
    N_VCopyToDevice_Cuda(y);
    N_VCopyToDevice_Cuda(abstol);

    /* Call CVodeCreate to create the solver memory and specify the
     * Backward Differentiation Formula */
    cvode_mem = CVodeCreate(CV_BDF);
    if (check_retval((void *)cvode_mem, "CVodeCreate", 0)) return(1);

    /* Call CVodeInit to initialize the integrator memory and specify the
     * user's right hand side function in y'=f(t,y), the inital time T0, and
     * the initial dependent variable vector y. */
    retval = CVodeInit(cvode_mem, calculate_rhs_primordial_cuda, T0, y);
    if (check_retval(&retval, "CVodeInit", 1)) return(1);

    /* Call CVodeSetUserData to attach the user data structure */
    retval = CVodeSetUserData(cvode_mem, &data);
    if (check_retval(&retval, "CVodeSetUserData", 1)) return(1);

    /* Call CVodeSVtolerances to specify the scalar relative tolerance
     * and vector absolute tolerances */
    retval = CVodeSVtolerances(cvode_mem, reltol, abstol);
    if (check_retval(&retval, "CVodeSVtolerances", 1)) return(1);

    /* Create sparse SUNMatrix for use in linear solves */
    A = SUNMatrix_cuSparse_NewBlockCSR(ngroups, nchem, nchem, nchem*nchem, cusp_handle);
    if(check_retval((void *)A, "SUNMatrix_cuSparse_NewBlockCSR", 0)) return(1);

    /* Set the sparsity pattern to be fixed so that the row pointers
     * and column indicies are not zeroed out by SUNMatZero */
    retval = SUNMatrix_cuSparse_SetFixedPattern(A, 1);

    /* Initialiize the Jacobian with its fixed sparsity pattern */
    blockJacInit(A);

    /* Create the SUNLinearSolver object for use by CVode */
    LS = SUNLinSol_cuSolverSp_batchQR(y, A, cusol_handle);
    if(check_retval((void *)LS, "SUNLinSol_cuSolverSp_batchQR", 0)) return(1);

    /* Call CVodeSetLinearSolver to attach the matrix and linear solver to CVode */
    retval = CVodeSetLinearSolver(cvode_mem, LS, A);
    if(check_retval(&retval, "CVodeSetLinearSolver", 1)) return(1);

    /* Set the user-supplied Jacobian routine Jac */
    retval = CVodeSetJacFn(cvode_mem, calculate_jacobian_primordial_cuda);
    if(check_retval(&retval, "CVodeSetJacFn", 1)) return(1);

    /* In loop, call CVode, print results, and test for error.
       Break out of loop when NOUT preset output times have been reached.  */
    printf(" \nGroup of independent 3-species kinetics problems\n\n");
    printf("number of groups = %d\n\n", ngroups);


    iout = 0;  tout = 1.0e13;

    retval = CVode(cvode_mem, tout, y, &t, CV_NORMAL);
    N_VCopyFromDevice_Cuda(y);
    for (groupj = 0; groupj < ngroups; groupj += 256) {
        printf("group %d: ", groupj);
        for (int i = 0; i < nchem; i++){
            printf("ydata[%d] = %0.5g\n", nchem*groupj+i, ydata[nchem*groupj+i]);
        }
        printf("\n");
    }

    /*
       if (check_retval(&retval, "CVode", 1)) break;
       if (retval == CV_SUCCESS) {
       iout++;
       tout *= TMULT;
       }

       if (iout == NOUT) break;
     */

    /* Print some final statistics */
    PrintFinalStats(cvode_mem, LS);

    /* Free y and abstol vectors */
    N_VDestroy(y);
    N_VDestroy(abstol);

    /* Free integrator memory */
    CVodeFree(&cvode_mem);

    /* Free the linear solver memory */
    SUNLinSolFree(LS);

    /* Free the matrix memory */
    SUNMatDestroy(A);

    /* Destroy the cuSOLVER and cuSPARSE handles */
    cusparseDestroy(cusp_handle);
    cusolverSpDestroy(cusol_handle);

    return(0);
}
void *setup_cvode_cuda_solver(CVRhsFn f, CVLsJacFn Jac, int NEQ, primordial_cuda_data *data, SUNLinearSolver LS, SUNMatrix A, N_Vector y, double reltol, N_Vector abstol, cusparseHandle_t *cusp_handle, cusolverSpHandle_t *cusol_handle)
{
    int retval, iout;
    void *cvode_mem;
    cvode_mem = NULL;
    /* Call CVodeCreate to create the solver memory and specify the
     * Backward Differentiation Formula */
    cvode_mem = CVodeCreate(CV_BDF);
    //if (check_retval((void *)cvode_mem, "CVodeCreate", 0)) return(1);

    /* Call CVodeInit to initialize the integrator memory and specify the
     * user's right hand side function in y'=f(t,y), the inital time T0, and
     * the initial dependent variable vector y. */
    retval = CVodeInit(cvode_mem, f, 0.0, y);
    //if (check_retval(&retval, "CVodeInit", 1)) return(1);

    /* Call CVodeSetUserData to attach the user data structure */
    retval = CVodeSetUserData(cvode_mem, data);
    //if (check_retval(&retval, "CVodeSetUserData", 1)) return(1);

    /* Call CVodeSVtolerances to specify the scalar relative tolerance
     * and vector absolute tolerances */
    retval = CVodeSVtolerances(cvode_mem, reltol, abstol);
    //if (check_retval(&retval, "CVodeSVtolerances", 1)) return(1);

    /* Create sparse SUNMatrix for use in linear solves */
    A = SUNMatrix_cuSparse_NewBlockCSR(BATCHSIZE, nchem, nchem, nchem*nchem, *cusp_handle);
    //if(check_retval((void *)A, "SUNMatrix_cuSparse_NewBlockCSR", 0)) return(1);

    /* Set the sparsity pattern to be fixed so that the row pointers
     * and column indicies are not zeroed out by SUNMatZero */
    retval = SUNMatrix_cuSparse_SetFixedPattern(A, 1);

    /* Initialiize the Jacobian with its fixed sparsity pattern */
    blockJacInit(A);

    /* Create the SUNLinearSolver object for use by CVode */
    LS = SUNLinSol_cuSolverSp_batchQR(y, A, *cusol_handle);
    //if(check_retval((void *)LS, "SUNLinSol_cuSolverSp_batchQR", 0)) return(1);

    /* Call CVodeSetLinearSolver to attach the matrix and linear solver to CVode */
    retval = CVodeSetLinearSolver(cvode_mem, LS, A);
    //if(check_retval(&retval, "CVodeSetLinearSolver", 1)) return(1);

    /* Set the user-supplied Jacobian routine Jac */
    retval = CVodeSetJacFn(cvode_mem, calculate_jacobian_primordial_cuda);
    //if(check_retval(&retval, "CVodeSetJacFn", 1)) return(1);

    return cvode_mem;
}

int grid_performance()
{
    unsigned long dims = 16384;
    double density, temperature, h2fraction, efraction;
    double d, T;
    clock_t start, end;
    double cpu_time_used;

    cudaEvent_t startCuda, stopCuda;
    cudaEventCreate(&startCuda);
    cudaEventCreate(&stopCuda);
    float milliseconds = 0;

    density = 1.0;
    temperature = 10.0;
    h2fraction = 1e-4;
    efraction = 1e-4;

    const int nd = 9;
    const int nT = 9;
    double *timelapsed = new double[nd*nT];
    double *darray      = new double[nd*nT];
    double *Tarray      = new double[nd*nT];
    double *output      = new double[nchem*nd*nT];

    for (int i = 0; i < nd; i++)
    {
        for (int j = 0; j < nT; j++)
        {
            d = density* pow(10., i);
            T = temperature* pow(2.2, j);

            start = clock();
            cudaEventRecord(startCuda);
            run_dengo_struct(d, T, h2fraction, efraction, dims, &output[(i*nT+j)*nchem]);
            
            cudaEventRecord(stopCuda);
            end = clock();
            cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

            printf("d=%0.5g; T=%0.5g;  %f seconds %lu cells; %0.5g percell \n", d, T, cpu_time_used, dims, cpu_time_used/ dims); 
            cudaEventSynchronize(stopCuda);
            cudaEventElapsedTime(&milliseconds, startCuda, stopCuda);

            //printf("took %f milliseconds to execute %lu; %0.5g percell \n", milliseconds, dims, milliseconds/ dims); 

            timelapsed[i*nT+j] = milliseconds*1e-3;
            darray    [i*nT+j] = d;
            Tarray    [i*nT+j] = T;


        }
    }

    // create a file
    hid_t file_id, dataset;
    hid_t datatype, dataspace;
    hsize_t dimsf[1], dimsO[1];
    herr_t status;
    dimsf[0] = nd*nT;
    dimsO[0] = nchem*nd*nT;

    file_id = H5Fcreate("performance.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    H5LTmake_dataset_double(file_id, "/time", 1, dimsf, timelapsed);
    H5LTmake_dataset_double(file_id, "/density", 1, dimsf, darray);
    H5LTmake_dataset_double(file_id, "/temperature", 1, dimsf, Tarray);
    H5LTmake_dataset_double(file_id, "/output", 1, dimsO, output);

    H5Fclose(file_id);

}

int test_scaling_dims()
{
    unsigned long dims = 4096;
    double density, temperature, h2fraction, efraction;
    clock_t start, end;
    double cpu_time_used;

    cudaEvent_t startCuda, stopCuda;
    cudaEventCreate(&startCuda);
    cudaEventCreate(&stopCuda);
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
        start = clock();
        cudaEventRecord(startCuda);

        run_dengo_struct(density, temperature, h2fraction, efraction, dims, &output[i*nchem]);
        
        dims *= 2;
        cudaEventRecord(stopCuda);
        end = clock();
        cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

        printf("took %f seconds to execute %lu \n", cpu_time_used, dims); 
        cudaEventSynchronize(stopCuda);

        cudaEventElapsedTime(&milliseconds, startCuda, stopCuda);

        printf("took %f milliseconds to execute %lu; %0.5g percell \n", milliseconds, dims, milliseconds/ dims); 
        // measured the time lapsed
        // and the cells needed
        ncells[i] = dims;
        timelapsed[i] = milliseconds*1e-3;
    }

    // create a file
    hid_t file_id, dataset;
    hid_t datatype, dataspace;
    hsize_t dimsf[1], dimsO[1];
    herr_t status;
    dimsf[0] = ntimes;
    dimsO[0] = ntimes*nchem;

    file_id = H5Fcreate("scaling.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    H5LTmake_dataset_double(file_id, "/time", 1, dimsf, timelapsed);
    H5LTmake_dataset_long(file_id, "/ncells", 1, dimsf, ncells);
    H5LTmake_dataset_double(file_id, "/output", 1, dimsO, output);
    H5Fclose(file_id);

}


int run_dengo_struct(double density, double temperature, double h2fraction, double efraction, unsigned long dims, double *output)
{
    // dengo solver interface
    DengoSolver s;

    double *input, *rtol, *atol; 
    input = new double[dims*nchem];
    rtol  = new double[1];
    atol  = new double[dims*nchem];

    initialize_long_ydata(input, dims, density, temperature, h2fraction, efraction);
    rtol[0] = 1.0e-5;

    for (int i = 0; i<nchem; i++) printf("input[%d] = %0.5g\n", i, input[i]);

    // main evolution;
    double dtf = pow(6.67e-8*mh*density, -0.5);
    double z   = 0.0;

    s.EvolveChemistry(dtf, z, input, rtol, atol, dims);

    for (int i = 0; i<nchem; i++) printf("output[%d] = %0.5g\n", i, input[i]);

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





int dengo_evolve_primordial_cuda (double dtf, double &dt, double z, double *input,
        double *rtol, double *atol, unsigned long dims, primordial_cuda_data *data, double *temp_array ){

    //-----------------------------------------------------
    // Function     : dengo_evolve_primordial_cuda
    // Description  : Main ODE solver function in dengo

    // Parameter    :   dtf     : Desired time to be reached by the solver
    //                  dt      : Pointer to the actual time reached by the solver
    //                  z       : Current redshift
    //                  input   : Array to store the initial value of each species, 
    //                            it will be updated with the value at time dt
    //                  rtol    : relative tolerance required for convergenece in each internal CVODE timesteps
    //                  atol    : absolute tolerance required for convergence in each interanl CVODE timesteps
    //                  dims    : dimension of the input array, i.e. no. of species * no. of cells
    //                  data    : primordial_cuda_data object that relay the reaction/cooling rates, and normalizations 
    //                  temp_array: temperature of each cell by the end of the evolution
    //                           
    //-----------------------------------------------------

    unsigned long i, j;
    int N = 10;
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
    for (i = 0; i < dims; i++) {
        j = i * N;
        input[j] /= 2.0; // H2_1
    j++;
        input[j] /= 2.0; // H2_2
    j++;
        input[j] /= 1.00794; // H_1
    j++;
        input[j] /= 1.00794; // H_2
    j++;
        input[j] /= 1.00794; // H_m0
    j++;
        input[j] /= 4.002602; // He_1
    j++;
        input[j] /= 4.002602; // He_2
    j++;
        input[j] /= 4.002602; // He_3
    j++;
        input[j] /= 1.0; // de
    j++;
    j++;
    }

    CVRhsFn f    = calculate_rhs_primordial_cuda;
    CVLsJacFn jf = calculate_jacobian_primordial_cuda;

    if (dt < 0) dt = dtf / 1e0;
    data->current_z = z;
    int niter = 0;
    int siter = 0;

    double floor_value = 1e-20;

    // prepare memory space on device
    realtype *ydata, *abstol_data;
    N_Vector y, abstol;
    SUNMatrix A;
    SUNLinearSolver LS;
    void *cvode_mem;
    int retval, iout;
    int neq, ngroups, groupj;
    y = abstol = NULL;
    A = NULL;
    LS = NULL;
    cvode_mem = NULL;

    neq = BATCHSIZE*nchem;

    // Initialize cuSolver and suSPARSE
    cusparseHandle_t   cusp_handle;
    cusolverSpHandle_t cusol_handle;
    cusparseCreate  (&cusp_handle);
    cusolverSpCreate(&cusol_handle);

    /////////////////////////////////////////////////////////////////
    // Initialize memory spaces for cvode solver //
    /////////////////////////////////////////////////////////////////
    // Create CUDA vector of length neq for I.C. and abstol
    y = N_VNew_Cuda(neq);
    if (check_retval ((void *)y, "N_VNew_Cuda", 0)) return(1);
    abstol = N_VNew_Cuda(neq);
    if (check_retval ((void *)abstol, "N_VNew_Cuda", 0)) return(1);

    ydata = N_VGetHostArrayPointer_Cuda(y);
    abstol_data = N_VGetHostArrayPointer_Cuda(abstol);

    // Initialize cuda vector
    for (unsigned long i = 0; i < neq; i++){
        ydata[i] = input[i];
        abstol_data[i] = input[i]*1.0e-6;
    }
    N_VCopyToDevice_Cuda(y);
    N_VCopyToDevice_Cuda(abstol);

    /* Create the SUNLinearSolver object for use by CVode */
    // LS = SUNLinSol_cuSolverSp_batchQR(y, A, cusol_handle);
    // if(check_retval((void *)LS, "SUNLinSol_cuSolverSp_batchQR", 0)) return(1);
    cvode_mem = setup_cvode_cuda_solver(f, jf, neq, data, LS, A, y, rtol[0], abstol, &cusp_handle, &cusol_handle);

    //////////////////////////////////////////////////////////////////
    // Main Evolution /////////
    /////////////////////////////////////////////////////////////////

    unsigned long ntimes = dims/ BATCHSIZE;
    unsigned long eqlast = dims % BATCHSIZE;
    unsigned long count = 0;

    realtype t;
    realtype reltol = rtol[0];
    // split the input by batchsize, 
    for (count = 0; count < ntimes; count++)
    {
        // update the yvector, and abstol
        for (int i = 0; i < neq; i++){
            ydata      [i] = input              [count*neq +i];
            abstol_data[i] = reltol*reltol*input[count*neq +i];
        }
        N_VCopyToDevice_Cuda(y);
        N_VCopyToDevice_Cuda(abstol);

        retval = CVodeReInit(cvode_mem, 0.0, y);
        retval = CVode(cvode_mem, dtf, y, &t, CV_NORMAL);
        N_VCopyFromDevice_Cuda(y);

        // copy the CVode output from cuda vector
        for (int i = 0; i < neq; i++){
            input[count*neq +i] = ydata[i];
        }
    }


    return 0;
}


void initialize_long_ydata(double *ydata, unsigned long NSYSTEM, double density, double temperature, double h2fraction, double efraction)
{
    double tiny_fraction = 1.0e-20;
    for (unsigned long i = 0; i < NSYSTEM; i++)
    {
        // H2I
        ydata[i*nchem]   = density*h2fraction*0.76 /2.;
        // H2II
        ydata[i*nchem+1] = density*tiny_fraction;
        // HI
        ydata[i*nchem+2] = density*0.76*(1-h2fraction);
        // HII
        ydata[i*nchem+3] = density*efraction;
        // H-
        ydata[i*nchem+4] = density*tiny_fraction;
        // HeI
        ydata[i*nchem+5] = density*0.24 / 4.;
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

int run_dengo_solver(double density, double temperature, double h2fraction, double efraction, unsigned long dims)
{
    primordial_cuda_data data = primordial_cuda_setup_data(NULL, NULL);
    primordial_cuda_read_cooling_tables( &data);
    primordial_cuda_read_rate_tables( &data);
    // set a final runtime
    double dtf = pow(6.67e-8*mh*density, -0.5);
    double z   = 0.0;

    // prepare input data
    double *input = (double *) malloc(sizeof(double)*dims*nchem);
    double *atol  = (double *) malloc(sizeof(double)*dims*nchem);
    double rtol  = 1.0e-5;
    double *temp;
    double dt;
    // initialize ydata

    initialize_long_ydata(input, dims, density, temperature, h2fraction, efraction);
    for (int i = 0; i<nchem; i++)
    {
        printf("input[%d] = %0.5g\n", i, input[i]);
    }
    dengo_evolve_primordial_cuda(dtf, dt, z, input, &rtol, atol, dims, &data, temp);

    for (int i = 0; i<nchem; i++)
    {
        printf("output[%d] = %0.5g\n", i, input[i]);
    }

    // supposingly this would return the output in "input" array
    double diff;
    for (unsigned long i = 0; i<dims; i++)
    {
        diff = (input[i] - input[i%nchem]) / input[i%nchem];
        if (fabs(diff) > rtol){
            printf("outputi %lu diff[%d] = %0.5g; ref = %0.5g; out = %0.5g\n", i,  i%nchem, diff, input[i], input[i%nchem]);
        }
    }

    return 0;


}

// helper function to lookfor the best configuration for kernels


// cuda performance tuner
// https://developer.nvidia.com/blog/cuda-pro-tip-occupancy-api-simplifies-launch-configuration/
// the idea now is to benchmark the best ratios for rhs_kernel and jacobian_kernel, and interpolation kernel

void launchInterpolationKernel(primordial_cuda_data *data)
{

    // initialize temperature;
    for (int i = 0; i < BATCHSIZE; i++)
    {
        data->Ts[i] = (double) 3000.0 * (i+10)/BATCHSIZE;
        data->logTs[i] = log(data->Ts[i]);
    }

    int blockSize;
    int minGridSize;
    int gridSize;

    // heuristically calculatate a block size that achieves the maximum multiprocessor level occupancy
    cudaOccupancyMaxPotentialBlockSize(&minGridSize, &blockSize, linear_interpolation_kernel, 0,0);

    // roundup according to BATCHSIZE

    gridSize = (BATCHSIZE + blockSize -1)/blockSize;
    printf("gridSize = %d; blockSize = %d\n", gridSize, blockSize);
    linear_interpolation_kernel<<<gridSize, blockSize>>>(*data);

    cudaDeviceSynchronize();
    // calculate theoretical occupancy
    int maxActiveBlocks;
    cudaOccupancyMaxActiveBlocksPerMultiprocessor( &maxActiveBlocks, 
            linear_interpolation_kernel, blockSize, 
            0);

    int device;
    cudaDeviceProp props;
    cudaGetDevice(&device);
    cudaGetDeviceProperties(&props, device);

    float occupancy = (maxActiveBlocks * blockSize / props.warpSize) / 
        (float)(props.maxThreadsPerMultiProcessor / 
                props.warpSize);

    printf("maxActiveBlocks: %d; blockSize: %d; props.warpSize: %d\n", maxActiveBlocks, blockSize, props.warpSize);
    printf("MaxThreadsperSM: %d\n", props.maxThreadsPerMultiProcessor);
    printf("Launched blocks of size %d. Theoretical occupancy: %f\n", 
            blockSize, occupancy);


    /*
       float time;
       cudaEvent_t start, stop;
       cudaEventCreate(&start);
       cudaEventCreate(&stop);
       cudaEventRecord(start, 0);
       for (int j = 0; j < 8192; j++){
       linear_interpolation_kernel<<<GRIDSIZE, BLOCKSIZE>>>(data);
       }
       cudaEventRecord(stop, 0);
       cudaEventSynchronize(stop);
       cudaEventElapsedTime(&time, start, stop);
       printf("Time to generate:  %3.1f ms \n", time);

       cudaDeviceSynchronize();
       cudaError_t cuerr = cudaGetLastError();
       if (cuerr != cudaSuccess) {
       fprintf(stderr,
       ">>> ERROR in interpolation_kernel: cudaGetLastError returned %s\n",
       cudaGetErrorName(cuerr));
       }
     */


}


void launchTemperatureKernel(primordial_cuda_data *data)
{

    int neq = BATCHSIZE*nchem;

    N_Vector y = N_VNew_Cuda(neq);
    double *ydata;
    ydata = N_VGetHostArrayPointer_Cuda(y);
    initialize_ydata(ydata, BATCHSIZE);
    N_VCopyToDevice_Cuda(y);


    ydata = N_VGetDeviceArrayPointer_Cuda(y);

    // launch temp kernel
    int blockSize;
    int minGridSize;
    int gridSize;

    // heuristically calculatate a block size that achieves the maximum multiprocessor level occupancy
    cudaOccupancyMaxPotentialBlockSize(&minGridSize, &blockSize, temperature_kernel, 0,0);

    // roundup according to BATCHSIZE
    gridSize = (BATCHSIZE + blockSize -1)/blockSize;
    printf("gridSize = %d; blockSize = %d\n", gridSize, blockSize);
    temperature_kernel<<<gridSize, blockSize>>>(ydata, *data);
    cudaDeviceSynchronize();
    // calculate theoretical occupancy
    int maxActiveBlocks;
    cudaOccupancyMaxActiveBlocksPerMultiprocessor( &maxActiveBlocks, 
            temperature_kernel, blockSize, 
            0);
    int device;
    cudaDeviceProp props;
    cudaGetDevice(&device);
    cudaGetDeviceProperties(&props, device);
    float occupancy = (maxActiveBlocks * blockSize / props.warpSize) / 
        (float)(props.maxThreadsPerMultiProcessor / 
                props.warpSize);

    printf("maxActiveBlocks: %d; blockSize: %d; props.warpSize: %d\n", maxActiveBlocks, blockSize, props.warpSize);
    printf("MaxThreadsperSM: %d\n", props.maxThreadsPerMultiProcessor);
    printf("Launched temperature Kernel blocks of size %d. Theoretical occupancy: %f\n", 
            blockSize, occupancy);

    // end launch temp kernel

}

void launchRhsKernel(primordial_cuda_data *data)
{

    double t = 1.0;
    int neq = BATCHSIZE*nchem;
    printf("done reading data %d\n", neq);

    N_Vector y    = N_VNew_Cuda(neq);
    N_Vector ydot = N_VNew_Cuda(neq);
    double *ydata, *ydotdata;
    ydata    = N_VGetHostArrayPointer_Cuda(y);
    ydotdata = N_VGetHostArrayPointer_Cuda(ydot);
    initialize_ydata(ydata, BATCHSIZE);
    initialize_ydata(ydotdata, BATCHSIZE);
    N_VCopyToDevice_Cuda(y);
    N_VCopyToDevice_Cuda(ydot);

    double *ydata_dev    = N_VGetDeviceArrayPointer_Cuda(y);
    double *ydotdata_dev = N_VGetDeviceArrayPointer_Cuda(ydot);
    // launch temp kernel
    int blockSize;
    int minGridSize;
    int gridSize;

    // heuristically calculatate a block size that achieves the maximum multiprocessor level occupancy
    cudaOccupancyMaxPotentialBlockSize(&minGridSize, &blockSize, rhs_kernel, 0,0);

    // roundup according to BATCHSIZE
    gridSize = (BATCHSIZE + blockSize -1)/blockSize;
    printf("gridSize = %d; blockSize = %d\n", gridSize, blockSize);
    rhs_kernel<<<gridSize, blockSize>>>(t,ydata_dev, ydotdata_dev, *data);
    cudaDeviceSynchronize();

    cudaError_t cuerr = cudaGetLastError();
    if (cuerr != cudaSuccess) {
        fprintf(stderr,
                ">>> ERROR in f: cudaGetLastError returned %s\n",
                cudaGetErrorName(cuerr));
    }
    // calculate theoretical occupancy
    int maxActiveBlocks;
    cudaOccupancyMaxActiveBlocksPerMultiprocessor( &maxActiveBlocks, 
            rhs_kernel, blockSize, 
            0);
    int device;
    cudaDeviceProp props;
    cudaGetDevice(&device);
    cudaGetDeviceProperties(&props, device);
    float occupancy = (maxActiveBlocks * blockSize / props.warpSize) / 
        (float)(props.maxThreadsPerMultiProcessor / 
                props.warpSize);
    printf("maxActiveBlocks: %d; blockSize: %d; props.warpSize: %d\n", maxActiveBlocks, blockSize, props.warpSize);
    printf("MaxThreadsperSM: %d\n", props.maxThreadsPerMultiProcessor);

    printf("Launched Rhs Kernel blocks of size %d. Theoretical occupancy: %f\n", 
            blockSize, occupancy);
    cudaDeviceSynchronize();
}


void launchJacobianKernel(primordial_cuda_data *data)

{

    int neq = BATCHSIZE*nchem;

    N_Vector y = N_VNew_Cuda(neq);
    double *ydata, *Jdata;
    ydata = N_VGetHostArrayPointer_Cuda(y);
    initialize_ydata(ydata, BATCHSIZE);
    N_VCopyToDevice_Cuda(y);

    ydata = N_VGetDeviceArrayPointer_Cuda(y);
    N_Vector ydot = N_VNew_Cuda(neq);

    printf("done preparing data\n");
    // also need to initialize jacobian data space

    /* Create sparse SUNMatrix for use in linear solves */
    SUNMatrix A;
    A = NULL;
    cusparseHandle_t cusp_handle;
    cusparseCreate(&cusp_handle);
    A = SUNMatrix_cuSparse_NewBlockCSR(BATCHSIZE, nchem, nchem, nchem*nchem, cusp_handle);
    /* Initialiize the Jacobian with its fixed sparsity pattern */
    blockJacInit(A);

    Jdata = SUNMatrix_cuSparse_Data(A);

    printf("just before launching\n");
    // launch temp kernel
    int blockSize;
    int minGridSize;
    int gridSize;

    // heuristically calculatate a block size that achieves the maximum multiprocessor level occupancy
    cudaOccupancyMaxPotentialBlockSize(&minGridSize, &blockSize, jacobian_kernel, 0,0);

    // roundup according to BATCHSIZE
    gridSize = (BATCHSIZE + blockSize -1)/blockSize;
    printf("gridSize = %d; blockSize = %d\n", gridSize, blockSize);
    jacobian_kernel<<<gridSize, blockSize>>>(ydata, Jdata, *data);
    cudaDeviceSynchronize();
    // calculate theoretical occupancy
    int maxActiveBlocks;
    cudaOccupancyMaxActiveBlocksPerMultiprocessor( &maxActiveBlocks, 
            jacobian_kernel, blockSize, 
            0);
    cudaError_t cuerr = cudaGetLastError();
    if (cuerr != cudaSuccess) {
        fprintf(stderr, ">>> ERROR in MaxActiveBlocks: cudaGetLastError returned %s\n",
                cudaGetErrorName(cuerr));
    }
    int device;
    cudaDeviceProp props;
    cudaGetDevice(&device);
    cudaGetDeviceProperties(&props, device);
    float occupancy = (maxActiveBlocks * blockSize / props.warpSize) / 
        (float)(props.maxThreadsPerMultiProcessor / 
                props.warpSize);
    printf("maxActiveBlocks: %d; blockSize: %d; props.warpSize: %d\n", maxActiveBlocks, blockSize, props.warpSize);
    printf("MaxThreadsperSM: %d\n", props.maxThreadsPerMultiProcessor);

    printf("Launched Jacobian Kernel blocks of size %d. Theoretical occupancy: %f\n", 
            blockSize, occupancy);
}
