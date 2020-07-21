#include "hdf5.h"
#include "hdf5_hl.h"

#include <stdio.h>
#include <cvode/cvode.h>                              /* prototypes for CVODE fcts., consts.          */
#include <nvector/nvector_cuda.h>                     /* access to cuda N_Vector                      */
#include <sunmatrix/sunmatrix_sparse.h>               /* access to sparse SUNMatrix                   */
#include <sunlinsol/sunlinsol_cusolversp_batchqr.h>   /* acess to cuSolverSp batch QR SUNLinearSolver */
#include <sundials/sundials_types.h>                  /* defs. of realtype, int              */
// copy the data thing to the device

typedef struct 
{
    double nbins;
    double lb;
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

    int neq;

} predator_prey_data;


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
    H5LTread_dataset_double(file_id, "/natural_death_predator", ratedata.r_natural_death_predator);
    H5Fclose(file_id);


    // temperature related pieces
    cudaMallocManaged(&ratedata.logTs, sizeof(double)*256);
    cudaMallocManaged(&ratedata.Tdef,  sizeof(double)*256);
    cudaMallocManaged(&ratedata.dTs_ge,  sizeof(double)*256);

    for (int i = 0; i< 256; i++)
    {
        ratedata.logTs[i] = (double )i/256*1.001;
    }
    return ratedata;
}


// use global memory
__global__
void interpolate_kernel(predator_prey_data data)
{
    int j = threadIdx.x + blockDim.x* blockIdx.x;

    int k;
    double Tdef, t1;
    double *rates_out1, *rates_out2, *rates_out3;
    double *rates_in1,  *rates_in2,  *rates_in3;

    if(j< 2)
    {
        k = __float2int_rz(data.idbin*data.logTs[j] - data.lb);
        t1 = data.lb + k*data.dbin;
        Tdef = (data.logTs[j] - t1) * data.idbin;

        // printf("Tdef[%d] = %0.5g; k = %d\n", j, data.logTs[j], k);
        
        rates_in1  = data.r_exp_growth_prey;
        rates_out1 = data.rs_exp_growth_prey;
        rates_out1[j] = Tdef*rates_in1[k+1] + (-rates_in1[k]*Tdef + rates_in1[k]);
        

        rates_in2  = data.r_predation;
        rates_out2 = data.rs_predation;
        rates_out2[j] = Tdef*rates_in2[k+1] + (-rates_in2[k]*Tdef + rates_in2[k]);

        rates_in3  = data.r_natural_death_predator;
        rates_out3 = data.rs_natural_death_predator;
        rates_out3[j] = Tdef*rates_in3[k+1] + (-rates_in3[k]*Tdef + rates_in3[k]);

        /*
        printf("rates_out3[%d] = %0.5g; @ Tdef = %0.5g\n", j, rates_out3[j], Tdef);
        printf("rates_out2[%d] = %0.5g; @ Tdef = %0.5g\n", j, rates_out2[j], Tdef);
        printf("rates_out1[%d] = %0.5g; @ Tdef = %0.5g\n", j, rates_out1[j], Tdef);
        */
    }
}

__global__
static void f_kernel(double y, double* ydata, double* ydotdata,
        predator_prey_data data, int neq, int ngroups)
{
    int i = blockIdx.x* blockDim.x + threadIdx.x;
    
    // GROUPSIZE: number of equations per group;
    int GROUPSIZE = 5;
    int groupj = i*GROUPSIZE;

    // get rate pointer
    double *exp_growth_prey        = data.rs_exp_growth_prey;
    double *natural_death_predator = data.rs_natural_death_predator;
    double *predation              = data.rs_predation;

    if (i < data.neq/5)
    {
        double dead_predator = ydata[groupj];
        double dead_prey     = ydata[groupj+1];
        double ge            = ydata[groupj+2];
        double predator      = ydata[groupj+3];
        double prey          = ydata[groupj+4];

        // species: dead_predator
        ydotdata[groupj] = natural_death_predator[i]*predator;
        // species: dead_Prey
        ydotdata[groupj+1] = predation[i]*prey*predator;
        // species: ge
        ydotdata[groupj+2] = 0.0;
        // species: predator
        ydotdata[groupj+3] = -natural_death_predator[i]*predator+ 0.75*predation[i]*predator*prey;
        // species: prey
        ydotdata[groupj+4] = exp_growth_prey[i]*prey - predation[i]*predator*prey;


        /*
        printf("dead_predator[%d] = %0.5g\n", groupj, ydata[groupj]);
        printf("deal_prey    [%d] = %0.5g\n", groupj+1, ydata[groupj+1]);
        printf("ge           [%d] = %0.5g\n", groupj+2, ydata[groupj+2]);
        printf("predator     [%d] = %0.5g\n", groupj+3, ydata[groupj+3]);
        printf("prey         [%d] = %0.5g\n", groupj+4, ydata[groupj+4]);
        */
    }
}
// Function called by the solver
static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
    predator_prey_data *udata = (predator_prey_data *) user_data;
    double *ydata    = N_VGetDeviceArrayPointer_Cuda(y);
    double *ydotdata = N_VGetDeviceArrayPointer_Cuda(ydot);

    int neq = udata->neq;

    unsigned block_size = 32;
    unsigned grid_size  = (neq+block_size-1)/ block_size;


    // should calculate temperature here
    // for now, lsts just use a constant 'temperature'
    for (int i = 0; i < 2*5; i++){
        udata->logTs[i] = (double) 2.0;
    }

    interpolate_kernel<<<grid_size, block_size>>>(*udata);

    f_kernel<<<grid_size, block_size>>>(t, ydata, ydotdata, *udata, neq/ 5, 5);

    return 0;
}

int Jacobian(realtype t, 
        N_Vector y, 
        N_Vector fy, 
        SUNMatrix J,
        void *user_data,
        N_Vector tmp1,
        N_Vector tmp2,
        N_Vector tmp3)
{
    predator_prey_data *data = (predator_prey_data*)user_data;

    /* Now We set up some temporaries */
    // CSR is what we choose
    sunindextype *rowptrs = SUNSparseMatrix_IndexPointers(J);
    sunindextype *colvals = SUNSparseMatrix_IndexValues(J);
    realtype *matrix_data = SUNSparseMatrix_Data(J);
    SUNMatZero(J);


    double *Tge = data->dTs_ge;
    double *exp_growth_prey = data->rs_exp_growth_prey;
    double *rexp_growth_prey= data->drs_exp_growth_prey;
    double *natural_death_predator = data->rs_natural_death_predator;
    double *rnatural_death_predator= data->drs_natural_death_predator;
    double *predation = data->rs_predation;
    double *rpredation= data->drs_predation;
    double dead_predator;
    double dead_prey;
    double ge;
    double predator;
    double prey;
    double z;
    double T;

    double *ydata = N_VGetDeviceArrayPointer_Cuda(y);

    int groupj;
    int ngroups = data->neq/5;
    int GROUPSIZE = 5;
    int nchem = 5;
    int j;
    int nnzper = GROUPSIZE*GROUPSIZE;

    rowptrs[0] = 0;
    rowptrs = &rowptrs[1];

    for (int i = 0; i < ngroups; i++)
    {
        dead_predator = ydata[GROUPSIZE*i];
        dead_prey     = ydata[GROUPSIZE*i+1];
        ge            = ydata[GROUPSIZE*i+2];
        predator      = ydata[GROUPSIZE*i+3];
        prey          = ydata[GROUPSIZE*i+4];
        
        j         = i*nnzper;
        rowptrs[i*GROUPSIZE+0] = j+ 5;
        rowptrs[i*GROUPSIZE+1] = j+10;
        rowptrs[i*GROUPSIZE+2] = j+15;
        rowptrs[i*GROUPSIZE+3] = j+20;
        rowptrs[i*GROUPSIZE+4] = j+25;

        // first row of block:dead predator
        matrix_data[nnzper*i] = 0.0;
        matrix_data[nnzper*i+1] = 0.0;
        matrix_data[nnzper*i+2] = 0.0;
        matrix_data[nnzper*i+3] = natural_death_predator[i];
        matrix_data[nnzper*i+4] = 0.0;

        colvals[nnzper*i] = GROUPSIZE*i;
        colvals[nnzper*i+1] = GROUPSIZE*i+1;
        colvals[nnzper*i+2] = GROUPSIZE*i+2;
        colvals[nnzper*i+3] = GROUPSIZE*i+3;
        colvals[nnzper*i+4] = GROUPSIZE*i+4;
        
        // second row of block: dead prey
        matrix_data[nnzper*i+5] = 0.0;
        matrix_data[nnzper*i+6] = 0.0;
        matrix_data[nnzper*i+7] = 0.0;
        matrix_data[nnzper*i+8] = predation[i]*predator;
        matrix_data[nnzper*i+9] = predation[i]*prey;

        colvals[nnzper*i+5] = GROUPSIZE*i;
        colvals[nnzper*i+6] = GROUPSIZE*i+1;
        colvals[nnzper*i+7] = GROUPSIZE*i+2;
        colvals[nnzper*i+8] = GROUPSIZE*i+3;
        colvals[nnzper*i+9] = GROUPSIZE*i+4;

        // third row of block: ge
        matrix_data[nnzper*i+10] = 0.0;
        matrix_data[nnzper*i+11] = 0.0;
        matrix_data[nnzper*i+12] = 0.0;
        matrix_data[nnzper*i+13] = 0.0;
        matrix_data[nnzper*i+14] = 0.0;

        colvals[nnzper*i+10] = GROUPSIZE*i;
        colvals[nnzper*i+11] = GROUPSIZE*i+1;
        colvals[nnzper*i+12] = GROUPSIZE*i+2;
        colvals[nnzper*i+13] = GROUPSIZE*i+3;
        colvals[nnzper*i+14] = GROUPSIZE*i+4;


        // forth row of block: predator
        matrix_data[nnzper*i+15] = 0.0;
        matrix_data[nnzper*i+16] = 0.0;
        matrix_data[nnzper*i+17] = 0.0;
        matrix_data[nnzper*i+18] = -natural_death_predator[i] + 0.75*predation[i]*prey;
        matrix_data[nnzper*i+19] = 0.75*predation[i]*prey;

        colvals[nnzper*i+15] = GROUPSIZE*i;
        colvals[nnzper*i+16] = GROUPSIZE*i+1;
        colvals[nnzper*i+17] = GROUPSIZE*i+2;
        colvals[nnzper*i+18] = GROUPSIZE*i+3;
        colvals[nnzper*i+19] = GROUPSIZE*i+4;
        
        // fifth row of block: predator
        matrix_data[nnzper*i+20] = 0.0;
        matrix_data[nnzper*i+21] = 0.0;
        matrix_data[nnzper*i+22] = 0.0;
        matrix_data[nnzper*i+23] =  -predation[i]*prey;
        matrix_data[nnzper*i+24] = exp_growth_prey[i] - predation[i]*predator;

        colvals[nnzper*i+20] = GROUPSIZE*i;
        colvals[nnzper*i+21] = GROUPSIZE*i+1;
        colvals[nnzper*i+22] = GROUPSIZE*i+2;
        colvals[nnzper*i+23] = GROUPSIZE*i+3;
        colvals[nnzper*i+24] = GROUPSIZE*i+4;
        /*
        // dead_predator by ge
        colvals[j]   = i*GROUPSIZE + 2;
        colvals[j+1] = i*GROUPSIZE+ 3;
        colvals[j+2] = i*GROUPSIZE + 2;
        colvals[j+3] = i*GROUPSIZE + 3;
        colvals[j+4] = i*GROUPSIZE + 4;
        colvals[j+5] = i*GROUPSIZE + 2;
        colvals[j+6] = i*GROUPSIZE + 3;
        colvals[j+7] = i*GROUPSIZE + 4;
        colvals[j+8] = i*GROUPSIZE + 2;
        colvals[j+9] = i*GROUPSIZE + 3;
        colvals[j+10]= i*GROUPSIZE + 4;


        rowptrs[i*GROUPSIZE+0] = j+0;
        rowptrs[i*GROUPSIZE+1] = j+2;
        rowptrs[i*GROUPSIZE+2] = j+5;
        rowptrs[i*GROUPSIZE+3] = j+5;
        rowptrs[i*GROUPSIZE+4] = j+8;

        matrix_data[j]   = predator*rnatural_death_predator[i];
        matrix_data[j+1] = natural_death_predator[i];
        matrix_data[j+2] = predator*prey*rpredation[i];
        matrix_data[j+3] = predation[i]*prey;
        matrix_data[j+4] = predation[i]*predator;
        matrix_data[j+5] = 0.75*predator*prey*rpredation[i] - predator*rnatural_death_predator[i];
        matrix_data[j+6] = -natural_death_predator[i] + 0.75*predation[i]*prey;
        matrix_data[j+7] = 0.75*predation[i]*predator;
        matrix_data[j+8] = -predator*prey*rpredation[i] + prey*rexp_growth_prey[i];
        matrix_data[j+9] = -predation[i]*prey;
        matrix_data[j+10]= exp_growth_prey[i] - predation[i]* predator;
*/

        
    }

    //SUNSparseMatrix_Print(J, stderr);

    return 0;
}

/* Private functions to output results */

static void PrintOutput(realtype t, realtype y1, realtype y2, realtype y3);

/* Private function to print final statistics */

static void PrintFinalStats(void *cvode_mem);

/* Private function to check function return values */

static int check_retval(void *returnvalue, const char *funcname, int opt);

/*
 *-------------------------------
 * Main Program
 *-------------------------------
 */

int main(int argc, char *argv[])
{
    cudaDeviceReset();

  realtype reltol, t, tout;
  realtype *ydata, *abstol_data;
  N_Vector y, abstol;
  SUNMatrix A;
  SUNLinearSolver LS;
  void *cvode_mem;
  int retval, iout, nnz;
  int neq, ngroups, groupj;

  y = abstol = NULL;
  A = NULL;
  LS = NULL;
  cvode_mem = NULL;

  /* Parse command line arguments */
  if (argc > 1) {
    ngroups = atoi(argv[1]);
  } else {
    ngroups = 2;
  }
  // Read our rate data
  predator_prey_data udata = init_rate_data();

  // GROUPSIZE: number of species
  int GROUPSIZE = 5;
  neq = ngroups * GROUPSIZE;

  udata.neq = neq;

  /* Create CUDA vector of length neq for I.C. and abstol */
  y = N_VNewManaged_Cuda(neq);
  if (check_retval((void *)y, "N_VNewManaged_Cuda", 0)) return(1);
  abstol = N_VNewManaged_Cuda(neq);
  if (check_retval((void *)abstol, "N_VNewManaged_Cuda", 0)) return(1);

  ydata = N_VGetHostArrayPointer_Cuda(y);
  abstol_data = N_VGetHostArrayPointer_Cuda(abstol);

  /* Initialize y */
  for (groupj = 0; groupj < neq; groupj += GROUPSIZE) {
    ydata[groupj]   = 1.0;
    ydata[groupj+1] = 2.0;
    ydata[groupj+2] = 3.0;
    ydata[groupj+3] = 5.0;
    ydata[groupj+4] = 10.0;
  }

  /* Set the scalar relative tolerance */
  reltol = 1.0e-4;

  /* Set the vector absolute tolerance */
  for (groupj = 0; groupj < neq; groupj += GROUPSIZE) {
    abstol_data[groupj]   = 1.0e-6;
    abstol_data[groupj+1] = 1.0e-6;
    abstol_data[groupj+2] = 1.0e-6;
    abstol_data[groupj+3] = 1.0e-6;
    abstol_data[groupj+4] = 1.0e-6;
  }

  /* Call CVodeCreate to create the solver memory and specify the
   * Backward Differentiation Formula */
  cvode_mem = CVodeCreate(CV_BDF);
  if (check_retval((void *)cvode_mem, "CVodeCreate", 0)) return(1);

  /* Call CVodeInit to initialize the integrator memory and specify the
   * user's right hand side function in y'=f(t,y), the inital time T0, and
   * the initial dependent variable vector y. */
  double T0 = 0.0;
  double T1 = 1.0;
  retval = CVodeInit(cvode_mem, f, T0, y);
  if (check_retval(&retval, "CVodeInit", 1)) return(1);

  /* Call CVodeSetUserData to attach the user data structure */
  retval = CVodeSetUserData(cvode_mem, &udata);
  if (check_retval(&retval, "CVodeSetUserData", 1)) return(1);

  /* Call CVodeSVtolerances to specify the scalar relative tolerance
   * and vector absolute tolerances */
  retval = CVodeSVtolerances(cvode_mem, reltol, abstol);
  if (check_retval(&retval, "CVodeSVtolerances", 1)) return(1);

  /* Create sparse SUNMatrix for use in linear solves */
  int NSPARSE = GROUPSIZE*GROUPSIZE;
  nnz = NSPARSE * ngroups;
  A = SUNSparseMatrix(neq, neq, nnz, CSR_MAT);
  if(check_retval((void *)A, "SUNSparseMatrix", 0)) return(1);

  /* Create the SUNLinearSolver object for use by CVode */
  LS = SUNLinSol_cuSolverSp_batchQR(y, A, ngroups, GROUPSIZE, NSPARSE);
  if(check_retval((void *)LS, "SUNLinSol_cuSolverSp_batchQR", 0)) return(1);

  /* Call CVodeSetLinearSolver to attach the matrix and linear solver to CVode */
  retval = CVodeSetLinearSolver(cvode_mem, LS, A);
  if(check_retval(&retval, "CVodeSetLinearSolver", 1)) return(1);

  /* Set the user-supplied Jacobian routine Jac */
  retval = CVodeSetJacFn(cvode_mem, Jacobian);
  if(check_retval(&retval, "CVodeSetJacFn", 1)) return(1);

  /* In loop, call CVode, print results, and test for error.
     Break out of loop when NOUT preset output times have been reached.  */
  int NOUT = 1;
  int TMULT = 1.5;
  
  printf(" \nGroup of independent 3-species kinetics problems\n\n");
  printf("number of groups = %d\n\n", ngroups);

  iout = 0;  tout = T1;
  while(1) {
    retval = CVode(cvode_mem, tout, y, &t, CV_NORMAL);

    for (groupj = 0; groupj < 1; groupj++) {
      printf("group %d: \n", groupj);
      printf("y[%d] = %0.5g\n", GROUPSIZE*groupj, ydata[GROUPSIZE*groupj]);
      printf("y[%d] = %0.5g\n", GROUPSIZE*groupj+1, ydata[GROUPSIZE*groupj+1]);
      printf("y[%d] = %0.5g\n", GROUPSIZE*groupj+2, ydata[GROUPSIZE*groupj+2]);
      printf("y[%d] = %0.5g\n", GROUPSIZE*groupj+3, ydata[GROUPSIZE*groupj+3]);
      printf("y[%d] = %0.5g\n", GROUPSIZE*groupj+4, ydata[GROUPSIZE*groupj+4]);
    }

    if (check_retval(&retval, "CVode", 1)) break;
    if (retval == CV_SUCCESS) {
      iout++;
      tout *= TMULT;
    }

    if (iout == NOUT) break;
  }

  /* Print some final statistics */
  PrintFinalStats(cvode_mem);

  /* Free y and abstol vectors */
  N_VDestroy(y);
  N_VDestroy(abstol);

  /* Free integrator memory */
  CVodeFree(&cvode_mem);

  /* Free the linear solver memory */
  SUNLinSolFree(LS);

  /* Free the matrix memory */
  SUNMatDestroy(A);

  return(0);
}


/*
 *-------------------------------
 * Private helper functions
 *-------------------------------
 */

static void PrintOutput(realtype t, realtype y1, realtype y2, realtype y3)
{
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("At t = %0.4Le      y =%14.6Le  %14.6Le  %14.6Le\n", t, y1, y2, y3);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("At t = %0.4e      y =%14.6e  %14.6e  %14.6e\n", t, y1, y2, y3);
#else
  printf("At t = %0.4e      y =%14.6e  %14.6e  %14.6e\n", t, y1, y2, y3);
#endif

  return;
}

/*
 * Get and print some final statistics
 */

static void PrintFinalStats(void *cvode_mem)
{
  long int nst, nfe, nsetups, nje, nni, ncfn, netf, nge;
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

  printf("\nFinal Statistics:\n");
  printf("nst = %-6ld nfe  = %-6ld nsetups = %-6ld nje = %ld\n",
	 nst, nfe, nsetups, nje);
  printf("nni = %-6ld ncfn = %-6ld netf = %-6ld    nge = %ld\n \n",
	 nni, ncfn, netf, nge);
}

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
