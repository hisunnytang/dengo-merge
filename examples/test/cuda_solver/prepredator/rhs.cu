
#include <stdio.h>

#include <cvode/cvode.h>                              /* prototypes for CVODE fcts., consts.          */
#include <nvector/nvector_cuda.h>                     /* access to cuda N_Vector                      */
#include <sunmatrix/sunmatrix_sparse.h>               /* access to sparse SUNMatrix                   */
#include <sunlinsol/sunlinsol_cusolversp_batchqr.h>   /* acess to cuSolverSp batch QR SUNLinearSolver */
#include <sundials/sundials_types.h>                  /* defs. of realtype, int              */


__global__
static void f_kernel(realtype t, realtype* ydata, realtype* ydotdata,
                     int neq, int ngroups)
{
  realtype y1, y2, y3, yd1, yd3;
  int i = blockIdx.x*blockDim.x + threadIdx.x;
  int groupj = i*GROUPSIZE;

  if (i < neq) {
    y1 = ydata[groupj]; y2 = ydata[groupj+1]; y3 = ydata[groupj+2];

    yd1 = ydotdata[groupj]   = RCONST(-0.04)*y1 + RCONST(1.0e4)*y2*y3;
    yd3 = ydotdata[groupj+2] = RCONST(3.0e7)*y2*y2;
          ydotdata[groupj+1] = -yd1 - yd3;
  }
}

__global__
static void prey_predator_f_kernel(realtype t, realtype* ydata, realtype* ydotdata,
        predator_prey_data *data,
        int neq, int ngroups)
{
    double dead_predator;
    double dead_prey;
    double ge;
    double predator;
    double prey;

    double *natural_death_predator = data->r_natural_death_predator;
    double *predation              = data->r_predation;
    double *exp_growth_prey        = data->r_exp_growth_prey;

    int i = blockIdx.x*blockDim.x + threadIdx.x;
    int groupj = i*GROUPSIZE;
    if (i < neq)
    {
        dead_predator = ydata[j];
        dead_prey     = ydata[j+1];
        ge            = ydata[j+2];
        predator      = ydata[j+3];
        prey          = ydata[j+4];

        // Species: dead_predator
        ydotdata[groupj] = natural_death_predator[j]*predator;
        // Species: dead_prey
        ydot_ptr[groupj+1] = predation[j]*predator*prey;
        // Species: ge
        ydot_ptr[groupj+2] = 0;
        // Species: predator
        ydot_ptr[groupj+3] = -natural_death_predator[j]*predator + 0.75*predation[j]*predator*prey;
        // Species: prey
        ydot_ptr[groupj+4] = exp_growth_prey[j]*prey - predation[j]*predator*prey;

        /*
        #ifdef SCALE_INPUT
        ydot_ptr[groupj]   *= inv_scale[groupj];
        ydot_ptr[groupj+1] *= inv_scale[groupj+1];
        ydot_ptr[groupj+2] *= inv_scale[groupj+2];
        ydot_ptr[groupj+3] *= inv_scale[groupj+3];
        ydot_ptr[groupj+5] *= inv_scale[groupj+4];
        #endif
        */
    }
}

__global__
static void prey_predator_temperature_kernel(realtype t, realtype* Ts, realtype *logTs,
        realtype* ydata)
{

    double dead_predator;
    double dead_prey;
    double ge;
    double predator;
    double prey;

    int i = blockIdx.x*blockDim.x + threadIdx.x;
    int groupj = i*GROUPSIZE;
    if (i < neq)
    {
        dead_predator = ydata[j];
        dead_prey     = ydata[j+1];
        ge            = ydata[j+2];
        predator      = ydata[j+3];
        prey          = ydata[j+4];
        // calculate temperature here
        // quasi equil cases we might need the rates data, ignore it for now.

        // sample
        // density =
        // density = 1.0*dead_predator + 1.0*dead_prey + 1.0*predator + 1.0*prey;
        // data->Ts[threadID][i] = density*ge*mh/(kb*(_gamma_m1*dead_predator + _gamma_m1*dead_prey + _gamma_m1*predator + _gamma_m1*prey));
        /* if (data->Ts[threadID][i] < data->bounds[0]) {
            data->Ts[threadID][i] = data->bounds[0];
        } else if (data->Ts[threadID][i] > data->bounds[1]) {
            data->Ts[threadID][i] = data->bounds[1];
        }
        data->logTs[threadID][i] = log(data->Ts[threadID][i]);
        data->invTs[threadID][i] = 1.0 / data->Ts[threadID][i];

        dge_dT = kb*(_gamma_m1*dead_predator + _gamma_m1*dead_prey + _gamma_m1*predator + _gamma_m1*prey)/(density*mh);
        data->dTs_ge[threadID][i] = 1.0 / dge_dT;
         */
        Ts[j] = 10.0;
        logTs[j] = 1.0;

    }
}

/*
 *-------------------------------
 * Functions called by the solver
 *-------------------------------
 */

/* Right hand side function. This just launches the CUDA kernel
   to do the actual computation. At the very least, doing this
   saves moving the vector data in y and ydot to/from the device
   every evaluation of f. */
static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  UserData *udata;
  realtype *ydata, *ydotdata;

  udata = (UserData*) user_data;
  ydata = N_VGetDeviceArrayPointer_Cuda(y);
  ydotdata = N_VGetDeviceArrayPointer_Cuda(ydot);

  unsigned block_size = 32;
  unsigned grid_size = (udata->neq + block_size - 1) / block_size;
  f_kernel<<<grid_size, block_size>>>(t, ydata, ydotdata, udata->neq, udata->ngroups);

prey_predator_temperature_kernel<<<grid_size, block_size>>>(t,Ts,logTs,ydata);
  cudaError_t cuerr = cudaGetLastError();
  if (cuerr != cudaSuccess) {
    fprintf(stderr,
            ">>> ERROR in f: cudaGetLastError returned %s\n",
            cudaGetErrorName(cuerr));
    return 1;
  }

  return(0);
}


// our modification
static int rhs_function(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  predator_prey_data *data = (predator_prey_data* ) user_data;
  realtype *ydata, *ydotdata;
  ydata = N_VGetDeviceArrayPointer_Cuda(y);
  ydotdata = N_VGetDeviceArrayPointer_Cuda(ydot);

  unsigned block_size = 32;
  unsigned grid_size = (udata->neq + block_size - 1) / block_size;

  // kernel to calculate temperature

  f_kernel<<<grid_size, block_size>>>(t, ydata, ydotdata, udata->neq, udata->ngroups);
  // kernel to update the rates

  // kerenel to update the f
  f_kernel<<<grid_size, block_size>>>(t, ydata, ydotdata, udata->neq, udata->ngroups);

  cudaError_t cuerr = cudaGetLastError();
  if (cuerr != cudaSuccess) {
    fprintf(stderr,
            ">>> ERROR in f: cudaGetLastError returned %s\n",
            cudaGetErrorName(cuerr));
    return 1;
  }

  return(0);
}

/*
 * ----------------------------
 * Calculate Sparse Jacobian
 *-----------------------------
 */
int calculate_sparse_jacobian_predator_prey( realtype t,
                                        N_Vector y, N_Vector fy,
                                        SUNMatrix J, void *user_data,
                                        N_Vector tmp1, N_Vector tmp2,
                                        N_Vector tmp3)
{
    /* We iterate over all of the rates */
    /* Calcuate temperature first */


    predator_prey_data *data = (predator_prey_data*)user_data;

    int nchem = 5;
    int nstrip = data->nstrip;
    int i, j;
    int NSPARSE = 11;

    #ifdef _OPENMP
    int threadID = omp_get_thread_num();
    #else
    int threadID = 0;
    #endif

    /* change N_Vector back to an array */
    double y_arr[ 5 * nstrip ];
    double *scale     = data->scale[threadID];
    double *inv_scale = data->inv_scale[threadID];

    //TODO: Here we assumed during the evaluation of jacobian
    // temperature is approximately constant,
    // i.e. close enough to the point evaluation of f(y)
    // such that the rates and temperature need not be interpolated or evalulated
    // again during the jacobian evaluation.
    // We havent really fully explored the effect of this `assumption`...
    // But it definitely boost the performance

    /*
    int flag;
    flag = predator_prey_calculate_temperature(data, y_arr , nstrip, nchem );
    if (flag > 0){
        // check if the temperature failed to converged
        return -1;
    }
    predator_prey_interpolate_rates(data, nstrip);
    */

    // predator_prey_calculate_temperature(data, y_arr, nstrip, nchem);
    // predator_prey_interpolate_rates(data, nstrip);

    /* Now We set up some temporaries */
    // CSR is what we choose
    sunindextype *rowptrs = SUNSparseMatrix_IndexPointers(J);
    sunindextype *colvals = SUNSparseMatrix_IndexValues(J);
    realtype *matrix_data = SUNSparseMatrix_Data(J);

    SUNMatZero(J);

    double *Tge = data->dTs_ge[threadID];
    double *exp_growth_prey = data->rs_exp_growth_prey[threadID];
    double *rexp_growth_prey= data->drs_exp_growth_prey[threadID];
    double *natural_death_predator = data->rs_natural_death_predator[threadID];
    double *rnatural_death_predator= data->drs_natural_death_predator[threadID];
    double *predation = data->rs_predation[threadID];
    double *rpredation= data->drs_predation[threadID];
    double dead_predator;
    double dead_prey;
    double ge;
    double predator;
    double prey;
    double z;
    double T;

    double mh = 1.66054e-24;
    double mdensity, inv_mdensity;

    double scale2, inv_scale1;

    j = 0;
    mdensity = 0.0;
    z = data->current_z;

    int k = 0;

    double *yvec_ptr = N_VGetArrayPointer(y);

    for ( i = 0; i < nstrip; i++ ){

        #ifdef SCALE_INPUT
        j = i * nchem;
        dead_predator = yvec_ptr[j]*scale[j];
        j++;

        dead_prey = yvec_ptr[j]*scale[j];
        j++;

        ge = yvec_ptr[j]*scale[j];
        j++;

        predator = yvec_ptr[j]*scale[j];
        j++;

        prey = yvec_ptr[j]*scale[j];
        j++;

        #else
        j = i * nchem;
        dead_predator = yvec_ptr[j];
        j++;

        dead_prey = yvec_ptr[j];
        j++;

        ge = yvec_ptr[j];
        j++;

        predator = yvec_ptr[j];
        j++;

        prey = yvec_ptr[j];
        j++;

        #endif

        mdensity = data->mdensity[threadID][i];
        inv_mdensity = 1.0 / mdensity;




        j = i * NSPARSE;
        // dead_predator by ge
        colvals[j + 0] = i * nchem + 2 ;
        matrix_data[ j + 0 ] = predator*rnatural_death_predator[i];


        matrix_data[ j + 0] *= Tge[i];
        // dead_predator by predator
        colvals[j + 1] = i * nchem + 3 ;
        matrix_data[ j + 1 ] = natural_death_predator[i];


        // dead_prey by ge
        colvals[j + 2] = i * nchem + 2 ;
        matrix_data[ j + 2 ] = predator*prey*rpredation[i];


        matrix_data[ j + 2] *= Tge[i];
        // dead_prey by predator
        colvals[j + 3] = i * nchem + 3 ;
        matrix_data[ j + 3 ] = predation[i]*prey;


        // dead_prey by prey
        colvals[j + 4] = i * nchem + 4 ;
        matrix_data[ j + 4 ] = predation[i]*predator;


        // predator by ge
        colvals[j + 5] = i * nchem + 2 ;
        matrix_data[ j + 5 ] = 0.75*predator*prey*rpredation[i] - predator*rnatural_death_predator[i];


        matrix_data[ j + 5] *= Tge[i];
        // predator by predator
        colvals[j + 6] = i * nchem + 3 ;
        matrix_data[ j + 6 ] = -natural_death_predator[i] + 0.75*predation[i]*prey;


        // predator by prey
        colvals[j + 7] = i * nchem + 4 ;
        matrix_data[ j + 7 ] = 0.75*predation[i]*predator;


        // prey by ge
        colvals[j + 8] = i * nchem + 2 ;
        matrix_data[ j + 8 ] = -predator*prey*rpredation[i] + prey*rexp_growth_prey[i];


        matrix_data[ j + 8] *= Tge[i];
        // prey by predator
        colvals[j + 9] = i * nchem + 3 ;
        matrix_data[ j + 9 ] = -predation[i]*prey;


        // prey by prey
        colvals[j + 10] = i * nchem + 4 ;
        matrix_data[ j + 10 ] = exp_growth_prey[i] - predation[i]*predator;


        rowptrs[ i * nchem +  0] = i * NSPARSE + 0;
        rowptrs[ i * nchem +  1] = i * NSPARSE + 2;
        rowptrs[ i * nchem +  2] = i * NSPARSE + 5;
        rowptrs[ i * nchem +  3] = i * NSPARSE + 5;
        rowptrs[ i * nchem +  4] = i * NSPARSE + 8;

        #ifdef SCALE_INPUT
        j = i * nchem;
        inv_scale1 = inv_scale[ j + 0 ];
        scale2     = scale    [ j + 2 ];
        matrix_data[ i * NSPARSE + 0]  *= inv_scale1*scale2;
        inv_scale1 = inv_scale[ j + 0 ];
        scale2     = scale    [ j + 3 ];
        matrix_data[ i * NSPARSE + 1]  *= inv_scale1*scale2;
        inv_scale1 = inv_scale[ j + 1 ];
        scale2     = scale    [ j + 2 ];
        matrix_data[ i * NSPARSE + 2]  *= inv_scale1*scale2;
        inv_scale1 = inv_scale[ j + 1 ];
        scale2     = scale    [ j + 3 ];
        matrix_data[ i * NSPARSE + 3]  *= inv_scale1*scale2;
        inv_scale1 = inv_scale[ j + 1 ];
        scale2     = scale    [ j + 4 ];
        matrix_data[ i * NSPARSE + 4]  *= inv_scale1*scale2;
        inv_scale1 = inv_scale[ j + 3 ];
        scale2     = scale    [ j + 2 ];
        matrix_data[ i * NSPARSE + 5]  *= inv_scale1*scale2;
        inv_scale1 = inv_scale[ j + 3 ];
        scale2     = scale    [ j + 3 ];
        matrix_data[ i * NSPARSE + 6]  *= inv_scale1*scale2;
        inv_scale1 = inv_scale[ j + 3 ];
        scale2     = scale    [ j + 4 ];
        matrix_data[ i * NSPARSE + 7]  *= inv_scale1*scale2;
        inv_scale1 = inv_scale[ j + 4 ];
        scale2     = scale    [ j + 2 ];
        matrix_data[ i * NSPARSE + 8]  *= inv_scale1*scale2;
        inv_scale1 = inv_scale[ j + 4 ];
        scale2     = scale    [ j + 3 ];
        matrix_data[ i * NSPARSE + 9]  *= inv_scale1*scale2;
        inv_scale1 = inv_scale[ j + 4 ];
        scale2     = scale    [ j + 4 ];
        matrix_data[ i * NSPARSE + 10]  *= inv_scale1*scale2;
        #endif

    }

    rowptrs[ i * nchem ] = i * NSPARSE ;
    return 0;
}



/*
 *-------------------------------
 * Main Program
 *-------------------------------
 */
int main(int argc, char *argv[])
{
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
    ngroups = 1000;
  }
  neq = ngroups * GROUPSIZE;

  const char * FileLocation = field_data->dengo_data_file;
  predator_prey_data *udata = predator_prey_setup_data( FileLocation, NULL, NULL);

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
    ydata[groupj+1] = 1.0;
    ydata[groupj+2] = 1.0;
    ydata[groupj+3] = 1.0;
    ydata[groupj+4] = 1.0;

  }

  /* Set the scalar relative tolerance */
  reltol = RTOL;

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
  retval = CVodeInit(cvode_mem, f, T0, y);
  if (check_retval(&retval, "CVodeInit", 1)) return(1);

  /* Call CVodeSetUserData to attach the user data structure */
  retval = CVodeSetUserData(cvode_mem, udata);
  if (check_retval(&retval, "CVodeSetUserData", 1)) return(1);

  /* Call CVodeSVtolerances to specify the scalar relative tolerance
   * and vector absolute tolerances */
  retval = CVodeSVtolerances(cvode_mem, reltol, abstol);
  if (check_retval(&retval, "CVodeSVtolerances", 1)) return(1);

  /* Create sparse SUNMatrix for use in linear solves */
  nnz = GROUPSIZE * GROUPSIZE * ngroups;
  A = SUNSparseMatrix(neq, neq, nnz, CSR_MAT);
  if(check_retval((void *)A, "SUNSparseMatrix", 0)) return(1);

  /* Create the SUNLinearSolver object for use by CVode */
  LS = SUNLinSol_cuSolverSp_batchQR(y, A, ngroups, GROUPSIZE, GROUPSIZE*GROUPSIZE);
  if(check_retval((void *)LS, "SUNLinSol_cuSolverSp_batchQR", 0)) return(1);

  /* Call CVodeSetLinearSolver to attach the matrix and linear solver to CVode */
  retval = CVodeSetLinearSolver(cvode_mem, LS, A);
  if(check_retval(&retval, "CVodeSetLinearSolver", 1)) return(1);

  // Jacobian defined;
  jac_f jf = calculate_sparse_jacobian_predator_prey;
  /* Set the user-supplied Jacobian routine Jac */
  retval = CVodeSetJacFn(cvode_mem, Jac);
  if(check_retval(&retval, "CVodeSetJacFn", 1)) return(1);

  /* In loop, call CVode, print results, and test for error.
     Break out of loop when NOUT preset output times have been reached.  */
  printf(" \nGroup of independent 3-species kinetics problems\n\n");
  printf("number of groups = %d\n\n", ngroups);

  iout = 0;  tout = T1;
  while(1) {
    retval = CVode(cvode_mem, tout, y, &t, CV_NORMAL);

    for (groupj = 0; groupj < 1; groupj++) {
      printf("group %d: ", groupj);
      PrintOutput(t, ydata[GROUPSIZE*groupj],
                     ydata[1+GROUPSIZE*groupj],
                     ydata[2+GROUPSIZE*groupj],
                     ydata[3+GROUPSIZE*groupj],
                     ydata[4+GROUPSIZE*groupj]);
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
