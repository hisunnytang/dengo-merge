#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#include <cvodes/cvodes.h>
#include <cvodes/cvodes_dense.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_math.h>

/* Accessor macros */

#define Ith(v,i) NV_Ith_S(v,i-1)
#define IJth(A,i,j) DENSE_ELEM(A,i-1,j-1)

/* User-defined vector and matrix accessor macros: Ith, IJth */

/* These macros are defined in order to write code which exactly matches
   the mathematical problem description given above.

   Ith(v,i) references the ith component of the vector v, where i is in
   the range [1..NEQ] and NEQ is defined below. The Ith macro is defined
   using the N_VIth macro in nvector.h. N_VIth numbers the components of
   a vector starting from 0.

   IJth(A,i,j) references the (i,j)th element of the dense matrix A, where
   i and j are in the range [1..NEQ]. The IJth macro is defined using the
   DENSE_ELEM macro in dense.h. DENSE_ELEM numbers rows and columns of a
   dense matrix starting from 0. */



/* Functions Called by the Solver */

static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);

static int Jac(long int N, realtype t,
               N_Vector y, N_Vector fy, DlsMat J, void *user_data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

/* Private functions to output results */

static void PrintOutput(realtype t, realtype y1, realtype y2, realtype y3);
static void PrintRootInfo(int root_f1, int root_f2);

/* Private function to print final statistics */

static void PrintFinalStats(void *cvode_mem);

/* Private function to check function return values */

static int check_flag(void *flagvalue, const char *funcname, int opt);



/*
 * -------------------------------------------------------
 *  MAIN PROGRAM
 * -------------------------------------------------------
 */

int main_solver( rhs_f f, jac_f Jac, 
                 double *input, double rtol, double atol,
                 int nchem, void *sdata, double T0, double T1)
{
    void *cvode_mem;
    int flag;
    N_Vector y, abstol;


    
    y = NULL;
    cvode_mem = NULL;
    
    /* Create serial vector of length NEQ for I.C. and abstol */
    y = N_VNew_Serial(nchem);
    if (check_flag((void *)y, "N_VNew_Serial" , 0)) return(1);

    abstol = N_VNew_Serial(nchem);
    if (check_flag((void *)abstol, "N_VNew_Serial", 0)) return(1);

    /* Initialize y */
    for (i=0; i<nchem; i++) { Ith(y,i+1)      = input[i]; }
    /* Set the vector absolute tolerance */
    for (i=0; i<nchem; i++) { Ith(abstol,i+1) = atol[i] ; }
    /* Set the scalar relative tolerance */
    reltol = rtol;


    /* Create CVODES object */
    cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
    if (check_flag((void *)cvode_mem, "CVodeCreate", 0)) return(1);
    
    /* Allocate space for CVODES */
    /* f: RHS function
     * y: Array of initial values
     * T0: Initial Time
     */

    /* Allocate space for CVODES */
    flag = CVodeInit(cvode_mem, f, T0, y);
    if (check_flag( &flag, "CVodeInit", 1)) return(1);

    
    /* Call CVodesSVtolerances to specify the scalar relative tolerance
     * and vector absolute tolerances */
    flag = CVodeSVtolerances(cvode_mem, reltol, abstol);
    if (check_flag(&flag, "CVodeSVtolerances", 1)) return(1);
    
    /* Attach User Data */
    flag = CVodeSetUserData(cvode_mem, sdata);
    if (check_flag(&flag, "CVodeSetUserData", 1)) return(1);

    /* Attach Linear Solver */
    flag = CVDense(cvode_mem, NEQ);
    if (check_flag(&flag,"CVDense", 1)) return(1);
    
    /* Attach the Jacobian Function */
    flag = CVDlsSetDenseJacFn(cvode_mem, Jac);
    if (check_flag(&flag, "CVDlsSetDenseJacFn", 1)) return(1);
    
    


    /* In Loop, call CVode, print results, and test for error.
     * Break out of loop when NOUT preset output times have been reached
     */
    for (iout = 1, tout=T1, iput <=NOUT, iout++, tout *=TMULT){

        flag = CVode( cvode_mem, tout, y, &t, CV_NORMAL);
        if (check_flag(&flag, "CVode", 1)) break;

        PrintOutput(cvode_mem, t, y);
    }

    /* Print final Statisitics */
    PrintFinalStats(cvode_mem, sensi);

    /* Free Memory */
    N_VDestroy_Serial(y);
    if (sensi) {
        N_VDestroyVectorArray_Serial(yS, NS);
    }
    CVodeFree(&cvode_mem);

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

static void PrintRootInfo(int root_f1, int root_f2)
{
  printf("    rootsfound[] = %3d %3d\n", root_f1, root_f2);

  return;
}

/* 
 * Get and print some final statistics
 */

static void PrintFinalStats(void *cvode_mem)
{
  long int nst, nfe, nsetups, nje, nfeLS, nni, ncfn, netf, nge;
  int flag;

  flag = CVodeGetNumSteps(cvode_mem, &nst);
  check_flag(&flag, "CVodeGetNumSteps", 1);
  flag = CVodeGetNumRhsEvals(cvode_mem, &nfe);
  check_flag(&flag, "CVodeGetNumRhsEvals", 1);
  flag = CVodeGetNumLinSolvSetups(cvode_mem, &nsetups);
  check_flag(&flag, "CVodeGetNumLinSolvSetups", 1);
  flag = CVodeGetNumErrTestFails(cvode_mem, &netf);
  check_flag(&flag, "CVodeGetNumErrTestFails", 1);
  flag = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
  check_flag(&flag, "CVodeGetNumNonlinSolvIters", 1);
  flag = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
  check_flag(&flag, "CVodeGetNumNonlinSolvConvFails", 1);

  flag = CVDlsGetNumJacEvals(cvode_mem, &nje);
  check_flag(&flag, "CVDlsGetNumJacEvals", 1);
  flag = CVDlsGetNumRhsEvals(cvode_mem, &nfeLS);
  check_flag(&flag, "CVDlsGetNumRhsEvals", 1);

  flag = CVodeGetNumGEvals(cvode_mem, &nge);
  check_flag(&flag, "CVodeGetNumGEvals", 1);

  printf("\nFinal Statistics:\n");
  printf("nst = %-6ld nfe  = %-6ld nsetups = %-6ld nfeLS = %-6ld nje = %ld\n",
	 nst, nfe, nsetups, nfeLS, nje);
  printf("nni = %-6ld ncfn = %-6ld netf = %-6ld nge = %ld\n \n",
	 nni, ncfn, netf, nge);
}

/*
 * Check function return value...
 *   opt == 0 means SUNDIALS function allocates memory so check if
 *            returned NULL pointer
 *   opt == 1 means SUNDIALS function returns a flag so check if
 *            flag >= 0
 *   opt == 2 means function allocates memory so check if returned
 *            NULL pointer 
 */

static int check_flag(void *flagvalue, const char *funcname, int opt)
{
  int *errflag;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && flagvalue == NULL) {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1); }

  /* Check if flag < 0 */
  else if (opt == 1) {
    errflag = (int *) flagvalue;
    if (*errflag < 0) {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
	      funcname, *errflag);
      return(1); }}

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && flagvalue == NULL) {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1); }

  return(0);
}
