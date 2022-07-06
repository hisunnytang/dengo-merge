#include<stdio.h>
#include<stdlib.h>
#include<string.h>

//#include <cvodes/cvodes.h>
//#include <cvodes/cvodes_dense.h>
//#include <nvector/nvector_serial.h>
//#include <sundials/sundials_types.h>
//#include <sundials/sundials_math.h>




#include <cvode/cvode.h>               /* prototypes for CVODE fcts., consts.  */
#include <nvector/nvector_serial.h>    /* access to serial N_Vector            */
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix            */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver      */
#include <cvode/cvode_direct.h>        /* access to CVDls interface            */
#include <sundials/sundials_types.h>   /* defs. of realtype, sunindextype      */
/* Define the data structures */
#include <cvdls_9species_solver.h>


/* Accessor macros */

#define Ith(v,i) NV_Ith_S(v,i-1)
#define IJth(A,i,j) SM_ELEMENT_D(A,i-1,j-1)

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

typedef int(*rhs_f)( realtype, N_Vector , N_Vector , void * );
typedef int(*jac_f)( realtype, N_Vector  , N_Vector , SUNMatrix , void *, N_Vector, N_Vector, N_Vector);


/* Functions Called by the Solver */

int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);

int Jac(realtype t, N_Vector y, N_Vector fy, SUNMatrix J,
        void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

/* Private functions to output results */

static void PrintOutput(void *cvode_mem, realtype t, N_Vector u);
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

int cvodes_main_solver( rhs_f f, jac_f Jac,
                 double *input , double *rtol, double *atol, int NEQ, void *sdata, double *dt_now)
{
    void *cvode_mem;
    int i, flag, flagr, iout;
    realtype t, reltol;
    N_Vector y, abstol;
    SUNLinearSolver LS;
    SUNMatrix A;

    y = abstol = NULL;
    cvode_mem = NULL;
    LS = NULL;
    A = NULL;

    cvdls_9species_data *data = (cvdls_9species_data*)sdata;


    /* Create serial vector of length NEQ for I.C. and abstol */
    y = N_VNew_Serial(NEQ);
    if (check_flag((void *)y, "N_VNew_Serial" , 0)) return(1);

    abstol = N_VNew_Serial(NEQ);
    if (check_flag((void *)abstol, "N_VNew_Serial", 0)) return(1);

    /* Initialize y
     * Rescale the input variables to unity
     */
    double scale;
    for (i=0; i<NEQ; i++) {
        data->scale[i] = input[i];
        scale = data->scale[i];
        Ith(y,i+1)      = input[i] / scale;
        fprintf(stderr, "scale[%d]: %0.5g \n", i, scale);
    }

    /* fixed the incoming abs tolerance */
    /* Set the vector absolute tolerance */
    for (i=0; i<NEQ; i++) {
        Ith(abstol,i+1) = 1.0e-9 ;
        }

    /* Set the scalar relative tolerance */
    reltol = 1.0e-9;


    /* Create CVODES object */
    cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
    if (check_flag((void *)cvode_mem, "CVodeCreate", 0)) return(1);

    /* Allocate space for CVODES */
    /* f: RHS function
     * y: Array of initial values
     * T0: Initial Time
     */

    /* Allocate space for CVODES */
    flag = CVodeInit(cvode_mem, f, 0.0, y);
    if (check_flag( &flag, "CVodeInit", 1)) return(1);

    flag = CVodeSetMaxNumSteps(cvode_mem, 500 );

    /* Call CVodesSVtolerances to specify the scalar relative tolerance
     * and vector absolute tolerances */
    flag = CVodeSVtolerances(cvode_mem, reltol, abstol);
    if (check_flag(&flag, "CVodeSVtolerances", 1)) return(1);

    /* Attach User Data */
    flag = CVodeSetUserData(cvode_mem, sdata);
    if (check_flag(&flag, "CVodeSetUserDatr", 1)) return(1);


    /*Create dense SUNMatrix for use in linear solves*/
    A = SUNDenseMatrix(NEQ, NEQ);
    /* Create dense SUNLinearSolver object for use by CVode */
    LS = SUNDenseLinearSolver(y, A);
    if(check_flag((void *)LS, "SUNDenseLinearSolver", 0)) return(1);

    /* Call CVDlsSetLinearSolver to attach the matrix and linear solver to CVode */
    flag = CVDlsSetLinearSolver(cvode_mem, LS, A);
    if(check_flag(&flag, "CVDlsSetLinearSolver", 1)) return(1);

    /* Set the user-supplied Jacobian routine Jac */
    flag = CVDlsSetJacFn(cvode_mem, Jac);
    if(check_flag(&flag, "CVDlsSetJacFn", 1)) return(1);

    // /* Attach Linear Solver */
    //flag = CVDense(cvode_mem, NEQ);
    //if (check_flag(&flag,"CVDense", 1)) return(1);

    /* Attach the Jacobian Function */
    //flag = CVDlsSetDenseJacFn(cvode_mem, Jac);
    //if (check_flag(&flag, "CVDlsSetDenseJacFn", 1)) return(1);



    t = 0;
    double tout = dt_now[0];

    // tout is the desired output time
    // evolve CVODE in one dt
    flag = CVode( cvode_mem, tout, y, &t, CV_NORMAL);
    /* t is the actual evolved time */
    dt_now[0]   = t;


    // rescale the input
    for (i=0; i<NEQ; i++) {
        scale = data->scale[i];
        input[i] = Ith(y,i+1)*scale;
        }

    /* Free Memory */
    N_VDestroy(y);
    N_VDestroy(abstol);
    CVodeFree(&cvode_mem);


    if (flag == CV_TOO_MUCH_WORK){
        /* The solver took mxstep internal steps but still could not reach tout */
        return 1;
    }

    if (flag == CV_CONV_FAILURE){
        /* Either convergence test failures occurred too many times
         * during one internal time step, or with |h| = hmin
         */
        return 1;
    }
    if (flag == CV_TOO_CLOSE){
        /* The initial time t0 and the final time
         * tout are too close to each other
         * and the user did not specify an
         * initial step size
         */
        return 1;
    }



    return 0;


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
