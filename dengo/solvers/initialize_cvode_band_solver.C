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
#include <sunmatrix/sunmatrix_band.h> /* access to band SUNMatrix            */
#include <sunlinsol/sunlinsol_band.h> /* access to band SUNLinearSolver      */
#include <cvode/cvode_direct.h>        /* access to CVDls interface            */
#include <sundials/sundials_types.h>   /* defs. of realtype, sunindextype      */
/* Define the data structures */
#include <{{solver_name}}_solver.h>
#include <cvode/cvode_spils.h>           /* access to CVSpils interface       */
#include <sunlinsol/sunlinsol_spgmr.h>   /* access to SPGMR SUNLinearSolver   */


typedef int(*rhs_f)( realtype, N_Vector , N_Vector , void * );
typedef int(*jac_f)( realtype, N_Vector  , N_Vector , SUNMatrix , void *, N_Vector, N_Vector, N_Vector);


/* Functions Called by the Solver */

int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);

int Jac(realtype t, N_Vector y, N_Vector fy, SUNMatrix J, 
        void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
static int check_flag(void *flagvalue, const char *funcname, int opt);

int cvode_solver( void *cvode_mem, double *output, int NEQ, double *dt, {{solver_name}}_data * data, N_Vector y){
    
    int flag, i;

    for (i=0; i<NEQ; i++) {
        NV_Ith_S(y,i)   = 1.0;
    }
    
    flag = CVodeReInit( cvode_mem, 0.0, y);
    //if (check_flag( &flag, "CVodeReInit", 1)) return(1);

    double tout = dt[0];

    flag = CVode( cvode_mem, tout, y, dt, CV_NORMAL);
    
    for ( i = 0; i < NEQ; i++){
       output[i] = NV_Ith_S(y, i);
    }

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

void *setup_cvode_band_solver( rhs_f f, jac_f Jac, double abstol, int NEQ, 
        {{solver_name}}_data *data, SUNLinearSolver LS, SUNMatrix A, N_Vector y){
    
    void *cvode_mem;
    cvode_mem = NULL;
    
    int i, flag;

    /* Create CVODES object */
    cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
    if (check_flag((void *)cvode_mem, "CVodeCreate", 0)) return(NULL);
    
    
    //
    // Initialize y 
    // Rescale the input variables to unity
    //
    double scale;
    for (i=0; i<NEQ; i++) {
        NV_Ith_S(y,i)   = 1.0;
    }
    /* Allocate space for CVODES */
    flag = CVodeInit(cvode_mem, f, 0.0, y);
    if (check_flag( &flag, "CVodeInit", 1)) return(NULL);

    flag = CVodeSetMaxNumSteps(cvode_mem, 5000 );

    /* Call CVodesSVtolerances to specify the scalar relative tolerance
     * and vector absolute tolerances */

    /* 
     * since inputs are scaled to 1 
     * abstol should therefore be the same as reltol 
     * 
     */
    flag = CVodeSStolerances(cvode_mem, abstol, abstol);
    if (check_flag(&flag, "CVodeSStolerances", 1)) return(NULL);
    
    /* Attach User Data */
    flag = CVodeSetUserData(cvode_mem, data);
    if (check_flag(&flag, "CVodeSetUserData", 1)) return(NULL);
        

    /* Call CVDlsSetLinearSolver to attach the matrix 
     * and linear solver to CVode */
    // flag = CVDlsSetLinearSolver(cvode_mem, LS, A);
    // if(check_flag(&flag, "CVDlsSetLinearSolver", 1)) return(NULL);
    
    /* Set the user-supplied Jacobian routine Jac */
    // flag = CVDlsSetJacFn(cvode_mem, Jac);
    // if(check_flag(&flag, "CVDlsSetJacFn", 1)) return(NULL);


    LS = SUNSPGMR(y, PREC_NONE, 0);
    if(check_flag(&flag, "SUNSPGMR", 1)) return(NULL);
    

    flag = CVSpilsSetLinearSolver(cvode_mem, LS);
    if(check_flag(&flag, "CVSpilsSetLinearSolver", 1)) return(NULL);
    /* Set the JAcobian-times-vector function */
    // flag = CVSpilsSetJacTimes(cvode_mem, NULL,Jac);
    // if(check_flag(&flag, "CVSpilsSetJacTimes", 1)) return(NULL);


    return cvode_mem;

}



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
