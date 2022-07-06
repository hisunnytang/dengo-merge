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
#include <primordial_solver.h>


#include <cvode/cvode_spils.h>           /* access to CVSpils interface       */
#include <sunlinsol/sunlinsol_spgmr.h>   /* access to SPGMR SUNLinearSolver   */


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


/* Functions Called by the Solver */

int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);

#ifndef CVSPILS
int Jac(realtype t, N_Vector y, N_Vector fy, SUNMatrix J,
        void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
#endif

#ifdef CVSPILS
int Jac(N_Vector v, N_Vector Jv, realtype t,
             N_Vector y, N_Vector fy,
             void *user_data, N_Vector tmp);
#endif

static int check_flag(void *flagvalue, const char *funcname, int opt);

int cvode_solver( void *cvode_mem, double *output, int NEQ, double *dt, primordial_data * data, N_Vector y , double reltol, N_Vector abstol){

    int flag, cvode_flag, i;

    flag = CVodeReInit( cvode_mem, 0.0, y);
    flag = CVodeSVtolerances(cvode_mem, reltol, abstol  );

    double tout = dt[0];
    cvode_flag = CVode( cvode_mem, tout, y, dt, CV_NORMAL);

    for ( i = 0; i < NEQ; i++){
       output[i] = NV_Ith_S(y, i);
    }

    #ifdef PRINT_CVODE_STATS
    long int nsteps, nfevals, nlinsetups, netfails;
    int qlast, qcur;
    realtype hinused, hlast, hcur, tcur;
    flag = CVodeGetIntegratorStats(cvode_mem, &nsteps, &nfevals, &nlinsetups, &netfails, &qlast, &qcur, &hinused, &hlast, &hcur, &tcur);

    fprintf(stderr, "-----Printing Integrator Stats------- \n");
    fprintf(stderr, "nsteps    : %ld \n", nsteps);
    fprintf(stderr, "nfevals   : %ld \n", nfevals);
    fprintf(stderr, "nlinsetups: %ld \n", nlinsetups);
    fprintf(stderr, "netfails  : %ld \n", netfails);

    N_Vector ele, eweight;
    ele = NULL;
    eweight = NULL;
    ele     = N_VNew_Serial( NEQ );
    eweight = N_VNew_Serial( NEQ );

    flag = CVodeGetEstLocalErrors(cvode_mem, ele);
    flag = CVodeGetErrWeights(cvode_mem, eweight);
    fprintf(stderr, "-----Printing Local Errors   ------- \n");
    for ( i = 0; i < NEQ; i++){
        fprintf(stderr, "Local Error[ %d ] = %0.5g \n", i, NV_Ith_S( ele, i) );
        fprintf(stderr, "Error Weight      = %0.5g \n", NV_Ith_S(eweight, i) );
        fprintf(stderr, "contributions to error test = %0.5g \n",  NV_Ith_S( ele, i) * NV_Ith_S(eweight, i) );
        fprintf(stderr, "-------------------------------\n");
    }
    #endif

    if (cvode_flag < -1){
        return 1;
    }

    return 0;

}

void *setup_cvode_solver( rhs_f f, jac_f Jac,  int NEQ,
        primordial_data *data, SUNLinearSolver LS, SUNMatrix A, N_Vector y, double reltol, N_Vector abstol){

    void *cvode_mem;
    cvode_mem = NULL;

    int i, flag;

    /* Create CVODES object */
    cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
    if (check_flag((void *)cvode_mem, "CVodeCreate", 0)) return(NULL);


    /* Allocate space for CVODES */
    flag = CVodeInit(cvode_mem, f, 0.0, y);
    if (check_flag( &flag, "CVodeInit", 1)) return(NULL);

    flag = CVodeSetMaxNumSteps(cvode_mem, 1000 );
    flag = CVodeSetStabLimDet(cvode_mem, SUNTRUE);

    flag = CVodeSVtolerances(cvode_mem, reltol, abstol);
    if (check_flag(&flag, "CVodeSStolerances", 1)) return(NULL);

    /* Attach User Data */
    flag = CVodeSetUserData(cvode_mem, data);
    if (check_flag(&flag, "CVodeSetUserData", 1)) return(NULL);

    #ifndef CVSPILS
    /* Call CVDlsSetLinearSolver to attach the matrix
     * and linear solver to CVode */
    flag = CVDlsSetLinearSolver(cvode_mem, LS, A);
    if(check_flag(&flag, "CVDlsSetLinearSolver", 1)) return(NULL);

    /* Set the user-supplied Jacobian routine Jac */
    flag = CVDlsSetJacFn(cvode_mem, Jac);
    if(check_flag(&flag, "CVDlsSetJacFn", 1)) return(NULL);
    #endif

    #ifdef CVSPILS
    LS = SUNSPGMR(y, PREC_NONE, 0);
    if(check_flag(&flag, "SUNSPGMR", 1)) return(NULL);


    flag = CVSpilsSetLinearSolver(cvode_mem, LS);
    if(check_flag(&flag, "CVSpilsSetLinearSolver", 1)) return(NULL);
    /* Set the JAcobian-times-vector function */
    // flag = CVSpilsSetJacTimes(cvode_mem, NULL,Jac);
    // if(check_flag(&flag, "CVSpilsSetJacTimes", 1)) return(NULL);

    #endif



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
