/* header files for CVODES/SUNDIALS */
#include <cvodes/cvodes.h>           /* prototypes for CVODE fcts. and consts. */
#include <nvector/nvector_serial.h>  /* serial N_Vector types, fcts., and macros */
#include <cvodes/cvodes_dense.h>     /* prototype for CVDense */
#include <sundials/sundials_dense.h> /* definitions DlsMat DENSE_ELEM */
#include <sundials/sundials_types.h> /* definition of type realtype */


int cvodes_main_solver( rhs_f f, jac_f jac, 
                 double *input, double *rtol, double *atol,
                 int nchem, void *sdata, double t0, double t1, double *t_now, void *cvode_mem);
