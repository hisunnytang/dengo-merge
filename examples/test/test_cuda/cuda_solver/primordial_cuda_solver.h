#include "hdf5.h"
#include "hdf5_hl.h"
#include <stdio.h>
#include <cvode/cvode.h>                  /* prototypes for CVODE fcts., consts.          */
#include <nvector/nvector_cuda.h>         /* access to cuda N_Vector                      */
#include <sunmatrix/sunmatrix_cusparse.h>             /* access to cusparse SUNMatrix                  */
#include <sunlinsol/sunlinsol_cusolversp_batchqr.h>   /* acess to cuSolverSp batch QR SUNLinearSolver */
#include <sundials/sundials_types.h>     /* defs. of realtype, int              */

#define ZERO    RCONST(0.0)
#define kb      RCONST(1.3806504e-16)
#define mh      RCONST(1.67e-24)
#define gamma   RCONST(5.0/3.0)
#define _gamma_m1 RCONST(1.0/ (gamma-1.0) )
#define nchem 10

#define BATCHSIZE 2048
#define BLOCKSIZE 512
#define GRIDSIZE (BATCHSIZE/BLOCKSIZE)

#define T0 RCONST(0.0)
#define T1 RCONST(1e10)
#define TMULT RCONST(10.0)
#define NOUT 12

typedef struct code_units
{
    int comoving_coordinates;
    double density_units;
    double length_units;
    double time_units;
    double velocity_units;
    double a_units;
    double a_value;

} code_units;

// define our datatype
typedef struct
{
    double nbins;
    double dbin;
    double idbin;
    double lb;
    double ub;

    double current_z;
    double *Ts;
    double *logTs;
    double *Tdef;
    double *dTs_ge;
    double *Tge;

    // cooling & chemical tables
    double *r_k01;
    double *rs_k01;
    double *drs_k01;
    double *r_k02;
    double *rs_k02;
    double *drs_k02;
    double *r_k03;
    double *rs_k03;
    double *drs_k03;
    double *r_k04;
    double *rs_k04;
    double *drs_k04;
    double *r_k05;
    double *rs_k05;
    double *drs_k05;
    double *r_k06;
    double *rs_k06;
    double *drs_k06;
    double *r_k07;
    double *rs_k07;
    double *drs_k07;
    double *r_k08;
    double *rs_k08;
    double *drs_k08;
    double *r_k09;
    double *rs_k09;
    double *drs_k09;
    double *r_k10;
    double *rs_k10;
    double *drs_k10;
    double *r_k11;
    double *rs_k11;
    double *drs_k11;
    double *r_k12;
    double *rs_k12;
    double *drs_k12;
    double *r_k13;
    double *rs_k13;
    double *drs_k13;
    double *r_k14;
    double *rs_k14;
    double *drs_k14;
    double *r_k15;
    double *rs_k15;
    double *drs_k15;
    double *r_k16;
    double *rs_k16;
    double *drs_k16;
    double *r_k17;
    double *rs_k17;
    double *drs_k17;
    double *r_k18;
    double *rs_k18;
    double *drs_k18;
    double *r_k19;
    double *rs_k19;
    double *drs_k19;
    double *r_k21;
    double *rs_k21;
    double *drs_k21;
    double *r_k22;
    double *rs_k22;
    double *drs_k22;
    double *r_k23;
    double *rs_k23;
    double *drs_k23;
    double *c_brem_brem;
    double *cs_brem_brem;
    double *dcs_brem_brem;
    double *c_ceHeI_ceHeI;
    double *cs_ceHeI_ceHeI;
    double *dcs_ceHeI_ceHeI;
    double *c_ceHeII_ceHeII;
    double *cs_ceHeII_ceHeII;
    double *dcs_ceHeII_ceHeII;
    double *c_ceHI_ceHI;
    double *cs_ceHI_ceHI;
    double *dcs_ceHI_ceHI;
    double *c_cie_cooling_cieco;
    double *cs_cie_cooling_cieco;
    double *dcs_cie_cooling_cieco;
    double *c_ciHeI_ciHeI;
    double *cs_ciHeI_ciHeI;
    double *dcs_ciHeI_ciHeI;
    double *c_ciHeII_ciHeII;
    double *cs_ciHeII_ciHeII;
    double *dcs_ciHeII_ciHeII;
    double *c_ciHeIS_ciHeIS;
    double *cs_ciHeIS_ciHeIS;
    double *dcs_ciHeIS_ciHeIS;
    double *c_ciHI_ciHI;
    double *cs_ciHI_ciHI;
    double *dcs_ciHI_ciHI;
    double *c_compton_comp_;
    double *cs_compton_comp_;
    double *dcs_compton_comp_;
    double *c_gammah_gammah;
    double *cs_gammah_gammah;
    double *dcs_gammah_gammah;
    double *c_gloverabel08_gael;
    double *cs_gloverabel08_gael;
    double *dcs_gloverabel08_gael;
    double *c_gloverabel08_gaH2;
    double *cs_gloverabel08_gaH2;
    double *dcs_gloverabel08_gaH2;
    double *c_gloverabel08_gaHe;
    double *cs_gloverabel08_gaHe;
    double *dcs_gloverabel08_gaHe;
    double *c_gloverabel08_gaHI;
    double *cs_gloverabel08_gaHI;
    double *dcs_gloverabel08_gaHI;
    double *c_gloverabel08_gaHp;
    double *cs_gloverabel08_gaHp;
    double *dcs_gloverabel08_gaHp;
    double *c_gloverabel08_gphdl;
    double *cs_gloverabel08_gphdl;
    double *dcs_gloverabel08_gphdl;
    double *c_gloverabel08_gpldl;
    double *cs_gloverabel08_gpldl;
    double *dcs_gloverabel08_gpldl;
    double *c_gloverabel08_h2lte;
    double *cs_gloverabel08_h2lte;
    double *dcs_gloverabel08_h2lte;
    double *c_h2formation_h2mcool;
    double *cs_h2formation_h2mcool;
    double *dcs_h2formation_h2mcool;
    double *c_h2formation_h2mheat;
    double *cs_h2formation_h2mheat;
    double *dcs_h2formation_h2mheat;
    double *c_h2formation_ncrd1;
    double *cs_h2formation_ncrd1;
    double *dcs_h2formation_ncrd1;
    double *c_h2formation_ncrd2;
    double *cs_h2formation_ncrd2;
    double *dcs_h2formation_ncrd2;
    double *c_h2formation_ncrn;
    double *cs_h2formation_ncrn;
    double *dcs_h2formation_ncrn;
    double *c_reHeII1_reHeII1;
    double *cs_reHeII1_reHeII1;
    double *dcs_reHeII1_reHeII1;
    double *c_reHeII2_reHeII2;
    double *cs_reHeII2_reHeII2;
    double *dcs_reHeII2_reHeII2;
    double *c_reHeIII_reHeIII;
    double *cs_reHeIII_reHeIII;
    double *dcs_reHeIII_reHeIII;
    double *c_reHII_reHII;
    double *cs_reHII_reHII;
    double *dcs_reHII_reHII;

    // for now, we ignore the temperature dependent Gamma
    /*
       double *g_gamma;
       double *g_dgamma_dT;
       double *gamma;
       double *dgamma_dT;
       double *_gamma_dT;
       double *g_gamma;
       double *g_dgamma_dT;
       double *gamma;
       double *dgamma_dT;
       double *_gamma_dT;
     */
    double *cie_optical_depth_approx;
    double *h2_optical_depth_approx;


} primordial_cuda_data;

typedef struct dengo_field_data
{
    unsigned long int nstrip;
    unsigned long int ncells;
    double reltol;
    double floor_value;

    double *density;
    double *H2_1_density;
    double *H2_2_density;
    double *H_1_density;
    double *H_2_density;
    double *H_m0_density;
    double *He_1_density;
    double *He_2_density;
    double *He_3_density;
    double *de_density;
    double *ge_density;

    double *CoolingTime;
    double *MolecularWeight;
    double *temperature;
    double *Gamma;

    int *grid_start;
    int *grid_end;
    int *grid_dimension;

    const char *dengo_data_file;
    code_units *units;
} dengo_field_data;


// read rate data
primordial_cuda_data primordial_cuda_setup_data(int *NumberOfFields, char ***FieldNames);

void primordial_cuda_read_rate_tables(primordial_cuda_data *data);

void primordial_cuda_read_cooling_tables(primordial_cuda_data *data);


// cuda kernels
__global__ void linear_interpolation_kernel(primordial_cuda_data data);

__global__ static void rhs_kernel(double y, double *ydata, double *ydotdata, primordial_cuda_data data);


__global__ void temperature_kernel(double* ydata, primordial_cuda_data data);

__global__ static void jacobian_kernel(realtype *ydata, realtype *Jdata, primordial_cuda_data data);


// Function Called by the solver
static int calculate_rhs_primordial_cuda(realtype t, N_Vector y, N_Vector ydot, void *user_data);
static int calculate_jacobian_primordial_cuda(realtype t, N_Vector y, N_Vector fy, SUNMatrix J, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
static int calculate_sparse_jacobian_primordial_cuda(realtype t, N_Vector y, N_Vector fy, SUNMatrix J,
        void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
static int blockJacInit(SUNMatrix J);
static int JacInit(SUNMatrix J);
static int blockSparseJacInit(SUNMatrix J);



// test functions

void test_rhs_function(primordial_cuda_data data);
void test_jacobian_function(primordial_cuda_data data);
int grid_performance();
int test_scaling_dims();
int run_dengo_struct(double density, double temperature, double h2fraction, double efraction, unsigned long dims, double *output);


// cvode helper function
static int check_retval(void *returnvalue, const char *funcname, int opt);
static void PrintFinalStats(void *cvode_mem, SUNLinearSolver LS);

// dengo run solver
int run_solver(int argc, char *argv[]);
void *setup_cvode_cuda_solver(CVRhsFn f, CVLsJacFn Jac, int NEQ, primordial_cuda_data *data, SUNLinearSolver LS, SUNMatrix A, N_Vector y, double reltol, N_Vector abstol, cusparseHandle_t *cusp_handle, cusolverSpHandle_t *cusol_handle);
int run_dengo_solver(double density, double temperature, double h2fraction, double efraction, unsigned long dims);

int dengo_evolve_primordial_cuda (double dtf, double &dt, double z, double *input,
        double *rtol, double *atol, unsigned long dims, primordial_cuda_data *data, double *temp_array );
void initialize_long_ydata(double *ydata, unsigned long NSYSTEM, double density, double temperature, double h2fraction, double efraction);
int run_dengo_solver(double density, double temperature, double h2fraction, double efraction, unsigned long dims);


void launchInterpolationKernel(primordial_cuda_data*);
void launchTemperatureKernel(primordial_cuda_data*);
void launchRhsKernel(primordial_cuda_data*);
void launchJacobianKernel(primordial_cuda_data *data);
