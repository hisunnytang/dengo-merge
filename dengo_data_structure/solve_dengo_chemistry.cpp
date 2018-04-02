#include "dengo.h"
#include "dengo_chemistry_data.h"
#include "dengo_types.h"


/* header files for CVODES/SUNDIALS */
#include <cvode/cvode.h>               /* prototypes for CVODE fcts., consts.  */
#include <nvector/nvector_serial.h>    /* access to serial N_Vector            */
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix            */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver      */
#include <cvode/cvode_direct.h>        /* access to CVDls interface            */
#include <sundials/sundials_types.h>   /* defs. of realtype, sunindextype      */
// Create a solver with dengo, maybe in python
//
//
// Initialize a dengo field data structures, 
// and code_units structures
// 
//
// Invoke the solver (solve_chemistry)
//



void Dengo_Init_FieldData( const int Che_NPG )
{
    // nothing to do if we are not using the original dengo
    if (DENGO_MODE != DENGO_MODE_ORI)  return;

    // allocate memory for storing Field data, like HI, HII
    Che_FieldData = new dengo_field_data;

    // initialization
    const int NDim = 3;

    // fields not evolving with time
    Che_FieldData->grid_rank        = NDim;
    Che_FieldData->grid_dimension   = new int [NDim];
    Che_FieldData->grid_start       = new int [NDim];

    // grid_dimension, grid_start, and grid_end are set by
    // GPU_GrackleSolver_Original() since the number of
    // patch groups advanced at a time is not a constant
    

    // fields set by Grackle_Prepare() during each time-step
    Che_FieldData->density                 = NULL;
    Che_FieldData->internal_energy         = NULL;
    Che_FieldData->grid_dx                 = NULL_REAL;

    // fields not supported yet
    Che_FieldData->x_velocity              = NULL;
    Che_FieldData->y_velocity              = NULL;
    Che_FieldData->z_velocity              = NULL;
    Che_FieldData->metal_density           = NULL;

    // fields not supported yet
    Che_FieldData->HI_density              = NULL;
    Che_FieldData->HII_density             = NULL;
    Che_FieldData->HM_density              = NULL;
    Che_FieldData->HeI_density             = NULL;
    Che_FieldData->HeII_density            = NULL;
    Che_FieldData->HeIII_density           = NULL;
    Che_FieldData->H2I_density             = NULL;
    Che_FieldData->H2II_density            = NULL;
    Che_FieldData->DI_density              = NULL;
    Che_FieldData->DII_density             = NULL;
    Che_FieldData->HDI_density             = NULL;
    Che_FieldData->e_density               = NULL;
    Che_FieldData->volumetric_heating_rate = NULL;
    Che_FieldData->specific_heating_rate   = NULL;
    Che_FieldData->RT_HI_ionization_rate   = NULL;
    Che_FieldData->RT_HeI_ionization_rate  = NULL;
    Che_FieldData->RT_HeII_ionization_rate = NULL;
    Che_FieldData->RT_H2_dissociation_rate = NULL;
    Che_FieldData->RT_heating_rate         = NULL;

} // FUNCTION : Grackle_Init_FieldData
    


void CPU_DengoSolver_Original(dengo_field_data *Che_FieldData, code_units Che_units, const int NPatchGroup, const real dt)
{
    // set grid_dimension, grid_start, and grid_end
   const int OptFac = 16;  // optimization factor
   if ( SQR(PS2)%OptFac != 0 )   Aux_Error( ERROR_INFO, "SQR(PS2) %% OptFac != 0 !!\n" );

   Che_FieldData->grid_dimension[0] = PS2*OptFac;
   Che_FieldData->grid_dimension[1] = 1;
   Che_FieldData->grid_dimension[2] = SQR(PS2)*NPatchGroup/OptFac;

   for (int d=0; d<3; d++)
   {
      Che_FieldData->grid_start[d] = 0;
      Che_FieldData->grid_end  [d] = Che_FieldData->grid_dimension[d] - 1;
   }

// invoke Grackle
// --> note that we use the OpenMP implementation in Grackle directly, which applies the parallelization to the first two
//     dimensiones of the input grid
// --> this approach is found to be much more efficient than parallelizing different patches or patch groups here
   if (  solve_chemistry( &Che_Units, Che_FieldData, dt ) == 0  )
      Aux_Error( ERROR_INFO, "Grackle solve_chemistry() failed !!\n" );

} // FUNCTION : CPU_GrackleSolver_Original


int solve_chemistry( code_units Che_Units, dengo_field_data *Che_FieldData, dt)
{   
    
    int status;
    int NTOT;

    int *rhs_f( realtype, N_Vector , N_Vector , void * );

    int *jac_f( long int, realtype, N_Vector , N_Vector , DlsMat , 
            void *, N_Vector, N_Vector, N_Vector);

    rhs_f f = calculate_rhs_cvdls_9species
    jac_f jf = calculate_jacobian_cvdls_9species
    
    NTOT = Che_FieldData->NSPECIES;

    cvdls_9species_data *rates_data = cvdls_9species_setup_data(NULL, NULL);
    
    // Loop over all grid_dimensions
    
    double *input = alloca(NTOT * sizeof(double));


    // Initialize initial temperature
    for (i = 0; i<NDIM; i++ ){
        data.Ts[i] = Che_FieldData->T[i]
    }

    for iter in range(niter):

        status = cvodes_main_solver( f, jf, input, rtol ,  atol, NSPECIES, <void *> data, dt)
        j = 0; 
        dt_local = dt[0];

        for i in range(dims):
            H2_1_int[i, iter] = input[j]
            j += 1
            H2_2_int[i, iter] = input[j]
            j += 1
            H_1_int[i, iter] = input[j]
            j += 1
            H_2_int[i, iter] = input[j]
            j += 1
            H_m0_int[i, iter] = input[j]
            j += 1
            He_1_int[i, iter] = input[j]
            j += 1
            He_2_int[i, iter] = input[j]
            j += 1
            He_3_int[i, iter] = input[j]
            j += 1
            de_int[i, iter] = input[j]
            j += 1
            ge_int[i, iter] = input[j]
            j += 1
            temp_int[i, iter] = data.Ts[i]

        if status == 0:
            result_int[iter] = 1
            ttot += dt_local
        elif status == 1:
            result_int[iter] = 0
            ttot += dt_local

        t_int[iter] = ttot
        dt_int[iter] = dt_local
        
        if status == 0:
            if iter % 100 == 0:
                print "Successful iteration[% 5i]: (%0.3e) %0.3e / %0.3e" % (iter, dt_local, ttot, tf)

            dt_local = 1.1*dt_local

            copy_array(input, prev, NTOT)
            # Reset the scaling array to match the new values
            copy_array(input, scale, NTOT)
            dt[0] = dt_local;
            if tf - ttot < dt_local:
                dt_local = tf - ttot
                dt[0] = dt_local;
        elif status == 1:
            dt[0] = dt_local/2.0;
            # copy_array(prev, input, NTOT)
            # Reset the scaling array to match the new values
            # copy_array(input, scale, NTOT)
            if dt[0] < 1e-50 * tf:
                print "dt too small (%0.3e / %0.3e) so breaking" % (dt[0], tf)
                break
            continue
        if ttot >= tf: break
    

    free(dt_arr)
    free(ttot_arr)
    free(success_arr)
    
    free(u0)
    free(s)
    free(gu)
    free(Ju)

}
