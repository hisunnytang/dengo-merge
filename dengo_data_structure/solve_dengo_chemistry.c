#include "dengo.h"
#include "dengo_chemistry_data.h"
#include "dengo_types.h"


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
    



