/***********************************************************************
/
/ Grackle variable types
/
/
/ Copyright (c) 2013, Enzo/Grackle Development Team.
/
/ Distributed under the terms of the Enzo Public Licence.
/
/ The full license is in the file LICENSE, distributed with this 
/ software.
************************************************************************/

#ifndef __GRACKLE_TYPES_H_
#define __GRACKLE_TYPES_H_
/***********************************************************************
/  
/ VARIABLE TYPES
/
************************************************************************/

typedef struct
{

  int grid_rank;
  int *grid_dimension;
  int *grid_start;
  int *grid_end;
  int ncells;
  double grid_dx;
    
  // This should be updated dynamically 
  // with dengo
  double *density;
  double *H_1_density;
  double *H_2_density;
  double *H_m0_density;
  double *He_1_density;
  double *He_2_density;
  double *He_3_density;
  double *H2_1_density;
  double *H2_2_density;
  double *DI_density;
  double *DII_density;
  double *HDI_density;
  double *de_density;
  double *ge_density;

  double *internal_energy;
  double *x_velocity;
  double *y_velocity;
  double *z_velocity;

  double *volumetric_heating_rate;
  double *specific_heating_rate;

  double *RT_heating_rate;
  double *RT_HI_ionization_rate;
  double *RT_HeI_ionization_rate;
  double *RT_HeII_ionization_rate;
  double *RT_H2_dissociation_rate;

} dengo_field_data;

typedef struct
{

  int comoving_coordinates;
  double density_units;
  double length_units;
  double time_units;
  double velocity_units;
  double a_units;
  double a_value;

} code_units;

#endif
