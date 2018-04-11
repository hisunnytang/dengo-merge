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

  double grid_dx;
    
  // This should be updated dynamically 
  // with dengo
  double *density;
  double *HI_density;
  double *HII_density;
  double *HM_density;
  double *HeI_density;
  double *HeII_density;
  double *HeIII_density;
  double *H2I_density;
  double *H2II_density;
  double *de_density;

  double *ge_density;


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
