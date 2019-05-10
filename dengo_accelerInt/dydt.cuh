#ifndef DYDT_HEAD
#define DYDT_HEAD

#include "header.cuh"

__device__ void dydt (const double, const double, const double * __restrict__, double * __restrict__, const mechanism_memory * __restrict__);

__device__ void interpolate_gamma( cvklu_data *, double , double*, double* );
__device__ void evaluate_temperature( double*, double*, const double *, const double, cvklu_data* );
__device__ void interpolate_reaction_rates( double*, double, cvklu_data* );
__device__ void interpolate_cooling_rates( double*, double, cvklu_data*);
__device__ void interpolate_dcrate_dT(double*, const double, cvklu_data* );
__device__ void interpolate_drrate_dT(double*, const double, cvklu_data*);

#endif
