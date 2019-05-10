#ifndef JACOB_HEAD
#define JACOB_HEAD

#include "header.cuh"
#include "dydt.cuh"
#include "gpu_memory.cuh"

__device__ void eval_jacob (const double, const double, const double * __restrict__, double * __restrict__, const mechanism_memory * __restrict__);
__device__ void eval_jacob_accurate (const double, const double, const double * __restrict__, double * __restrict__ , const mechanism_memory *);
#endif
