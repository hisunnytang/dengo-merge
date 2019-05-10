#include "gpu_memory.cuh"
/*
size_t required_mechanism_size() {
  //returns the total required size for the mechanism per thread
  size_t mech_size = 0;
  //y
  mech_size += NSP;
  //dy
  mech_size += NSP;
  //conc
  mech_size += NSP;
  //fwd_rates
  mech_size += FWD_RATES;
  //rev_rates
  mech_size += REV_RATES;
  //spec_rates
  mech_size += NSP;
  //cp
  mech_size += NSP;
  //h
  mech_size += NSP;
  //dBdT
  mech_size += NSP;
  //jac
  mech_size += NSP * NSP;
  //var
  mech_size += 1;
  //pres_mod
  mech_size += PRES_MOD_RATES;
  //y_device
  mech_size += NSP;
  //pres_device
  mech_size += 1;
  return mech_size * sizeof(double);
}
*/

size_t required_mechanism_size() {
  //returns the total required size for the mechanism per thread
  size_t mech_size = 0;
  //y (9 species + ge)
  mech_size += NSP;
  //dy
  mech_size += NSP;
 //reaction_rates
  mech_size += REACTION_RATES;
  //cooling_rates
  mech_size += COOLING_RATES;
  // drate_dT
  mech_size += REACTION_RATES;
  // dcrate_dT
  mech_size += COOLING_RATES;
  //jac; which is not used here YET!
  mech_size += NSP * NSP;
  //temperature
  mech_size += 1;
  // density
  mech_size += 1;
  //y_device
  mech_size += NSP;
  // y scale
  mech_size += NSP;
  // y inv_scale
  mech_size += NSP;
  // temp_array
  mech_size += NSP;
  // work1
  mech_size += NSP;
  // work2
  mech_size += NSP;
  return mech_size * sizeof(double);
}


void initialize_gpu_memory(int padded, mechanism_memory** h_mem, mechanism_memory** d_mem, cvklu_data** h_chem_data)
{
  // Allocate storage for the device struct
  cudaErrorCheck( cudaMalloc(d_mem, sizeof(mechanism_memory)) );

  // Allocate the device arrays on the host pointer
  cudaErrorCheck( cudaMalloc(&((*h_mem)->y),  NSP * padded * sizeof(double)) );
  cudaErrorCheck( cudaMalloc(&((*h_mem)->dy), NSP * padded * sizeof(double)) );
  cudaErrorCheck( cudaMalloc(&((*h_mem)->reaction_rates), REACTION_RATES * padded * sizeof(double)) );
  cudaErrorCheck( cudaMalloc(&((*h_mem)->cooling_rates),  COOLING_RATES  * padded * sizeof(double)) );
  cudaErrorCheck( cudaMalloc(&((*h_mem)->temperature),    1 * padded * sizeof(double)) );
  cudaErrorCheck( cudaMalloc(&((*h_mem)->density),    1 * padded * sizeof(double)) );
  cudaErrorCheck( cudaMalloc(&((*h_mem)->scale ), NSP * padded * sizeof(double) ));
  cudaErrorCheck( cudaMalloc(&((*h_mem)->inv_scale ), NSP * padded* sizeof(double) ));
  cudaErrorCheck( cudaMalloc(&((*h_mem)->temp_array ), NSP * padded* sizeof(double) ));
cudaErrorCheck( cudaMalloc(&((*h_mem)->work1), NSP * padded* sizeof(double) ));
cudaErrorCheck( cudaMalloc(&((*h_mem)->work2), NSP * padded* sizeof(double) ));




  cudaErrorCheck( cudaMalloc(&((*h_mem)->var), 1 * padded * sizeof(double)) );
  cudaErrorCheck( cudaMalloc(&((*h_mem)->chemistry_data), sizeof( cvklu_data ) )); 
  //jacobian is not needed yet: 
  cudaErrorCheck( cudaMalloc(&((*h_mem)->jac), NSP * NSP * padded * sizeof(double)) );
  cudaErrorCheck( cudaMalloc(&((*h_mem)->drrate_dT), REACTION_RATES * padded * sizeof(double)) );
  cudaErrorCheck( cudaMalloc(&((*h_mem)->dcrate_dT), COOLING_RATES  * padded * sizeof(double)) );
  cudaErrorCheck( cudaMalloc(&((*h_mem)->dTs_ge), padded * sizeof(double)) );
  cudaErrorCheck( cudaMalloc(&((*h_mem)->h2_optical_depth_approx), padded * sizeof(double)) );

  // calls to rhs func and jacobian
  cudaErrorCheck( cudaMalloc(&((*h_mem)->rhs_call), sizeof(int)));
  cudaErrorCheck( cudaMalloc(&((*h_mem)->jac_call), sizeof(int)));

 
  // initialize the arrays with values
  cudaErrorCheck( cudaMemset((*h_mem)->y,1.0, NSP * padded * sizeof(double)) );
  cudaErrorCheck( cudaMemset((*h_mem)->dy, 0, NSP * padded * sizeof(double)) );

  cudaErrorCheck( cudaMemset((*h_mem)->drrate_dT,0.0, REACTION_RATES * padded * sizeof(double)) );
  cudaErrorCheck( cudaMemset((*h_mem)->dcrate_dT,0.0, COOLING_RATES * padded * sizeof(double)) );


  cudaErrorCheck( cudaMemset((*h_mem)->rhs_call, 0, sizeof(int) )) ;
  cudaErrorCheck( cudaMemset((*h_mem)->jac_call, 0, sizeof(int) )) ;

  cudaErrorCheck( cudaMemset((*h_mem)->temperature, 0, 1 * padded * sizeof(double)) );
  cudaErrorCheck( cudaMemcpy((*h_mem)->chemistry_data, *h_chem_data, sizeof(cvklu_data), cudaMemcpyHostToDevice ));
 
 
  cudaErrorCheck( cudaMemcpy(*d_mem, *h_mem, sizeof(mechanism_memory), cudaMemcpyHostToDevice) );


}



void free_gpu_memory(mechanism_memory** h_mem, mechanism_memory** d_mem)
{
  cudaErrorCheck(cudaFree((*h_mem)->y));
  cudaErrorCheck(cudaFree((*h_mem)->dy));
  cudaErrorCheck(cudaFree((*h_mem)->reaction_rates));
  cudaErrorCheck(cudaFree((*h_mem)->cooling_rates));
  cudaErrorCheck(cudaFree((*h_mem)->temperature));
  cudaErrorCheck(cudaFree((*h_mem)->scale));
  cudaErrorCheck(cudaFree((*h_mem)->inv_scale));
  cudaErrorCheck(cudaFree(*d_mem));
}
