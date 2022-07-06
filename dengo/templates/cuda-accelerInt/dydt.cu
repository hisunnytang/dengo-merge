#include "header.cuh"
#include "gpu_memory.cuh"

extern __shared__ double y_shared[];

__device__ void {{solver_name}}_interpolate_gamma( {{solver_name}}_data *rate_data, double T
        {%- for sp in network.interpolate_gamma_species_name  -%}
        ,double *gamma{{sp}}, double *dgamma{{sp}}_dT
        {%- endfor} )
{
  int tid, bin_id, zbin_id;
  double t1, t2;
  double Tdef, log_temp_out;
  int no_photo = 0;
  double lb = log(rate_data->bounds[0]);

  log_temp_out = log(T);
  bin_id = (int) ( rate_data->idbin * ( log_temp_out -  lb ) );
  if ( bin_id <= 0) {
    bin_id = 0;
  } else if ( bin_id >= rate_data->nbins) {
    bin_id = rate_data->nbins - 1;
  }

  //printf( "bin_id = %d; temp_out = %0.5g \n", bin_id, temp_out[tid]);
  t1 = (lb + (bin_id    ) * rate_data->dbin);
  t2 = (lb + (bin_id + 1) * rate_data->dbin);
  Tdef = (log_temp_out - t1)/(t2 - t1);

  {%- for sp in network.interpolate_gamma_species_name%}
  *gamma{{sp}}  = rate_data->g_gamma{{sp}}[bin_id] +
        Tdef * (rate_data->g_gamma{{sp}}[bin_id+1] - rate_data->g_gamma{{sp}}[bin_id]);
  *dgamma{{sp}}_dT = data->g_dgamma{{sp}}_dT[bin_id] +
        Tdef * (rate_data->g_dgamma{{sp}}_dT[bin_id+1] - rate_data->g_dgamma{{sp}}_dT[bin_id]);
  {% endfor %}

}


__device__ void evaluate_temperature( double* T, double* dTs_ge, const double *y, const double mdensity, cvklu_data *rate_data )
{
  // note that y array is loaded into the shared memory
  // the indexing of y has to be access with S_INDEX instead

  // iterate temperature to convergence
  double t, tnew, tdiff;
  double dge, dge_dT;
  double gammaH2, dgammaH2_dT, _gammaH2_m1;

  int count = 0;
  int MAX_ITERATION = 100;
  double gamma     = 5./3.;
  double _gamma_m1 = 1.0 / (gamma - 1.0);
  double kb = 1.3806504e-16; // Boltzamann constant [erg/K]
  // prepare t, tnew for the newton's iteration;

  {% for sp in network.interpolate_gamma_species_name | sort %}
  double gamma{{sp}};
  double dgamma{{sp}}_dT;
  double _gamma{{sp}}_m1;
  {% endfor %}

  t     = *T;
  if (t != t) t = 1000.0;
  tnew  = 1.1*t;
  tdiff = tnew - t;

  while ( tdiff/ tnew > 1.0e-8 ){
    // We do Newton's Iteration to calculate the temperature
    // Since gammaH2 is dependent on the temperature too!
    interpolate_gamma( rate_data, t {%- for sp in network.interpolate_gamma_species_name  -%}
                                    ,double *gamma{{sp}}, double *dgamma{{sp}}_dT
                                    {%- endfor});

    {{ solver_name }}_interpolate_gamma(data, i);
    {%- for sp in network.interpolate_gamma_species_name | sort %}
    gamma{{sp}} = data->gamma{{sp}}[i];
    dgamma{{sp}}_dT = data->dgamma{{sp}}_dT[i];
    _gamma{{sp}}_m1 = 1.0 / (gamma{{sp}} - 1.0);
    // fprintf(stderr, ":gamma{{sp}} %0.5g , dgamma{{sp}}_dT: %.5g \n", gamma{{sp}}, dgamma{{sp}}_dT  );
    {% endfor %}


    // update gammaH2
    // The derivatives of  sum (nkT/(gamma - 1)/mh/density) - ge
    // This is the function we want to minimize
    // which should only be dependent on the first part
    dge_dT = {{network.temperature_calculation(derivative_dge_dT=True, replace_by_array = True,
               array_name = "y", array_index = "S_INDEX" )}};

    //This is the change in ge for each iteration
    dge = {{network.temperature_calculation(get_dge=True, replace_by_array = True,
              array_name = "y", array_index = "S_INDEX" )}};


    //This is the change in ge for each iteration
    tnew = t - dge/dge_dT;
    count += 1;

    tdiff = fabs(t - tnew);
    t     = tnew;
    if (count > MAX_ITERATION){
      printf("T[tid = %d] failed to converge (iteration: %d); at T = %0.3g \n", T_ID, count, tnew );
    }
    if ( t!= t && T_ID == 0){
      printf("T[tid = %d] is %0.5g, count = %d; ge = %0.5g, gamma_H2 = %0.5g \n", T_ID, t, count, y[S_INDEX(9)], gammaH2);
      t = 1000.0;
      for (int i = 0; i < 10; i++){
          printf("y[S_INDEX(%d)] = %0.5g \n", i, y[S_INDEX(i)]);
      }
      break;
    }

  }
  // update the temperature;
  *T = t;

  // update the dge_dT with the converged temp
  dge_dT = {{network.temperature_calculation(derivative_dge_dT=True, replace_by_array = True,
               array_name = "y", array_index = "S_INDEX" )}};

  *dTs_ge = 1.0 / dge_dT;

  // printf("T[tid = %d] is %0.5g, count = %d; ge = %0.5g, gamma_H2 = %0.5g \n", tid, t, count, y[INDEX(9)], gammaH2);

}


__device__ void interpolate_reaction_rates( double *reaction_rates_out, double temp_out, cvklu_data *rate_data)
{

    int tid, bin_id, zbin_id;
    double t1, t2;
    double Tdef, dT, invTs, log_temp_out;
    int no_photo = 0;
    double lb = log(rate_data->bounds[0]);

    tid = threadIdx.x + blockDim.x * blockIdx.x;

    log_temp_out = log(temp_out);
    bin_id = (int) ( rate_data->idbin * ( log_temp_out -  lb ) );
    if ( bin_id <= 0) {
        bin_id = 0;
    } else if ( bin_id >= rate_data->nbins) {
        bin_id = rate_data->nbins - 1;
    }

    //printf( "bin_id = %d; temp_out = %0.5g \n", bin_id, temp_out[tid]);
    t1 = (lb + (bin_id    ) * rate_data->dbin);
    t2 = (lb + (bin_id + 1) * rate_data->dbin);
    Tdef = (log_temp_out - t1)/(t2 - t1);
    dT    = t2 - t1;
    invTs = 1.0 / temp_out;

    // rate_out is a long 1D array
    // NRATE is the number of rate required by the solver network
    {%- for name, rate in network.reactions | dictsort %}
    reaction_rates_out[INDEX({{loop.index0}})] = rate_data->r_{{name}}[bin_id] +
            Tdef * (rate_data->r_{{name}}[bin_id+1] - rate_data->r_{{name}}[bin_id]);
    {% endfor %}

}

__device__ void interpolate_cooling_rates( double *cooling_rates_out, double temp_out, cvklu_data *rate_data)
{

    int tid, bin_id, zbin_id;
    double t1, t2;
    double Tdef, log_temp_out;
    int no_photo = 0;
    double lb = log(rate_data->bounds[0]);

    tid = threadIdx.x + blockDim.x * blockIdx.x;

    log_temp_out = log(temp_out);
    bin_id = (int) ( rate_data->idbin * ( log_temp_out -  lb ) );

/*
    if (T_ID == 0){
    printf( "bin_id = %d; temp_out = %0.5g \n", bin_id, temp_out);
    }
*/

    if ( bin_id <= 0) {
        bin_id = 0;
    } else if ( bin_id >= rate_data->nbins) {
        bin_id = rate_data->nbins - 1;
    }
    t1 = (lb + (bin_id    ) * rate_data->dbin);
    t2 = (lb + (bin_id + 1) * rate_data->dbin);
    Tdef = (log_temp_out - t1)/(t2 - t1);

    // rate_out is a long 1D array
    // NRATE is the number of rate required by the solver network
    {%- for name, rate in network.cooling_actions | dictsort %}
    {%- for name2 in rate.tables | sort %}
    cooling_rates_out[INDEX({{loop.index0}})] = rate_data->c_{{name}}_{{name2}}[bin_id] +
            Tdef * (rate_data->c_{{name}}_{{name2}}[bin_id+1] - rate_data->c_{{name}}_{{name2}}[bin_id]);
    {% endfor %}
    {% endfor%}
}

__device__ void interpolate_dcrate_dT(double *dcr_dT, const double temp_out, cvklu_data *rate_data ){
    int tid, bin_id, zbin_id;
    double t1, t2;
    double Tdef, dT, inv_Ts, log_temp_out;
    int no_photo = 0;
    double lb = log(rate_data->bounds[0]);

    tid = threadIdx.x + blockDim.x * blockIdx.x;

    log_temp_out = log(temp_out);
    bin_id = (int) ( rate_data->idbin * ( log_temp_out -  lb ) );
    if ( bin_id <= 0) {
        bin_id = 0;
    } else if ( bin_id >= rate_data->nbins) {
        bin_id = rate_data->nbins - 1;
    }

    //printf( "bin_id = %d; temp_out = %0.5g \n", bin_id, temp_out[tid]);
    t1 = (lb + (bin_id    ) * rate_data->dbin);
    t2 = (lb + (bin_id + 1) * rate_data->dbin);
    Tdef = (log_temp_out - t1)/(t2 - t1);
    dT    = t2 - t1;
    inv_Ts = temp_out;

    double _dT = 1.0 / dT / temp_out ;
    {%- for name, rate in network.cooling_actions | dictsort %}
    {%- for name2 in rate.tables | sort %}
    dcr_dT[INDEX({{loop.index0}})] = (rate_data->c_{{name}}_{{name2}}[bin_id+1] - rate_data->c_{{name}}_{{name2}}[bin_id]) *_dT ;
    {% endfor %}
    {% endfor %}

    //cie_optical_depth_approx: 26
    dcr_dT[INDEX(26)] = 0.0;

}

__device__ void interpolate_drrate_dT(double *drr_dT, const double temp_out, cvklu_data *rate_data ){
    int tid, bin_id, zbin_id;
    double t1, t2;
    double Tdef, dT, inv_Ts, log_temp_out;
    int no_photo = 0;
    double lb = log(rate_data->bounds[0]);

    tid = threadIdx.x + blockDim.x * blockIdx.x;

    log_temp_out = log(temp_out);
    bin_id = (int) ( rate_data->idbin * ( log_temp_out -  lb ) );
    if ( bin_id <= 0) {
        bin_id = 0;
    } else if ( bin_id >= rate_data->nbins) {
        bin_id = rate_data->nbins - 1;
    }

    //printf( "bin_id = %d; temp_out = %0.5g \n", bin_id, temp_out[tid]);
    t1 = (lb + (bin_id    ) * rate_data->dbin);
    t2 = (lb + (bin_id + 1) * rate_data->dbin);
    Tdef = (log_temp_out - t1)/(t2 - t1);
    dT    = t2 - t1;
    inv_Ts = temp_out;

    double _dT = 1.0 / dT / temp_out;

    {%- for name, rate in network.reactions | dictsort %}
    drr_dT[INDEX({{loop.index0}})] = (rate_data->r_{{name}}[bin_id+1] - rate_data->r_{{name}}[bin_id]) *_dT ;
    {% endfor %}

}


__device__ void dydt (const double t, const double pres, const double * __restrict__ y_in, double * __restrict__ dy, const mechanism_memory * d_mem) {


  int tid = threadIdx.x + blockDim.x * blockIdx.x;
//  int NSPECIES = 10;

  double * local_reaction_rates = d_mem->reaction_rates;
  double * local_cooling_rates  = d_mem->cooling_rates ;

  // scale related piece
  // double * y = d_mem->temp_array; // working space for scaling the variable back;

  cvklu_data *rate_data = d_mem->chemistry_data;

  // these should be retreieved from d_mem object
  double T_local  = d_mem->temperature[T_ID];
  double Tge      = d_mem->dTs_ge[T_ID];

  const double mdensity = d_mem->density[T_ID];
  const double inv_mdensity = 1.0 / mdensity;
  const double h2_optical_depth_approx = d_mem->h2_optical_depth_approx[T_ID];


  // scaling the input vector back to cgs units
  #ifdef SCALE_INPUT
  const double * __restrict__ scale = d_mem->scale;
  const double * __restrict__ inv_scale = d_mem->inv_scale;
  #pragma unroll
  for (int i = 0; i < 10; i++){
    y_shared[S_INDEX(i)] = y_in[INDEX(i)]*scale[INDEX(i)];
    // printf( "y_in[%d] = %0.5g; scale[%d] = %0.5g\n", i, y_in[INDEX(i)], i, scale[INDEX(i)] );
  }
  #else
  #pragma unroll
  for (int i = 0; i < 10; i++){
    y_shared[S_INDEX(i)] = y_in[INDEX(i)];
  }
  #endif

  evaluate_temperature ( &T_local, &Tge, y_shared, mdensity, rate_data );
  interpolate_reaction_rates( local_reaction_rates, T_local, rate_data);
  interpolate_cooling_rates ( local_cooling_rates , T_local, rate_data);

  {%- for sp in network.required_species %}
  // {{sp.name}}
  {{ network.print_ccode( sp, assign_to = "dy[INDEX({{loop.index0}})]",
          replace_species = True, species_name = "y_shared", species_index = "S_INDEX",
          replace_reaction = True, reaction_name = "local_reaction_rates", reaction_index = "INDEX",
          replace_cooling  = True, cooling_name = "local_cooling_rates", cooling_index = "INDEX") }}
  {% if species.name == "ge" %}
  dy[INDEX({{loop.index0}})] *= inv_mdensity;
  {% endif %}

  {% endfor %}

  #ifdef SCALE_INPUT
  // scaling the dydt vector back to code untis
  #pragma unroll
  for (int i = 0; i< 10; i++){
    dy[INDEX(i)] *= inv_scale[INDEX(i)];
  }
  #endif

/*
  if ( T_ID == 0 ){
    *d_mem->rhs_call += 1;
    printf("t = %0.5g; rhs_call = %d\n", t, *d_mem->rhs_call );
  }
*/

/*
  if ( T_ID == 0 ){
    printf("time = %0.5g, at temp = %0.5g\n", t, T_local);
    for (int i = 0; i< 10; i++){
      printf("from tid[%d]: dy[%d] = %0.5g, y = %0.5g at t = %0.5g \n", T_ID, i, dy[INDEX(i)], y_in[INDEX(i)], t);
    }
  }
*/

//  printf(" \n");
//  }
}
