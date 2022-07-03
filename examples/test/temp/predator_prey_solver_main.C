#include "predator_prey_solver.h"

int solve_chemistry(int argc, char **argv);
int solve_chemistry_enzo(int argc, char **argv);

// Handy Function to compare identical array
int compare_array(double *array, double rtol, unsigned long N){
    /* Check if the array element are the same within tolerance level */
    double ref, diff;
    ref = array[0];
    for (unsigned long i = 0; i < N; i++){
        diff = fabs(ref - array[i])/ref;
        if (diff > rtol){
            printf("Exceeded tolerance level (%0.5g); ref = %0.5g; array[%lu] = %0.5g\n", rtol, ref, i, array[i]);
            return 1;
        }
    }
    return 0;
}

int main(int argc, char **argv) {
    /*
    if (argc > 1){
       predator_prey_main(argc, argv);
       return 0;
    }
    */
    // solve_chemistry(argc, argv);
    solve_chemistry_enzo(argc, argv);
}

// Sample on how to use dengo with primordial chemistry
// with predator_prey_solve_chemistry_enzo
int solve_chemistry_enzo(int argc, char **argv) {
    dengo_field_data *field_data = (dengo_field_data *) malloc(sizeof(dengo_field_data));
    int Nx = 32;
    int Ny = 32;
    int Nz = 32;

    int N = Nx*Ny*Nz;
    field_data->ncells = N; 
    double density = 1.0e0; // in cm^-3
    double T = 1000.0; // in K
    double mH, k, tiny;

    mH = 1.67e-24;
    k  = 1.380e-16;
    tiny = 1.0e-20; 
    density *= mH;
    double G = 6.67259e-8;
    code_units *units = (code_units *) malloc(sizeof(code_units));
    units->density_units = 1.0;
    units->length_units = 1.0;
    units->time_units = 1.0;
    units->velocity_units = 1.0;
    double *dead_predator_density = (double*) malloc(N * sizeof(double));
    double *dead_prey_density = (double*) malloc(N * sizeof(double));
    double *ge_density = (double*) malloc(N * sizeof(double));
    double *predator_density = (double*) malloc(N * sizeof(double));
    double *prey_density = (double*) malloc(N * sizeof(double));
    double *cooling_time = (double *) malloc( N * sizeof(double) );
    double *gamma = (double * ) malloc( N * sizeof(double) );
    double *temperature = (double *) malloc( N * sizeof(double) );
    double *mean_molecular_weight = (double *) malloc( N * sizeof(double) );
    double *density_arr = (double *) malloc( N * sizeof(double) );

    double reltol;
    double *abstol = (double*) malloc(sizeof(double)* N * 5);
    reltol = 1.0e-5;

    for ( int i = 0; i < N; i++){
        predator_density[i] = tiny*density;
        dead_prey_density[i] = tiny*density;
        prey_density[i] = tiny*density;
        dead_predator_density[i] = tiny*density;
        H2_1_density[i] = 1.0e-5 * density;
        H_1_density[i]  = 0.76   * density;
        He_1_density[i] = 0.24   * density;
	    de_density[i]    = 1.0e-5 * density;
	    H_2_density[i]   = 1.0e-5 * density;

	    // ge ~ nkT / (gamma - 1)/ rho; gamaa ~ 5/3
        ge_density[i]   = 3.0/2.0 * k * T / mH;
	    density_arr[i] = (1+2.0*1e-5)*density;
    }

    for ( int i = 0; i < N * 5; i++ ){
    	abstol[i] = tiny * reltol;
    }
    field_data->predator_density = predator_density;
    field_data->dead_prey_density = dead_prey_density;
    field_data->ge_density = ge_density;
    field_data->prey_density = prey_density;
    field_data->dead_predator_density = dead_predator_density;

    field_data->density         = density_arr;
    field_data->CoolingTime     = cooling_time;
    field_data->Gamma           = gamma;
    field_data->temperature     = temperature;
    field_data->MolecularWeight = mean_molecular_weight;

    int gstart[3];
    int gend[3];
    int gd[3];

    gstart[0] = 0;
    gstart[1] = 0;
    gstart[2] = 0;

    gend[0] = Nx - 1;
    gend[1] = Ny - 1;
    gend[2] = Nz - 1;

    gd[0] = Nx;
    gd[1] = Ny;
    gd[2] = Nz;

    field_data->grid_start     = &gstart[0];
    field_data->grid_end       = &gend[0];
    field_data->grid_dimension = &gd[0];


    const char *fileloc = "/home/kwoksun2/dengo_install/predator_prey_tables.h5";
    field_data->dengo_data_file = fileloc;
    field_data->reltol = reltol;

    units->a_value = 1.0;
    units->a_units = 1.0;

    double dt = 1.0 / sqrt(G * density) ;
    fprintf(stderr, "MAX_NCELLS = %d \n", MAX_NCELLS);
    predator_prey_solve_chemistry_enzo( units, field_data, dt );
    dengo_estimate_cooling_time_enzo( units, field_data);
    fprintf(stderr, "dead_predator = %0.5g\n", field_data->dead_predator_density[0] / mH );
    fprintf(stderr, "dead_prey = %0.5g\n", field_data->dead_prey_density[0] / mH );
    fprintf(stderr, "ge = %0.5g\n", field_data->ge_density[0] );
    fprintf(stderr, "predator = %0.5g\n", field_data->predator_density[0] / mH );
    fprintf(stderr, "prey = %0.5g\n", field_data->prey_density[0] / mH );
    fprintf(stderr, "CoolingTime = %0.5g\n", field_data->CoolingTime[0]);

    /*
    double *Pressure = (double *) malloc(sizeof(double)*N);
    double *Temperature = (double *) malloc(sizeof(double)*N);
    double *Gamma  = (double *) malloc(sizeof(double)*N);
    field_data->Gamma       = Gamma;
    field_data->Pressure    = Pressure;
    field_data->temperature = Temperature;

    dengo_calculate_pressure_enzo   (units, field_data);
    dengo_calculate_temperature_enzo(units, field_data);
    dengo_calculate_gamma_enzo      (units, field_data);


    for (int i = 0; i < 1; i++){
    fprintf(stderr, "Gamma    = %0.5g ", field_data->Gamma[i]);
    fprintf(stderr, "Pressure = %0.5g ", field_data->Pressure[i]);
    fprintf(stderr, "Temperature = %0.5g\n", field_data->temperature[i]);
    }

    fprintf(stderr, "\nGamma    = %0.5g ", field_data->Gamma[0]);
    compare_array(field_data->Gamma, 1.0e-4, N);
    fprintf(stderr, "\nPressure = %0.5g ", field_data->Pressure[i]);
    compare_array(field_data->Pressure, 1.0e-4, N);
    fprintf(stderr, "\n Temperature = %0.5g\n", field_data->temperature[i]);
    compare_array(field_data->temperature, 1.0e-4, N);
    */

    unsigned long d;
    // lets just compare everything!!!!!!
    double ref0, frac;
    if (compare_array(field_data->predator_density, reltol, N) == 1)
        printf("predator is not consistent\n");
    if (compare_array(field_data->dead_prey_density, reltol, N) == 1)
        printf("dead_prey is not consistent\n");
    if (compare_array(field_data->ge_density, reltol, N) == 1)
        printf("ge is not consistent\n");
    if (compare_array(field_data->prey_density, reltol, N) == 1)
        printf("prey is not consistent\n");
    if (compare_array(field_data->dead_predator_density, reltol, N) == 1)
        printf("dead_predator is not consistent\n");if (compare_array(field_data->CoolingTime, reltol, N) == 1)
        printf("CoolingTime is not consistent\n");

    free(field_data);
    free(dead_predator_density);
    free(dead_prey_density);
    free(ge_density);
    free(predator_density);
    free(prey_density);
    free(cooling_time);
    free(gamma);
    free(temperature);
    free(mean_molecular_weight);
    free(abstol);
    free(units);
    free(density_arr);
}
