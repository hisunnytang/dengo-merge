#include "cvklu_solver.h"

int main(int argc, char **argv) {
    //cvklu_main(argc, argv);

    
    dengo_field_data *field_data = (dengo_field_data *) malloc(sizeof(dengo_field_data));
    
    int N = 64*64*64;
    field_data->nstrip = N; 
    double density = 1.0e3; // in cm^-3
    double T = 3000.0; // in K
    double mH, k, tiny;

    mH = 1.67e-24;
    k  = 1.3806488e-16;
    tiny = 1.0e-20; 
    
//    density *= mH;
    
    code_units *units = (code_units *) malloc(sizeof(code_units));
    units->density_units = 1.0;
    units->length_units = 1.0;
    units->time_units = 1.0;
    units->velocity_units = 1.0;

    double *H2_1_density  = (double *) malloc(N * sizeof(double));
    double *H2_2_density= (double *) malloc(N * sizeof(double));
    double *H_1_density= (double *) malloc(N * sizeof(double));
    double *H_2_density= (double *) malloc(N * sizeof(double));
    double *H_m0_density= (double *) malloc(N * sizeof(double));
    double *He_1_density= (double *) malloc(N * sizeof(double));
    double *He_2_density= (double *) malloc(N * sizeof(double));
    double *He_3_density= (double *) malloc(N * sizeof(double));
    double *de_density= (double *) malloc(N * sizeof(double));
    double *ge_density= (double *) malloc(N * sizeof(double));
    
    for ( int i = 0; i < N; i++){
    H2_1_density[i] = tiny * density;
    H2_2_density[i] = tiny * density;
    H_1_density[i] = 0.24 * density;
    H_2_density[i] = tiny * density;
    H_m0_density[i] = tiny * density;
    He_1_density[i] = 0.76 * density;
    He_2_density[i] = tiny * density;
    He_3_density[i] = tiny * density;
    de_density[i] = tiny * density;
    ge_density[i] = 3.0 / 2.0 * k * T / mH;   
    }

    field_data->H2_1_density = H2_1_density;
    field_data->H2_2_density = H2_2_density;
    field_data->H_1_density = H_1_density;
    field_data->H_2_density = H_2_density;
    field_data->H_m0_density = H_m0_density;
    field_data->He_1_density = He_1_density;
    field_data->He_2_density = He_2_density;
    field_data->He_3_density = He_3_density;
    field_data->de_density = de_density;
    field_data->ge_density = ge_density;

    double dt = 299204.0 * 1.0e10 / sqrt(density) ;
    fprintf(stderr, "MAX_NCELLS = %d \n", MAX_NCELLS);
    cvklu_solve_chemistry_dt( units, field_data, dt );


    

}
