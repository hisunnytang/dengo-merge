#include "cvklu_solver.h"

int main(int argc, char **argv) {
    if (argc > 1){
       cvklu_main(argc, argv);
       return 0;
    }

    
    
    dengo_field_data *field_data = (dengo_field_data *) malloc(sizeof(dengo_field_data));
    
    int N = 64*64;
    field_data->ncells = N; 
    double density = 1.0e2; // in cm^-3
    double T = 2000.0; // in K
    double mH, k, tiny;

    mH = 1.67e-24;
    k  = 1.3806488e-16;
    tiny = 1.0e-20; 
    
    density *= mH;
    
    code_units *units = (code_units *) malloc(sizeof(code_units));
    units->density_units = 1.0;
    units->length_units = 1.0;
    units->time_units = 1.0;
    units->velocity_units = 1.0;
    double *H2_1_density = (double*) malloc(N * sizeof(double));
    
    double *H2_2_density = (double*) malloc(N * sizeof(double));
    
    double *ge_density = (double*) malloc(N * sizeof(double));
    
    double *He_1_density = (double*) malloc(N * sizeof(double));
    
    double *H_m0_density = (double*) malloc(N * sizeof(double));
    
    double *He_3_density = (double*) malloc(N * sizeof(double));
    
    double *He_2_density = (double*) malloc(N * sizeof(double));
    
    double *H_1_density = (double*) malloc(N * sizeof(double));
    
    double *de_density = (double*) malloc(N * sizeof(double));
    
    double *H_2_density = (double*) malloc(N * sizeof(double));
    
    double *cooling_time = (double *) malloc( N * sizeof(double) );
    double *gamma = (double * ) malloc( N * sizeof(double) );
    double *temperature = (double *) malloc( N * sizeof(double) );
    double *mean_molecular_weight = (double *) malloc( N * sizeof(double) );

    for ( int i = 0; i < N; i++){
        
        
        
        H2_2_density[i] = tiny*density;
        
        
        
        
        
        
        
        H_m0_density[i] = tiny*density;
        
        
        
        He_3_density[i] = tiny*density;
        
        
        
        He_2_density[i] = tiny*density;
        
        
        
        
        
        de_density[i] = tiny*density;
        
        
        
        H_2_density[i] = tiny*density;
        
        
        H2_1_density[i] = 1.0e-3 * density;
        H_1_density[i]  = 0.76   * density;
        He_1_density[i] = 0.24   * density;
        ge_density[i]   = 3.0/2.0 * k * T / mH;
    }
    field_data->H2_1_density = H2_1_density;
    
    field_data->H2_2_density = H2_2_density;
    
    field_data->ge_density = ge_density;
    
    field_data->He_1_density = He_1_density;
    
    field_data->H_m0_density = H_m0_density;
    
    field_data->He_3_density = He_3_density;
    
    field_data->He_2_density = He_2_density;
    
    field_data->H_1_density = H_1_density;
    
    field_data->de_density = de_density;
    
    field_data->H_2_density = H_2_density;
    
    field_data->CoolingTime     = cooling_time;
    field_data->Gamma           = gamma;
    field_data->temperature     = temperature;
    field_data->MolecularWeight = mean_molecular_weight;

    const char *fileloc = "/home/kwoksun2/dengo_install/cvklu_tables.h5";
    field_data->dengo_data_file = fileloc;

    double dt = 299204.0 * 1.0e10 / sqrt(density / mH ) ;
    fprintf(stderr, "MAX_NCELLS = %d \n", MAX_NCELLS);
    cvklu_solve_chemistry_dt( units, field_data, dt );
    
    free(field_data);
    free(H2_1_density);
    
    free(H2_2_density);
    
    free(ge_density);
    
    free(He_1_density);
    
    free(H_m0_density);
    
    free(He_3_density);
    
    free(He_2_density);
    
    free(H_1_density);
    
    free(de_density);
    
    free(H_2_density);
    
    free(cooling_time);
    free(gamma);
    free(temperature);
    free(mean_molecular_weight);

}