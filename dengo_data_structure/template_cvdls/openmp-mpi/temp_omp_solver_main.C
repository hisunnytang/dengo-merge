#include "dengo_types.h"
#include "temp_omp_solver.h"

temp_omp_data *temp_omp_setup_data(
    int *NumberOfFields, char ***FieldNames);

int write_dengo_data_to_file( dengo_field_data *);    
int read_init_data_to_dengo(dengo_field_data * );
int temp_omp_solve_chemistry_dt( dengo_field_data *, 
                                       temp_omp_data*, 
                                       double);
double dengo_evolve_temp_omp_field_data (double dtf, double &dt, double z, dengo_field_data *field_data,
            double *rtol, double *atol, long long dims, temp_omp_data *data);

int main(int argc, char** argv){
    
    // Test the data object created in dengo
    //
    // create initial field data
    //
     
    dengo_field_data *field_data = (dengo_field_data *) malloc(sizeof(dengo_field_data));

    int flag = read_init_data_to_dengo( field_data);
   
    fprintf(stderr, "data are now read in \n");

    double idt = 1e10;

    
    temp_omp_data *rate_data = temp_omp_setup_data(NULL, NULL);
    
    // try to make sure it runs, and reads from the dengo_field_data
    fprintf(stderr,"rates too!");

    double tfinal;
    int ncells = field_data->ncells;
    double z = 0.0;
    double dt[1];
    double *atol, *rtol; 
    dt[0] = idt;

    tfinal = dengo_evolve_temp_omp_field_data(idt, dt[0], z, field_data, rtol, atol, ncells,rate_data );
    
    free(rate_data);
    flag = write_dengo_data_to_file( field_data);
    

    // temp_omp_main(argc, argv);
}
