#include "dengo_types.h"
#include "omp_test_cvdls_solver.h"

omp_test_cvdls_data *omp_test_cvdls_setup_data(
    int *NumberOfFields, char ***FieldNames);

int write_dengo_data_to_file( dengo_field_data *);    
int read_init_data_to_dengo(dengo_field_data * );
int omp_test_cvdls_solve_chemistry_dt( dengo_field_data *, 
                                       omp_test_cvdls_data*, 
                                       double);
double dengo_evolve_omp_test_cvdls_field_data (double dtf, double &dt, double z, dengo_field_data *field_data,
            double rtol, double atol, long long dims, omp_test_cvdls_data *data);


int main(int argc, char** argv){
    
    // Test the data object created in dengo
    //
    // create initial field data
    //
    
    dengo_field_data *field_data = (dengo_field_data *) malloc(sizeof(dengo_field_data));

    int flag = read_init_data_to_dengo( field_data);
   
    fprintf(stderr, "data are now read in \n");

    double idt = 1e10;

    
    omp_test_cvdls_data *rate_data = omp_test_cvdls_setup_data(NULL, NULL);
    
    // try to make sure it runs, and reads from the dengo_field_data
    fprintf(stderr,"rates too!");

    double tfinal;
    int ncells = field_data->ncells;
    double z = 0.0;
    double dt[1];
    double atol, rtol; 
    atol = 1e-4;
    dt[0] = idt;
    tfinal = dengo_evolve_omp_test_cvdls_field_data(idt, dt[0], z, field_data, rtol, atol, ncells,rate_data );
    fprintf(stderr, "tfinal now: %0.5g \n", tfinal);

    flag = write_dengo_data_to_file( field_data);
    
    // omp_test_cvdls_main(argc, argv);
}
