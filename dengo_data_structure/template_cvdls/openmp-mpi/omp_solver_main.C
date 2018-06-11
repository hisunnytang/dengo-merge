#include "dengo_types.h"
#include "omp_solver.h"

omp_data *omp_setup_data(
    int *NumberOfFields, char ***FieldNames);

int write_dengo_data_to_file( dengo_field_data *);    
int read_init_data_to_dengo(dengo_field_data * );
int omp_solve_chemistry_dt( dengo_field_data *, 
                                       omp_data*, 
                                       double);

int main(int argc, char** argv){
    
    
    omp_main(argc, argv);

    return 1;
}

/*
int _main(int argc, char** argv){
    
    // Test the data object created in dengo
    //
    // create initial field data
    //

    dengo_field_data *field_data = (dengo_field_data *) malloc(sizeof(dengo_field_data));

    int flag = read_init_data_to_dengo( field_data);
   
    fprintf(stderr, "data are now read in \n");

    double idt = 1e10;

    
    omp_data *rate_data = omp_setup_data(NULL, NULL);
    
    // try to make sure it runs, and reads from the dengo_field_data
    fprintf(stderr,"rates too!");
    for (int i = 0; i<5; i++){
    flag = omp_solve_chemistry_dt( field_data, rate_data, idt );
    }

    fprintf(stderr,"writing data!!!!!\n");
    flag = write_dengo_data_to_file( field_data);
}
*/
