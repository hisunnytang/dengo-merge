#include "dengo_types.h"
#include "cvdls_9species_solver.h"

cvdls_9species_data *cvdls_9species_setup_data(
    int *NumberOfFields, char ***FieldNames);

int write_dengo_data_to_file( dengo_field_data *);    
int read_init_data_to_dengo(dengo_field_data * );
int cvdls_9species_solve_chemistry_dt( dengo_field_data *, 
                                       cvdls_9species_data*, 
                                       double);


int main(int argc, char** argv){
    
    // Test the data object created in dengo
    //
    // create initial field data
    //

    dengo_field_data *field_data = (dengo_field_data *) malloc(sizeof(dengo_field_data));

    int flag = read_init_data_to_dengo( field_data);
   
    fprintf(stderr, "data are now read in \n");

    double idt = 1e10;

    
    cvdls_9species_data *rate_data = cvdls_9species_setup_data(NULL, NULL);
    
    fprintf(stderr,"rates too!");
    for (int i = 0; i<5; i++){
    flag = cvdls_9species_solve_chemistry_dt( field_data, rate_data, idt );
    }

    fprintf(stderr,"writing data!!!!!\n");
    flag = write_dengo_data_to_file( field_data);
}
