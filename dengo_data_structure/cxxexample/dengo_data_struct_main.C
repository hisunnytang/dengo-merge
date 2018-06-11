#include "dengo_types.h"
#include "cvdls_9species_solver.h"
#include <time.h>

//cvdls_9species_data *cvdls_9species_setup_data(
//    int *NumberOfFields, char ***FieldNames);

int write_dengo_data_to_file( dengo_field_data *);    
int read_init_data_to_dengo(dengo_field_data * );
int cvdls_9species_solve_chemistry_dt( dengo_field_data *, cvdls_9species_data *, double);


    

dengo_field_data *field_data = (dengo_field_data *) malloc(sizeof(dengo_field_data));

cvdls_9species_data *rate_data = cvdls_9species_setup_data(NULL, NULL);

int flag = read_init_data_to_dengo( field_data);

int main(int argc, char** argv){
    
    // Test the data object created in dengo
    //
    // create initial field data
    //
    fprintf(stderr, "answer: %0.6g\n", field_data->H2I_density[0]);
   

    double idt = 1.0e10;
    
    

    struct timeval start, stop;
    double secs = 0;
    gettimeofday(&start, NULL);

    flag = cvdls_9species_solve_chemistry_dt( field_data, rate_data, idt ); 
    
    gettimeofday(&stop, NULL);
    secs = (double)(stop.tv_usec - start.tv_usec) / 1000000 + (double)(stop.tv_sec - start.tv_sec);
    printf("time taken %f\n",secs);

    //flag = write_dengo_data_to_file( field_data);
fprintf(stderr, "answer: %0.6g\n", field_data->H2I_density[0]);
    write_dengo_data_to_file( field_data );
    
    free(field_data->H2I_density);
    free(field_data->H2II_density);
    free(field_data->HI_density);
    free(field_data->HII_density);
    free(field_data->HM_density);
    free(field_data->HeI_density);
    free(field_data->HeII_density);
    free(field_data->HeIII_density);
    free(field_data->de_density);
    free(field_data->ge_density);

    free(field_data);
    free(rate_data);

}
