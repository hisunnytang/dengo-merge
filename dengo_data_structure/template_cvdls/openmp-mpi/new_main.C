#include "omp.h"
/* stdlib, hdf5, local includes */
#include "dengo_types.h"
#include "stdio.h"
#include "stdlib.h"
#include "_temp_omp_solver.h"


_temp_omp_data *_temp_omp_setup_data(
    int *NumberOfFields, char ***FieldNames);

double dengo_evolve__temp_omp_field_data (double dtf, double &dt, double z, dengo_field_data *field_data,
            double rtol, double atol, long long dims, _temp_omp_data *data);

int write_dengo_data_to_file( dengo_field_data *);    
int read_init_data_to_dengo(dengo_field_data * );



int main(){
    
    dengo_field_data *field_data = (dengo_field_data *)malloc(sizeof(dengo_field_data));
    
    _temp_omp_data *rate_data = _temp_omp_setup_data(NULL,NULL);

     int flag = read_init_data_to_dengo( field_data);
   
    fprintf(stderr, "data are now read in \n");

    double idt = 1e10;

    double tfinal;
    int ncells = field_data->ncells;
    double z = 0.0;
    double dt[1];
    double atol, rtol; 
    atol = 1e-6;
    dt[0] = idt;
    tfinal = dengo_evolve__temp_omp_field_data(idt, dt[0], z, field_data, rtol, atol, ncells,rate_data );
    fprintf(stderr, "tfinal now: %0.5g \n", tfinal);

    fprintf(stderr,"writing data!!!!!\n");
    flag = write_dengo_data_to_file( field_data);

}
