#include "test_cvdls_solver.h"

int main(int argc, char **argv) {
    double t_final;

    char filename[50];
    strcpy(filename, argv[1] );
    
    int status;
    status = test_cvdls_main(filename);
    if ( (status) < 0.9999999 ){
        return 1;
    }
    return 0;
}   
