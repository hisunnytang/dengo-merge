#include "testing_solver.h"

int main(int argc, char **argv) {
    clock_t start = clock();
    testing_main(argc, argv);

    clock_t stop = clock();
    double elapsed = ( (double)(stop - start) ) / CLOCKS_PER_SEC;
    // fprintf(stderr ,"at density: %0.5g and temp: %0.5g \n", density, temperature);
    fprintf(stderr, "elapsed time: %0.9g\n", elapsed);
}
