#include "dls_solver.h"

int main(int argc, char **argv) {

    ProfilerStart("prof.out");
    dls_main(argc, argv);
    ProfilerFlush();
    ProfilerStop();
}
