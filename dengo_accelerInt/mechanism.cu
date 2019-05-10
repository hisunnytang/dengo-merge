#include <stdio.h>
#include "mechanism.cuh"
#include "gpu_memory.cuh"

/*
    //apply masking of ICs for cache optimized mechanisms
    void apply_mask(double* y_specs) {
        double temp [NSP];
        memcpy(temp, y_specs, NSP * sizeof(double));
        y_specs[0] = temp[0];
        y_specs[1] = temp[1];
        y_specs[2] = temp[2];
        y_specs[3] = temp[3];
        y_specs[4] = temp[4];
        y_specs[5] = temp[5];
        y_specs[6] = temp[6];
        y_specs[7] = temp[7];
        y_specs[8] = temp[9];
        y_specs[9] = temp[10];
        y_specs[10] = temp[11];
        y_specs[11] = temp[12];
        y_specs[12] = temp[8];
    }
    //reverse masking of ICs for cache optimized mechanisms
    void apply_reverse_mask(double* y_specs) {
        double temp [NSP];
        memcpy(temp, y_specs, NSP * sizeof(double));
        y_specs[0] = temp[0];
        y_specs[1] = temp[1];
        y_specs[2] = temp[2];
        y_specs[3] = temp[3];
        y_specs[4] = temp[4];
        y_specs[5] = temp[5];
        y_specs[6] = temp[6];
        y_specs[7] = temp[7];
        y_specs[8] = temp[12];
        y_specs[9] = temp[8];
        y_specs[10] = temp[9];
        y_specs[11] = temp[10];
        y_specs[12] = temp[11];
    }
*/

void set_same_initial_conditions(int NUM, double** y_host, double** var_host) 
{

/*
{'H2_1': array([10000.]),
 'H2_2': array([10000.]),
 'H_1': array([7.6e+13]),
 'H_2': array([10000.]),
 'H_m0': array([10000.]),
 'He_1': array([2.4e+13]),
 'He_2': array([10000.]),
 'He_3': array([10000.]),
 'T': array([3000.]),
 'de': array([12495.12442156]),
 'density': array([1.e+14]),
 'ge': array([3.02681398e+11])}
*/

    // set intial temperature, units [K]
    double T0 = 2000.0;
    double density = 1.0e13;
    double mH = 1.67e-24;
    double k  = 1.3806488e-16;
    double tiny = 1.0e-20;


    (*y_host) = (double*)malloc(NUM * NSP * sizeof(double));
    (*var_host) = (double*)malloc(NUM * sizeof(double));
    //load temperature and mass fractions for all threads (cells)
    printf("NUM = %d; NSP = %d \n", NUM, NSP ); 
    int j = 1;
    for (int i = 0; i < NUM; ++i) {
        (*y_host)[i] = T0;
        //loop through species
	j = 0;
    	// H2I
	(*y_host)[i + NUM * j] = density * tiny / 2.0;
        j += 1;
        // H2II
	(*y_host)[i + NUM * j] = density * tiny / 2.0;
	j += 1;
        // HI
	(*y_host)[i + NUM * j] = 0.76 * density / 1.00794;
	j += 1;
        // HII
	(*y_host)[i + NUM * j] = density * tiny / 1.00794;
	j += 1;
        // H-
	(*y_host)[i + NUM * j] = density * tiny / 1.00794;
	j += 1;
        // HeI
	(*y_host)[i + NUM * j] = 0.24 * density / 4.002602;
	j += 1;
        // HeII
	(*y_host)[i + NUM * j] = tiny * density / 4.002602;
	j += 1;
        // HeIII
	(*y_host)[i + NUM * j] = tiny * density / 4.002602;
	j += 1;
        // electron (de)
	(*y_host)[i + NUM * j] = tiny * density;
	j += 1;
        // internal energy (ge)
	(*y_host)[i + NUM * j] = 3.0 / 2.0 * k * T0 / mH;
	j += 1;
    }

}

