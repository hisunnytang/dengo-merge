#include "primordial_solver.h"

int solve_chemistry(int argc, char **argv);
int solve_chemistry_enzo(int argc, char **argv);

int main(int argc, char **argv) {
    if (argc > 1){
       primordial_main(argc, argv);
       return 0;
    }
    // solve_chemistry(argc, argv);
    // solve_chemistry_enzo(argc, argv);
}

/*
// Sample on how to use dengo with primordial chemistry
// with primordial_solve_chemistry
int solve_chemistry(int argc, char **argv) {
    dengo_field_data *field_data = (dengo_field_data *) malloc(sizeof(dengo_field_data));
    int N = 64*64;
    field_data->ncells = N; 
    double density = 1.0e8; // in cm^-3
    double T = 1000.0; // in K
    double mH, k, tiny;

    mH = 1.67e-24;
    k  = 1.380e-16;
    tiny = 1.0e-20; 
    density *= mH;
    double G = 6.67259e-8;
    code_units *units = (code_units *) malloc(sizeof(code_units));
    units->density_units = 1.0;
    units->length_units = 1.0;
    units->time_units = 1.0;
    units->velocity_units = 1.0;
    double *H2_1_density = (double*) malloc(N * sizeof(double));
    double *H2_2_density = (double*) malloc(N * sizeof(double));
    double *H_1_density = (double*) malloc(N * sizeof(double));
    double *H_2_density = (double*) malloc(N * sizeof(double));
    double *H_m0_density = (double*) malloc(N * sizeof(double));
    double *He_1_density = (double*) malloc(N * sizeof(double));
    double *He_2_density = (double*) malloc(N * sizeof(double));
    double *He_3_density = (double*) malloc(N * sizeof(double));
    double *de_density = (double*) malloc(N * sizeof(double));
    double *ge_density = (double*) malloc(N * sizeof(double));
    double *cooling_time = (double *) malloc( N * sizeof(double) );
    double *gamma = (double * ) malloc( N * sizeof(double) );
    double *temperature = (double *) malloc( N * sizeof(double) );
    double *mean_molecular_weight = (double *) malloc( N * sizeof(double) );

    double *reltol = (double*) malloc(sizeof(double));
    double *abstol = (double*) malloc(sizeof(double)* N * 10);
    reltol[0] = 1.0e-5;

    for ( int i = 0; i < N; i++){
        He_2_density[i] = tiny*density;
        H_m0_density[i] = tiny*density;
        H2_2_density[i] = tiny*density;
        He_3_density[i] = tiny*density;
        H2_1_density[i] = 1.0e-5 * density;
        H_1_density[i]  = 0.76   * density;
        He_1_density[i] = 0.24   * density;
	de_density[i]    = 1.0e-5 * density;
	H_2_density[i]   = 1.0e-5 * density;

	// ge ~ nkT / (gamma - 1)/ rho; gamaa ~ 5/3
        ge_density[i]   = 3.0/2.0 * k * T / mH;
	//density[i] = (0.76+0.24+1.0e-5+1.0e-5)*density;
    }

    for ( int i = 0; i < N * 10; i++ ){
    	abstol[i] = tiny * reltol[0];
    }
    field_data->H_1_density = H_1_density;
    field_data->He_2_density = He_2_density;
    field_data->He_1_density = He_1_density;
    field_data->ge_density = ge_density;
    field_data->H_m0_density = H_m0_density;
    field_data->H2_2_density = H2_2_density;
    field_data->de_density = de_density;
    field_data->He_3_density = He_3_density;
    field_data->H_2_density = H_2_density;
    field_data->H2_1_density = H2_1_density;
    //field_data->density         = density;
    field_data->CoolingTime     = cooling_time;
    field_data->Gamma           = gamma;
    field_data->temperature     = temperature;
    field_data->MolecularWeight = mean_molecular_weight;

    const char *fileloc = "/home/kwoksun2/dengo_install/primordial_tables.h5";
    field_data->dengo_data_file = fileloc;
    field_data->reltol = reltol[0];

    double dt = 1.0 / sqrt(G * density) ;
    fprintf(stderr, "MAX_NCELLS = %d \n", MAX_NCELLS);
    primordial_solve_chemistry( units, field_data, dt );
    dengo_estimate_cooling_time( units, field_data);
    fprintf(stderr, "H2_1 = %0.5g\n", field_data->H2_1_density[0] / mH );
    fprintf(stderr, "H2_2 = %0.5g\n", field_data->H2_2_density[0] / mH );
    fprintf(stderr, "H_1 = %0.5g\n", field_data->H_1_density[0] / mH );
    fprintf(stderr, "H_2 = %0.5g\n", field_data->H_2_density[0] / mH );
    fprintf(stderr, "H_m0 = %0.5g\n", field_data->H_m0_density[0] / mH );
    fprintf(stderr, "He_1 = %0.5g\n", field_data->He_1_density[0] / mH );
    fprintf(stderr, "He_2 = %0.5g\n", field_data->He_2_density[0] / mH );
    fprintf(stderr, "He_3 = %0.5g\n", field_data->He_3_density[0] / mH );
    fprintf(stderr, "de = %0.5g\n", field_data->de_density[0] / mH );
    fprintf(stderr, "ge = %0.5g\n", field_data->ge_density[0] );
    fprintf(stderr, "CoolingTime = %0.5g\n", field_data->CoolingTime[0]);
    
    // at the low density limit,
    // Hm is important catalyst for H2 formation
    // we compare the relative difference of of the solution
    // between different cells
    // since we starts with the same initial conditions,
    // it should give the same result, at the worst, 
    // results within the relative tolerance level
    
    unsigned long d;
    // lets just compare everything!!!!!!
    double ref0, frac;
    double H_1;
    ref0 = field_data->H_1_density[0];
    for (d = 1; d < N; d++){
	H_1 = field_data->H_1_density[d];
	frac = fabs(H_1-ref0)/ref0;
	if (frac > reltol[0]){
	    fprintf(stderr, "H_1[%lu] = %0.5g; diff = %0.5g\n", d, H_1, frac);
	}
    }
    double He_2;
    ref0 = field_data->He_2_density[0];
    for (d = 1; d < N; d++){
	He_2 = field_data->He_2_density[d];
	frac = fabs(He_2-ref0)/ref0;
	if (frac > reltol[0]){
	    fprintf(stderr, "He_2[%lu] = %0.5g; diff = %0.5g\n", d, He_2, frac);
	}
    }
    double He_1;
    ref0 = field_data->He_1_density[0];
    for (d = 1; d < N; d++){
	He_1 = field_data->He_1_density[d];
	frac = fabs(He_1-ref0)/ref0;
	if (frac > reltol[0]){
	    fprintf(stderr, "He_1[%lu] = %0.5g; diff = %0.5g\n", d, He_1, frac);
	}
    }
    double ge;
    ref0 = field_data->ge_density[0];
    for (d = 1; d < N; d++){
	ge = field_data->ge_density[d];
	frac = fabs(ge-ref0)/ref0;
	if (frac > reltol[0]){
	    fprintf(stderr, "ge[%lu] = %0.5g; diff = %0.5g\n", d, ge, frac);
	}
    }
    double H_m0;
    ref0 = field_data->H_m0_density[0];
    for (d = 1; d < N; d++){
	H_m0 = field_data->H_m0_density[d];
	frac = fabs(H_m0-ref0)/ref0;
	if (frac > reltol[0]){
	    fprintf(stderr, "H_m0[%lu] = %0.5g; diff = %0.5g\n", d, H_m0, frac);
	}
    }
    double H2_2;
    ref0 = field_data->H2_2_density[0];
    for (d = 1; d < N; d++){
	H2_2 = field_data->H2_2_density[d];
	frac = fabs(H2_2-ref0)/ref0;
	if (frac > reltol[0]){
	    fprintf(stderr, "H2_2[%lu] = %0.5g; diff = %0.5g\n", d, H2_2, frac);
	}
    }
    double de;
    ref0 = field_data->de_density[0];
    for (d = 1; d < N; d++){
	de = field_data->de_density[d];
	frac = fabs(de-ref0)/ref0;
	if (frac > reltol[0]){
	    fprintf(stderr, "de[%lu] = %0.5g; diff = %0.5g\n", d, de, frac);
	}
    }
    double He_3;
    ref0 = field_data->He_3_density[0];
    for (d = 1; d < N; d++){
	He_3 = field_data->He_3_density[d];
	frac = fabs(He_3-ref0)/ref0;
	if (frac > reltol[0]){
	    fprintf(stderr, "He_3[%lu] = %0.5g; diff = %0.5g\n", d, He_3, frac);
	}
    }
    double H_2;
    ref0 = field_data->H_2_density[0];
    for (d = 1; d < N; d++){
	H_2 = field_data->H_2_density[d];
	frac = fabs(H_2-ref0)/ref0;
	if (frac > reltol[0]){
	    fprintf(stderr, "H_2[%lu] = %0.5g; diff = %0.5g\n", d, H_2, frac);
	}
    }
    double H2_1;
    ref0 = field_data->H2_1_density[0];
    for (d = 1; d < N; d++){
	H2_1 = field_data->H2_1_density[d];
	frac = fabs(H2_1-ref0)/ref0;
	if (frac > reltol[0]){
	    fprintf(stderr, "H2_1[%lu] = %0.5g; diff = %0.5g\n", d, H2_1, frac);
	}
    }

    double ct;
    ref0 = field_data->CoolingTime[0];
    for (d = 1; d < N; d++){
	ct = field_data->CoolingTime[d];
	frac = fabs(ct-ref0)/ref0;
	if (frac > reltol[0]){
	    fprintf(stderr, "CoolingTime[%lu] = %0.5g; diff = %0.5g\n",d,ct, frac);
	}
    }

    free(field_data);
    free(H2_1_density);
    free(H2_2_density);
    free(H_1_density);
    free(H_2_density);
    free(H_m0_density);
    free(He_1_density);
    free(He_2_density);
    free(He_3_density);
    free(de_density);
    free(ge_density);
    free(cooling_time);
    free(gamma);
    free(temperature);
    free(mean_molecular_weight);
    free(abstol);
    free(reltol);
}


// Sample on how to use dengo with primordial chemistry
// with primordial_solve_chemistry_enzo
int solve_chemistry_enzo(int argc, char **argv) {
    dengo_field_data *field_data = (dengo_field_data *) malloc(sizeof(dengo_field_data));
    int N = 16*16*16;
    field_data->ncells = N; 
    double density = 1.0e8; // in cm^-3
    double T = 1000.0; // in K
    double mH, k, tiny;

    mH = 1.67e-24;
    k  = 1.380e-16;
    tiny = 1.0e-20; 
    density *= mH;
    double G = 6.67259e-8;
    code_units *units = (code_units *) malloc(sizeof(code_units));
    units->density_units = 1.0;
    units->length_units = 1.0;
    units->time_units = 1.0;
    units->velocity_units = 1.0;
    double *H2_1_density = (double*) malloc(N * sizeof(double));
    double *H2_2_density = (double*) malloc(N * sizeof(double));
    double *H_1_density = (double*) malloc(N * sizeof(double));
    double *H_2_density = (double*) malloc(N * sizeof(double));
    double *H_m0_density = (double*) malloc(N * sizeof(double));
    double *He_1_density = (double*) malloc(N * sizeof(double));
    double *He_2_density = (double*) malloc(N * sizeof(double));
    double *He_3_density = (double*) malloc(N * sizeof(double));
    double *de_density = (double*) malloc(N * sizeof(double));
    double *ge_density = (double*) malloc(N * sizeof(double));
    double *cooling_time = (double *) malloc( N * sizeof(double) );
    double *gamma = (double * ) malloc( N * sizeof(double) );
    double *temperature = (double *) malloc( N * sizeof(double) );
    double *mean_molecular_weight = (double *) malloc( N * sizeof(double) );
    double *density_arr = (double *) malloc( N * sizeof(double) );

    double *reltol = (double*) malloc(sizeof(double));
    double *abstol = (double*) malloc(sizeof(double)* N * 10);
    reltol[0] = 1.0e-5;

    for ( int i = 0; i < N; i++){
        He_2_density[i] = tiny*density;
        H_m0_density[i] = tiny*density;
        H2_2_density[i] = tiny*density;
        He_3_density[i] = tiny*density;
        H2_1_density[i] = 1.0e-5 * density;
        H_1_density[i]  = 0.76   * density;
        He_1_density[i] = 0.24   * density;
	de_density[i]    = 1.0e-5 * density;
	H_2_density[i]   = 1.0e-5 * density;

	// ge ~ nkT / (gamma - 1)/ rho; gamaa ~ 5/3
        ge_density[i]   = 3.0/2.0 * k * T / mH;
	density_arr[i] = (1+2.0*1e-5)*density;
    }

    for ( int i = 0; i < N * 10; i++ ){
    	abstol[i] = tiny * reltol[0];
    }
    field_data->H_1_density = H_1_density;
    field_data->He_2_density = He_2_density;
    field_data->He_1_density = He_1_density;
    field_data->ge_density = ge_density;
    field_data->H_m0_density = H_m0_density;
    field_data->H2_2_density = H2_2_density;
    field_data->de_density = de_density;
    field_data->He_3_density = He_3_density;
    field_data->H_2_density = H_2_density;
    field_data->H2_1_density = H2_1_density;

    field_data->density         = density_arr;
    field_data->CoolingTime     = cooling_time;
    field_data->Gamma           = gamma;
    field_data->temperature     = temperature;
    field_data->MolecularWeight = mean_molecular_weight;

    int gstart[3];
    int gend[3];
    int gd[3];

    gstart[0] = 0;
    gstart[1] = 0;
    gstart[2] = 0;

    gend[0] = 15;
    gend[1] = 15;
    gend[2] = 15;

    gd[0] = 16;
    gd[1] = 16;
    gd[2] = 16;

    field_data->grid_start = &gstart[0];
    field_data->grid_end   = &gend[0];
    field_data->grid_dimension = &gd[0];


    const char *fileloc = "/home/kwoksun2/dengo_install/primordial_tables.h5";
    field_data->dengo_data_file = fileloc;
    field_data->reltol = reltol[0];

    units->a_value = 1.0;
    units->a_units = 1.0;

    double dt = 1.0 / sqrt(G * density) ;
    fprintf(stderr, "MAX_NCELLS = %d \n", MAX_NCELLS);
    primordial_solve_chemistry_enzo( units, field_data, dt );
    dengo_estimate_cooling_time_enzo( units, field_data);
    fprintf(stderr, "H2_1 = %0.5g\n", field_data->H2_1_density[0] / mH );
    fprintf(stderr, "H2_2 = %0.5g\n", field_data->H2_2_density[0] / mH );
    fprintf(stderr, "H_1 = %0.5g\n", field_data->H_1_density[0] / mH );
    fprintf(stderr, "H_2 = %0.5g\n", field_data->H_2_density[0] / mH );
    fprintf(stderr, "H_m0 = %0.5g\n", field_data->H_m0_density[0] / mH );
    fprintf(stderr, "He_1 = %0.5g\n", field_data->He_1_density[0] / mH );
    fprintf(stderr, "He_2 = %0.5g\n", field_data->He_2_density[0] / mH );
    fprintf(stderr, "He_3 = %0.5g\n", field_data->He_3_density[0] / mH );
    fprintf(stderr, "de = %0.5g\n", field_data->de_density[0] / mH );
    fprintf(stderr, "ge = %0.5g\n", field_data->ge_density[0] );
    fprintf(stderr, "CoolingTime = %0.5g\n", field_data->CoolingTime[0]);

    unsigned long d;
    // lets just compare everything!!!!!!
    double ref0, frac;
    double H_1;
    ref0 = field_data->H_1_density[0];
    for (d = 1; d < N; d++){
	H_1 = field_data->H_1_density[d];
	frac = fabs(H_1-ref0)/ref0;
	if (frac > reltol[0]){
	    fprintf(stderr, "H_1[%lu] = %0.5g; diff = %0.5g\n", d, H_1, frac);
	}
    }
    double He_2;
    ref0 = field_data->He_2_density[0];
    for (d = 1; d < N; d++){
	He_2 = field_data->He_2_density[d];
	frac = fabs(He_2-ref0)/ref0;
	if (frac > reltol[0]){
	    fprintf(stderr, "He_2[%lu] = %0.5g; diff = %0.5g\n", d, He_2, frac);
	}
    }
    double He_1;
    ref0 = field_data->He_1_density[0];
    for (d = 1; d < N; d++){
	He_1 = field_data->He_1_density[d];
	frac = fabs(He_1-ref0)/ref0;
	if (frac > reltol[0]){
	    fprintf(stderr, "He_1[%lu] = %0.5g; diff = %0.5g\n", d, He_1, frac);
	}
    }
    double ge;
    ref0 = field_data->ge_density[0];
    for (d = 1; d < N; d++){
	ge = field_data->ge_density[d];
	frac = fabs(ge-ref0)/ref0;
	if (frac > reltol[0]){
	    fprintf(stderr, "ge[%lu] = %0.5g; diff = %0.5g\n", d, ge, frac);
	}
    }
    double H_m0;
    ref0 = field_data->H_m0_density[0];
    for (d = 1; d < N; d++){
	H_m0 = field_data->H_m0_density[d];
	frac = fabs(H_m0-ref0)/ref0;
	if (frac > reltol[0]){
	    fprintf(stderr, "H_m0[%lu] = %0.5g; diff = %0.5g\n", d, H_m0, frac);
	}
    }
    double H2_2;
    ref0 = field_data->H2_2_density[0];
    for (d = 1; d < N; d++){
	H2_2 = field_data->H2_2_density[d];
	frac = fabs(H2_2-ref0)/ref0;
	if (frac > reltol[0]){
	    fprintf(stderr, "H2_2[%lu] = %0.5g; diff = %0.5g\n", d, H2_2, frac);
	}
    }
    double de;
    ref0 = field_data->de_density[0];
    for (d = 1; d < N; d++){
	de = field_data->de_density[d];
	frac = fabs(de-ref0)/ref0;
	if (frac > reltol[0]){
	    fprintf(stderr, "de[%lu] = %0.5g; diff = %0.5g\n", d, de, frac);
	}
    }
    double He_3;
    ref0 = field_data->He_3_density[0];
    for (d = 1; d < N; d++){
	He_3 = field_data->He_3_density[d];
	frac = fabs(He_3-ref0)/ref0;
	if (frac > reltol[0]){
	    fprintf(stderr, "He_3[%lu] = %0.5g; diff = %0.5g\n", d, He_3, frac);
	}
    }
    double H_2;
    ref0 = field_data->H_2_density[0];
    for (d = 1; d < N; d++){
	H_2 = field_data->H_2_density[d];
	frac = fabs(H_2-ref0)/ref0;
	if (frac > reltol[0]){
	    fprintf(stderr, "H_2[%lu] = %0.5g; diff = %0.5g\n", d, H_2, frac);
	}
    }
    double H2_1;
    ref0 = field_data->H2_1_density[0];
    for (d = 1; d < N; d++){
	H2_1 = field_data->H2_1_density[d];
	frac = fabs(H2_1-ref0)/ref0;
	if (frac > reltol[0]){
	    fprintf(stderr, "H2_1[%lu] = %0.5g; diff = %0.5g\n", d, H2_1, frac);
	}
    }

    double ct;
    ref0 = field_data->CoolingTime[0];
    for (d = 1; d < N; d++){
	ct = field_data->CoolingTime[d];
	frac = fabs(ct-ref0)/ref0;
	if (frac > reltol[0]){
	    fprintf(stderr, "CoolingTime[%lu] = %0.5g; diff = %0.5g\n",d,ct, frac);
	}
    }

    free(field_data);
    free(H2_1_density);
    free(H2_2_density);
    free(H_1_density);
    free(H_2_density);
    free(H_m0_density);
    free(He_1_density);
    free(He_2_density);
    free(He_3_density);
    free(de_density);
    free(ge_density);
    free(cooling_time);
    free(gamma);
    free(temperature);
    free(mean_molecular_weight);
    free(abstol);
    free(reltol);
}
*/