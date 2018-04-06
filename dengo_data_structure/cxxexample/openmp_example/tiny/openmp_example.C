#include "omp.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>

void Test2( double *n, int m )
{
  printf( "<T:%d> - %0.5g, %d\n", omp_get_thread_num(), n[0], m );
  *n = 100 + m;
}

int main()
{

    double happy[100];
    for (int i = 0; i<100; i++){
        happy[i] = i+0.1;
    }

#pragma omp parallel for private(i) 
for( int d = 0; d < 100; d ++ ){
    Test2( &happy[d], d );
    int N = 100;
    for (int i = 222; i < 229; i++){
        fprintf(stderr, "i: %d; d:%d; dN+i:%d\n", i, d, d*N+i);
        //fprintf( stderr, "happy[%d] = %0.5g \n", d*N + i, happy[d] );
    //fprintf(stderr, "max_thread: %d",omp_get_max_threads());
}
}
}
