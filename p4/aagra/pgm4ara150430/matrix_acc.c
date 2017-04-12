/* Author - Ayush Agrawal
 File - matrix_acc.c*/
 
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#define N 1000

float A[N][N], B[N][N], C[N][N]; // declaring matrices of NxN size
int main ()
{

	int i, j, m; // indices for matrix multiplication
	struct timeval start_time, stop_time, elapsed_time;  // timers


	for(i=0;i<N;i++) {
		for(j=0;j<N;j++) {
			A[i][j]= 1;
			B[i][j]= 1;
		}
	}


	gettimeofday(&start_time,NULL); // Unix timer
	
	// Compute matrix multiplication.
	#pragma acc kernels copyin(A,B) copy(C)
	for (i = 0; i < N; ++i) {
		for (j = 0; j < N; ++j) {
			for (m = 0; m < N; ++m) {
				C[i][j] += A[i][m] * B[m][j];
			}
		}
	}
		
    gettimeofday(&stop_time,NULL);
    timersub(&stop_time, &start_time, &elapsed_time); // Unix time subtract routine
    printf("Total time was %f seconds.\n", elapsed_time.tv_sec+elapsed_time.tv_usec/1000000.0);
	
	/*for(i=0;i<N;i++) {
		for(j=0;j<N;j++) {
					printf("%f\n",C[i][j]);
				}
			}*/

	return 0;
}
