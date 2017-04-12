/******************************************************************************
* FILE: mm_acc.c  
*   OpenACC matrix multiplication
* Name: David
* Program 4
******************************************************************************/
//#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

//#define N 10        /* number of rows in matrix A */
					/* number of columns in matrix A */
					/* number of columns in matrix B */
#define N 1000

int main (int argc, char *argv[]) 
{
int	tid, nthreads, i, j, k, chunk;
float	a[N][N],           /* matrix A to be multiplied */
		b[N][N],           /* matrix B to be multiplied */
		c[N][N];           /* result matrix C */

double start,stop;		
struct timeval start_time, stop_time, elapsed_time;
		
chunk = 10;                    /* set loop iteration chunk size */

/*** Spawn a parallel region explicitly scoping all variables ***/
   #pragma acc data copyin (a[0: N * N ], b[0: N * N]) copyout (a[0: N * N ], b[0: N * N],c[0: N * N ])
  {
  
  /*** Initialize matrices ***/
  #pragma acc kernels
  for (i=0; i<N; i++)
    for (j=0; j<N; j++)
      a[i][j]= j-i;
  #pragma acc kernels
  for (i=0; i<N; i++)
    for (j=0; j<N; j++)
      b[i][j]= N-j+i;
  #pragma acc kernels
  for (i=0; i<N; i++)
    for (j=0; j<N; j++)
      c[i][j]= 0;

  /*** Do matrix multiply sharing iterations on outer loop ***/
  /*** Display who does which iterations for demonstration purposes ***/
  
  //start = omp_get_wtime();
  gettimeofday(&start_time,NULL);
  #pragma acc region
  {
  #pragma acc loop independent vector(16)
  for (i=0; i<N; i++)    
  {
	#pragma acc loop independent vector(16)
    for(j=0; j<N; j++)       
      for (k=0; k<N; k++)
        c[i][j] += a[i][k] * b[k][j];
  }
  }
  }   /*** End of acc parallel region ***/

  gettimeofday(&stop_time,NULL);
  timersub(&stop_time, &start_time, &elapsed_time); // Unix time subtract routine

  
/*** Print results ***/
printf("******************************************************\n");
printf("Total time was %f seconds.\n", elapsed_time.tv_sec+elapsed_time.tv_usec/1000000.0);
//printf("Total time was %f seconds.\n", stop-start);
//double t = stop-start;
double t = elapsed_time.tv_sec+elapsed_time.tv_usec/1000000.0;

double flop = 2.0*N*N*N/(t*1000000000.0);
printf("The program's GFlops:  %f\n",flop);
	/*
for (i=0; i<N; i++)
  {
  for (j=0; j<N; j++) 
    printf("%6.2f   ", a[i][j]);
  printf("\n"); 
  }

printf("\n");
printf("\n");
  
for (i=0; i<N; i++)
  {
  for (j=0; j<N; j++) 
    printf("%6.2f   ", b[i][j]);
  printf("\n"); 
  }
	
printf("Result Matrix:\n");
for (i=0; i<N; i++)
  {
  for (j=0; j<N; j++) 
    printf("%6.2f   ", c[i][j]);
  printf("\n"); 
  }*/
printf("******************************************************\n");
printf ("Done.\n");

}