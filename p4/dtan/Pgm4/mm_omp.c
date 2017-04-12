/******************************************************************************
* FILE: mm_omp.c  
*   OpenMp matrix multiplication
* Name: David
* Program 4
******************************************************************************/
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

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
		
chunk = 10;                    /* set loop iteration chunk size */

/*** Spawn a parallel region explicitly scoping all variables ***/
#pragma omp parallel shared(a,b,c,nthreads,chunk) private(tid,i,j,k)
  {
  tid = omp_get_thread_num();
  if (tid == 0)
    {
    nthreads = omp_get_num_threads();
    printf("Starting matrix multiple example with %d threads\n",nthreads);
    printf("Initializing matrices...\n");
    }
  /*** Initialize matrices ***/
  #pragma omp for schedule (static, chunk) 
  for (i=0; i<N; i++)
    for (j=0; j<N; j++)
      a[i][j]= j-i;
  #pragma omp for schedule (static, chunk)
  for (i=0; i<N; i++)
    for (j=0; j<N; j++)
      b[i][j]= N-j+i;
  #pragma omp for schedule (static, chunk)
  for (i=0; i<N; i++)
    for (j=0; j<N; j++)
      c[i][j]= 0;

  /*** Do matrix multiply sharing iterations on outer loop ***/
  /*** Display who does which iterations for demonstration purposes ***/
  printf("Thread %d starting matrix multiply...\n",tid);
  
  start = omp_get_wtime();
  #pragma omp for schedule (static, chunk)
  for (i=0; i<N; i++)    
    {
    printf("Thread=%d did row=%d\n",tid,i);
    for(j=0; j<N; j++)       
      for (k=0; k<N; k++)
        c[i][j] += a[i][k] * b[k][j];
    }
  }   /*** End of parallel region ***/

  stop = omp_get_wtime();
  
/*** Print results ***/
printf("******************************************************\n");

printf("Total time was %f seconds.\n", stop-start);
double t = stop-start;
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