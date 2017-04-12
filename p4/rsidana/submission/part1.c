// File Name :    OpenMp.c
// Course:        CS 6376 Spring 2017
// Author:        Rachna Sidana
// Assignment:    Program 4 - Part 1
// Description:   This program applies openMp directives to the Matrix Multiplication Algorithm and 
//                calculates the execution time of code that has been parallelized.

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

void InitFirstMatrix(int** matA, int n);
void InitSecondMatrix(int** matB, int n);
void Print_matrix(int** mat, int n);
void mult(int** matA, int** matB, int** matC, int n);
int verifyMult(int** matA, int** matB, int **matC, int n);
int main()
{
  int m, n, p, q, c, d, k, sum = 0, verification=1;
  struct timeval start_time, stop_time, elapsed_time;  // timers
  double  numFlops;
  float gflops;
  printf("Enter the value of 'N' ?\n");
  scanf("%d",&n);
  int *matA[n], *matB[n], *matC[n];
  for (int i=0; i<n; i++){
         matA[i] = (int *)malloc(n * sizeof(int));
         matB[i] = (int *)malloc(n * sizeof(int));
         matC[i] = (int *)malloc(n * sizeof(int));
       }
  InitFirstMatrix(matA,n);
  InitSecondMatrix(matB,n);
  gettimeofday(&start_time,NULL);
  mult(matA, matB, matC, n);
  gettimeofday(&stop_time,NULL);
  timersub(&stop_time, &start_time, &elapsed_time); // Unix time subtract routine
  printf("Total time was %f seconds.\n", elapsed_time.tv_sec+elapsed_time.tv_usec/1000000.0);
  numFlops = (2.0f*n*n*n-n*n)/1000000000.0f;
  gflops = numFlops/(elapsed_time.tv_sec+elapsed_time.tv_usec/1000000.0);
//printf("Flops :  %f .\n",numFlops); 
 printf("GFlops :  %f .\n",gflops);	
  verification =  verifyMult(matA,matB,matC,n);
/*printf("The input matrix A is :\n");
for (c = 0; c < n; c++) {
      for (d = 0; d < n; d++)
        printf("%d\t", matA[c][d]);

      printf("\n");
    }

printf("\nThe input matrix B is :\n");
for (c = 0; c < n; c++) {
      for (d = 0; d < n; d++)
        printf("%d\t", matB[c][d]);

      printf("\n");
    }
printf("\nThe output matrix C is :\n");
for (c = 0; c < n; c++) {
      for (d = 0; d < n; d++)
        printf("%d\t", matC[c][d]);
 
      printf("\n");
    } */


  if(verification==1)
printf("\nThe 'verifyMult() function' returned: 1\n");
    printf("\nResult is correct !\n");
  return 0;
}


/*-------------------------------------------------------------------
 * Function:  Init_matrix
 * Purpose:   initialize first matrix
 * In arg:    n
 * Out arg:   matA
 --------------------------------------------------------------------*/
void InitFirstMatrix(int** matA, int n) {
   int i, j;
   for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++)
         matA[i][j]=j-i;
   }
 }
/*-------------------------------------------------------------------
 * Function:  Init_matrix
 * Purpose:   initialize second matrix
 * In arg:    n
 * Out arg:   matB
 -----------------------------------------------------------------------*/
void InitSecondMatrix(int** matB, int n) {
   int i, j;
   for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++)
         matB[i][j]=n-j+i;
   }
 }
 /*-------------------------------------------------------------------
 * Function:  Mult
 * Purpose:   multiply first and second matrix
 * In arg:    matA, matB, n
 * Out arg:   matC
 --------------------------------------------------------------------*/

void mult(int** matA, int** matB, int** matC, int n){
  int c,d,sum=0,k;
  //the innermost loop is not parallelizable
#pragma omp parallel for private(d,k) reduction(+:sum)
  for (c = 0; c < n; c++) {
    for (d = 0; d < n; d++) {
        for (k = 0; k < n; k++) {
          sum = sum + matA[c][k]*matB[k][d];
        }
 
        matC[c][d] = sum;
        sum = 0;
      }
    }
}  
/*-------------------------------------------------------------------
 * Function:  verifyMult
 * Purpose:   verify the multiplication of first and second matrix
 * In arg:    matA, matB, matC n
 * Out arg:   int 1 if matrix multiplication is correct, 0 otherwise
 --------------------------------------------------------------------*/

int verifyMult(int** matA, int** matB, int **matC, int n){
  int i,j,k,sum=0;
int *matD[n];
  for (int i=0; i<n; i++){
         matD[i] = (int *)malloc(n * sizeof(int));
       }

  for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) {
        for (k = 0; k < n; k++) {
          sum = sum + (k-i)*(n-j+k);
        }
        matD[i][j] = sum;
        if(matC[i][j]!=matD[i][j])
          return 0;
        sum = 0;
      }
    }
    return 1;
}
