// --------------------------------------------------
// smm.c
// Jay Butera
//
// rec_matmul - recursive matrix multiplication
// --------------------------------------------------
#include "smm.h"
#include "common.h"
#include <stdio.h>
#include <omp.h>

// rec_matmul
// Jay Butera
//
// Description: Recursive matrix multiply for global variables A * B = C
//
// crow/ccol - I/P - int - row and column start index for matrix C
// arow/acol - I/P - int - row and column start index for matrix A
// brow/bcol - I/P - int - row and column start index for matrix B
// l/m/n     - I/P - int - matrix dimensions. A is (lxm) B is (mxn) C is (lxn)
void rec_matmul (int crow, int ccol,
                 int arow, int acol,
                 int brow, int bcol,
                 int l, int m, int n)
                 //float** a, float** b, float** c, int N)
{
    int lhalf[3], mhalf[3], nhalf[3];
    int i,j,k;
    DTYPE *aptr, *bptr, *cptr;

    // THRESHOLD usually corresponds to cache size
    if (m * n > THRESHOLD) {
        lhalf[0] = 0; lhalf[1] = l/2; lhalf[2] = l - l/2;
        mhalf[0] = 0; mhalf[1] = m/2; mhalf[2] = m - m/2;
        nhalf[0] = 0; nhalf[1] = n/2; nhalf[2] = n - n/2;

        // Recursively divide matrices A and B
        for (i = 0; i < 2; i++)
            for (j = 0; j < 2; j++)
                for (k = 0; k < 2; k++) {
                    rec_matmul( crow+lhalf[i], ccol+nhalf[j],
                        arow+lhalf[i], acol+mhalf[k],
                        brow+mhalf[k], bcol+nhalf[j],
                        lhalf[i+1], mhalf[k+1], nhalf[j+1]);
                        //a,b,c,N);
                }
    }
    else {
#ifdef USE_ACC // OpenACC
        #pragma acc data copyin (A[arow:arow+l][acol:acol+m], B[brow:brow+m][bcol:bcol+n], C[crow:crow+l][ccol:ccol+n]) copyout (C[crow:crow+l][ccol:ccol+n])
#endif
#ifdef USE_OMP // OpenMP
        #pragma omp parallel shared(A,B,C) private(j,k)
#endif
        {
#ifdef USE_ACC
            #pragma acc loop independent
#endif
#ifdef USE_OMP
        #pragma omp parallel for schedule(static)
#endif
            for (i = 0; i < l; i++) {
                //printf("thread %d did row %d\n", tid, i);
#ifdef USE_ACC
                #pragma acc loop independent
#endif
                for (j = 0; j < n; j++) {
                    cptr = &C[crow+i][ccol+j];
                    aptr = &A[arow+i][acol];
                    bptr = &B[brow][bcol+j];

                    for (k = 0; k < m; k++) {
                        *cptr += *(aptr++) * *bptr;
                        bptr += SIZE;
                    }
                }
            }
        }
    }
}
