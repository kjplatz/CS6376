/** 
 * @file    floyd.c
 * @author  Kevin Gillespie
 * @brief   floyd's algorithm to calculate the shortest path between all pairs
 *          of a give graph. We're assuming symmetric matricies. 
 *
 */

/* Includes */
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include "mpi.h"

/* Defines */
// #define ARRAY_DIM 2048  
#define ARRAY_DIM   8192 
#define MAX_VAL     ARRAY_DIM

#define MIN(a,b) (((a)<(b))?(a):(b))

/* Globals */
 

/**
 * @name     print_array
 * @brief    print the given 2d array
 * @param 
 *       @name   array
 *       @dir    I
 *       @type   float
 *       @brief  2d array to print
 *
 ******************************************************************************/
void print_array(float* array[])
{
    unsigned i, j;

    // Print the j labels
    printf(" j:");
    for(j = 0; j < ARRAY_DIM; j++) {
        printf("%d  ", j);
    }
    printf("\ni ");
    for(j = 0; j < ARRAY_DIM; j++) {
        printf("---");
    }

    // Print the array and i labels
    for(i = 0; i < ARRAY_DIM; i++) {
        printf("\n%d| ", i);
        for(j = 0; j < ARRAY_DIM; j++) {
            printf("%02d ", (unsigned)array[i][j]);
        }
    }
    printf("\n\n");
}

/**
 * @name     main
 * @brief    main function for floyd.c
 * @param 
 *       @name   argc
 *       @dir    I
 *       @type   int 
 *       @brief  Number of arguments in argv.
 * @param 
 *       @name   argv
 *       @dir    I
 *       @type   char*[]
 *       @brief  Command line arguments.
 *
 * @returns 0 for success, error status otherwise
 *
 ******************************************************************************/
int main(int argc, char *argv[])
{
    unsigned i,j,k;
    struct timeval start_time, stop_time, elapsed_time;
    double etime, flops;

    int id, p, i_offset, root;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    i_offset = id*(ARRAY_DIM/p);

    // Initialize the array
    float* array[ARRAY_DIM/p];
    float* tmp = (float*)malloc(ARRAY_DIM*sizeof(float)); 
    for(i = 0; i < (ARRAY_DIM/p); i++) {
        array[i] = (float*)malloc(ARRAY_DIM * sizeof(float));
        for(j = 0; j < ARRAY_DIM; j++) {

            // Offset the value of i 
            int i_offset_temp = i+i_offset; 
            if(i_offset_temp == j) {
                array[i][j] = 0;
            } else if((i_offset_temp-j) == 1) {
                array[i][j] = 1;
            } else if((j-i_offset_temp) == 1) {
                array[i][j] = 1;
            } else {
                array[i][j] = MAX_VAL;
            }
        }
    }

    if(ARRAY_DIM < 32) {
        printf("Initial array\n");
        print_array(array);
    }
 
    MPI_Barrier(MPI_COMM_WORLD);

    gettimeofday(&start_time, NULL);

    // Execute the algorithm
    for(k = 0; k < ARRAY_DIM; k++) {

        // Have to copy out the kth row 
        root = k/(ARRAY_DIM/p);
        if(root == id) {
            int offset_temp = k - i_offset;
            #pragma omp parallel for private(i) schedule(static)
            for(i = 0; i < ARRAY_DIM; i++) {
                tmp[i] = array[offset_temp][i];
            }
        }

        MPI_Bcast(tmp, ARRAY_DIM, MPI_FLOAT, root, MPI_COMM_WORLD);
        #pragma omp parallel for private(i,j) schedule(static)
        for(i = 0; i < (ARRAY_DIM/p); i++) {
            for(j = 0; j < ARRAY_DIM; j++) {
                array[i][j] = MIN(array[i][j], (array[i][k] + tmp[j]));
            }
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);   

    gettimeofday(&stop_time, NULL);

    timersub(&stop_time, &start_time, &elapsed_time);    
    etime = elapsed_time.tv_sec+elapsed_time.tv_usec/1000000.0;

    if(ARRAY_DIM < 32) {
        printf("\nFinal array\n");
        print_array(array);
    }

    // Check the output of the matrix
    for(i = 0; i < ARRAY_DIM/p; i++) {
        for(j = 0; j < ARRAY_DIM; j++) {
            int i_offset_temp = i+i_offset; 
            if(array[i][j] != (float)abs(i_offset_temp-j)) {
                printf("Array error! i = %d j= %d array[i][j] = %d\n", i_offset_temp, j, (unsigned)array[i][j]);
                return 1;
            }
        }
    }


    if(0 == id) {
        flops = ((double)2 * (double)ARRAY_DIM * (double)ARRAY_DIM * (double)ARRAY_DIM)/etime;
        printf("%d, %f, %f, %d\n", ARRAY_DIM, etime, flops, omp_get_max_threads());
    }
 
    MPI_Finalize();
    return 0;
}
