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
#include "openacc.h"

/* Defines */
#define ARRAY_DIM 16
// #define ARRAY_DIM  8192 
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
void print_array(float array[ARRAY_DIM][ARRAY_DIM])
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
    float* restrict tmp = (float*)malloc(ARRAY_DIM*sizeof(float)); 

    int id, p, i_offset, root;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    i_offset = id*(ARRAY_DIM/p);


    printf("p = %d id = %d i_offset = %d\n", p, id, i_offset);

    // Initialize the array
    // float* restrict  array[ARRAY_DIM/p];
    float array[ARRAY_DIM/p][ARRAY_DIM];
    for(i = 0; i < (ARRAY_DIM/p); i++) {
        // array[i] = (float*)malloc(ARRAY_DIM * sizeof(float));
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

    memcpy(tmp, array[0], ARRAY_DIM*sizeof(float));

    if(ARRAY_DIM < 32) {
        printf("Initial array\n");
        print_array(array);
    }

    
    int start_offset, end_offset;
    int numgpus = acc_get_num_devices(acc_device_nvidia);
    int gpunum;
    
    printf("numgpus = %d\n", numgpus);

    MPI_Barrier(MPI_COMM_WORLD);
     
    gettimeofday(&start_time, NULL);

    // Execute the algorithm    

    for(gpunum = 0; gpunum < numgpus; ++gpunum) {
        acc_set_device_num(gpunum, acc_device_nvidia);
        
        start_offset = gpunum*((ARRAY_DIM/p)/numgpus);
        end_offset = start_offset + ((ARRAY_DIM/p)/numgpus);

        #pragma acc enter data create(array[start_offset:end_offset], tmp[0:ARRAY_DIM]) async(gpunum)
        #pragma acc update device(array[start_offset:end_offset], tmp[0:ARRAY_DIM]) async(gpunum)

    }

    acc_async_wait_all();

    for(k = 0; k < ARRAY_DIM; k++) {

        for(gpunum = 0; gpunum < numgpus; ++gpunum) {
            acc_set_device_num(gpunum, acc_device_nvidia);
        
            start_offset = gpunum*((ARRAY_DIM/p)/numgpus);
            end_offset = start_offset + ((ARRAY_DIM/p)/numgpus);

            #pragma acc kernels async(gpunum)
            for(i = start_offset; i < end_offset; i++) {
                for(j = 0; j < ARRAY_DIM; j++) {
                    array[i][j] = MIN(array[i][j], (array[i][k] + tmp[j]));
                }
            }
        }

    
        root = k/(ARRAY_DIM/p);

        for(gpunum = 0; gpunum < numgpus; ++gpunum) {
            acc_set_device_num(gpunum, acc_device_nvidia);
        
            start_offset = gpunum*((ARRAY_DIM/p)/numgpus);
            end_offset = start_offset + ((ARRAY_DIM/p)/numgpus);

            // Have to copy out the kth row 
            int i_offset_temp = k-i_offset;
            if(root == id && (i_offset_temp >= start_offset) && (i_offset_temp < end_offset))  {     
                #pragma acc kernels async(gpunum)
                for(i = 0; i < ARRAY_DIM; i++) {
                    tmp[i] = array[i_offset_temp][i];
                }
                #pragma acc update host(tmp[0:ARRAY_DIM]) async(gpunum)
            }
        }

        acc_async_wait_all();
        MPI_Bcast(tmp, ARRAY_DIM, MPI_FLOAT, root, MPI_COMM_WORLD);

        for(gpunum = 0; gpunum < numgpus; ++gpunum) {
            acc_set_device_num(gpunum, acc_device_nvidia);
        

            #pragma acc update device(tmp[0:ARRAY_DIM])
        }
    }

    acc_async_wait_all();

    for(gpunum = 0; gpunum < numgpus; ++gpunum) {
        acc_set_device_num(gpunum, acc_device_nvidia);
        
        start_offset = gpunum*((ARRAY_DIM/p)/numgpus);
        end_offset = start_offset + ((ARRAY_DIM/p)/numgpus);

        #pragma update host(array[start_offset:end_offset])
    }

    acc_async_wait_all();
    MPI_Barrier(MPI_COMM_WORLD);
    
    gettimeofday(&stop_time, NULL);
    timersub(&stop_time, &start_time, &elapsed_time);    
    etime = elapsed_time.tv_sec+elapsed_time.tv_usec/1000000.0;

    if(ARRAY_DIM < 32) {
        printf("\nFinal array\n");
        print_array(array);
    }

    // Check the output of the matrix
    int error = 0;
    for(i = 0; (i < ARRAY_DIM/p && !error); i++) {
        for(j = 0; (j < ARRAY_DIM && !error); j++) {
            int i_offset_temp = i+i_offset;
            if(array[i][j] != (float)abs(i_offset_temp-j)) {
                printf("Array error! i = %d j= %d array[i][j] = %d\n", i_offset_temp, j, (unsigned)array[i][j]);
                error = 1;
            }
        }
    }

    flops = ((double)2 * (double)ARRAY_DIM * (double)ARRAY_DIM * (double)ARRAY_DIM)/etime;
    if(0 == id) {
        printf("%d, %f, %f, %d\n", ARRAY_DIM, etime, flops, omp_get_max_threads());
    }

    MPI_Finalize();
    
    return 0;
}
