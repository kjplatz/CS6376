/*************************************************
 * Laplace Serial C Version
 *
 * Temperature is initially 0.0
 * Boundaries are as follows:
 *
 *      0         T         0
 *   0  +-------------------+  0
 *      |                   |
 *      |                   |
 *      |                   |
 *   T  |                   |  T
 *      |                   |
 *      |                   |
 *      |                   |
 *   0  +-------------------+ 100
 *      0         T        100
 *
 *  John Urbanic, PSC 2014
 *
 ************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <openacc.h>
#include "openacc.h"
#include <sys/time.h>

// size of plate
#define COLUMNS    4096
#define ROWS       4096

// largest permitted change in temp (This value takes about 3400 steps)
#define MAX_TEMP_ERROR 0.01

float Temperature[ROWS+2][COLUMNS+2];      // temperature grid
float Temperature_last[ROWS+2][COLUMNS+2]; // temperature grid from last iteration

//   helper routines
void initialize();
void track_progress(int iter);

int main(int argc, char *argv[]) {

    int i, j;                                            // grid indexes
    int max_iterations;                                  // number of iterations
    int iteration=1;                                     // current iteration
    float dt=100;                                       // largest change in t
    struct timeval start_time, stop_time, elapsed_time;  // timers

    printf("Maximum iterations [100-4000]?\n");
    scanf("%d", &max_iterations);


    initialize();                   // initialize Temp_last including boundary conditions
    gettimeofday(&start_time,NULL); // Unix timer

    // do until error is minimal or until max steps
    #pragma acc data copyin(Temperature_last) copy(Temperature)
    while ( dt > MAX_TEMP_ERROR && iteration <= max_iterations ) {
        void *dev_Temperature, *dev_Temperature_last;
        dev_Temperature = acc_deviceptr(Temperature);
        dev_Temperature_last = acc_deviceptr(Temperature_last);
        dt = 0.0; // reset largest temperature change

        // main calculation: average my four neighbors
        #pragma acc kernels loop reduction(max:dt)
        for(i = 1; i <= ROWS; i++) {
            for(j = 1; j <= COLUMNS; j++) {
                Temperature[i][j] = 0.25 * (Temperature_last[i+1][j] + Temperature_last[i-1][j] +
                                            Temperature_last[i][j+1] + Temperature_last[i][j-1]);
                dt = fmax(dt, fabs(Temperature[i][j] - Temperature_last[i][j]));
            }
        }
        

#ifndef USE_ACC_MEMCPY
        // copy grid to old grid for next iteration and find latest dt
        #pragma acc kernels loop 
        for(i = 1; i <= ROWS; i++){
            for(j = 1; j <= COLUMNS; j++){
	      Temperature_last[i][j] = Temperature[i][j];
            }
        }
#else
        acc_memcpy_device(dev_Temperature_last, dev_Temperature, sizeof(Temperature));
#endif

        // periodically print test values
        if((iteration % 100) == 0) {
            printf( "Iteration: %d DeltaT: %.2f\n", iteration, dt );
            // #pragma acc data update host(Temperature)
            // acc_update_host(Temperature, sizeof(Temperature));
 	    // track_progress(iteration);
        }

	    iteration++;
    }

    gettimeofday(&stop_time,NULL);
    timersub(&stop_time, &start_time, &elapsed_time); // Unix time subtract routine

    printf("\nMax error at iteration %d was %f\n", iteration-1, dt);
    printf("Total time was %f seconds.\n", elapsed_time.tv_sec+elapsed_time.tv_usec/1000000.0);

    float gf = ROWS * COLUMNS * 5;
    gf /= (elapsed_time.tv_sec + elapsed_time.tv_usec/1000000.0);
    gf /= (1024 * 1024);
    printf( "Performance = %5.2f Gflops\n", gf );

}


// initialize plate and boundary conditions
// Temp_last is used to to start first iteration
void initialize(){

    int i,j;

    memset( (void*)Temperature, 0, sizeof(Temperature));
    memset( (void*)Temperature_last, 0, sizeof(Temperature_last));

    // these boundary conditions never change throughout run
    // set left side to 0 and right to a linear increase
    for(i = 0; i <= ROWS+1; i++) {
        Temperature[i][COLUMNS+1] = Temperature_last[i][COLUMNS+1] = (100.0/ROWS)*i;
    }
    
    // set top to 0 and bottom to linear increase
    for(j = 0; j <= COLUMNS+1; j++) {
        Temperature[ROWS+1][j] = Temperature_last[ROWS+1][j] = (100.0/COLUMNS)*j;
    }

    track_progress(0);
}


// print diagonal in bottom right corner where most action is
void track_progress(int iteration) {

    int i;

    printf("---------- Iteration number: %d ------------\n", iteration);
    printf( "Temp   :" );
    for(i = ROWS-5; i <= ROWS+1; i++) {
        printf("[%d,%d]: %5.2f  ", i, i, Temperature[i][i]);
    }
    printf("\n");
}
