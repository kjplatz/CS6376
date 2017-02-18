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
#include <string.h>
#include <math.h>
#include <sys/time.h>

#include "config.h"

double Temperature[ROWS+2][COLUMNS+2];      // temperature grid
double Temperature_last[ROWS+2][COLUMNS+2]; // temperature grid from last iteration

//   helper routines
void initialize();
void track_progress(int iter, double dt, double*, double* );

double jacobi_loop( int row, double *restrict Temp, double *restrict Temp_last );

int main(int argc, char *argv[]) {

    int i, j;                                            // grid indexes
    int max_iterations;                                  // number of iterations
    int iteration=1;                                     // current iteration
    double dt=100;                                       // largest change in t
    struct timeval start_time, stop_time, elapsed_time;  // timers

    printf("Maximum iterations [100-4000]?\n");
    scanf("%d", &max_iterations);

    initialize();                   // initialize Temp_last including boundary conditions
    gettimeofday(&start_time,NULL); // Unix timer

    double* Temp = (double*)Temperature;
    double* Temp_last = (double*)Temperature_last;

    // do until error is minimal or until max steps
    while ( dt > MAX_TEMP_ERROR && iteration <= max_iterations ) {

        dt = 0.0; // reset largest temperature change
        // main calculation: average my four neighbors
        #pragma omp parallel for reduction(max:dt)
        for(i = 1; i <= ROWS; i++) {
            dt = fmax( dt, jacobi_loop( i, Temp, Temp_last ) );
        }
        
        double* tmp = Temp;
        Temp = Temp_last;
        Temp_last = tmp;

        // periodically print test values
        if((iteration % 100) == 0) {
 	    track_progress(iteration, dt, Temp, Temp_last);
        }

	iteration++;
    }

    gettimeofday(&stop_time,NULL);
    timersub(&stop_time, &start_time, &elapsed_time); // Unix time subtract routine

    printf("\nMax error at iteration %d was %f\n", iteration-1, dt);
    printf("Total time was %f seconds.\n", elapsed_time.tv_sec+elapsed_time.tv_usec/1000000.0);

}


// initialize plate and boundary conditions
// Temp_last is used to to start first iteration
void initialize(){

    int i,j;

    memset( Temperature_last, 0, sizeof(Temperature_last) );
    memset( Temperature, 0, sizeof(Temperature) );
    for(i = 0; i <= ROWS+1; i++){
        for (j = 0; j <= COLUMNS+1; j++){
            Temperature_last[i][j] = 0.0;
        }
    }

    // these boundary conditions never change throughout run

    // set left side to 0 and right to a linear increase
    for(i = 0; i <= ROWS+1; i++) {
        Temperature_last[i][0] = 0.0;
        Temperature[i][COLUMNS+1] = Temperature_last[i][COLUMNS+1] = (100.0/ROWS)*i;
    }
    
    // set top to 0 and bottom to linear increase
    for(j = 0; j <= COLUMNS+1; j++) {
        Temperature_last[0][j] = 0.0;
        Temperature[ROWS+1][j] = Temperature_last[ROWS+1][j] = (100.0/COLUMNS)*j;
    }
}


// print diagonal in bottom right corner where most action is
void track_progress(int iteration, double dt, double *restrict Temperature, double *restrict Temperature_last) {

    int i;

    printf("---------- Iteration number: %d ------------\n", iteration);
    printf( "[%d,%d]: %5.2f  ", 250, 900, Temperature[250*(COLUMNS+2)+900] );
    for(i = ROWS-5; i <= ROWS; i++) {
        printf("[%d,%d]: %5.2f  ", i, i, Temperature[i*(COLUMNS+2)+i]);
    }
    printf( "  max error=%7.4f  ", dt );
    printf("\n");
}
