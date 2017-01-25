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
#include <sys/time.h>

// size of plate
#define COLUMNS    1000
#define ROWS       1000

// largest permitted change in temp (This value takes about 3400 steps)
#define MAX_TEMP_ERROR 0.01

double TempA[ROWS+2][COLUMNS+2];      // temperature grid
double TempB[ROWS+2][COLUMNS+2];      // temperature grid from last iteration

double *Temperature=&TempA[0][0];
double *Temperature_last=&TempB[0][0];

//   helper routines
void initialize();
void track_progress(int iter);


int main(int argc, char *argv[]) {

    int i, j;                                            // grid indexes
    int max_iterations=100;                              // number of iterations
    int iteration=1;                                     // current iteration
    double dt=100;                                       // largest change in t
    struct timeval start_time, stop_time, elapsed_time;  // timers

    if ( argc > 1 ) max_iterations = atoi( argv[1] );

    printf("Maximum iterations: %d\n", max_iterations );


    initialize();                   // initialize Temp_last including boundary conditions
    gettimeofday(&start_time,NULL); // Unix timer

    char doagain = 1;
    // do until error is minimal or until max steps
    while ( doagain && iteration <= max_iterations ) {
        doagain = 0;

        double *LastA, *LastB, *LastC, *CurrTemp;
        // main calculation: average my four neighbors
        #pragma omp parallel for private(LastA, LastB, LastC, CurrTemp, j)
        for(i = 1; i <= ROWS; i++) {
            LastA = Temperature_last+( (i-1) * (COLUMNS + 2) );
            LastB = Temperature_last+( i * (COLUMNS + 2) );
            LastC = Temperature_last+( (i+1) * (COLUMNS + 2) );

            CurrTemp = Temperature+( i * (COLUMNS + 2) );

            __builtin_prefetch( LastA+8 );
            __builtin_prefetch( LastA+16 );
            __builtin_prefetch( LastB+8 );
            __builtin_prefetch( LastB+16 );
            __builtin_prefetch( LastC+8 );
            __builtin_prefetch( LastC+16 );
            __builtin_prefetch( CurrTemp+8 );
            __builtin_prefetch( CurrTemp+16 );
            
            for(j = 1; j <= COLUMNS; j++) {
                if ( !(j % 8) ) {
                    __builtin_prefetch( LastA+j+16 );
                    __builtin_prefetch( LastB+j+16 );
                    __builtin_prefetch( LastC+j+16 );
                    __builtin_prefetch( CurrTemp+j+16 );
                }
                CurrTemp[j] = 0.25 * (LastA[j] + LastB[j-1] + LastB[j+1] + LastC[j]);
                /*Temperature[i][j] = 0.25 * (Temperature_last[i+1][j] + Temperature_last[i-1][j] +
                                            Temperature_last[i][j+1] + Temperature_last[i][j-1]); */
            }

            if ( !doagain && fabs( LastB[j] - CurrTemp[j] ) > MAX_TEMP_ERROR ) {
                #pragma omp atomic write
                doagain = 1;
            }
        }
        
        double* Temp = Temperature;
        Temperature = Temperature_last;
        Temperature_last = Temp;

        // periodically print test values
        if((iteration % 100) == 0) {
 	    track_progress(iteration);
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

    for(i = 0; i <= ROWS+1; i++){
        for (j = 0; j <= COLUMNS+1; j++){
            Temperature_last[i*(COLUMNS+2)+j] = 0.0;
        }
    }

    // these boundary conditions never change throughout run

    // set left side to 0 and right to a linear increase
    for(i = 0; i <= ROWS+1; i++) {
        Temperature_last[i*(COLUMNS+2)] = 0.0;
        Temperature_last[i*(COLUMNS+2)+COLUMNS+1] = (100.0/ROWS)*i;
    }
    
    // set top to 0 and bottom to linear increase
    for(j = 0; j <= COLUMNS+1; j++) {
        Temperature_last[j] = 0.0;
        Temperature_last[(ROWS+1)*(COLUMNS*2)+j] = (100.0/COLUMNS)*j;
    }
}


// print diagonal in bottom right corner where most action is
void track_progress(int iteration) {

    int i;

    printf("---------- Iteration number: %d ------------\n", iteration);
    for(i = ROWS-5; i <= ROWS; i++) {
        printf("[%d,%d]: %5.2f  ", i, i, Temperature[i*(COLUMNS+2)+i]);
    }
    printf("\n");
}
