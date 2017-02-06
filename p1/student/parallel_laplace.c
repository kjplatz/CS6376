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
#include <omp.h>

// size of plate
#define COLUMNS    1000
#define ROWS       1000
#define CHUNK (ROWS / 10)

// largest permitted change in temp (This value takes about 3400 steps)
#define MAX_TEMP_ERROR 0.01

double Temperature[ROWS+2][COLUMNS+2];      // temperature grid
double Temperature_last[ROWS+2][COLUMNS+2]; // temperature grid from last iteration

//   helper routines
void initialize();
void track_progress(int iter);
void sim_block(int, int);


int main(int argc, char *argv[]) {

    int i, j;                                            // grid indexes
    int max_iterations;                                  // number of iterations
    int iteration=1;                                     // current iteration
    double dt=100;                                       // largest change in t
    struct timeval start_time, stop_time, elapsed_time;  // timers

    int thread_count = omp_get_num_threads();

    printf("Maximum iterations [100-4000]?\n");
    //scanf("%d", &max_iterations);
    max_iterations=10000;
	//omp_set_nuthreads(4);

    initialize();                   // initialize Temp_last including boundary conditions
	printf("Num threads: %d\n", omp_get_num_threads());
	printf("Num procs: %d\n", omp_get_num_procs());
	omp_set_num_threads( omp_get_num_procs() );
	printf("Num threads: %d\n", omp_get_num_threads());

    double j_avg = 0, u_avg = 0;

    gettimeofday(&start_time,NULL); // Unix timer
    double start = omp_get_wtime( );

    // do until error is minimal or until max steps
    while ( dt > MAX_TEMP_ERROR && iteration <= max_iterations ) {
        /*********************************
         * BLOCK LAPLACE
        *********************************/

	/*
        # pragma omp parallel for //private(j)// schedule(runtime)
        for(i = 1; i <= ROWS-(CHUNK-1); i+=CHUNK) {
            # pragma omp parallel for //private(j)// schedule(runtime)
            for(j = 1; j <= COLUMNS-(CHUNK-1); j+=CHUNK) {
                sim_block(i,j);
            }
        }

        dt = 0.0; // reset largest temperature change

        // copy grid to old grid for next iteration and find latest dt
        //clock_t start1 = clock(), diff1;
        # pragma omp parallel for reduction(max:dt) private (j) schedule(static)
        for(i = 1; i <= ROWS; i++){
            for(j = 1; j <= COLUMNS; j++){
                dt = fmax( fabs(Temperature[i][j]-Temperature_last[i][j]), dt);
                Temperature_last[i][j] = Temperature[i][j];
            }
        }
	*/
        /*********************************/

        // main calculation: average my four neighbors
        # pragma omp parallel for private(j)// schedule(runtime)
        for(i = 1; i <= ROWS; i++) {
            for(j = 1; j <= COLUMNS; j++) {
                Temperature[i][j] = 0.25 * (Temperature_last[i+1][j] + Temperature_last[i-1][j] +
                                            Temperature_last[i][j+1] + Temperature_last[i][j-1]);
            }
        }

        dt = 0.0; // reset largest temperature change

        // copy grid to old grid for next iteration and find latest dt
        # pragma omp parallel for reduction(max:dt) private (j) schedule(static)
        for(i = 1; i <= ROWS; i++){
            for(j = 1; j <= COLUMNS; j++){
                dt = fmax( fabs(Temperature[i][j]-Temperature_last[i][j]), dt);
                Temperature_last[i][j] = Temperature[i][j];
            }
        }

        // periodically print test values
        if((iteration % 100) == 0) {
         track_progress(iteration);
         printf("dt: %f, T[1000][1000]: %f\n", dt, Temperature_last[1000][1000]);
         //printf("Jacobi update running avg: %d milliseconds\n", j_avg);
         //printf("Max and matrix copy running avg: %d milliseconds\n", u_avg);
        }

    iteration++;
    }

    gettimeofday(&stop_time,NULL);

    double end = omp_get_wtime();
    double etime = end-start;

    double GFLOPS = (double)(iteration-1 )*5.f*ROWS*COLUMNS / (1000000000.f*(end-start));
    //printf("iters: %d\n", iteration-1);
    printf("GFLOPS: %f\n", GFLOPS);

    timersub(&stop_time, &start_time, &elapsed_time); // Unix time subtract routine

    printf("T[%d][%d] = %f\n", 250,900,Temperature[250][900]);
    printf("\nMax error at iteration %d was %f\n", iteration-1, dt);
    printf("Total time was %f seconds.\n", elapsed_time.tv_sec+elapsed_time.tv_usec/1000000.0);

}

void sim_block(int i, int j) {
    int row_end = i+(CHUNK-1);
    int col_end = j+(CHUNK-1);
    int col_start = j;
    for(; i <= row_end; i++)
        for (j=col_start; j <= col_end; j++) {
            Temperature[i][j] = 0.25 * (Temperature_last[i+1][j] + Temperature_last[i-1][j] +
                                        Temperature_last[i][j+1] + Temperature_last[i][j-1]);
        }
}

// initialize plate and boundary conditions
// Temp_last is used to to start first iteration
void initialize(){

    int i,j;

    for(i = 0; i <= ROWS+1; i++){
        for (j = 0; j <= COLUMNS+1; j++){
            Temperature_last[i][j] = 0.0;
        }
    }

    // these boundary conditions never change throughout run

    // set left side to 0 and right to a linear increase
    for(i = 0; i <= ROWS+1; i++) {
        Temperature_last[i][0] = 0.0;
        Temperature_last[i][COLUMNS+1] = (100.0/ROWS)*i;
    }
    
    // set top to 0 and bottom to linear increase
    for(j = 0; j <= COLUMNS+1; j++) {
        Temperature_last[0][j] = 0.0;
        Temperature_last[ROWS+1][j] = (100.0/COLUMNS)*j;
    }
}


// print diagonal in bottom right corner where most action is
void track_progress(int iteration) {

    int i,j;

    printf("---------- Iteration number: %d ------------\n", iteration);
    /*
    for(i = 0; i <= ROWS; i+= 100) {
        for (j = 0; j <= COLUMNS; j+=100)
            printf("%5.2f  ", Temperature[i][i]);
        printf("| %5.2f \n", Temperature_last[i][1001]);
    }
    printf("---------------------------------------------------------------------------\n");
    for(i = 0; i <= ROWS; i+= 100)
            printf("%5.2f  ", Temperature_last[1001][i]);
        printf("\n");
            //printf("[%d,%d]: %5.2f  ", i, i, Temperature[i][i]);
    */
    printf("Test point [%d,%d]: %5.2f  ", 250, 900, Temperature[250][900]);
    printf("\n");
}
