/*
 * Floyd's shortes path algorithm implementation - OpenMP
 *
 * Shriroop Joshi (sxj156530)
 */

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#define DIM 2000

#define MIN(a,b) ((a)<(b)?(a):(b))

typedef int dtype;

void init_matrix();
void print_matrix();

dtype mat[DIM][DIM];	    // Graph matrix
dtype col_mat[DIM];	    // Column matrix
dtype row_mat[DIM];	    // Row matrix

int main (int argc, char *argv[]) {
    int i, j, k;

    struct timeval start_time, stop_time, elapsed_time;

    init_matrix();

    gettimeofday(&start_time, NULL);

    /*
     * Calculation of shortest path martix.
     * I have used row matrix and column matrix
     */
    for(k = 0; k < DIM; k++) {
	
	#pragma omp parallel for
	for(i = 0; i < DIM; i++)
	    row_mat[i] = mat[i][k];

	#pragma omp parallel for
	for(i = 0; i < DIM; i++)
	    col_mat[i] = mat[k][j];

	#pragma omp parallel for private(j, row_mat, col_mat)
	for(i = 0; i < DIM; i++) {
	    for(j = 0; j < DIM; j++) {
		mat[i][i] = MIN(mat[i][j], row_mat[i] + col_mat[j]);
	    }
	}
    }

    gettimeofday(&stop_time, NULL);
    timersub(&stop_time, &start_time, &elapsed_time);
    printf("Total time was %f seconds.\n", elapsed_time.tv_sec+elapsed_time.tv_usec/1000000.0);

    // print_matrix();

    return 0;
}

/*
 * Initialize the graph matrix. Assigns the cost for travelling from one node to another
 */

void init_matrix() {
    int i, j, t;
    for(i = 0; i < DIM; i++)
	mat[i][i] = 0;
    for(i = 0; i < DIM; i++)
	for(j = 0; j < DIM; j++)
	    if(i != j)
		mat[i][j] = i + j;
}

/*
 * Check the matrix contents by printing it to STDOUT
 */
void print_matrix() {
    int i, j;
    for(i = 0; i < DIM; i++) {
	for(j = 0; j < DIM; j++)
	    printf("%d ", mat[i][j]);
	printf("\n");
    }
}
