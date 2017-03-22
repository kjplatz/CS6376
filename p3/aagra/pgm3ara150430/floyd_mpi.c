/* file - floyd_mpi.c 
Author - Ayush Agrawal
Description - Implements Floyd alogorithm using MPI
*/

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <mpi.h>

const int INFINITY = 1000000;

void generateMatrix(int mat[], int n);
void printMatrix(int mat[], int n);
int min(int x, int y);
void floyd(int p, int n, int local_mat[], int my_rank);


int main(int argc, char* argv[])
{
    int p;
    int my_rank;
    int n = 2000; // number of nodes
    int* mat;
    int* local_mat;
    int* temp_mat;
	struct timeval start_time, stop_time, elapsed_time;  // timers
	
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

   
    generateMatrix( mat, n);
  
    /* Buffer allocation for local rows */
    local_mat = malloc(n * (n/p) * sizeof(int));

    /* Buffer allocation for revised matrix */
    temp_mat = malloc(n * n * sizeof(int));

    /* Broadcasts n (number of "cities") to each processor */
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

    /* Distributes matrix among the processors */
    MPI_Scatter(mat, n * (n/p), MPI_INT,
            local_mat, n * (n/p), MPI_INT, 0, MPI_COMM_WORLD);

	gettimeofday(&start_time,NULL); // Unix  start timer
	
    /* Uses Floyd's algo to compute shortest paths */
    floyd(p, n, local_mat, my_rank);

	gettimeofday(&stop_time,NULL);// Unix end timer
    timersub(&stop_time, &start_time, &elapsed_time); // Unix time subtract routine
    printf("Total time was %f seconds.\n", elapsed_time.tv_sec+elapsed_time.tv_usec/1000000.0);  
	
    /* Gathers the data */
    MPI_Gather(local_mat, n * (n/p), MPI_INT, temp_mat,
            n * (n/p), MPI_INT, 0, MPI_COMM_WORLD);

    /* Prints the matrix */
    /* if (my_rank == 0)
    {
        printf("The solution is:\n");
        printf("\n");
        printMatrix(temp_mat, n);
        printf("\n");
    } */

    MPI_Finalize();
    return(0);
} 

/* Creates adjacency matrix */
void generateMatrix(int mat[], int n)
{
    int i, j;

    for (i = 0; i < n; i++)
       for (j = 0; j < n; j++)
	        if (i==j)
                mat[i * n + j]=0;
            else
                mat[i * n + j]= (int)( 3.0 * rand() / ( RAND_MAX + 1.0 ) );

} 

/* Prints adjacency matrix */
void printMatrix(int mat[], int n)
{
    int i, j;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            if (mat[i * n + j] == INFINITY)
                printf("i ");
            else
                printf("%d ", mat[i * n + j]);
        }
        printf("\n");
    }
}

/* Returns minimum of two ints--used in Floyd's algorithm */
int min(int x, int y)
{
    if (x < y)
        return x;
    else
        return y;
} 

/* Floyd's Algorithm */
void floyd(int p, int n, int local_mat[], int my_rank)
{
    int* tmp;
    int offset;
    int root;

    tmp = malloc(n * sizeof(int));

    int k;
    for (k = 0; k < n; k++) {
        root = k / (n / p);
        if (my_rank == root) {
            offset = k % (n / p);
            int j;
            for (j = 0; j < n; j++)
                tmp[j] = local_mat[offset * n + j];
        }
        MPI_Bcast(tmp, n, MPI_INT, root, MPI_COMM_WORLD);
        int i;
        for (i = 0; i < n / p; i++)
        {
            int j;
            for (j = 0; j < n; j++)
            {
                local_mat[i * n + j] = min(local_mat[i * n + j], local_mat[i*n + k] + tmp[j]);
            }
        }
    }
    free(tmp);
} 

