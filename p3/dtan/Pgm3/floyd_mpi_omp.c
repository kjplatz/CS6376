/* Program3      floyd_mpi.c
 * Name: David Tan
 * 
 * n is evenly divisible by p.
 * 0 on the diagonal.
 * The matrix is distributed by block rows.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h> 
#include <omp.h>
#include <mpi.h>
#include <math.h>
#include <sys/time.h>


void Read_matrix(int local_mat[], int n, int my_rank, int p, MPI_Comm comm);

void Print_matrix(int local_mat[], int n, int my_rank, int p, MPI_Comm comm);

void Floyd(int local_mat[], int n, int my_rank, int p, MPI_Comm comm);
int Owner(int k, int p, int n);
void Copy_row(int local_mat[], int n, int p, int row_k[], int k);
void Print_row(int local_mat[], int n, int my_rank, int i);

int main(int argc, char* argv[]) {
	int  n;
    int* local_mat;
    MPI_Comm comm;
    int p, my_rank;
	struct timeval start_time, stop_time, elapsed_time;  // timers


    MPI_Init(&argc, &argv);
    comm = MPI_COMM_WORLD;
    MPI_Comm_size(comm, &p);
    MPI_Comm_rank(comm, &my_rank);

    if (my_rank == 0) {
		printf("How many vertices?\n");
        scanf("%d", &n);
    }
    MPI_Bcast(&n, 1, MPI_INT, 0, comm);
    local_mat = malloc(n*n/p*sizeof(int));  //rectange block row

    Read_matrix(local_mat, n, my_rank, p, comm);
	
	if(my_rank == 0) gettimeofday(&start_time,NULL); // Unix timer

    Floyd(local_mat, n, my_rank, p, comm);

	if(my_rank == 0){
		gettimeofday(&stop_time,NULL);
		timersub(&stop_time, &start_time, &elapsed_time); // Unix time subtract routine
		printf("Total time was %f seconds.\n", elapsed_time.tv_sec+elapsed_time.tv_usec/1000000.0);
		double t = elapsed_time.tv_sec+elapsed_time.tv_usec/1000000.0;
		
		double temp = n/1000.0;
		//double flop = 2*n*n*n/(t*1000000000);
		double flop = 2*temp*temp*temp/t;
		printf("The program's GFlops:  %f\n",flop);
	
	}
	
    //if (my_rank == 0) printf("Result matirx:\n");
    //Print_matrix(local_mat, n, my_rank, p, comm);

    free(local_mat);
    MPI_Finalize();

    return 0;
}

/*---------------------------------------------------------------------
 * Read in matrix on process 0 and scatter it by row distribution
 */
void Read_matrix(int local_mat[], int n, int my_rank, int p, MPI_Comm comm) { 

	int i, j;
	int* temp_mat = NULL;

    if (my_rank == 0) {
		temp_mat = malloc(n*n*sizeof(int));  //square temporary matrix
		
		//build test matrix
		#pragma omp parallel for private(i,j)
	 	for (i = 0; i < n; i++){
			for (j = 0; j < n; j++){
				//temp_mat[i*n+j] = rand()%100;
				//if(temp_mat[i*n+j]==0)temp_mat[i*n+j] = 99999;
				temp_mat[i*n+j] = 99999;
				if( (i-j==1) || (j-i==1)) temp_mat[i*n+j]=1;
				if( (i==0)&&(j==(n-1))) temp_mat[i*n+j]=1;
				if( (j==0)&&(i==(n-1))) temp_mat[i*n+j]=1; 
				if(i == j) temp_mat[i*n+j] = 0;
			}
		}
		// MPI_Scatter( void* send_data, int send_count, MPI_Datatype send_datatype,
		//              void* recv_data, int recv_count, MPI_Datatype recv_datatype,
		//			    int root, MPI_Comm communicator)
        MPI_Scatter(temp_mat, n*n/p, MPI_INT, local_mat, n*n/p, MPI_INT, 0, comm);
        free(temp_mat);
    } else {
        MPI_Scatter(temp_mat, n*n/p, MPI_INT, 
                  local_mat, n*n/p, MPI_INT, 0, comm);
    }

}

/*---------------------------------------------------------------------
 * distributed floyd's algorithm for
 * block rows
 */
void Floyd(int local_mat[], int n, int my_rank, int p, MPI_Comm comm) {
    int global_k, local_i, global_j, temp;
    int root;
    int* row_k = malloc(n*sizeof(int)); //temp storage
	#pragma omp parallel for private(local_i,global_j,global_k)  schedule(static,16)
    for (global_k = 0; global_k < n; global_k++) {
		root = global_k/(n/p);
        if (my_rank == root){
			int j;
			int local_k = k % (n/p);

			for (j = 0; j < n; j++)
				row_k[j] = local_mat[local_k*n + j];
		}
        MPI_Bcast(row_k, n, MPI_INT, root, comm);
        for (local_i = 0; local_i < n/p; local_i++)
            for (global_j = 0; global_j < n; global_j++) {
                temp = local_mat[local_i*n + global_k] + row_k[global_j];
                if (temp < local_mat[local_i*n+global_j])
                    local_mat[local_i*n + global_j] = temp;
            }
    }
    free(row_k);
}

