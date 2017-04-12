/* CS 6376.001 - Parallel Processing
 * Program 4 - Cannon's Algorithm
 * Tarunesh Verma (txv140330)
 * tarunesh.verma@utd.edu

 * Implementation of Cannon's Matrix Multiplication Algorithm for square matrices.
 * A, B, and C are matrices with dimensions N x N. The operation is A X B = C.
 *
 * Matrices A and B are divided across processes in blocks of sizes N/sqrt(P),
 * where N is the dimension of the matrices and P is the number of processes.
 *
 * In this program, the size of the block matrices is chosen to be 16K x 16K;
 * thus the dimension N of the three matrices turns out as N = 16K * sqrt(P).
 * Based on the above values, each process (node) generates its own blocks of
 * A and B. These blocks are multiplied to create C, which stores the partial
 * result. Blocks A and B are then rotated using MPI commands, i.e. each node 
 * sends its blocks to another node and receives a new pair to multiply. This
 * result is added to the node's running total ofthe block C. Thus, each node
 * is responsible for computing a part of C. 
 * The final matrix, along with the running time is then output.
 * 
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

/* Datatype of the matrix elements. */
#ifdef dbl
	typedef double DTYPE;		
	#define MPI_TYPE MPI_DOUBLE
#else
	typedef int DTYPE;
	#define MPI_TYPE MPI_INT
#endif

/* Specify whether matrices should be printed. */
#ifdef prt
	#define PRINT_MATRICES 1
#else
	#define PRINT_MATRICES 0
#endif

#define THRESHOLD 512           /* Max block size allowable in a process (node). */
#define PROC_SIZE 4096          /* Block size in a process. This value equals N/sqrt(P), where N is the dimension of the three matrices. */
#define ROOT_PROC 0             /* Variable to denote Process 0. */

#define BLOCK_LOW(id,p,n)  ((id)*(n)/(p))
#define BLOCK_HIGH(id,p,n) (BLOCK_LOW((id)+1,p,n)-1)
#define BLOCK_SIZE(id,p,n) (BLOCK_HIGH(id,p,n)-BLOCK_LOW(id,p,n)+1)

/* Function prototypes. */
DTYPE** create_matrix();
void mesh_setup(int *, int, int *, MPI_Comm *, int *, int *);
void multiply_matrices(int, int, int, int, int, int, int, int, int, DTYPE**, DTYPE**, DTYPE**);
void print_checkerboard_matrix(void **a, MPI_Datatype, int, int, MPI_Comm);

int main(int argc, char** argv)
{
	int num_procs, rank;                             /* Number of processes and rank of each process. */
	int mesh_size[2];                                /* Array to store the dimensions of the block matrix. */
	int mesh_coords[2];                              /* Array to store the coordinates of the block matrix in the process mesh. */
	int periodic[2];                                 /* Array to store the periodicity of the block matrix. */
	int mesh_rank;                                   /* Array to store the rank of the block matrix in the process mesh. */
	int shift_coords[2][2] = { { 0, 0 },{ 0, 0 } };  /* Array to store the coordinates for the rearrangement and rotation of the matrices. */
	double total_time;                               /* Time taken across all nodes. */
	double process_time;                             /* Time taken in one node. */
	MPI_Comm mesh_comm;                              /* Communicator for the process mesh. */
	MPI_Status status;                               /* MPI Status variable for MPI_Sendrecv_replace. */

	MPI_Init(&argc, &argv);                          /* Initialize MPI. */
	MPI_Comm_size(MPI_COMM_WORLD, &num_procs);       /* Get the numbe of processes. */
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);            /* Get the rank of each process. */

	int N = PROC_SIZE * sqrt(num_procs);             /* Determine the size of the three matrices. */

	/* If the number of processes is not a perfect square, exit. */
	if (sqrt(num_procs) * sqrt(num_procs) != num_procs)
	{
		MPI_Finalize();
		if (rank == 0)
			printf("Number of Processes (Nodes) should be a perfect square. Exiting...\n");
		exit(0);
	}
	
	/* Create a process mesh. */
	mesh_setup(mesh_size, num_procs, periodic, &mesh_comm, &mesh_rank, mesh_coords);

	process_time = - MPI_Wtime();                    /* Start measuring time for performance evaluation. */

	/* Create the three matrices A, B, and C. */
	DTYPE** A = create_matrix();
	DTYPE** C = create_matrix();
	DTYPE** B = create_matrix();

	/* Create values in block matrices across nodes. */
	for (int i = 0; i < PROC_SIZE; i++)
		for (int j = 0; j < PROC_SIZE; j++)
		{
			C[i][j] = 0;
			A[i][j] = (mesh_coords[1]*PROC_SIZE + j) - (mesh_coords[0]*PROC_SIZE + i);
			B[i][j] = N - (mesh_coords[1]*PROC_SIZE + j) + (mesh_coords[0]*PROC_SIZE + i);
		}

	if(rank == ROOT_PROC)
		printf("Matrices A, B, and C created. \n");
	if(PRINT_MATRICES && (rank == ROOT_PROC))
	{
		printf("Matrix A: \n-----------------------");
		print_checkerboard_matrix((void **) A, MPI_TYPE, N, N, mesh_comm);
		printf("Matrix B: \n-----------------------");
		print_checkerboard_matrix((void **) B, MPI_TYPE, N, N, mesh_comm);
	}

	/* Determine Coordinates for the initial rearrangement of A and B. */
	MPI_Cart_shift(mesh_comm, 1, mesh_coords[0], &shift_coords[0][1], &shift_coords[0][0]);         /* Blocks of Matrix A are shifted i places left in the process mesh. */
	MPI_Cart_shift(mesh_comm, 0, mesh_coords[1], &shift_coords[1][1], &shift_coords[1][0]);         /* Blocks of Matrix B are shifted j places above in the process mesh. */

	/* Rearrange matrices A and B (by replacing blocks with the appropriate source and destination). */
	MPI_Sendrecv_replace(&(A[0][0]), PROC_SIZE*PROC_SIZE, MPI_TYPE, shift_coords[0][1], 0, shift_coords[0][0], 0, mesh_comm, &status);
	MPI_Sendrecv_replace(&(B[0][0]), PROC_SIZE*PROC_SIZE, MPI_TYPE, shift_coords[1][1], 0, shift_coords[1][0], 0, mesh_comm, &status);

	if (rank == ROOT_PROC)
	{
		printf("Matrices A and B rearranged for multiplication. \n");
		printf("Multiplying A and B... \n");
	}

	/* Multiply and rotate. */
	#pragma omp parallel for
	for (int m = 0; m < sqrt(num_procs); m++)
	{
		multiply_matrices(0, 0, 0, 0, 0, 0, PROC_SIZE, PROC_SIZE, PROC_SIZE, A, B, C);
		
		MPI_Cart_shift(mesh_comm, 1, 1, &shift_coords[0][1], &shift_coords[0][0]);               /* Blocks of Matrix A are shifted 1 place left in the process mesh. */
		MPI_Sendrecv_replace(&(A[0][0]), PROC_SIZE*PROC_SIZE, MPI_TYPE, shift_coords[0][1], 0, shift_coords[0][0], 0, mesh_comm, &status);

		MPI_Cart_shift(mesh_comm, 0, 1, &shift_coords[1][1], &shift_coords[1][0]);               /* Blocks of Matrix B are shifted 1 place above in the process mesh. */
		MPI_Sendrecv_replace(&(B[0][0]), PROC_SIZE*PROC_SIZE, MPI_TYPE, shift_coords[1][1], 0, shift_coords[1][0], 0, mesh_comm, &status);
	}

	MPI_Barrier(MPI_COMM_WORLD);                                                                     /* Wait for multiplication to finish across all processes. */

	if(PRINT_MATRICES && (rank == ROOT_PROC))
	{
		printf("Result Matrix: \n---------------------------\n");
		print_checkerboard_matrix((void **) C, MPI_TYPE, N, N, mesh_comm);
	}

	MPI_Barrier(MPI_COMM_WORLD);                                                                     /* Wait for printing to finish. */

	process_time += MPI_Wtime();                                                                     /* Find time taken by the process. */
	MPI_Reduce(&process_time, &total_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);               /* Use reduction to find the process which takes the maximum time. */
	if (rank == ROOT_PROC)                                                                           /* Print results. */
	{
		printf("\n-------------------------------------------------------------------------------------------------------------\n");
		printf("N = %d, Processes = %d, Time = %8.3f seconds.\n", N, num_procs, total_time);
		printf("Performance = %6.3f GFLOPS. \n", (2*N*N*N)/(1000000000.0*total_time));
		printf("-------------------------------------------------------------------------------------------------------------\n");
	}

	/* Deallocate matrices and MPI variables. */
	free(A);
	free(B);
	free(C);
	MPI_Finalize();
	
	return 0;
}

/* Dynamically allocate memory for 2D matrices of dimensions PROC_SIZE x PROC_SIZE. */
DTYPE** create_matrix()
{
	DTYPE* array = (DTYPE *)malloc(sizeof(DTYPE)*PROC_SIZE*PROC_SIZE);
	DTYPE** matrix = (DTYPE **)malloc(sizeof(DTYPE*)*PROC_SIZE);
	for (int i = 0; i < PROC_SIZE; i++)
		matrix[i] = &(array[PROC_SIZE*i]);

	return matrix;
}

/* Create a process mesh to handle the separate blocks of the entire matrix. */
void mesh_setup(int *mesh_size, int num_procs, int *periodic, MPI_Comm *mesh_comm, int *mesh_rank, int *mesh_coords)
{
	mesh_size[0] = mesh_size[1] = 0;                                        /* Initialize mesh_size array. */
	MPI_Dims_create(num_procs, 2, mesh_size);                               /* Create a balanced distribution of processes (nodes). */
	periodic[0] = periodic[1] = 1;                                          /* Set periodicity to 1 (the matrices wrap around the edges). */
	MPI_Cart_create(MPI_COMM_WORLD, 2, mesh_size, periodic, 0, mesh_comm);  /* Create a process mesh and return a communicator for it (mesh_comm). */
	MPI_Comm_rank(*mesh_comm, mesh_rank);                                   /* Determine the ranks within the process mesh. */
	MPI_Cart_coords(*mesh_comm, *mesh_rank, 2, mesh_coords);                /* Determine the process coordinates using the rank in the process mesh. */
	
	MPI_Barrier(MPI_COMM_WORLD);                                            /* Wait until all process reach this routine. */
}

/* Multiply Block Matrices.*/
void multiply_matrices(
	int c_row, int c_col,                                                                   /* Corner of Block C. */
	int a_row, int a_col,                                                                   /* Corner of Block A. */
	int b_row, int b_col,                                                                   /* Corner of Block B. */
	int l, int m, int n,                                                                    /* Dimensions of matrices for division. */
	DTYPE** A, DTYPE** B, DTYPE** C								/* Matrices A, B, and C. */
)
{
	int a_quad[3], b_quad[3], c_quad[3];                                                    /* Quadrants for matrices A, B, and C. */
	DTYPE* aptr, * bptr, * cptr;                                                            /* Variables for adding. */

	if (m * n > THRESHOLD)                                                                  /* If the matrix doesn't fit in the dimension divide recursively.*/
	{
		a_quad[0] = 0; a_quad[1] = l / 2; a_quad[2] = l - l / 2;                        /* Determine quadrants half the current size. */
		b_quad[0] = 0; b_quad[1] = m / 2; b_quad[2] = m - m / 2;
		c_quad[0] = 0; c_quad[1] = n / 2; c_quad[2] = n - n / 2;

		for (int i = 0; i < 2; i++)                                                     /* Recursively divide the matrix. */
			for (int j = 0; j < 2; j++)
				for (int k = 0; k < 2; k++)
					multiply_matrices(
						c_row + a_quad[i], c_col + c_quad[j],
						a_row + a_quad[i], a_col + b_quad[k],
						b_row + b_quad[k], b_col + c_quad[j],
						a_quad[i + 1], b_quad[k + 1], c_quad[j + 1],
						A, B, C
					);
	}
	else                                                                                    /* Otherwise, multiply the blocks. */
		#pragma omp parallel for private(j, cptr, apt, bptr)
		for (int i = 0; i < l; i++)
			for (int j = 0; j < n; j++)
			{
				cptr = &C[c_row + i][c_col + j];
				aptr = &A[a_row + i][a_col];
				bptr = &B[b_row][b_col + j];
				for (int k = 0; k < m; k++)
				{
					*cptr += *(aptr++) * *bptr;                             /* Add the address in C. */
					bptr += PROC_SIZE;                                      /* Move to the next block for B. */
				}
			}

}

void print_subvector(
	void        *a,       /* IN - Array pointer */
	MPI_Datatype dtype,   /* IN - Array type */
	int          n)       /* IN - Array size */
{
	int i;

	for (i = 0; i < n; i++) {
		if (dtype == MPI_DOUBLE)
			printf("%6.3f ", ((double *)a)[i]);
		else {
			if (dtype == MPI_FLOAT)
				printf("%6.3f ", ((float *)a)[i]);
			else if (dtype == MPI_INT)
				printf("%6d ", ((int *)a)[i]);
		}
	}
}

void print_checkerboard_matrix(
	void       **a,            /* IN -2D matrix */
	MPI_Datatype dtype,        /* IN -Matrix element type */
	int          m,            /* IN -Matrix rows */
	int          n,            /* IN -Matrix columns */
	MPI_Comm     grid_comm)    /* IN - Communicator */
{
	void      *buffer;         /* Room to hold 1 matrix row */
	int        coords[2];      /* Grid coords of process
							   sending elements */
	int        datum_size;     /* Bytes per matrix element */
	int        els;            /* Elements received */
	int        grid_coords[2]; /* Coords of this process */
	int        grid_id;        /* Process rank in grid */
	int        grid_period[2]; /* Wraparound */
	int        grid_size[2];   /* Dims of process grid */
	int        i, j, k;
	void      *laddr;          /* Where to put subrow */
	int        local_cols;     /* Matrix cols on this proc */
	int        p;              /* Number of processes */
	int        src;            /* ID of proc with subrow */
	MPI_Status status;         /* Result of receive */

	MPI_Comm_rank(grid_comm, &grid_id);
	MPI_Comm_size(grid_comm, &p);
	datum_size = sizeof(int);

	MPI_Cart_get(grid_comm, 2, grid_size, grid_period,
		grid_coords);
	local_cols = BLOCK_SIZE(grid_coords[1], grid_size[1], n);

	buffer = malloc((size_t)n*datum_size);

	/* For each row of the process grid */
	for (i = 0; i < grid_size[0]; i++) {
		coords[0] = i;

		/* For each matrix row controlled by the process row */
		for (j = 0; j < BLOCK_SIZE(i, grid_size[0], m); j++) {

			/* Collect the matrix row on grid process 0 and
			print it */
			if (!grid_id) {
				for (k = 0; k < grid_size[1]; k++) {
					coords[1] = k;
					MPI_Cart_rank(grid_comm, coords, &src);
					els = BLOCK_SIZE(k, grid_size[1], n);
					laddr = buffer + (BLOCK_LOW(k, grid_size[1], n) * datum_size);
					if (src == 0) {
						memcpy(laddr, a[j], els * datum_size);
					}
					else {
						MPI_Recv(laddr, els, dtype, src, 0,
							grid_comm, &status);
					}
				}
				print_subvector(buffer, dtype, n);
				putchar('\n');
			}
			else if (grid_coords[0] == i) {
				MPI_Send(a[j], local_cols, dtype, 0, 0, grid_comm);
			}
		}
	}
	if (!grid_id) {
		free(buffer);
		putchar('\n');
	}
}
