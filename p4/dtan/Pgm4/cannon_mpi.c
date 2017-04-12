/******************************************************************************
* FILE: mm_mpi.c  
*   Cannon MPI matrix multiplication
* Name: David
* Program 4
******************************************************************************/


#include<stdio.h>
#include<mpi.h>
#include<math.h>
#include<stdlib.h>
#include<time.h>

#define N 1600

void cannon(int n, float *a, float *b, float *c, float *new_c, MPI_Comm comm)
{
	int i;
	int nlocal;
	int npes, dims[2], periods[2];
	int myrank, my2drank, mycoords[2];
	int uprank, downrank, leftrank, rightrank, coords[2];
	int shiftsource, shiftdest;
	MPI_Status status;
	MPI_Comm comm_2d;

	/* communicator info */
	MPI_Comm_size(comm, &npes);
	MPI_Comm_rank(comm, &myrank);

	/* prepare Cartesian */
	dims[0] = dims[1] = sqrt(npes); /* 4 4 */

	/* set periods of wraparound connections */
	periods[0] = periods[1] = 1;

	/* set up the Cartesian topology, with rank reordering */
	MPI_Cart_create(comm, 2, dims, periods, 1, &comm_2d);
	/* Get the rank and coordinates with respect to the new topology */
	MPI_Comm_rank(comm_2d, &my2drank);
	MPI_Cart_coords(comm_2d, my2drank, 2, mycoords);

	/* Compute ranks of the up and left shifts */
	MPI_Cart_shift(comm_2d, 1, -1, &rightrank, &leftrank);
	MPI_Cart_shift(comm_2d, 0, -1, &downrank, &uprank);

	/* Determine the dimension of the local matrix block */
	nlocal = n/dims[0];

	/* Perform the initial matrix alignment. First for A and then for B */
	MPI_Cart_shift(comm_2d, 1, -mycoords[0], &shiftsource, &shiftdest);
	MPI_Sendrecv_replace(a, nlocal*nlocal, MPI_FLOAT, shiftdest, 1, shiftsource, 1, comm_2d, &status);
	MPI_Cart_shift(comm_2d, 0, -mycoords[1], &shiftsource, &shiftdest);
	MPI_Sendrecv_replace(b, nlocal*nlocal, MPI_FLOAT, shiftdest, 1, shiftsource, 1, comm_2d, &status);

	/* main computation loop */
	for (i=0; i<dims[0]; i++)
	{
		int x, y, k;
		for (x=0; x<nlocal; x++)
			for (y=0; y<nlocal; y++)
				for (k=0; k<nlocal; k++)
					c[x*nlocal+y] += a[x*nlocal+k]*b[k*nlocal+y];
			
		/* Shift matrix a left by one */
		MPI_Sendrecv_replace(a, nlocal*nlocal, MPI_FLOAT, leftrank, 1, rightrank, 1, comm_2d, &status);
		/* Shift matrix b up by one */
		MPI_Sendrecv_replace(b, nlocal*nlocal, MPI_FLOAT, uprank, 1, downrank, 1, comm_2d, &status);
	}

	/* Change back to the original distribution of a and b */
	MPI_Cart_shift(comm_2d, 1, +mycoords[0], &shiftsource, &shiftdest);
	MPI_Sendrecv_replace(a, nlocal*nlocal, MPI_FLOAT, shiftdest, 1, shiftsource, 1, comm_2d, &status);
	MPI_Cart_shift(comm_2d, 0, +mycoords[1], &shiftsource, &shiftdest);
	MPI_Sendrecv_replace(b, nlocal*nlocal, MPI_FLOAT, shiftdest, 1, shiftsource, 1, comm_2d, &status);
	MPI_Comm_free(&comm_2d); /* Free communicator */
    
    int idx;
	int new_n = n;
	int idx_col;
	for(idx = 0; idx < nlocal; idx++)
	{
		for(idx_col = 0; idx_col < nlocal; idx_col++)
		{
			int new_row = mycoords[0]*nlocal+idx;
			int new_col = mycoords[1]*nlocal+idx_col;
			new_c[new_row*new_n+new_col] = c[idx*nlocal+idx_col];
		}
	}

	if(myrank != 0){
		MPI_Reduce(new_c, new_c, new_n*new_n, MPI_FLOAT, MPI_SUM, 0, comm);
	}else{
		MPI_Reduce(MPI_IN_PLACE, new_c, new_n*new_n, MPI_FLOAT, MPI_SUM, 0, comm);
	}
} /* end of cannon */

int main(int argc, char *argv[])
{
	/* locate memory for A B C matrix */
    float *a = calloc(N * N, sizeof(float));
	float *b = calloc(N * N, sizeof(float));
	float *c = calloc(N * N, sizeof(float));
	
	int i, j;
	/* initialize matrix A and B */
	for(i = 0; i < N; i++){
		for(j = 0; j < N; j++){
			a[i*N + j] = j - i;
		}
    }
	for(i = 0; i < N; i++){
		for(j = 0; j < N; j++){
			b[i*N + j] = N - j + i;
		}
    }
	/*MPI initialize*/
    MPI_Init(&argc,&argv);
	int num_proc, myrank;     /*number of process,  proc id*/
	MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	
	int num_block = N/sqrt(num_proc); /* block size in a row(col) 
									ex. 16/sqrt(16) = 4   */
	int row_block;
	int col_block;
	int row;
	int col;

	int n_block = N/num_block;  /* ex. 16/4 = 4  number of blocks */

	float block_a[num_proc][num_block*num_block]; /* ex.  [16][4*4] */
	float block_b[num_proc][num_block*num_block];
	float block_c[num_proc][num_block*num_block];

	col_block = myrank % n_block;
	row_block = (myrank - col_block) / n_block;

	int index_block = 0;
	for(row = row_block * num_block; row < row_block * num_block + num_block; row++)
	{
		for(col = col_block * num_block; col < col_block * num_block + num_block; col++)
		{
			block_a[myrank][index_block] = a[row * N + col];
			block_b[myrank][index_block] = b[row * N + col];
			block_c[myrank][index_block] = 0;
			index_block++;
		}
	}

    MPI_Barrier(MPI_COMM_WORLD);
    double start = MPI_Wtime();
	/* start of cannon */
    cannon(N, block_a[myrank], block_b[myrank], block_c[myrank], &c[0], MPI_COMM_WORLD);
	/* end of cannon */
    MPI_Barrier(MPI_COMM_WORLD);
    double end = MPI_Wtime();
    
    if(myrank == 0){
        /* Serial multiplication */
        /* only for testing */
		/*
		float *d = calloc(N * N, sizeof(float));
		
		int k;
		for (i=0; i<N; i++)
			for (j=0; j<N; j++)
				for (k=0; k<N; k++)
					d[i*N+j] += a[i*N+k]*b[k*N+j];
		

        printf("\n\n");
        for (i=0; i<N; i++){
			for (j=0; j<N; j++)
				printf("%8.2f  ", c[i*N + j]);
			printf("\n");
		}
			
			
		/* compare result with serial mm */	
		/*
		int equal = 1;
		int i = 0;
		for(i = 0; i < N*N; i++){
			if(d[i] != c[i]){
				equal = 0;
			} 
		}
        if(equal){
            printf("\n Equal with serial mm! \n");
        }
		*/
        printf("\nTime: %.4f secounds\n", (end - start));
		double t = end - start;
		double flop = 2.0*N*N*N/(t*1000000000.0);
		printf("Gfliops: %6.2f\n\n", flop);
    }

	MPI_Finalize();
	return 0;
}

