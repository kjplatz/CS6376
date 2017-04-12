/*Author:Anay Sonawane
File Name: cannon_openmpi.c
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <mpi.h>

#define BLOCK_SIZE N/p

void print_matrix(int *A, int N){

int i,j;

for(i=0;i<N;i++){
	for(j=0; j<N;j++)
	{
		printf("%d ",A[i*N+j]);
	}

	printf("\n");
}

}

//initially fill the matrix
void fill_matrix(int* A, int*B, int N){

int i,j;

for(i=0;i<N;i++){
	for(j=0; j<N;j++)
	{
		A[i*N+j]= j-i;
		B[i*N+j]= N-j+i;
	}
}

}

//Matrix multiplication
void matrix_mult(int n, int *a, int *b, int *c)
{
	int i, j, k;
	#pragma omp parallel for private(j,k)
	for (i=0; i<n; i++){
		for (j=0; j<n; j++)
			for (k=0; k<n; k++)
				c[i*n+j] += a[i*n+k]*b[k*n+j];
	}
}

void main(int argc, char *argv[])
{

int i,j,k, block_row, block_column;
int *A, *B, *C;
int N,threads,p;

int dims[2],periods[2];
int id, my2drank, mycoords[2];
int uprank, downrank, leftrank, rightrank;
int shiftsource, shiftdest;
double elapsed_time, total_time, perf;
FILE *fp;

N = atoi(argv[1]);

A = (int*)malloc(N*N*sizeof(int));
B = (int*)malloc(N*N*sizeof(int));
C = (int*)calloc(N*N, sizeof(int));

fill_matrix(A,B,N);

MPI_Init(&argc, &argv);
MPI_Comm_size(MPI_COMM_WORLD, &threads);
p=sqrt(threads);

int A_Block[threads][BLOCK_SIZE*BLOCK_SIZE];
int B_Block[threads][BLOCK_SIZE*BLOCK_SIZE];
int C_Block[threads][BLOCK_SIZE*BLOCK_SIZE];

MPI_Comm_rank(MPI_COMM_WORLD, &id);

if(id==0)fp = fopen("mpi.csv","a+");

//Calculate the co ordinates of the sub matrices acccording to the rank of MPI thread
block_row = (id/p) * BLOCK_SIZE;
block_column = (id%p) * BLOCK_SIZE;

k=0;

//Divide the matrices in to submatrices of size BLOCK_SIZE
for(i=block_row;i<block_row+BLOCK_SIZE;i++){
	for(j=block_column;j<block_column+BLOCK_SIZE;j++){
		A_Block[id][k]=A[i*N+j];
		B_Block[id][k]=B[i*N+j];
 		C_Block[id][k]=0;
		k++;
	}

}
	
MPI_Barrier(MPI_COMM_WORLD);

elapsed_time = -MPI_Wtime();

MPI_Status status;
MPI_Comm comm_2d;

dims[0] = dims[1] = p;
periods[0] = periods[1] = 1;

//Create the cartesian topolgy
MPI_Cart_create(MPI_COMM_WORLD,2,dims,periods,1,&comm_2d);
MPI_Comm_rank(comm_2d, &my2drank);
MPI_Cart_coords(comm_2d, my2drank, 2, mycoords);
MPI_Cart_shift(comm_2d, 1, -1, &rightrank, &leftrank);
MPI_Cart_shift(comm_2d, 0, -1, &downrank, &uprank);

int nlocal = N/dims[0];

//initial allignement of matrix
MPI_Cart_shift(comm_2d, 1, -mycoords[0], &shiftsource,&shiftdest);
MPI_Sendrecv_replace(A_Block[id], nlocal*nlocal, MPI_INT, shiftdest, 1,shiftsource, 1, comm_2d, &status);
MPI_Cart_shift(comm_2d, 0, -mycoords[1], &shiftsource,&shiftdest);
MPI_Sendrecv_replace(B_Block[id], nlocal*nlocal, MPI_INT, shiftdest, 1,shiftsource, 1, comm_2d, &status);

for(i=0; i<dims[0]; i++){
	//product of matrix
	matrix_mult(nlocal,A_Block[id],B_Block[id],C_Block[id]); 
	//Left shift the matrxi A by nlocal size
	MPI_Sendrecv_replace(A_Block[id], nlocal*nlocal, MPI_INT, leftrank,1, rightrank, 1, comm_2d, &status);
	//Up shift the matrix B by nlocal size
	MPI_Sendrecv_replace(B_Block[id], nlocal*nlocal, MPI_INT, uprank, 1, downrank, 1, comm_2d, &status);
}

//Get the original distribution of matrices
MPI_Cart_shift(comm_2d, 1, +mycoords[0], &shiftsource,&shiftdest);
MPI_Sendrecv_replace(A_Block[id], nlocal*nlocal, MPI_INT, shiftdest, 1,shiftsource, 1, comm_2d, &status);
MPI_Cart_shift(comm_2d, 0, +mycoords[1], &shiftsource,&shiftdest);
MPI_Sendrecv_replace(B_Block[id], nlocal*nlocal, MPI_INT, shiftdest, 1,shiftsource, 1, comm_2d, &status);
MPI_Comm_free(&comm_2d);


//Save the partial output in the matrix C
for(block_row = 0; block_row < nlocal; block_row++)
{
	for(block_column= 0; block_column < nlocal; block_column++)
	{
		i = mycoords[0]*nlocal+block_row;
		j = mycoords[1]*nlocal+block_column;
		C[i*N+j] = C_Block[id][block_row*nlocal+block_column];
	}
}

if(id != 0)
{
	MPI_Reduce(C, C, N*N, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
}
else
{
	MPI_Reduce(MPI_IN_PLACE, C, N*N, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
}

MPI_Barrier(MPI_COMM_WORLD);

elapsed_time += MPI_Wtime();

if(id==0){
	//print_matrix(C,N);
	perf = 2*pow(N,3)/(elapsed_time*pow(10,9));
	printf ("Processes=%d Total_time=%10.6f Peformance=%lf\n",
			threads, elapsed_time, perf);
        fprintf(fp,"%d,%d,%lf,%lf\n", threads,N,elapsed_time, perf);
}
MPI_Finalize();

}

