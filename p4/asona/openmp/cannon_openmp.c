/*Author:Anay Sonawane
  File Name: cannon_openmp.c
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <omp.h>


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

void shift_matrix_left(int* A, int N, int block_size, int flag)
{	
int i,j,k,l,tmp,shift=0;

//#pragma omp parallel for private(i,j,k)
for(i=0;i<N;i+=block_size){
	if(flag==1)shift=block_size;
	for(l=i;l<i+block_size;l++){
		k=0;
		while(k<shift){
			for(j=0;j<N;j++){
				if(j==0)tmp=A[l*N];
				if(j==N-1)A[l*N+j]=tmp;
				else A[l*N+j]=A[l*N+j+1];
			}
		k++;

		}		
	}
	if(flag==0)shift = shift + block_size;
}
}

void shift_matrix_up(int* A, int N, int block_size, int flag)
{	
int i,j,k,l,tmp,shift=0;

//#pragma omp parallel for private(i,j,k)
for(j=0;j<N;j+=block_size){
	if(flag==1)shift=block_size;
	for(l=j;l<j+block_size;l++){
		k=0;
		while(k<shift){
			for(i=0;i<N;i++){
				if(i==0)tmp=A[i*N+l];
				if(i==N-1)A[i*N+l]=tmp;
				else A[i*N+l]=A[(i+1)*N+l];
			}
			k++;
		}
	}
	if(flag==0)shift = shift + block_size;   	
		
}
}

void matrix_mult(int*A, int*B, int*C, int N, int block_size, int p){

int id,i,j,k,i_block, j_block, block_row, block_column;
int *A_block, *B_block, *C_block;

int threads= p*p;

#pragma omp parallel default(none) private(id,i,j,k,i_block,j_block,block_row,block_column,A_block,B_block,C_block) shared(A,B,C,N,p, block_size) num_threads(threads) 
	{
	id = omp_get_thread_num();
//	printf("%d\n",id);
	block_row = (id/p) * block_size;
	block_column = (id%p) * block_size;

	A_block = (int*)malloc(block_size*block_size*sizeof(int));
	B_block = (int*)malloc(block_size*block_size*sizeof(int));
	C_block = (int*)malloc(block_size*block_size*sizeof(int));

	//Allocate blocks according threads number
	for(i=block_row,i_block=0;i<block_row+block_size;i++,i_block++){
		for(j=block_column, j_block=0;j<block_column+block_size;j++,j_block++){
			A_block[i_block*block_size+j_block]=A[i*N+j];
			B_block[i_block*block_size+j_block]=B[i*N+j];
 			C_block[i_block*block_size+j_block]=C[i*N+j];

		}

	}
	
	//block multiply
	for(i=0;i<block_size;i++){
		for(j=0;j<block_size;j++){
			for(k=0;k<block_size;k++){
				C_block[i*block_size+j] += A_block[i*block_size+k]*B_block[k*block_size+j]; 
			}
		}

	}
	
	//Save output back to matrix C
	for(i=block_row,i_block=0;i<block_row+block_size;i++,i_block++){
		for(j=block_column, j_block=0;j<block_column+block_size;j++,j_block++){
 			C[i*N+j]=C_block[i_block*block_size+j_block];

		}

	}

	
	}
}


void main(int argc, char *argv[])
{

int *A, *B, *C;
int N,i,threads,p;
struct timeval start_time, stop_time, elapsed_time;
double total_time,perf;
int block_size;

FILE *fp=fopen("openmp.csv","a+");

N = atoi(argv[1]);
threads = omp_get_max_threads();
p=sqrt(threads);
block_size = N/p;

A = (int*)malloc(N*N*sizeof(int));
B = (int*)malloc(N*N*sizeof(int));
C = (int*)calloc(N*N, sizeof(int));

fill_matrix(A,B,N);
//print_matrix(A,N);
//print_matrix(B,N);

//do the first alignment of matrices
shift_matrix_left(A,N,block_size,0);
shift_matrix_up(B,N,block_size,0);
//print_matrix(A,N);

gettimeofday(&start_time, NULL);

for(i=0;i<p;i++)
{
	//multiplication
	matrix_mult(A,B,C,N,block_size,p);
	//left shift Matrix A
	shift_matrix_left(A,N,block_size,1);
	//Right shift matrix B
	shift_matrix_up(B,N,block_size,1);
}


gettimeofday(&stop_time, NULL);

timersub(&stop_time, &start_time, &elapsed_time);

total_time =(double)( elapsed_time.tv_sec+elapsed_time.tv_usec/1000000.0);

perf = 2*pow(N,3)/(total_time*pow(10,9));
printf("Threads=%d Total_time=%lf Performance=%lf \n",threads,total_time,perf);
fprintf(fp,"%d,%d,%lf,%lf\n",threads,N,total_time,perf);

fclose(fp);
//print_matrix(C,N);

}

