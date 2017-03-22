//Author: Anay Sonawane
//File Name: floyd_mpi.c

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include "mpi.h"

#define MIN(a,b) ((a)<(b)?(a):(b))
#define BLOCK_LOW(id,p,n)  ((id)*(n)/(p))
#define BLOCK_HIGH(id,p,n) (BLOCK_LOW((id)+1,p,n)-1)
#define BLOCK_SIZE(id,p,n) \
                     (BLOCK_HIGH(id,p,n)-BLOCK_LOW(id,p,n)+1)
#define BLOCK_OWNER(j,p,n) (((p)*((j)+1)-1)/(n))


void compute_shortest_paths(int,int,int*,int);

int main(int argc, char *argv[])
{

int n, i, j;
int *arr;
double elapsed_time, total_time, perf;
int id,p;
int *local_arr;
int* temp_arr;
FILE *fp;

MPI_Init(&argc, &argv);
MPI_Comm_rank(MPI_COMM_WORLD, &id);
MPI_Comm_size(MPI_COMM_WORLD, &p);

if(id==0)
{
	//printf("please enter vertices:");
	//scanf("%d",&n);
     tasks = omp_get_max_threads();
	//get the matrix dimensions
	n = atoi(argv[1]);
	arr = (int*)malloc(n*n*sizeof(int));
    fp = fopen("mpi_job.csv","a+");
	//FILL THE ADJACENCY MATRIX
	for (i=0;i<n;i++){
	
		for(j=0;j<n;j++){
			if(i==j){
				arr[i*n + j]=0;
			}
			else if((i==j-1) || (j==i-1)){
				arr[i*n + j] = 1;
			}
			else{
				arr[i*n + j]=n;
			}			
		}		
	}
       //broadcast value n to all the nodes      
       MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);	
}
else{
MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
}

/*if(id == 0){

	printf("input:\n");

	for(i=0;i<n;i++)
	{		
		for(j=0;j<n;j++)
		{
			printf("%d ", arr[i*n+j]);
		}
		printf("\n");
	}
}*/

/* allocate buffer for the local rows*/
local_arr = malloc(n * (n/p) * sizeof(int));
/* allocate buffer for the temp array*/
temp_arr = malloc(n * n * sizeof(int));
/* Divide and send the matrix to other processes*/
MPI_Scatter(arr, n * (n/p), MPI_INT,local_arr, n * (n/p), MPI_INT, 0, MPI_COMM_WORLD);

MPI_Barrier(MPI_COMM_WORLD);

elapsed_time = -MPI_Wtime();

//floyds' algorithm
compute_shortest_paths(id,p,local_arr,n);

elapsed_time += MPI_Wtime(); 

//to calculate the maximum time required for the processing at any node
MPI_Reduce (&elapsed_time, &total_time, 1, MPI_DOUBLE, MPI_MAX, 0,
      MPI_COMM_WORLD);

if (id==0){
         perf = 2*pow(n,3)/(total_time*pow(10,9));
	 printf ("Processes=%d Tasks_per_process=%d Total_time=%10.6f Peformance=%lf\n",
			p, tasks, total_time, perf);
         fprintf(fp,"%d,%d,%lf,%lf\n", p, tasks,  total_time, perf);
}

//collect output from every node
MPI_Gather(local_arr, n * (n/p), MPI_INT, temp_arr,
            n * (n/p), MPI_INT, 0, MPI_COMM_WORLD);

//display the output
/*if(id == 0){
	printf("output: \n");
	for(i=0;i<n;i++)
	{
		for(j=0;j<n;j++)
		{	
			printf("%d ", temp_arr[i*n+j]);
		}
		printf("\n");
	}
}*/
MPI_Finalize();

}


void compute_shortest_paths (int id, int p, int a[], int n)
{
   int		i, j, k;
   int		offset;
   int		root;		
   int		*tmp;	

   tmp = (int *) malloc (n * sizeof(int));

   for (k = 0; k < n; k++)
   {  root = BLOCK_OWNER(k, p, n);
      if (root == id)
      {  offset = k - BLOCK_LOW(id, p, n);
         for (j = 0; j < n; j++)
            tmp[j] = a[offset*n+j];
      }

      MPI_Bcast (tmp, n, MPI_INT, root, MPI_COMM_WORLD);

      for (i = 0; i < n/p; i++)
         for (j = 0; j < n; j++)
            a[i*n+j] = MIN(a[i*n + j], a[i*n+k] + tmp[j]);
   }
   free (tmp);
}
