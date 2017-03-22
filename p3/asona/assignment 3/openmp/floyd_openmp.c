//Author: ANAY Sonawane
//File Name: floyd_openmp.c

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <omp.h>

#define MIN(a,b) ((a)<(b)?(a):(b))

void compute_shortest_paths(int*,int);

int main(int argc, char *argv[])
{

int i, j, n, threads;
int *arr;
double total_time,perf;
struct timeval start_time, stop_time, elapsed_time;

FILE *fp=fopen("floyd_openmp.csv","a+");

//Get the dimensions/vertices of graph
printf("Enter the dimensions of matrix: \n");
scanf("%d",&n);

threads = omp_get_max_threads();

arr = (int*)malloc(n*n*sizeof(int));

for (i=0;i<n;i++){

	for(j=0;j<n;j++){
		if(i==j){
			arr[i*n + j]=0;
		}
		else if((i==j-1) || (j==i -1)){
			arr[i*n + j] = 1;
		}
		else{
			arr[i*n + j]=n;
		}
	}
}


//Print the array
/*printf("input:\n");
for(i=0;i<n;i++)
{
	for(j=0;j<n;j++)
	{
		printf("%d ", arr[i*n+j]);
	}
	printf("\n");
}*/

gettimeofday(&start_time, NULL);

compute_shortest_paths(arr,n);

gettimeofday(&stop_time, NULL);

timersub(&stop_time, &start_time, &elapsed_time);

total_time =(double)( elapsed_time.tv_sec+elapsed_time.tv_usec/1000000.0);
perf = 2*pow(n,3)/(total_time*pow(10,9));
printf("Threads=%d Total_time=%lf Performance=%lf \n",threads,total_time,perf);
fprintf(fp,"%d,%lf,%lf\n",threads,total_time,perf);
fclose(fp);
/*printf("output: \n");
for(i=0;i<n;i++)
{
	for(j=0;j<n;j++)
	{
		printf("%d ", arr[i*n+j]);
	}
	printf("\n");
}
*/


}


void compute_shortest_paths (int *a, int n)
{
   int	i, j, k, dummy;
  
   #pragma omp parallel for private(k)
   for (k = 0; k < n; k++)
   {  
      #pragma omp parallel for private(i,j) 
      for (i = 0; i < n; i++)
         for (j = 0; j < n; j++)
         	a[i*n + j]= MIN(a[i*n + j], a[i*n + k] + a[k*n + j]);
   }
}
