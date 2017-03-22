//Name: Anay Sonawane
//Author: floyd_openacc.c

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>

#define MIN(a,b) ((a)<(b)?(a):(b))

void compute_shortest_paths(int*,int);

int main(int argc, char *argv[])
{

int n, i, j;
int* arr;
struct timeval start_time, stop_time, elapsed_time;
double total_time, perf;

printf("Enter the dimensions of matrix: \n");
scanf("%d",&n);

arr = (int*)malloc(n*n*sizeof(int));

for (i=0;i<n;i++){

	for(j=0;j<n;j++){
		if(i==j){
			arr[i*n+j]=0;
		}
		else if((i==j-1) || (j==i-1)){
			arr[i*n+j]=1;
		}
		else{
			arr[i*n+j]=n;
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

printf("Total_Time=%f Performance=%lf\n",total_time,perf);

printf("output: \n");
for(i=0;i<n;i++)
{
	for(j=0;j<n;j++)
	{
		printf("%d ", arr[i*n+j]);
	}
	printf("\n");
}


}


void compute_shortest_paths (int* a, int n)
{
   int	i, j, k;
  #pragma acc data copy(a[0:n*n])
  {
   for (k = 0; k < n; k++)
   {
      #pragma acc parallel
      {
      #pragma acc loop independent
      for (i = 0; i < n; i++){
         #pragma acc loop independent
         for (j = 0; j < n; j++)
         	a[i*n+j]= MIN(a[i*n+j], a[i*n+k] + a[k*n+j]);
      } 
     }
   }
 }	
}
