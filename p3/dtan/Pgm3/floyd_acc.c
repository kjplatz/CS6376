/**************
*  File:  floyd_acc.c   program3
*  Name: David Tan
*
 */
#include <stdio.h>
#include <stdlib.h>
//#include <omp.h>
#include <math.h>
#include <sys/time.h>
int main(int argc, char* argv[]) {
	int  n = 5000;
	double start,stop;
	int dist[n][n];
	int i,j;
	struct timeval start_time, stop_time, elapsed_time;  // timers
	srand(time(NULL));
//	#pragma acc kernels
	for (i = 0; i < n; i++){
		for (j = 0; j < n; j++){
			dist[i][j] = 99999;
			if( (i-j==1) || (j-i==1)) dist[i][j]=1;
			if( (i==0)&&(j==(n-1))) dist[i][j]=1;
			if( (j==0)&&(i==(n-1))) dist[i][j]=1; 
			if(i == j) dist[i][j] = 0;
		}
	}
	
	//for (i = 0; i < n; i++) {
        //    for (j = 0; j < n; j++)
        //        printf("%d ", dist[i][j]);
        //    printf("\n");
       // }
	//Floyd
	int k, temp;
	//start = omp_get_wtime();
	gettimeofday(&start_time,NULL); // Unix timer
	//#pragma acc kernels
	#pragma acc data copy(dist)
	{
		#pragma acc parallel private(dist)
	        	#pragma acc loop 
           		for (k = 0; k < n; k++) {
	        	//#pragma omp parallel for private(i,j)
	        		#pragma acc loop
                		for (i = 0; i < n; i++){
	            			//#pragma acc loop vector(128)
                    			for (j = 0; j < n; j++) {
                        //temp = dist[i][k] + dist[k][j];
                        //if (temp < dist[i][j])
                        //    dist[i][j] = temp;
                        			dist[i][j] = fmin(dist[i][k] + dist[k][j], dist[i][j]);
                   			 }
	        		}   
            		}
	}
	
	gettimeofday(&stop_time,NULL);
	timersub(&stop_time, &start_time, &elapsed_time); // Unix time subtract routine
	
	//printf("The solution is:\n");
	//for (i = 0; i < n; i++) {
        //    for (j = 0; j < n; j++)
        //        printf("%d ", dist[i][j]);
        //    printf("\n");
        //}
	
	//stop = omp_get_wtime();
	double t = elapsed_time.tv_sec+elapsed_time.tv_usec/1000000.0;
	printf("Total time was %f seconds.\n", t);
	//double t = stop-start;
	double flop = 2.0*n*n*n/(t*1000000000.0);
	printf("The program's GFlops:  %f\n",flop);
	

    return 0;
}  
	
	
	
