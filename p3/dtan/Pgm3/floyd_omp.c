/**************
*  File:  floyd_omp.c   program3
*  Name: David Tan
*
*/
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h>
int main(int argc, char* argv[]) {
	int  n = 10;
	double start,stop;
	int dist[n][n];
	int i,j;
	srand(time(NULL));
	#pragma omp parallel for private(i,j)
	for (i = 0; i < n; i++){
		for (j = 0; j < n; j++){
			dist[i][j] = 99999;
			if( (i-j==1) || (j-i==1)) dist[i][j]=1;
			if( (i==0)&&(j==(n-1))) dist[i][j]=1;
			if( (j==0)&&(i==(n-1))) dist[i][j]=1; 
			if(i == j) dist[i][j] = 0;
		}
	}
        

	
	//Floyd
	int k, temp;
	start = omp_get_wtime();
	#pragma omp parallel for private(i,j,k)  schedule(static,16) num_threads(28)
        for (k = 0; k < n; k++) {
	    //#pragma omp parallel for private(i,j)
            for (i = 0; i < n; i++){
                for (j = 0; j < n; j++) {
                    //temp = dist[i][k] + dist[k][j];
                    //if (temp < dist[i][j])
                    //    dist[i][j] = temp;
                    dist[i][j] = fmin(dist[i][k] + dist[k][j], dist[i][j]);
                }
	    }   
        }
	stop = omp_get_wtime();

	//printf("The solution is:\n");
	//for (i = 0; i < n; i++) {
        //    for (j = 0; j < n; j++)
        //        printf("%d ", dist[i][j]);
        //    printf("\n");
        //}
	
	printf("Total time was %f seconds.\n", stop-start);
	double t = stop-start;
	double flop = 2.0*n*n*n/(t*1000000000.0);
	printf("The program's GFlops:  %f\n",flop);
	

    return 0;
}  
	
	
	
