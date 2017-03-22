/* file - floyd_omp.c 
Author - Ayush Agrawal
Description - Implements Floyd alogorithm using OpenMP
*/

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#define n 2000 /* Then number of nodes in the graph */
int dist[n][n]; /* dist[i][j] is the length of the edge between i and j if
            it exists, or 0 if it does not */

/* Prints the adjacency matrix dist*/
void printDist() {
    int i, j;
    printf("    ");
    for (i = 0; i < n; ++i)
        printf("%4c", 'A' + i);
    printf("\n");
    for (i = 0; i < n; ++i) {
        printf("%4c", 'A' + i);
        for (j = 0; j < n; ++j)
            printf("%4d", dist[i][j]);
        printf("\n");
    }
    printf("\n");
}

/* Implements the floyds algorithm for shortest path on adjacency matrix dist*/ 
void floyd() {
    int i, j, k;
 
    for (k = 0; k < n; ++k)
        #pragma omp parallel for private(i,j)
        for (i = 0; i < n; ++i)
            for (j = 0; j < n; ++j)
                    if ((dist[i][k] + dist[k][j] < dist[i][j]) || (dist[i][j] == 0))// checks for the min distance
                        dist[i][j] = dist[i][k] + dist[k][j];
}
 
int main(int argc, char *argv[]) {
	struct timeval start_time, stop_time, elapsed_time;  // timers
    int i, j;
	/* initialize the dist matrix */
    for (i = 0; i < n; ++i)
        for (j = 0; j < n; ++j)
            if (i==j)
                dist[i][j]=0;
            else
                dist[i][j]= (int)( 3.0 * rand() / ( RAND_MAX + 1.0 ) ); // assigns a number between 0 - 2 for an edge
                                                                           
                                                                                 
    gettimeofday(&start_time,NULL); // Unix  start timer
    floyd();
    gettimeofday(&stop_time,NULL);// Unix end timer
    timersub(&stop_time, &start_time, &elapsed_time); // Unix time subtract routine
    printf("Total time was %f seconds.\n", elapsed_time.tv_sec+elapsed_time.tv_usec/1000000.0);                                                                                                         
    return 0;
}
