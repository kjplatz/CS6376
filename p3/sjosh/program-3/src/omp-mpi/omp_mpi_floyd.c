/*
 *   Floyd's all-pairs shortest path - MPI and OpenMP implementation
 *
 *   Shriroop Joshi (sxj156530)
 *   Modified program by Michael J. Quinn
 */

#include <stdio.h>
#include <mpi.h>
#include "my-mpi.h"

typedef int dtype;
#define MPI_TYPE MPI_INT

#define DIM 2000

void compute_shortest_path(int, int, int **, int);

int main (int argc, char *argv[]) {
   dtype** a;         /* Doubly-subscripted array */
   dtype*  storage;   /* Local portion of array elements */
   int     i, j, k;
   int     id;        /* Process rank */
   int     m;         /* Rows in matrix */
   int     n;         /* Columns in matrix */
   int     p;         /* Number of processes */
   double  time, max_time;

   void compute_shortest_paths (int, int, int**, int);

   MPI_Init (&argc, &argv);
   MPI_Comm_rank (MPI_COMM_WORLD, &id);
   MPI_Comm_size (MPI_COMM_WORLD, &p);

    m = n = DIM;

    // Generation of matrix
   storage = (dtype *) malloc(m * n * sizeof(dtype));
   for(i = 0; i < m; i++) {
	for(j = 0; j < n; j++) {
	    if(i == j)
		storage[i * m + j] = 0;
	    else
		storage[i * m + j] = i + j;
	}
   }
   a = (dtype **) malloc(m * sizeof(int));
   for(i = 0; i < m; i++)
	a[i] = &storage[i * n];

   if (m != n) terminate (id, "Matrix must be square\n");

/*
   print_row_striped_matrix ((void **) a, MPI_TYPE, m, n,
      MPI_COMM_WORLD);
*/
   MPI_Barrier (MPI_COMM_WORLD);
   time = -MPI_Wtime();
   compute_shortest_paths (id, p, (dtype **) a, n);
   time += MPI_Wtime();
   MPI_Reduce (&time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0,
      MPI_COMM_WORLD);
   if (!id) printf ("Floyd, matrix size %d, %d processes: %6.2f seconds\n",
      n, p, max_time);
/*
   print_row_striped_matrix ((void **) a, MPI_TYPE, m, n,
      MPI_COMM_WORLD);
*/
   MPI_Finalize();
}


/*
 * Algorithm - uses Column matrix implementation
 */
void compute_shortest_paths (int id, int p, dtype **a, int n)
{
   int  i, j, k;
   int  offset;   /* Local index of broadcast row */
   int  root;     /* Process controlling row to be bcast */
   int* tmp;      /* Holds the broadcast row */

   tmp = (dtype *) malloc (n * sizeof(dtype));
   for (k = 0; k < n; k++) {
      root = BLOCK_OWNER(k,p,n);
      if (root == id) {
         offset = k - BLOCK_LOW(id,p,n);
	#pragma omp parallel for
         for (j = 0; j < n; j++)
            tmp[j] = a[offset][j];
      }
      MPI_Bcast (tmp, n, MPI_TYPE, root, MPI_COMM_WORLD);
	#pragma omp parallel for private(j, tmp)
      for (i = 0; i < BLOCK_SIZE(id,p,n); i++)
         for (j = 0; j < n; j++)
            a[i][j] = MIN(a[i][j],a[i][k]+tmp[j]);
   }
   free (tmp);
}
