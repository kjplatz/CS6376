#include<stdio.h>
#define MIN(a,b)           ((a)<(b)?(a):(b))
#include <sys/time.h>

/****************************************************************************/
int main (int argc, char *argv[])
{
   int **restrict a ;/* Doubly-subscripted array */

   int *storage;/* Local portion of array elements */
   int i, j, k;
   int m;/* Rows in matrix */
   int n;/* Columns in matrix */
   struct timeval start_time, stop_time, elapsed_time;  // timers

   void read_matrix(char*, void***, void**, int*, int*);
   void print_matrix(int**, int, int);

   read_matrix (argv[1], (void *) &a, (void *) &storage, &m, &n);

   if (m != n)
     printf("Error: Matrix must be square\n");

   //print_matrix(a, m, n);

   gettimeofday(&start_time,NULL); // Unix timer
  
  for (k = 0; k < n; k++){
    #pragma acc parallel loop
    for (i = 0; i < n; i++){
      for (j = 0; j < n; j++){
             a[i][j] = MIN(a[i][j], a[i][k] + a[k][j]);
      }
    }
  }

   gettimeofday(&stop_time,NULL);
   timersub(&stop_time, &start_time, &elapsed_time); // Unix time subtract routine
   printf("Total time was %f seconds.\n", elapsed_time.tv_sec+elapsed_time.tv_usec/1000000.0);
   printf("Performance: %.2f\n", 2.00 * n * n * n /(elapsed_time.tv_sec + elapsed_time.tv_usec / 1000000.0) / 1000000000.0);

   //        print_matrix(a, m, n);
}

/*
 *   Process p-1 opens a file and inputs a two-dimensional
 *   matrix, reading and distributing blocks of rows to the
 *   other processes.
 */

void read_matrix (
		  char        *s,/* IN - File name */
		 int     ***subs,     /* OUT - 2D submatrix indices */
		 int       **storage,  /* OUT - Submatrix stored here */
		  int         *m,        /* OUT - Matrix rows */
		  int         *n         /* OUT - Matrix cols */
		  )
{
  int          i;
  FILE        *infileptr;    /* Input file pointer */
  void       **lptr;         /* Pointer into 'subs' */
  void        *rptr;         /* Pointer into 'storage' */
  
 int          x;            /* Result of read */
 int data[100000000];
 int* rows[10000];
 int          temp;

  infileptr = fopen (s, "r");
  if (infileptr == NULL){
    *m = 0;
  }
  else {
    fread (m, sizeof(int), 1, infileptr);
    fread (n, sizeof(int), 1, infileptr);
  }

  *storage = data;
  *subs = rows;

  temp =*m * *n * sizeof(int);
  x = fread (*storage, sizeof(int), *m * *n, infileptr);

  //construct 2D array
  lptr = *subs;
  rptr = (void *)*storage;

  for (i = 0; i < *m; i++) {
    *(lptr++)= (void *) rptr;
    rptr += *n * sizeof(int);
  }

  fclose (infileptr);
}

/*
 *   Print elements of a doubly-subscripted array.
 */

void print_matrix (
		   int       **a,       /* IN - Doubly-subscripted array */
		   int          rows,    /* IN - Matrix rows */
		   int          cols)    /* IN - Matrix cols */
{
  int i, j;

  for (i = 0; i < rows; i++) {
    for (j = 0; j < cols; j++) {
      printf ("%6d ", a[i][j]);
    }
    putchar ('\n');
  }
}

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/
