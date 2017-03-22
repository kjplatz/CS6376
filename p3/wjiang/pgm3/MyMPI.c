/*
 *   MyMPI.c -- A library of matrix/vector
 *   input/output/redistribution functions
 *
 *   Programmed by Michael J. Quinn
 *
 *   Last modification: 4 September 2002
 */

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "MyMPI.h"


/***************** MISCELLANEOUS FUNCTIONS *****************/

/*
 *   Given MPI_Datatype 't', function 'get_size' returns the
 *   size of a single datum of that data type.
 */

int get_size (MPI_Datatype t) {
  if (t == MPI_BYTE) return sizeof(char);
  if (t == MPI_DOUBLE) return sizeof(double);
  if (t == MPI_FLOAT) return sizeof(float);
  if (t == MPI_INT) return sizeof(int);
  printf ("Error: Unrecognized argument to 'get_size'\n");
  fflush (stdout);
  MPI_Abort (MPI_COMM_WORLD, TYPE_ERROR);
}


/*
 *   Function 'my_malloc' is called when a process wants
 *   to allocate some space from the heap. If the memory
 *   allocation fails, the process prints an error message
 *   and then aborts execution of the program.
 */

void *my_malloc (
		 int id,     /* IN - Process rank */
		 int bytes)  /* IN - Bytes to allocate */
{
  void *buffer;
  if ((buffer = malloc ((size_t) bytes)) == NULL) {
    printf ("Error: Malloc failed for process %d\n", id);
    fflush (stdout);
    MPI_Abort (MPI_COMM_WORLD, MALLOC_ERROR);
  }
  return buffer;
}


/*
 *   Function 'terminate' is called when the program should
 *   not continue execution, due to an error condition that
 *   all of the processes are aware of. Process 0 prints the
 *   error message passed as an argument to the function.
 *
 *   All processes must invoke this function together!
 */

void terminate (
		int   id,            /* IN - Process rank */
		char *error_message) /* IN - Message to print */
{
  if (!id) {
    printf ("Error: %s\n", error_message);
    fflush (stdout);
  }
  MPI_Finalize();
  exit (-1);
}

/*
 *   Process p-1 opens a file and inputs a two-dimensional
 *   matrix, reading and distributing blocks of rows to the
 *   other processes.
 */

void read_row_striped_matrix (
			      char        *s,        /* IN - File name */
			      void      ***subs,     /* OUT - 2D submatrix indices */
			      void       **storage,  /* OUT - Submatrix stored here */
			      MPI_Datatype dtype,    /* IN - Matrix element type */
			      int         *m,        /* OUT - Matrix rows */
			      int         *n,        /* OUT - Matrix cols */
			      MPI_Comm     comm)     /* IN - Communicator */
{
  int          datum_size;   /* Size of matrix element */
  int          i;
  int          id;           /* Process rank */
  FILE        *infileptr;    /* Input file pointer */
  int          local_rows;   /* Rows on this proc */
  void       **lptr;         /* Pointer into 'subs' */
  int          p;            /* Number of processes */
  void        *rptr;         /* Pointer into 'storage' */
  MPI_Status   status;       /* Result of receive */
  int          x;            /* Result of read */

  MPI_Comm_size (comm, &p);
  MPI_Comm_rank (comm, &id);
  datum_size = get_size (dtype);

  /* Process p-1 opens file, reads size of matrix,
     and broadcasts matrix dimensions to other procs */

  if (id == (p-1)) {
    infileptr = fopen (s, "r");
    if (infileptr == NULL) *m = 0;
    else {
      fread (m, sizeof(int), 1, infileptr);
      fread (n, sizeof(int), 1, infileptr);
    }
  }
  MPI_Bcast (m, 1, MPI_INT, p-1, comm);

  if (!(*m)) MPI_Abort (MPI_COMM_WORLD, OPEN_FILE_ERROR);

  MPI_Bcast (n, 1, MPI_INT, p-1, comm);

  local_rows = BLOCK_SIZE(id,p,*m);

  /* Dynamically allocate matrix. Allow double subscripting
     through 'a'. */

  *storage = (void *) my_malloc (id,
				 local_rows * *n * datum_size);
  //printf("size: %d\n", local_rows * *n * datum_size);
  //printf("*storage: %d\n", *storage);
  *subs = (void **) my_malloc (id, local_rows * PTR_SIZE);

  lptr = (void *) &(*subs[0]);
  rptr = (void *) *storage;
  for (i = 0; i < local_rows; i++) {
    *(lptr++)= (void *) rptr;
    rptr += *n * datum_size;
  }

  /* Process p-1 reads blocks of rows from file and
           sends each block to the correct destination process.
	   The last block it keeps. */

  if (id == (p-1)) {
    for (i = 0; i < p-1; i++) {
      x = fread (*storage, datum_size,
		 BLOCK_SIZE(i,p,*m) * *n, infileptr);
      MPI_Send (*storage, BLOCK_SIZE(i,p,*m) * *n, dtype,
		i, DATA_MSG, comm);
    }
    x = fread (*storage, datum_size, local_rows * *n,
	       infileptr);
    fclose (infileptr);
  } else
    MPI_Recv (*storage, local_rows * *n, dtype, p-1,
	      DATA_MSG, comm, &status);
}

/******************** OUTPUT FUNCTIONS ********************/

/*
 *   Print elements of a doubly-subscripted array.
 */

void print_submatrix (
		      void       **a,       /* OUT - Doubly-subscripted array */
		      MPI_Datatype dtype,   /* OUT - Type of array elements */
		      int          rows,    /* OUT - Matrix rows */
		      int          cols)    /* OUT - Matrix cols */
{
  int i, j;

  for (i = 0; i < rows; i++) {
    for (j = 0; j < cols; j++) {
      if (dtype == MPI_DOUBLE)
	printf ("%6.3f ", ((double **)a)[i][j]);
      else {
	if (dtype == MPI_FLOAT)
	  printf ("%6.3f ", ((float **)a)[i][j]);
	else if (dtype == MPI_INT)
	  printf ("%6d ", ((int **)a)[i][j]);
      }
    }
    putchar ('\n');
  }
}


/*
 *   Print elements of a singly-subscripted array.
 */

void print_subvector (
		      void        *a,       /* IN - Array pointer */
		      MPI_Datatype dtype,   /* IN - Array type */
		      int          n)       /* IN - Array size */
{
  int i;

  for (i = 0; i < n; i++) {
    if (dtype == MPI_DOUBLE)
      printf ("%6.3f ", ((double *)a)[i]);
    else {
      if (dtype == MPI_FLOAT)
	printf ("%6.3f ", ((float *)a)[i]);
      else if (dtype == MPI_INT)
	printf ("%6d ", ((int *)a)[i]);
    }
  }
}

/*
 *   Print a matrix that is distributed in row-striped
 *   fashion among the processes in a communicator.
 */

void print_row_striped_matrix (
			       void **a,            /* IN - 2D array */
			       MPI_Datatype dtype,  /* IN - Matrix element type */
			       int m,               /* IN - Matrix rows */
			       int n,               /* IN - Matrix cols */
			       MPI_Comm comm)       /* IN - Communicator */
{
  MPI_Status  status;          /* Result of receive */
  void       *bstorage;        /* Elements received from
				  another process */
  void      **b;               /* 2D array indexing into
				  'bstorage' */
  int         datum_size;      /* Bytes per element */
  int         i;
  int         id;              /* Process rank */
  int         local_rows;      /* This proc's rows */
  int         max_block_size;  /* Most matrix rows held by
				  any process */
  int         prompt;          /* Dummy variable */
  int         p;               /* Number of processes */

  MPI_Comm_rank (comm, &id);
  MPI_Comm_size (comm, &p);
  local_rows = BLOCK_SIZE(id,p,m);
  if (!id) {
    print_submatrix (a, dtype, local_rows, n);
    if (p > 1) {
      datum_size = get_size (dtype);
      max_block_size = BLOCK_SIZE(p-1,p,m);
      bstorage = my_malloc (id,
			    max_block_size * n * datum_size);
      b = (void **) my_malloc (id,
			       max_block_size * datum_size);
      b[0] = bstorage;
      for (i = 1; i < max_block_size; i++) {
	b[i] = b[i-1] + n * datum_size;
      }
      for (i = 1; i < p; i++) {
	MPI_Send (&prompt, 1, MPI_INT, i, PROMPT_MSG,
		  MPI_COMM_WORLD);
	MPI_Recv (bstorage, BLOCK_SIZE(i,p,m)*n, dtype,
		  i, RESPONSE_MSG, MPI_COMM_WORLD, &status);
	print_submatrix (b, dtype, BLOCK_SIZE(i,p,m), n);
      }
      free (b);
      free (bstorage);
    }
    putchar ('\n');
  } else {
    MPI_Recv (&prompt, 1, MPI_INT, 0, PROMPT_MSG,
	      MPI_COMM_WORLD, &status);
    MPI_Send (*a, local_rows * n, dtype, 0, RESPONSE_MSG,
	      MPI_COMM_WORLD);
  }
}

