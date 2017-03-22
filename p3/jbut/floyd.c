/*****************
 * floyd.c
 * Author: Jay Butera
 *
 * FUNCTIONS
 *
 * file_to_mat(char*, int, int)
 * | Read binary file of float32s representing the adjacency matrix of a graph.
 * | In the case of MPI, each node only reads its contiguous block
 *
 * compute_shortest_paths(int, int, float**)
 * | Implements floyds algorithm.
 * | In MPI, each n^2 iteration of k is parallelized into each node's block.
 * | The node containing the kth row broadcasts to all other nodes.
 * | Using OpenACC and OMP further parallellize the block on each node.
 *
 * print_matrix(int, int, int)
 * | Prints matrix A
 * | With MPI, each node sends its block to node p-1 which prints the matrix as
 * | a whole.
 *****************/

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <omp.h>
#include <sys/time.h>

#define ROWS 10
#define COLS ROWS

#define BLOCK_LOW(id,p) (id * ROWS) / p
#define BLOCK_OWNER(k,p,n) (((p)*(k+1)-1)/n)
#define MIN(x,y) (x > y) ? y : x

float** A;
float* Astorage;

/*
 * file_to_mat
 * Jay Butera
 *
 * Description
 * Read binary file of float32s representing the adjacency matrix of a graph.
 * In the case of MPI, each node only reads its contiguous block
 *
 * Arguments
 * filename - I/P - char* - file name string
 * id - I/P - int - node id
 * p - I/P - int - proccessor count
 */
void file_to_mat(char* filename, int id, int p) {
    FILE* f = fopen(filename, "rb");

    int start_row = (id * ROWS) / p;
    int end_row   = ((id+1) * ROWS) / p - 1;

    printf("Proc %d reading lines [%d-%d]\n", id, start_row, end_row);

    fseek(f, start_row * COLS * sizeof(float), SEEK_SET);  // Jump to the end of the file
    long offset = (end_row - start_row + 1) * COLS * sizeof(float);

    fread((void*)(Astorage), offset, 1, f);

    fclose(f);
}

/*
 * compute_shortest_paths
 * Jay Butera
 *
 * Description
 * Implements floyds algorithm.
 * In MPI, each n^2 iteration of k is parallelized into each node's block.
 * The node containing the kth row broadcasts to all other nodes.
 * Using OpenACC and OMP further parallellize the block on each node.
 *
 * Arguments
 * id - I/P - int - node id
 * p - I/P - int - proccessor count
 * a - I/O - float** - adj. matrix to shortest paths
 */
void compute_shortest_paths (int id, int p, float** a) {
    int i, j, k;
    int offset;
    int root;
    float* tmp;

    int start_row = (id * ROWS) / p;
    int end_row   = ((id+1) * ROWS) / p - 1;
    int local_rows = (end_row - start_row)+1;

    tmp = (float*) malloc (COLS * sizeof(float));
    //#pragma acc data, copy(a)
    for (k = 0; k < ROWS; k++) {

        root = BLOCK_OWNER(k,p,ROWS);

        if (k == ROWS-1) {
            root = p - 1;
        }

        if (root == id) {
            offset = k - BLOCK_LOW(id,p);
            #pragma omp parallel for
            for (j = 0; j < ROWS; j++)
                tmp[j] = a[offset][j];
        }

        MPI_Bcast(tmp, ROWS, MPI_FLOAT, root, MPI_COMM_WORLD);
        //#pragma acc kernels
        #pragma omp parallel for private(j)
        for (i = 0; i < local_rows; i++)
            for (j = 0; j < ROWS; j++)
                a[i][j] = MIN(a[i][j], a[i][k]+tmp[j]);

    }

    free( tmp );
}

/*
 * print_matrix
 * Jay Butera
 *
 * Description
 * Prints matrix A
 * With MPI, each node sends its block to node p-1 which prints the matrix as
 * a whole.
 *
 * Arguments
 * id - I/P - int - node id
 * num_procs - I/P - int - proccessor count
 * local_rows - I/P - int - number of rows assigned to the block associated
 *                          with the node running the program
 */
print_matrix(int id, int num_procs, int local_rows) {
    // Print matrix
    if (id == num_procs-1) {
        // Tmp matrix to store nodes' blocks
        //----------------------------
        float* Astorage_tmp = (float*) malloc(local_rows * COLS * sizeof(float));
        if (Astorage_tmp == NULL) {
            printf("Astorage mem could not allocate\n");
            exit(0);
        }

        float** A_tmp = (float**) malloc(local_rows * sizeof(float*));
        if (A_tmp == NULL) {
            printf("A mem could not allocate\n");
            exit(0);
        }

        int i;
        for (i = 0; i < local_rows; i++) {
            A_tmp[i] = &Astorage_tmp[i * COLS];
        }
        //----------------------------

        // Recieve each block from p-1 nodes
        MPI_Status status;
        for (id = 0; id < num_procs-1; id++) {
            MPI_Recv(
                     Astorage_tmp,
                     (ROWS / num_procs) * COLS,
                     MPI_FLOAT,
                     id,
                     0,
                     MPI_COMM_WORLD,
                     &status);

            // Print sub matrix
            int j;
            for (i = 0; i < (ROWS / num_procs ); i++) {
                for (j = 0; j < COLS; j++)
                    printf("%10.2f ", A_tmp[i][j]);
                printf("\n");
            }
        }

        // Print sub matrix
        int j;
        for (i = 0; i < local_rows; i++) {
            for (j = 0; j < COLS; j++)
                printf("%10.2f ", A[i][j]);
            printf("\n");
        }

        free(Astorage_tmp);
        free(A_tmp);
    }
    else {
        printf("Send %d elements from proc %d\n", local_rows * COLS, id);
        MPI_Send(Astorage, local_rows * COLS, MPI_FLOAT, num_procs-1, 0, MPI_COMM_WORLD);
    }
}

int main (int argc, char** argv) {
    int id; // Process id
    int num_procs; // Number of processors

    // Time record
    struct timeval start_time, stop_time, elapsed_time;

    // Init MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    // Indices for contiguous block assignment
    int start_row = (id * ROWS) / num_procs;
    int end_row   = ((id+1) * ROWS) / num_procs - 1;
    int local_rows = (end_row - start_row)+1;

    // Allocate space
    Astorage = (float*) malloc(local_rows * COLS * sizeof(float)); if (Astorage == NULL) {
        printf("Astorage mem could not allocate\n");
        exit(0);
    }

    A = (float**) malloc(local_rows * sizeof(float*));
    if (A == NULL) {
        printf("A mem could not allocate\n");
        exit(0);
    }

    int i;
    for (i = 0; i < local_rows; i++) {
        A[i] = &Astorage[i * COLS];
    }

    // Read in block from file
    file_to_mat("mp_mat", id, num_procs);


    //-----------
    // Start time
    //-----------
    gettimeofday(&start_time,NULL);

    compute_shortest_paths(id, num_procs, A);
    MPI_Barrier(MPI_COMM_WORLD);

    //---------
    // End time
    //---------
    gettimeofday(&stop_time,NULL);
    timersub(&stop_time, &start_time, &elapsed_time); // Unix time subtract routine

    // 2N^3/T
    if (!id) {
        float GFLOPS = (float)(2.f*ROWS*ROWS*ROWS) / (1000000000.f*(elapsed_time.tv_sec+elapsed_time.tv_usec/1000000.0));

        printf("elapsed time (s): %f\n", ((elapsed_time.tv_sec+elapsed_time.tv_usec/1000000.0)));
        printf("GFLOPS: %f\n", GFLOPS);
    }

    // Print full matrix
    //print_matrix(id, num_procs, local_rows);

    // Dealloc
    free(A);
    free(Astorage);

    MPI_Finalize();
    return 0;
}
