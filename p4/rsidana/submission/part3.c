// File Name :    part3.c
// Course:        CS 6376 Spring 2017
// Author:        Rachna Sidana
// Assignment:    Program 4 - Part 3
// Description:   This program applies MPI directives to Cannon's Algorithm and 
//                calculates the execution time of code that has been parallelized.

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include <sys/time.h>
#define COLS  1000 // no. of columns 
#define ROWS  1000 // no. of rows

int main(int argc, char **argv) {

    MPI_Init(&argc, &argv);
    int p, rank, gsum;
    MPI_Comm_size(MPI_COMM_WORLD, &p); // no. of processes
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); 
    char i; double  numFlops;
    float gflops;
    MPI_Request request;
    struct timeval start_time, stop_time, elapsed_time;  // timers
    char a[ROWS*COLS] , a1[ROWS*COLS], c[ROWS][COLS], calc_c[ROWS][COLS];
    const int NPROWS=sqrt(p);  // number of rows in _decomposition_ 
    const int NPCOLS=sqrt(p);  // number of cols in _decomposition_ 
    const int BLOCKROWS = ROWS/sqrt(p);  // number of rows in _block_ 
    const int BLOCKCOLS = COLS/sqrt(p); // number of cols in _block_ 
    char b[BLOCKROWS*BLOCKCOLS], b1[BLOCKROWS*BLOCKCOLS],local_c[BLOCKROWS][BLOCKCOLS];
    if (rank == 0) 
    {
        for (int ii=0; ii<ROWS; ii++) {
            for (int jj=0; jj<COLS; jj++) {
                a[ii*COLS+jj]=jj-ii;
	       }	
        }
        for (int ii=0; ii<ROWS; ii++) {
            for (int jj=0; jj<COLS; jj++) {
             a1[ii*COLS+jj]=ROWS-jj+ii;
	       }
	   }
       for (int i = 0; i < ROWS; i++) {
            for (int j = 0; j < COLS; j++) {
                for (int k = 0; k < COLS; k++) {
                    gsum = gsum + (k-i)*(ROWS-j+k);
                    }
                c[i][j] = gsum;
                printf("%3d ",(int)c[i][j]);
                gsum = 0;
            }
        printf("\n");
        }
    }

    if (p != NPROWS*NPCOLS) {
        fprintf(stderr,"Error: number of PEs %d != %d x %d\n", p, NPROWS, NPCOLS);
        MPI_Finalize();
        exit(-1);
    }
    
    for (int ii=0; ii<BLOCKROWS*BLOCKCOLS; ii++) b[ii] = 0;
    for (int ii=0; ii<BLOCKROWS*BLOCKCOLS; ii++) b1[ii] = 0;
    for (int ii=0; ii<NPROWS; ii++) 
        for (int jj=0; jj<NPCOLS; jj++) local_c[ii][jj]=0; 	
    MPI_Datatype blocktype;
    MPI_Datatype blocktype2;

    //The size of the block that is to be transferred to processes is nearly square and its 
    // size is n/sqrt(p)*n/sqrt(p)
    MPI_Type_vector(BLOCKROWS, BLOCKCOLS, COLS, MPI_CHAR, &blocktype2);
    MPI_Type_create_resized( blocktype2, 0, sizeof(char), &blocktype);
    MPI_Type_commit(&blocktype);

    int disps[NPROWS*NPCOLS];
    int counts[NPROWS*NPCOLS];
    for (int ii=0; ii<NPROWS; ii++) {
        for (int jj=0; jj<NPCOLS; jj++) {
            disps[ii*NPCOLS+jj] = ii*COLS*BLOCKROWS+jj*BLOCKCOLS;
            counts [ii*NPCOLS+jj] = 1;
       		 }
            }
    int coord[2];
    int periods[2] = {1,1};
    int dims[2] ;
    dims[0]=NPROWS;
    dims[1]=NPCOLS;
    int sqrt_p = sqrt(p);
    int receivesFrom_rank, sendTo_rank;
    int b_receivesFrom_rank, b_sendTo_rank;
    MPI_Status status;
    MPI_Comm cartcomm; 
    //creating the cartesian topology
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 1, &cartcomm); 
    //Scatter the matrices on the different processes. The matrices are square/ square-like.
    MPI_Scatterv(a, counts, disps, blocktype, b, BLOCKROWS*BLOCKCOLS, MPI_CHAR, 0, cartcomm);
    MPI_Scatterv(a1, counts, disps, blocktype, b1, BLOCKROWS*BLOCKCOLS, MPI_CHAR, 0, cartcomm);

	for (int proc=0; proc<p; proc++) {
        if (proc == rank) { 
            //assign co-ordinates to the created cartesian topology
           MPI_Cart_coords(cartcomm,rank,2,coord);
	   }
    }
    //Skew the input matrices.
    // Move row i by i positions left on Matrix A
	MPI_Cart_shift(cartcomm,1,-1*coord[0],&receivesFrom_rank,&sendTo_rank);
	// Move column j by j positions up on Matrix B
	MPI_Cart_shift(cartcomm,0,-1*coord[1],&b_receivesFrom_rank,&b_sendTo_rank);

    gettimeofday(&start_time,NULL);
    MPI_Isend(&b, BLOCKROWS*BLOCKCOLS, MPI_CHAR,sendTo_rank,0,cartcomm,&request);
    MPI_Isend(&b1, BLOCKROWS*BLOCKCOLS, MPI_CHAR,b_sendTo_rank,0,cartcomm,&request);
    MPI_Recv(&b, BLOCKROWS*BLOCKCOLS, MPI_CHAR,receivesFrom_rank,0,cartcomm,&status);
    MPI_Recv(&b1, BLOCKROWS*BLOCKCOLS, MPI_CHAR,b_receivesFrom_rank,0,cartcomm,&status);
    
	for(int iter=0;iter<sqrt(p);iter++)
    {
	int sum =0;
	for (int c = 0; c < BLOCKROWS; c++) {
    	for (int d = 0; d < BLOCKCOLS; d++) {
       	for (int k = 0; k < BLOCKCOLS; k++) {
          sum += b[c*BLOCKCOLS+k]*b1[k*BLOCKCOLS+d];
	}
        local_c[c][d] = local_c[c][d]+sum;
        sum = 0;
      }
    }

    // shift row by 1 to left and column by 1 to up.
    int current_pos_x,current_pos_y, leftRank,upRank;
    MPI_Cart_shift(cartcomm,1,-1,&current_pos_x,&leftRank);
    MPI_Isend(&b, BLOCKROWS*BLOCKCOLS, MPI_CHAR,leftRank,0,cartcomm,&request);
    MPI_Cart_shift(cartcomm,0,-1,&current_pos_y,&upRank);
    MPI_Isend(&b1, BLOCKROWS*BLOCKCOLS, MPI_CHAR,upRank,0,cartcomm,&request);
    MPI_Recv(&b, BLOCKROWS*BLOCKCOLS, MPI_CHAR,current_pos_x,0,cartcomm,&status);
    MPI_Recv(&b1, BLOCKROWS*BLOCKCOLS, MPI_CHAR,current_pos_y,0,cartcomm,&status);
    }
    MPI_Barrier(cartcomm);
    if (0 == rank) {
    gettimeofday(&stop_time,NULL);
    timersub(&stop_time, &start_time, &elapsed_time); // Unix time subtract routine
    printf("Total time was %f seconds.\n", elapsed_time.tv_sec+elapsed_time.tv_usec/1000000.0);
    numFlops = (2.0f*ROWS*ROWS*ROWS-ROWS*ROWS)/1000000000.0f;
    gflops = numFlops/(elapsed_time.tv_sec+elapsed_time.tv_usec/1000000.0);
    printf("GFlops :  %f .\n",gflops); 
    }

    // uncomment this if you want to see the local matrices on various processes:

/*    printf("Final Local Matrix C on Process: %d\n", rank);
     for (int ii=0; ii<BLOCKROWS; ii++) {
       for (int jj=0; jj<BLOCKCOLS; jj++) {
           printf("%3d ",(int)local_c[ii][jj]);
           }
       printf("\n");
    }*/

/*-------------------------------------------------------------------
 * VERIFICATION OF CORRECT RESULTS :
 * Purpose:   verify the multiplication of first and second matrix
 * Uncomment the below code to gather matrices on P0 and verify the result of matrix multiplication
 * int verification = 1, if matrix multiplication is corect, 0 otherwise
 --------------------------------------------------------------------*/



/*
MPI_Gatherv(&local_c,BLOCKROWS*BLOCKCOLS, MPI_CHAR, calc_c,counts,disps , blocktype, 0, MPI_COMM_WORLD);

if(rank==0)
{

printf("------ ----------   Gathered matrix on P: %d\n",rank);
for (int ii=0; ii<ROWS; ii++) {
            for (int jj=0; jj<COLS; jj++) {
             printf("%3d ",calc_c[ii][jj]);
               }
printf("\n");
           }
printf("\n");
int verification=1;
for (int ii=0; ii<ROWS; ii++) {
            for (int jj=0; jj<COLS; jj++) {
                if(calc_c[ii][jj]!=c[ii][jj])
                {
                        printf(" WRONG !! ");
                        verification=0;
                }

        }
}

printf("Verification = %d\n",verification);
if(verification==1)
printf("\n******** Correct Result *********\n ");

}*/



    MPI_Finalize();
    return 0;
}

