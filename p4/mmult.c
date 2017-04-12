#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifndef SIZE
#define SIZE 8
#endif

typedef float value_type;

void init_A( value_type X[SIZE][SIZE] );
void init_B( value_type X[SIZE][SIZE] );
void print_array( value_type X[SIZE][SIZE]);

void matrix_mult( value_type A[SIZE][SIZE], value_type B[SIZE][SIZE], 
                  value_type C[SIZE][SIZE] ) {
    for( int i=0; i<SIZE; ++i ) {
        for( int j=0; j<SIZE; ++j ) {
            C[i][j] = 0.0f;
            for( int k=0; k<SIZE; ++k ) {
                C[i][j] += A[i][k] * B[k][j];    
            }
        }
    }
}

void verify( value_type A[SIZE][SIZE],
             value_type B[SIZE][SIZE],
             value_type C[SIZE][SIZE] );

int main( int argc, char** argv ) {
    int size = 0;
    value_type A[SIZE][SIZE];
    value_type B[SIZE][SIZE];
    value_type C[SIZE][SIZE];

    init_A( A );
    init_B( B );

    matrix_mult( A, B, C );

    printf( "Array A\n" );
    print_array( A );
    printf( "Array B\n" );
    print_array( B );
    printf( "Array C\n" );
    print_array( C );

    verify( A, B, C );
}

void print_array( value_type X[SIZE][SIZE] ) {
    for( int i=0; i<SIZE; ++i ) {
        printf( "Row %d", i );
        for( int j=0; j<SIZE; ++j ) {
            printf( " %5.2f", X[i][j] );
        }
        printf( "\n" );
    }
}

void init_A( value_type A[SIZE][SIZE] ) {
    srand(time(NULL));
    for( int i=0; i<SIZE; ++i ) {
        printf( "A Row %d", i );
        for( int j=0; j<SIZE; ++j ) {
            value_type x = rand() % 100;
            A[i][j] = x;
            printf( " %5.2f",  x);
        }
        printf( "\n" );
    }
}

void init_B( value_type B[SIZE][SIZE] ) {
    memset( B, 0, SIZE*SIZE*sizeof(value_type));
    for( int i=0; i<SIZE; ++i ) {
        B[SIZE-i-1][i] = 1;
    }
}

void verify( value_type A[SIZE][SIZE],
             value_type B[SIZE][SIZE],
             value_type C[SIZE][SIZE] ) {
    for( int i=0; i<SIZE; ++i ) {
        for( int j=0; j<SIZE; ++j ) {
            if ( A[i][j] != C[i][SIZE-j-1] ) {
                fprintf( stderr, "Incorrect output detected at [%d][%d]!\n", i, j );
                fprintf( stderr, "A[%d][%d] = %f, C[%d][%d] = %f\n", i, j, A[i][j], i, (SIZE-j-1), C[i][SIZE-j-1] );
                printf( "Array A\n" );
                print_array( A );
                printf( "Array B\n" );
                print_array( B );
                printf( "Array C\n" );
                print_array( C );
                abort();
            }
        }
    }
    
    fprintf( stderr, "Output verified\n" );
}
