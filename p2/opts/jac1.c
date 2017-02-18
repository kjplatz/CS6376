/*************************************************
 * Laplace Serial C Version
 *
 * Temperature is initially 0.0
 * Boundaries are as follows:
 *
 *      0         T         0
 *   0  +-------------------+  0
 *      |                   |
 *      |                   |
 *      |                   |
 *   T  |                   |  T
 *      |                   |
 *      |                   |
 *      |                   |
 *   0  +-------------------+ 100
 *      0         T        100
 *
 *  John Urbanic, PSC 2014
 *
 ************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <sys/time.h>
#include "config.h"

double jacobi_loop( double* Temp, double* Temp_last ) {
    double dt = 0.0;
    int i, j;
    #pragma openmp parallel for private(j,row0,row1,row2) reduction(max:dt)
    for(i = 1; i <= ROWS; i++) {
        int row0 = (i-1) * (COLUMNS+2);
        int row1 = row0 + (COLUMNS+2);
        int row2 = row0 + 2*(COLUMNS+2);

        for(j = 1; j <= COLUMNS; j++) {
            Temp[row1+j] = 0.25 * ( Temp_last[row0+j] + Temp_last[row2+j] +
                                    Temp_last[row1+j-1] + Temp_last[row1+j+1] );
            dt = fmax( dt, fabs( Temp[row1+j] - Temp_last[row1+j] ));
        }
    }
    return dt;
}
