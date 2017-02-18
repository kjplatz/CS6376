#include <math.h>
#include "config.h"
#include <stdio.h>

extern double Temperature[ROWS+2][COLUMNS+2];
extern double Temperature_last[ROWS+2][COLUMNS+2];

double jacobi_loop( int row, double* Temp, double* Temp_last ) {
    double dt = 0.0;
    if ( row == 1000 ) printf( "Temperature[%d][1000] = %5.2f, Temperature_last[%d][1000] = %5.2f\n", row, Temperature[row][1000], row, Temperature_last[row][1000] );
    for( int j=1; j <= COLUMNS; j++ ) {
        Temperature[row][j] = 0.25 * ( Temperature_last[row][j-1] + Temperature_last[row][j+1] +
                                       Temperature_last[row-1][j] + Temperature_last[row+1][j] );

        // dt = fmax( dt, fabs( Temperature[row][j] - Temperature_last[row][j] ) );
    }

    return dt;
}
