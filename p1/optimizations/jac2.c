#include <math.h>
#include <stdio.h>
#include "config.h"

extern double Temperature[ROWS+2][COLUMNS+2];
extern double Temperature_last[ROWS+2][COLUMNS+2];

double jacobi_loop( int row, double* Temp, double* Temp_last ) {
    double dt = 0.0;
    int off1 = row * (COLUMNS+2);
    int off0 = off1 - (COLUMNS+2);
    int off2 = off1 + (COLUMNS+2);
    if ( row == 1000 ) printf( "Temperature[%d][1000] = %5.2f, Temperature_old[%d][1000] = %5.2f\n", row, Temp[off1+1000], row, Temp_last[off1+1000] );

    for( int j=1; j <= COLUMNS; j++ ) {
        Temp[off1+j] = 0.25 * ( Temp_last[off1+j-1] + Temp_last[off1+j+1] +
                                Temp_last[off0+j] + Temp_last[off2+j] );

        dt = fmax( dt, fabs( Temp[off1+j] - Temp_last[off1+j] ) );
    }
    return dt;
}
