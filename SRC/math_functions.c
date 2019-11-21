#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "fixed_parameters.h"
#include "math_functions.h"


void laplacian_filter2D (float *lp, float *inputMatrix) {
    
    int i, k;
    
    
    #pragma omp parallel for private(i,k)  schedule(runtime)
    for (i = HALFORDERDF; i < Nx-HALFORDERDF; i++) {
        #pragma omp simd
        #pragma vector aligned
        for (k = HALFORDERDF; k < Nz-HALFORDERDF; k++) {
            
            lp[k + i*Nz] = 2.0f*AC08*( inputMatrix[(k) + (i)*Nz] )
                                                                + AC18 * ( inputMatrix[(k) + (i-1)*Nz] + inputMatrix[(k) + (i+1)*Nz] + inputMatrix[(k-1) + (i)*Nz] + inputMatrix[(k+1) + (i)*Nz] )
                                                                + AC28 * ( inputMatrix[(k) + (i-2)*Nz] + inputMatrix[(k) + (i+2)*Nz] + inputMatrix[(k-2) + (i)*Nz] + inputMatrix[(k+2) + (i)*Nz] )
                                                                + AC38 * ( inputMatrix[(k) + (i-3)*Nz] + inputMatrix[(k) + (i+3)*Nz] + inputMatrix[(k-3) + (i)*Nz] + inputMatrix[(k+3) + (i)*Nz] )
                                                                + AC48 * ( inputMatrix[(k) + (i-4)*Nz] + inputMatrix[(k) + (i+4)*Nz] + inputMatrix[(k-4) + (i)*Nz] + inputMatrix[(k+4) + (i)*Nz] );
        }
    }
}