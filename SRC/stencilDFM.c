#include <omp.h>

#include "fixed_parameters.h"
#include "stencilDFM.h"


void iso_acoustic_wave_equation(float *pond, float *waveField2, float *waveField3) {         // Ok!
  
  int i, k;
  
  // Redefining Dimensions to considering attenuation...
  int dim1 = thicknessLeft + Nx + thicknessRight;
  int dim2 = thicknessUp   + Nz + thicknessDown;
  
  int halforderdf = HALFORDERDF;
  
  #pragma omp parallel for private(i,k)  schedule(runtime)
        for (i = halforderdf; i < dim1-halforderdf; i++) {
                #pragma omp simd
                #pragma vector aligned
                for (k = halforderdf; k < dim2-halforderdf; k++) {

                        waveField3[k + i*dim2] = pond[k + i*dim2] * ( 2.0f*AC08*( waveField2[(k) + (i)*dim2] )
                                                                         + AC18 * ( waveField2[(k) + (i-1)*dim2] + waveField2[(k) + (i+1)*dim2] + waveField2[(k-1) + (i)*dim2] + waveField2[(k+1) + (i)*dim2] )
                                                                         + AC28 * ( waveField2[(k) + (i-2)*dim2] + waveField2[(k) + (i+2)*dim2] + waveField2[(k-2) + (i)*dim2] + waveField2[(k+2) + (i)*dim2] )
                                                                         + AC38 * ( waveField2[(k) + (i-3)*dim2] + waveField2[(k) + (i+3)*dim2] + waveField2[(k-3) + (i)*dim2] + waveField2[(k+3) + (i)*dim2] )
                                                                         + AC48 * ( waveField2[(k) + (i-4)*dim2] + waveField2[(k) + (i+4)*dim2] + waveField2[(k-4) + (i)*dim2] + waveField2[(k+4) + (i)*dim2] ))
								         + 2.0f * ( waveField2[(k) + (i)*dim2] ) - waveField3[(k) + (i)*dim2];
                }
        }
}