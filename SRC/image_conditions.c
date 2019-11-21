#include <stdio.h>
#include <stdlib.h>

#include "fixed_parameters.h"
#include "image_conditions.h"


void cross_correlation(float *imageCondition, float *sourceWavefield, float *receiverWavefield) {
    
    int i, k;
    
    for (i = 0; i < Nx; i++) {
      for (k = 0; k < Nz; k++) {
	
	    imageCondition[i*Nz + k] = imageCondition[i*Nz + k] + sourceWavefield[(i+thicknessLeft)*(thicknessUp+Nz+thicknessDown) + (k+thicknessUp)] * receiverWavefield[(i+thicknessLeft)*(thicknessUp+Nz+thicknessDown) + (k+thicknessUp)]; 
      }
    }
}

void normalized_cross_correlation(float *imageCondition, float *selfImageCondition, float *sourceWavefield, float *receiverWavefield) {
    
    int i, k;
    float epsilon = 1.0f;
    
    for (i = 0; i < Nx; i++) {
      for (k = 0; k < Nz; k++) {
	
	    imageCondition[i*Nz + k] = imageCondition[i*Nz + k] + ( sourceWavefield[(i+thicknessLeft)*(thicknessUp+Nz+thicknessDown) + (k+thicknessUp)] * receiverWavefield[(i+thicknessLeft)*(thicknessUp+Nz+thicknessDown) + (k+thicknessUp)] );
            
            selfImageCondition[i*Nz + k] = selfImageCondition[i*Nz + k] + ( sourceWavefield[(i+thicknessLeft)*(thicknessUp+Nz+thicknessDown) + (k+thicknessUp)]*sourceWavefield[(i+thicknessLeft)*(thicknessUp+Nz+thicknessDown) + (k+thicknessUp)] ); 
      }
    }
}