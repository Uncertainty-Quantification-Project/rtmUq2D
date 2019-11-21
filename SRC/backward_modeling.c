#include <stdio.h>
#include <stdlib.h>

#include "fixed_parameters.h"
#include "backward_modeling.h"
#include "absorbingBoundaryCondition.h"
#include "stencilDFM.h"
#include "image_conditions.h"
#include "pi_file.h"


int adjointModeling(int simId, int *shot, float *velocity, float *wavelet, float *seismogram, float *imageCondition, float *selfImageCondition) {
    
    FILE *filePointerWavefield2;
    
    char sourceWavefieldFile_02[VEC_LEN];
    
    int nR, nW;
    int i, k;
    int dimX, dimZ;
    int time, tau, irec, pos_rec;
    int x_shot, z_shot;
    int z_receiver;
    int controler;
    
    // 2D Fields, Models, and Thickness...
    float *coef, *wavefield01, *wavefield02, *sourceWavefield;
    float *point_aux;
    
    float *a_dampingProfileLeft, *a_dampingProfileRight, *a_dampingProfileTop, *a_dampingProfileBottom;
    float *b_dampingProfileLeft, *b_dampingProfileRight, *b_dampingProfileTop, *b_dampingProfileBottom;
    
    float *psi_borderLeft, *psi_borderRight, *psi_borderTop, *psi_borderBottom;
    float *csi_borderLeft, *csi_borderRight, *csi_borderTop, *csi_borderBottom;
    
    
    // Redefining Dimensions to consider attenuation...
    dimX = thicknessLeft + Nx + thicknessRight;
    dimZ = thicknessUp   + Nz + thicknessDown;
     
    
    // +++ Allocating memory for 2D Fields and Models +++
    sourceWavefield = (float *)calloc(dimX*dimZ, sizeof(float));
    
    wavefield01     = (float *)calloc(dimX*dimZ, sizeof(float));
    wavefield02     = (float *)calloc(dimX*dimZ, sizeof(float));
    coef            = (float *)calloc(dimX*dimZ, sizeof(float));
    
    // +++ Allocating memory for Damping profiles of the CPML thickness +++
    a_dampingProfileLeft  = (float *) calloc(thicknessLeft, sizeof(float));
    a_dampingProfileRight = (float *) calloc(thicknessRight, sizeof(float));
    a_dampingProfileTop = (float *) calloc(thicknessUp, sizeof(float));
    a_dampingProfileBottom  = (float *) calloc(thicknessDown, sizeof(float));

    b_dampingProfileLeft  = (float *) calloc(thicknessLeft, sizeof(float));
    b_dampingProfileRight = (float *) calloc(thicknessRight, sizeof(float));
    b_dampingProfileTop = (float *) calloc(thicknessUp, sizeof(float));
    b_dampingProfileBottom  = (float *) calloc(thicknessDown, sizeof(float));
    
    // ............................. Constructing Psi's & Csi's	.......................... //
    psi_borderLeft = (float *) calloc(thicknessLeft*dimZ, sizeof(float));
    psi_borderRight = (float *) calloc(thicknessRight*dimZ, sizeof(float));
    psi_borderTop = (float *) calloc(dimX*thicknessUp, sizeof(float));
    psi_borderBottom = (float *) calloc(dimX*thicknessDown, sizeof(float));

    csi_borderLeft = (float *) calloc(thicknessLeft*dimZ, sizeof(float));
    csi_borderRight = (float *) calloc(thicknessRight*dimZ, sizeof(float));
    csi_borderTop = (float *) calloc(dimX*thicknessUp, sizeof(float));
    csi_borderBottom = (float *) calloc(dimX*thicknessDown, sizeof(float));
        
    
    // Calculating ...
    float aux_sampGrid = sampling/spacingGrid;
    for (i = 0; i < dimX; i++) {
        for (k = 0; k < dimZ; k++) {
            
            coef[i*dimZ + k] = (velocity[k + i*dimZ]*aux_sampGrid) * (velocity[k + i*dimZ]*aux_sampGrid);
        }
    }
    
    
    // ............................. Constructing CPML Attenuation Profile .......................... //
    attenuation_profile_cpml(a_dampingProfileLeft, a_dampingProfileRight, a_dampingProfileTop, a_dampingProfileBottom, b_dampingProfileLeft, b_dampingProfileRight, b_dampingProfileTop, b_dampingProfileBottom);
    
    
    
    // Considering absorbing layer to the shot locations...
    x_shot      = shot[0] + thicknessLeft;
    z_shot      = shot[1] + thicknessUp;
    
    z_receiver  = thicknessUp   + depthReceiver;
    // ..................
    
    sprintf(sourceWavefieldFile_02, sourceWavefieldFile2, simId);
    filePointerWavefield2 = fopen (sourceWavefieldFile_02,"rb"); //Abrir arquivo para escrever em binario "wb"-> "write binary"
    controler = (numberOfTimeStep / TIME_SLICES) - 1;
    for (time = 0; time < numberOfTimeStep; time++) {
        
        tau = numberOfTimeStep - time - 1;
        
        iso_acoustic_wave_equation(coef, wavefield01, wavefield02);
        
	aiso_attenuation_cpml(wavefield01, wavefield02, coef, a_dampingProfileLeft, b_dampingProfileLeft, a_dampingProfileRight, b_dampingProfileRight, a_dampingProfileTop, b_dampingProfileTop, a_dampingProfileBottom, b_dampingProfileBottom, psi_borderLeft, psi_borderRight, psi_borderTop, psi_borderBottom, csi_borderLeft, csi_borderRight, csi_borderTop, csi_borderBottom);
	
	
        for (irec = 0; irec < numberOfReceivers; irec++) {
            
            pos_rec = firstLocationReceiver + (spacingReceiver - 1) * irec;
            wavefield02 [(thicknessLeft+pos_rec)*dimZ + z_receiver] = seismogram[irec*numberOfTimeStep + tau];
        }
        
        if( (tau) % TIME_SLICES == 0 ) {
            
            fseek (filePointerWavefield2 , sizeof(float)*controler*(dimX*dimZ), SEEK_SET);
            fread (sourceWavefield, sizeof(float), dimX*dimZ, filePointerWavefield2);
            
            normalized_cross_correlation(imageCondition, selfImageCondition, sourceWavefield, wavefield02);
	    
            controler = controler - 1;
        }
        
        point_aux = &wavefield02[0];
	wavefield02 = &wavefield01[0];
	wavefield01 = point_aux;
        point_aux = NULL;
    }
    fclose(filePointerWavefield2);
    filePointerWavefield2 = NULL;
    point_aux = NULL;
    
    nR = remove(sourceWavefieldFile_02);
    if( nR == 0 ) {
    } else {
        perror("Error Removing source wavefield file!");
        return 1;
    }
    
    
    //	+++ Free memory +++
    free(coef);                     coef = NULL;
    free(wavefield01);              wavefield01 = NULL;
    free(wavefield02);              wavefield02 = NULL;
    free(sourceWavefield);	    sourceWavefield = NULL;
    
    free(a_dampingProfileLeft);	    a_dampingProfileLeft = NULL;
    free(a_dampingProfileRight);    a_dampingProfileRight = NULL;
    free(a_dampingProfileTop);      a_dampingProfileTop = NULL;
    free(a_dampingProfileBottom);   a_dampingProfileBottom = NULL;
    free(b_dampingProfileLeft);     b_dampingProfileLeft = NULL;
    free(b_dampingProfileRight);    b_dampingProfileRight = NULL;
    free(b_dampingProfileTop);      b_dampingProfileTop = NULL;
    free(b_dampingProfileBottom);   b_dampingProfileBottom = NULL;
    
    free(psi_borderLeft);   psi_borderLeft = NULL;
    free(psi_borderRight);   psi_borderRight = NULL;
    free(psi_borderTop);   psi_borderTop = NULL;
    free(psi_borderBottom);   psi_borderBottom = NULL;

    free(csi_borderLeft);   csi_borderLeft = NULL;
    free(csi_borderRight);   csi_borderRight = NULL;
    free(csi_borderTop);   csi_borderTop = NULL;
    free(csi_borderBottom);   csi_borderBottom = NULL;
    
    return 1;
}
