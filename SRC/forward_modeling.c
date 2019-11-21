#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

#include "fixed_parameters.h"
#include "forward_modeling.h"
#include "absorbingBoundaryCondition.h"
#include "stencilDFM.h"
#include "pi_file.h"


void isotropicAcousticModeling(int simId, int *shot, float *velocity, float *wavelet) {
    
    FILE *filePointerWavefield = NULL;
    
    char sourceWavefieldFile_02[VEC_LEN];
    
    int nS, nW;
    int i, k;
    int dimX, dimZ;
    int time, irec;
    int x_shot, z_shot;
    int z_receiver;
    int controler;
    
    // 2D Fields, Models, and Thickness...
    float *coef, *wavefield01, *wavefield02;
    float *point_aux;
    
    float *a_dampingProfileLeft, *a_dampingProfileRight, *a_dampingProfileTop, *a_dampingProfileBottom;
    float *b_dampingProfileLeft, *b_dampingProfileRight, *b_dampingProfileTop, *b_dampingProfileBottom;
    
    float *psi_borderLeft, *psi_borderRight, *psi_borderTop, *psi_borderBottom;
    float *csi_borderLeft, *csi_borderRight, *csi_borderTop, *csi_borderBottom;
    
    
    // Redefining Dimensions to considering attenuation...
    dimX = thicknessLeft + Nx + thicknessRight;
    dimZ = thicknessUp   + Nz + thicknessDown;
     
    
    // +++ Allocating memory for 2D Fields and Models +++
    wavefield01   = (float *)calloc(dimX*dimZ, sizeof(float));
    wavefield02   = (float *)calloc(dimX*dimZ, sizeof(float));
    coef          = (float *)calloc(dimX*dimZ, sizeof(float));
    
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
    
    #pragma omp parallel for private(i,k)  schedule(runtime)
    for (i = 0; i < dimX; i++) {
        #pragma omp simd
        #pragma vector aligned
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
    
    sprintf(sourceWavefieldFile_02, sourceWavefieldFile, simId);
    filePointerWavefield = fopen (sourceWavefieldFile_02,"wb"); //Abrir arquivo para escrever em binario "wb"-> "write binary"
    controler = 0;
    for (time = 0; time < numberOfTimeStep; time++) {
        
        iso_acoustic_wave_equation(coef, wavefield01, wavefield02);
	   
	aiso_attenuation_cpml(wavefield01, wavefield02, coef, a_dampingProfileLeft, b_dampingProfileLeft, a_dampingProfileRight, b_dampingProfileRight, a_dampingProfileTop, b_dampingProfileTop, a_dampingProfileBottom, b_dampingProfileBottom, psi_borderLeft, psi_borderRight, psi_borderTop, psi_borderBottom, csi_borderLeft, csi_borderRight, csi_borderTop, csi_borderBottom);
        
        
        if (time < sampligWavelet) {
            
            // right hand side forcing term...
            wavefield02[z_shot + x_shot*dimZ] = wavefield02[z_shot + x_shot*dimZ] + wavelet[time] * coef[x_shot*dimZ + z_shot];
        }
        
        if( (time) % TIME_SLICES == 0 ) {
                
                nS = fseek (filePointerWavefield , sizeof(float)*controler*(dimX*dimZ), SEEK_SET);
                
                // Writing Source Wavefield ---------------------------------------------------------------------------
                nW = fwrite (wavefield02, sizeof(float), dimX*dimZ, filePointerWavefield);
                if (nW != dimX*dimZ) {
                    printf("\nError reading file (fread in writing source wavefield in forward modeling)\n");
                    exit(1);
                }
                // ----------------------------------------------------------------------------------------------------
                
                controler = controler + 1;
        }
        
        point_aux = &wavefield02[0];
	wavefield02 = &wavefield01[0];
	wavefield01 = point_aux;
        
    }
    fclose(filePointerWavefield);
    filePointerWavefield = NULL;
    point_aux = NULL;
    
    
    //	+++ Free memory +++
    free(coef);                     coef = NULL;
    free(wavefield01);              wavefield01 = NULL;
    free(wavefield02);              wavefield02 = NULL;
    
    free(a_dampingProfileLeft);	    a_dampingProfileLeft = NULL;
    free(a_dampingProfileRight);    a_dampingProfileRight = NULL;
    free(a_dampingProfileTop);      a_dampingProfileTop = NULL;
    free(a_dampingProfileBottom);   a_dampingProfileBottom = NULL;
    free(b_dampingProfileLeft);     b_dampingProfileLeft = NULL;
    free(b_dampingProfileRight);    b_dampingProfileRight = NULL;
    free(b_dampingProfileTop);      b_dampingProfileTop = NULL;
    free(b_dampingProfileBottom);   b_dampingProfileBottom = NULL;
    
    free(psi_borderLeft);           psi_borderLeft = NULL;
    free(psi_borderRight);          psi_borderRight = NULL;
    free(psi_borderTop);            psi_borderTop = NULL;
    free(psi_borderBottom);         psi_borderBottom = NULL;

    free(csi_borderLeft);           csi_borderLeft = NULL;
    free(csi_borderRight);          csi_borderRight = NULL;
    free(csi_borderTop);            csi_borderTop = NULL;
    free(csi_borderBottom);         csi_borderBottom = NULL;
}
