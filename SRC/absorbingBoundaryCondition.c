#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <omp.h>

#include "fixed_parameters.h"
#include "absorbingBoundaryCondition.h"



void attenuation_profile_cpml(float *a_Left, float *a_Right, float *a_Top, float *a_Bottom, float *b_Left, float *b_Right, float *b_Top, float *b_Bottom) {

        int i, ii, k, kk;
        
        int dim1, dim2;

        float dominantFrequency,

              thicknessLeft_meters, thicknessRight_meters, thicknessTop_meters, thicknessBottom_meters,

              maximumSigma_x, maximumSigma_z,

              xPosition_meters, zPosition_meters,

              abscissa_meters,
              abscissa_normalized,

              modelOrigin_x, modelEnd_x,
              modelOrigin_z, modelEnd_z,

              MAXIMUMALPHA;


        float *alpha_left = NULL, *alpha_right = NULL,
              *alpha_top = NULL, *alpha_bottom = NULL,

              *sigma_left = NULL, *sigma_right = NULL,
              *sigma_top = NULL, *sigma_bottom = NULL;

        const float pi 			= M_PI;
	const float attenuationCoef 	= 3.0f*pi;
	const float maximumVelocity 	= 4000.0f;
	const float cutFrequency 	= 60.0;


        alpha_left   = (float *) calloc (thicknessLeft, sizeof(float));
        alpha_right  = (float *) calloc (thicknessRight, sizeof(float));
        alpha_top    = (float *) calloc (thicknessUp, sizeof(float));
        alpha_bottom = (float *) calloc (thicknessDown, sizeof(float));

        sigma_left   = (float *) calloc (thicknessLeft, sizeof(float));
        sigma_right  = (float *) calloc (thicknessRight, sizeof(float));
        sigma_top    = (float *) calloc (thicknessUp, sizeof(float));
        sigma_bottom = (float *) calloc (thicknessDown, sizeof(float));

        dominantFrequency = cutFrequency/(3.0f*sqrt(pi));


        thicknessLeft_meters   = (thicknessLeft - ORDERDF) * spacingGrid;
        thicknessRight_meters  = (thicknessRight - ORDERDF) * spacingGrid;
        thicknessTop_meters    = (thicknessUp - ORDERDF) * spacingGrid;
        thicknessBottom_meters = (thicknessDown - ORDERDF) * spacingGrid;

	
        maximumSigma_x = (attenuationCoef) * (((NPOWER + 1.0f) * maximumVelocity)/(thicknessLeft_meters));

        MAXIMUMALPHA = pi*dominantFrequency;
        
        // Redefining Dimensions to considering attenuation...
        dim1 = thicknessLeft + Nx + thicknessRight;
        dim2 = thicknessUp   + Nz + thicknessDown;


        // ............................. Attenuation Profile Left Border ............................. //
        for (i = HALFORDERDF; i < thicknessLeft-HALFORDERDF; i++) {

                xPosition_meters = spacingGrid * (i - HALFORDERDF);
                abscissa_meters = thicknessLeft_meters - xPosition_meters;
                abscissa_normalized = abscissa_meters / thicknessLeft_meters;

                sigma_left[i] = maximumSigma_x * powf(abscissa_normalized, NPOWER);
                alpha_left [i] = MAXIMUMALPHA * (1.0f - abscissa_normalized) + 0.1f * MAXIMUMALPHA;

                b_Left[i] = powf(M_E, -(sigma_left[i] + alpha_left[i]) * sampling);
                a_Left[i] = sigma_left[i] * (b_Left[i] - 1.0f) / (sigma_left[i] + alpha_left[i]);
        }


        // ............................. Attenuation Profile Right Border ............................. //
        ii = 0;
        for (i = dim1 - (thicknessRight-HALFORDERDF); i < dim1-HALFORDERDF; i++) {

                abscissa_meters = spacingGrid + spacingGrid * (ii);
                abscissa_normalized = abscissa_meters / thicknessRight_meters;

                sigma_right[ii] = maximumSigma_x * powf(abscissa_normalized, NPOWER);
                alpha_right[ii] = MAXIMUMALPHA * (1.0f - abscissa_normalized) + 0.1f * MAXIMUMALPHA;

                b_Right[ii + HALFORDERDF] = powf(M_E, -(sigma_right[ii] + alpha_right[ii]) * sampling);
                a_Right[ii + HALFORDERDF] = sigma_right[ii] * (b_Right[ii + HALFORDERDF] - 1.0f) / (sigma_right[ii] + alpha_right[ii]);

                ii++;
        }


        // ............................. Attenuation Profile Top Border ............................. //
        for (k = HALFORDERDF; k < thicknessUp-HALFORDERDF; k++) {
                
                xPosition_meters = spacingGrid * (k - HALFORDERDF);
                abscissa_meters = thicknessTop_meters - xPosition_meters;
                abscissa_normalized = abscissa_meters / thicknessTop_meters;

                sigma_top[k] = maximumSigma_x * powf(abscissa_normalized, NPOWER);
                alpha_top [k] = MAXIMUMALPHA * (1.0f - abscissa_normalized) + 0.1f * MAXIMUMALPHA;

                b_Top[k] = powf(M_E, -(sigma_top[k] + alpha_top[k]) * sampling);
                a_Top[k] = sigma_top[k] * (b_Top[k] - 1.0f) / (sigma_top[k] + alpha_top[k]);
        }


        // ............................. Attenuation Profile Bottom Border ............................. //;
        kk = 0;
        for (k = dim2 - (thicknessDown-HALFORDERDF); k < dim2-HALFORDERDF; k++) {

                abscissa_meters = spacingGrid + spacingGrid * (kk);
                abscissa_normalized = abscissa_meters / thicknessBottom_meters;

                sigma_bottom[kk] = maximumSigma_x * powf(abscissa_normalized, NPOWER);
                alpha_bottom[kk] = MAXIMUMALPHA * (1.0f - abscissa_normalized) + 0.1f * MAXIMUMALPHA;

                b_Bottom[kk + HALFORDERDF] = powf(M_E, -(sigma_bottom[kk] + alpha_bottom[kk]) * sampling);
                a_Bottom[kk + HALFORDERDF] = sigma_bottom[kk] * (b_Bottom[kk + HALFORDERDF] - 1.0f) / (sigma_bottom[kk] + alpha_bottom[kk]);

                kk++;
        }
	
	
	free(alpha_left);	alpha_left = NULL;
        free(alpha_right);	alpha_right = NULL;
        free(alpha_top);	alpha_top = NULL;
        free(alpha_bottom);	alpha_bottom = NULL;

        free(sigma_left);	sigma_left = NULL;
        free(sigma_right);	sigma_right = NULL;
        free(sigma_top);	sigma_top = NULL;
        free(sigma_bottom);	sigma_bottom = NULL;
	
}

void aiso_attenuation_cpml(float *waveField1, float *waveField2, float *Coef, float *a_left, float *b_left, float *a_right, float *b_right, float *a_top, float *b_top, float *a_bottom, float *b_bottom, float *psi_borderLeft, float *psi_borderRight, float *psi_borderTop, float *psi_borderBottom, float *csi_borderLeft, float *csi_borderRight, float *csi_borderTop, float *csi_borderBottom){

        int i, k, ii, kk;
        int dim1, dim2;
        
        float h1, h2, h3;

        float dpsi_borderLeft, dpsi_borderRight, dpsi_borderTop, dpsi_borderBottom;

        h1 = 1.0f/spacingGrid;
        h2 = spacingGrid*spacingGrid;
        h3 = 1.0f/(spacingGrid*spacingGrid);
        
        // Redefining Dimensions to considering attenuation...
        dim1 = thicknessLeft + Nx + thicknessRight;
        dim2 = thicknessUp   + Nz + thicknessDown;
	
	//Auxiliary variables...
	int dimThickR = dim1-thicknessRight;
	int dimThickD = dim2-thicknessDown;

        // ............................. Evaluating the auxiliary variable of the Left Border ............................... //
	#pragma omp parallel
	{
	#pragma omp for private(i,k)  schedule(runtime)
        for (i = HALFORDERDF; i < thicknessLeft-HALFORDERDF; i++) {
		#pragma omp simd
                #pragma vector aligned
                for (k = HALFORDERDF; k < dim2-HALFORDERDF; k++) {

                        psi_borderLeft [(dim2)*(i) + (k)] =
                                b_left [i] * psi_borderLeft [(dim2)*(i) + (k)] +
                                a_left [i] * (
                                        B48*(waveField1 [(dim2)*(i+4) + (k)] - waveField1 [(dim2)*(i-4) + (k)]) +
                                        B38*(waveField1 [(dim2)*(i+3) + (k)] - waveField1 [(dim2)*(i-3) + (k)]) +
                                        B28*(waveField1 [(dim2)*(i+2) + (k)] - waveField1 [(dim2)*(i-2) + (k)]) +
                                        B18*(waveField1 [(dim2)*(i+1) + (k)] - waveField1 [(dim2)*(i-1) + (k)])
                                        ) * h1;
                }
        }
        
        // ............................. Addition of the corrective terms of the Left Border ................................ //
        #pragma omp for private(i, k, dpsi_borderLeft)  schedule(runtime)
        for (i = HALFORDERDF; i < thicknessLeft-HALFORDERDF; i++) {
		#pragma omp simd
                #pragma vector aligned
                for (k = HALFORDERDF; k < dim2-HALFORDERDF; k++) {
                        
                        dpsi_borderLeft = (
                                B48*(psi_borderLeft[(dim2)*(i+4) + (k)] - psi_borderLeft[(dim2)*(i-4) + (k)]) +
                                B38*(psi_borderLeft[(dim2)*(i+3) + (k)] - psi_borderLeft[(dim2)*(i-3) + (k)]) +
                                B28*(psi_borderLeft[(dim2)*(i+2) + (k)] - psi_borderLeft[(dim2)*(i-2) + (k)]) +
                                B18*(psi_borderLeft[(dim2)*(i+1) + (k)] - psi_borderLeft[(dim2)*(i-1) + (k)])
                                ) * h1;

                        csi_borderLeft[(dim2)*(i) + (k)] =
                                b_left[i]*csi_borderLeft[(dim2)*(i) + (k)] +
                                a_left[i]*((
                                                   AC48*(waveField1[(dim2)*(i+4) + (k)] + waveField1[(dim2)*(i-4) + (k)]) +
                                                   AC38*(waveField1[(dim2)*(i+3) + (k)] + waveField1[(dim2)*(i-3) + (k)]) +
                                                   AC28*(waveField1[(dim2)*(i+2) + (k)] + waveField1[(dim2)*(i-2) + (k)]) +
                                                   AC18*(waveField1[(dim2)*(i+1) + (k)] + waveField1[(dim2)*(i-1) + (k)]) +
                                                   AC08*waveField1[(dim2)*(i) + (k)]
                                                   ) * h3 + dpsi_borderLeft);

                        waveField2[(dim2)*(i) + (k)] = waveField2[(dim2)*(i) + (k)] + (Coef[(dim2)*(i) + (k)]*h2)*(dpsi_borderLeft + csi_borderLeft[(dim2)*(i) + (k)]);
                }
        }
        
        // ............................. Evaluating the auxiliary variable of the Top Border ................................ //
        #pragma omp for private(i,k)  schedule(runtime)
        for (i = HALFORDERDF; i < dim1-HALFORDERDF; i++) {
	        #pragma omp simd
                #pragma vector aligned
                for (k = HALFORDERDF; k < thicknessUp-HALFORDERDF; k++) {

                        psi_borderTop[(thicknessUp)*(i) + (k)] =
                                b_top[k] * psi_borderTop[(thicknessUp)*(i) + (k)] +
                                a_top[k] * (
                                        B48*(waveField1[(dim2)*(i) + (k+4)] - waveField1[(dim2)*(i) + (k-4)]) +
                                        B38*(waveField1[(dim2)*(i) + (k+3)] - waveField1[(dim2)*(i) + (k-3)]) +
                                        B28*(waveField1[(dim2)*(i) + (k+2)] - waveField1[(dim2)*(i) + (k-2)]) +
                                        B18*(waveField1[(dim2)*(i) + (k+1)] - waveField1[(dim2)*(i) + (k-1)])
                                        ) * h1;
                }
        }
        
        // ............................. Addition of the corrective terms of the Top Border ................................. //
        #pragma omp for private(i, k, dpsi_borderTop)  schedule(runtime)
        for (i = 0; i < dim1; i++) {
		#pragma omp simd
                #pragma vector aligned
                for (k = HALFORDERDF; k < thicknessUp-HALFORDERDF; k++) {

                        dpsi_borderTop = (
                                B48*(psi_borderTop[(thicknessUp)*(i) + (k+4)] - psi_borderTop[(thicknessUp)*(i) + (k-4)]) +
                                B38*(psi_borderTop[(thicknessUp)*(i) + (k+3)] - psi_borderTop[(thicknessUp)*(i) + (k-3)]) +
                                B28*(psi_borderTop[(thicknessUp)*(i) + (k+2)] - psi_borderTop[(thicknessUp)*(i) + (k-2)]) +
                                B18*(psi_borderTop[(thicknessUp)*(i) + (k+1)] - psi_borderTop[(thicknessUp)*(i) + (k-1)])
                                ) * h1;

                        csi_borderTop[(thicknessUp)*(i) + (k)] =
                                b_top[k] * csi_borderTop[(thicknessUp)*(i) + (k)] +
                                a_top[k] * ((
                                                    AC48*(waveField1[(dim2)* (i) + (k+4)] + waveField1[(dim2)*(i) + (k-4)]) +
                                                    AC38*(waveField1[(dim2)* (i) + (k+3)] + waveField1[(dim2)*(i) + (k-3)]) +
                                                    AC28*(waveField1[(dim2)* (i) + (k+2)] + waveField1[(dim2)*(i) + (k-2)]) +
                                                    AC18*(waveField1[(dim2)* (i) + (k+1)] + waveField1[(dim2)*(i) + (k-1)]) +
                                                    AC08*waveField1[(dim2)*(i) + (k)]
                                                    ) * h3 + dpsi_borderTop);

                        waveField2[(dim2)*(i) + (k)] = waveField2[(dim2)*(i) + (k)] + (Coef[(dim2)*(i) + (k)]*h2)*(dpsi_borderTop + csi_borderTop[(thicknessUp)*(i) + (k)]);
                }
        }

        // ............................. Evaluating the auxiliary variable of the Right Border .............................. //
        #pragma omp for private(i,k)  schedule(runtime)
        for (i = dim1-(thicknessRight-HALFORDERDF); i < dim1-HALFORDERDF; i++) {
		#pragma omp simd
                #pragma vector aligned
                for (k = HALFORDERDF; k < dim2-HALFORDERDF; k++) {

                        psi_borderRight [(dim2)*(i-(dimThickR)) + (k)] =
                                b_right [i-(dimThickR)] * psi_borderRight [(dim2)*(i-(dimThickR)) + (k)] +
                                a_right [i-(dimThickR)] * (
                                        B48*(waveField1 [(dim2)*(i+4) + (k)] - waveField1 [(dim2)*(i-4) + (k)]) +
                                        B38*(waveField1 [(dim2)*(i+3) + (k)] - waveField1 [(dim2)*(i-3) + (k)]) +
                                        B28*(waveField1 [(dim2)*(i+2) + (k)] - waveField1 [(dim2)*(i-2) + (k)]) +
                                        B18*(waveField1 [(dim2)*(i+1) + (k)] - waveField1 [(dim2)*(i-1) + (k)])
                                        ) * h1;
                }
        }
        
        // ............................. Addition of the corrective terms of the Right Border ............................... //
	#pragma omp for private(i, k, dpsi_borderRight)  schedule(runtime)
        for (i = dim1 - (thicknessRight-HALFORDERDF); i < dim1-HALFORDERDF; i++) {
		#pragma omp simd
                #pragma vector aligned
                for (k = HALFORDERDF; k < dim2-HALFORDERDF; k++) {

                        dpsi_borderRight = (
                                B48*(psi_borderRight[(dim2)*(i-(dimThickR)+4) + (k)] - psi_borderRight[(dim2)*(i-(dimThickR)-4) + (k)]) +
                                B38*(psi_borderRight[(dim2)*(i-(dimThickR)+3) + (k)] - psi_borderRight[(dim2)*(i-(dimThickR)-3) + (k)]) +
                                B28*(psi_borderRight[(dim2)*(i-(dimThickR)+2) + (k)] - psi_borderRight[(dim2)*(i-(dimThickR)-2) + (k)]) +
                                B18*(psi_borderRight[(dim2)*(i-(dimThickR)+1) + (k)] - psi_borderRight[(dim2)*(i-(dimThickR)-1) + (k)])
                                ) * h1;

                        csi_borderRight[(dim2)*(i-(dimThickR)) + (k)] =
                                b_right[i-(dimThickR)] * csi_borderRight[(dim2)*(i-(dimThickR)) + (k)] +
                                a_right[i-(dimThickR)] * ((
                                                        AC48*(waveField1[(dim2)* (i+4) + (k)] + waveField1[(dim2)*(i-4) + (k)]) +
                                                        AC38*(waveField1[(dim2)* (i+3) + (k)] + waveField1[(dim2)*(i-3) + (k)]) +
                                                        AC28*(waveField1[(dim2)* (i+2) + (k)] + waveField1[(dim2)*(i-2) + (k)]) +
                                                        AC18*(waveField1[(dim2)* (i+1) + (k)] + waveField1[(dim2)*(i-1) + (k)]) +
                                                        AC08*waveField1[(dim2)*(i) + (k)]
                                                        ) * h3 + dpsi_borderRight);

                        waveField2[(dim2)*(i) + (k)] = waveField2[(dim2)*(i) + (k)] + (Coef[(dim2)*(i) + (k)]*h2)*(dpsi_borderRight + csi_borderRight[(dim2)*(i-(dimThickR)) + (k)]);

                }
        }

        // ............................. Evaluating the auxiliary variable of the Bottom Border ............................. //
        #pragma omp for private(i,k)  schedule(runtime)
        for (i = HALFORDERDF; i < dim1-HALFORDERDF; i++) {
		#pragma omp simd
                #pragma vector aligned
                for (k = dim2-(thicknessDown-HALFORDERDF); k < dim2-HALFORDERDF; k++) {

                        psi_borderBottom[(thicknessDown)*(i) + (k-(dimThickD))] =
                                b_bottom[k-(dimThickD)] * psi_borderBottom[(thicknessDown)*(i) + (k-(dimThickD))] +
                                a_bottom[k-(dimThickD)] * (
                                        B48*(waveField1[(dim2)*(i) + (k+4)] - waveField1[(dim2)*(i) + (k-4)]) +
                                        B38*(waveField1[(dim2)*(i) + (k+3)] - waveField1[(dim2)*(i) + (k-3)]) +
                                        B28*(waveField1[(dim2)*(i) + (k+2)] - waveField1[(dim2)*(i) + (k-2)]) +
                                        B18*(waveField1[(dim2)*(i) + (k+1)] - waveField1[(dim2)*(i) + (k-1)])
                                        ) * h1;
                }
        }


        // ............................. Addition of the corrective terms of the Bottom Border .............................. //
        #pragma omp for private(i, k, dpsi_borderBottom)  schedule(runtime)
        for (i = HALFORDERDF; i < dim1-HALFORDERDF; i++) {
		#pragma omp simd
                #pragma vector aligned
                for (k = dim2-(thicknessDown-HALFORDERDF); k < dim2-HALFORDERDF; k++) {

                        dpsi_borderBottom = (
                                B48*(psi_borderBottom[(thicknessDown)*(i) + (k-(dimThickD)+4)] - psi_borderBottom[(thicknessDown)*(i) + (k-(dimThickD)-4)]) +
                                B38*(psi_borderBottom[(thicknessDown)*(i) + (k-(dimThickD)+3)] - psi_borderBottom[(thicknessDown)*(i) + (k-(dimThickD)-3)]) +
                                B28*(psi_borderBottom[(thicknessDown)*(i) + (k-(dimThickD)+2)] - psi_borderBottom[(thicknessDown)*(i) + (k-(dimThickD)-2)]) +
                                B18*(psi_borderBottom[(thicknessDown)*(i) + (k-(dimThickD)+1)] - psi_borderBottom[(thicknessDown)*(i) + (k-(dimThickD)-1)])
                                ) * h1;

                        csi_borderBottom[(thicknessDown)*(i) + (k-(dimThickD))] =
                                b_bottom[k-(dimThickD)] * csi_borderBottom[(thicknessDown)*(i) + (k-(dimThickD))] +
                                a_bottom[k-(dimThickD)] * ((
                                                        AC48*(waveField1[(dim2)* (i) + (k+4)] + waveField1[(dim2)*(i) + (k-4)]) +
                                                        AC38*(waveField1[(dim2)* (i) + (k+3)] + waveField1[(dim2)*(i) + (k-3)]) +
                                                        AC28*(waveField1[(dim2)* (i) + (k+2)] + waveField1[(dim2)*(i) + (k-2)]) +
                                                        AC18*(waveField1[(dim2)* (i) + (k+1)] + waveField1[(dim2)*(i) + (k-1)]) +
                                                        AC08*waveField1[(dim2)*(i) + (k)]
                                                        ) * h3 + dpsi_borderBottom);

                        waveField2[(dim2)*(i) + (k)] = waveField2[(dim2)*(i) + (k)] + (Coef[(dim2)*(i) + (k)]*h2)*(dpsi_borderBottom + csi_borderBottom[(thicknessDown)*(i) + (k-(dimThickD))]);
                }
        }

	}
	
}