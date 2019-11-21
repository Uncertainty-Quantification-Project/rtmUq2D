#include <stdio.h>
#include "fixed_parameters.h"
#include "generic_functions.h"

void extend_matrix(float *extendedModel, float *model) {

	int i, k;
	
// ......................................  ............................................//
	for (i = 0; i < Nx; i++) {
		for (k = 0; k < Nz; k++) {

			extendedModel[(k+thicknessUp) + (i+thicknessLeft)*(Nz+thicknessUp+thicknessDown)] = model[k + i*Nz];
		}
	}

// ............................. Thickness Left ............................................ //
	for (i = 0; i < thicknessLeft; i++) {
		for (k = thicknessUp; k < (Nz+thicknessUp); k++) {

			extendedModel[k + i*(Nz+thicknessUp+thicknessDown)] = extendedModel[k + thicknessLeft*(Nz+thicknessUp+thicknessDown)];
		}
	}

// ............................. Thickness Right ............................................ //
	for (i = (Nx+thicknessLeft); i < (Nx+thicknessLeft+thicknessRight); i++) {
		for (k = thicknessUp; k < (Nz+thicknessDown); k++) {

			extendedModel[k + i*(Nz+thicknessUp+thicknessDown)] = extendedModel[k + (Nx+thicknessLeft-1)*(Nz+thicknessUp+thicknessDown)];
		}
        }
        
// ............................. Thickness Up ........................................... //
	  for (i = 0; i < (Nx+thicknessLeft+thicknessRight); i++) {
		for (k = 0; k < thicknessUp; k++) {

			extendedModel[k + i*(Nz+thicknessUp+thicknessDown)] = extendedModel[thicknessUp + i*(Nz+thicknessUp+thicknessDown)];
		}
	}

// ............................. Thickness Down ........................................... //
	for (i = 0; i < (Nx+thicknessLeft+thicknessRight); i++) {
		for (k = (Nz+thicknessUp); k < Nz+thicknessUp+thicknessDown; k++) {

			extendedModel[k + i*(Nz+thicknessUp+thicknessDown)] = extendedModel[(Nz+thicknessUp-1) + i*(Nz+thicknessUp+thicknessDown)];
		}
	}
}


void cleaningImage(int waterLayer, float *stack) {
	
	int i, k;
	
	for (i = 0; i < Nx; i++) {
		for (k = 0; k < waterLayer; k++) {
			
			stack[i*Nz + k] = 0.0;
		}
	}
}