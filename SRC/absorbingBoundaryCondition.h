#ifndef PI_ABSORBINGBOUNDARYCONDITION_H
#define PI_ABSORBINGBOUNDARYCONDITION_H


#define NPOWER 2.0f

// Discretization order for FDM (order/2)
#define ORDERDF 8

// Half of discretization order for FDM (order/2)
#define HALFORDERDF 4

// Second Convolutional's derivative coefficients 8 order
#define AC08 -2.97399944f
#define AC18 +1.70507669f
#define AC28 -0.25861812f
#define AC38 +0.04577745f
#define AC48 -0.00523630f

// Frist Taylor's derivative coefficients 8 order
#define B18 +0.8f
#define B28 -0.2f
#define B38 +0.025396825f
#define B48 -0.001785714f


void attenuation_profile_cpml(float *a_Left, float *a_Right, float *a_Top, float *a_Bottom, float *b_Left, float *b_Right, float *b_Top, float *b_Bottom);

void aiso_attenuation_cpml(float *waveField1, float *waveField2, float *Coef, float *a_left, float *b_left, float *a_right, float *b_right, float *a_top, float *b_top, float *a_bottom, float *b_bottom, float *psi_borderLeft, float *psi_borderRight, float *psi_borderTop, float *psi_borderBottom, float *csi_borderLeft, float *csi_borderRight, float *csi_borderTop, float *csi_borderBottom);

#endif
