#ifndef FORWARD_MODELING_H
#define FORWARD_MODELING_H

#define TIME_SLICES 20
#define VEC_LEN 1000

// -------- Temporary FILES...
static char sourceWavefieldFile[] = "./TEMPORARY_FILES/sourceWavefield_%05d.bin";
// ----------------------------------

void isotropicAcousticModeling(int simId, int *shot, float *velocity, float *wavelet);

#endif