#ifndef BACKWARD_MODELING_H
#define BACKWARD_MODELING_H

#define TIME_SLICES 20
#define VEC_LEN 1000

// -------- Temporary FILES...
static char sourceWavefieldFile2[] = "./TEMPORARY_FILES/sourceWavefield_%05d.bin";
// ----------------------------------

int adjointModeling(int simId, int *shot, float *velocity, float *wavelet, float *seismogram, float *imageCondition, float *selfImageCondition);

#endif