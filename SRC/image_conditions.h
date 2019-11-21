#ifndef IMAGE_CONDITIONS_H
#define IMAGE_CONDITIONS_H

void cross_correlation(float *imageCondition, float *sourceWavefield, float *receiverWavefield);

void normalized_cross_correlation(float *imageCondition, float *selfImageCondition, float *sourceWavefield, float *receiverWavefield);

#endif