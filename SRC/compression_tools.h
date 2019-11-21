#ifndef COMPRESSION_TOOLS_H
#define COMPRESSION_TOOLS_H

/* Compress Array with ZFP Library*/
int zfp2d_lossyCompression(void** compression, float* array, int dimension1, int dimension2, float tolerance);


/* Decompress Array with ZFP Library*/
int zfp2d_lossyDecompression(char *zfp_data_02, float* array, int dimension1, int dimension2, float tolerance);

#endif