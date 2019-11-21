#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#include "zfp.h"


/* Compress Array */
int zfp2d_lossyCompression(void** compression, float* array, int dimension1, int dimension2, float tolerance) {
    
    FILE *pFile;
    
    
    uint num_threads_zfp;
    uint set_num_threads_zfp = 0;
    int OMP_SUCCESS, POLICY_SUPPORTED;
    
//     zfp_exec_policy exec = zfp_exec_serial;
    zfp_exec_policy exec = zfp_exec_omp;
    
    int status = 0;    /* return value: 0 = success */
    zfp_type type;     /* array scalar type */
    zfp_field* field;  /* array meta data */
    zfp_stream* zfp;   /* compressed stream */
    void* buffer;      /* storage for compressed stream */
    size_t bufsize;    /* byte size of compressed buffer */
    bitstream* stream; /* bit stream to write to or read from */
    size_t zfpsize;    /* byte size of compressed stream */

    /* allocate meta data for the 3D array a[nz][ny][nx] */
    type = zfp_type_float;
    field = zfp_field_2d(&array[0], type, dimension1, dimension2);

    /* allocate meta data for a compressed stream */
    zfp = zfp_stream_open(NULL);
    
    /* set compression mode and parameters via one of three functions */
    /*  zfp_stream_set_rate(zfp, rate, type, 3, 0); */
    /*  zfp_stream_set_precision(zfp, precision); */
    zfp_stream_set_accuracy(zfp, tolerance);
    
    /* allocate buffer for compressed data */
    bufsize = zfp_stream_maximum_size(zfp, field);
    buffer = malloc(bufsize);
    
    /* associate bit stream with allocated buffer */
    stream = stream_open(buffer, bufsize);
    zfp_stream_set_bit_stream(zfp, stream);
    zfp_stream_rewind(zfp);
    
    
    if ( zfp_stream_set_execution(zfp,  exec) ) {
        
        OMP_SUCCESS     = zfp_stream_set_omp_threads(zfp, set_num_threads_zfp);
        
//         num_threads_zfp = zfp_stream_omp_threads(zfp);
//         printf("Thread ID = %d\n", omp_get_thread_num());
//         printf("\n%d  |  %u\n", OMP_SUCCESS, num_threads_zfp);
            
        zfpsize = zfp_compress(zfp, field);
    
    }
        
    *compression = buffer;

    /* clean up */
    zfp_field_free(field);
    zfp_stream_close(zfp);
    stream_close(stream);
    buffer = NULL;

    return zfpsize;
}


/* Decompress Array */
int zfp2d_lossyDecompression(char *zfp_data_02, float* array, int dimension1, int dimension2, float tolerance) {
    
    FILE *pFile;
    
    int nW;
    int status = 0;    /* return value: 0 = success */
    zfp_type type;     /* array scalar type */
    zfp_field* field;  /* array meta data */
    zfp_stream* zfp;   /* compressed stream */
    void* buffer;      /* storage for compressed stream */
    size_t bufsize;    /* byte size of compressed buffer */
    bitstream* stream; /* bit stream to write to or read from */
    size_t zfpsize;    /* byte size of compressed stream */

    /* allocate meta data for the 3D array a[nz][ny][nx] */
    type = zfp_type_float;
    field = zfp_field_2d(&array[0], type, dimension1, dimension2);

    /* allocate meta data for a compressed stream */
    zfp = zfp_stream_open(NULL);
    
    /* set compression mode and parameters via one of three functions */
    /*  zfp_stream_set_rate(zfp, rate, type, 3, 0); */
    /*  zfp_stream_set_precision(zfp, precision); */
    zfp_stream_set_accuracy(zfp, tolerance);

    /* allocate buffer for compressed data */
    bufsize = zfp_stream_maximum_size(zfp, field);
    buffer = malloc(bufsize);

    /* associate bit stream with allocated buffer */
    stream = stream_open(buffer, bufsize);
    zfp_stream_set_bit_stream(zfp, stream);
    zfp_stream_rewind(zfp);
    
    /* read compressed stream and decompress array */
//     zfpsize = fread(buffer, 1, bufsize, stdin);
    pFile = fopen (zfp_data_02, "rb");
    fread(buffer, 1, bufsize, pFile);
    fclose (pFile);
    
    zfpsize = zfp_decompress(zfp, field);
//     printf("\n\n %lu", zfpsize);
    
    if (!zfpsize) {
      fprintf(stderr, "Decompression Failed\n");
      status = 1;
    }
    

    /* clean up */
    zfp_field_free(field);
    zfp_stream_close(zfp);
    stream_close(stream);
    free(buffer);

  return status;
}