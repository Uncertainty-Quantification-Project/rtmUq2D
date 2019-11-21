#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "pi_file.h"



FILE *file;

int reading_file(char *path, float *matrix, int dim1, int dim2) {

	int nR;

	file = fopen (path,"rb");

	if (file == NULL) {
		printf("\nError oppening file (fopen in reading_file function)\n");
                exit(1);
	}


	nR = fread (matrix, sizeof(float), dim1*dim2, file);

	if (nR != dim1*dim2) {
		printf("\nError reading file (fread in reading_file function)\n");
		exit(1);
	}

	fclose(file);
	file = NULL;


	return 0;
}

int recording_file(char *path, float *matrix, int dim1, int dim2) {

	int nW;

	file = fopen (path,"wb");

	if (file == NULL) {
		printf("\nError oppening file (fopen in recording_file function)\n");
		exit(1);
	}


	nW = fwrite (matrix, sizeof(float), dim1*dim2, file);

        if (nW != dim1*dim2) {
		printf("\nError reading file (fread in recording_file function)\n");
		exit(1);
	}

	fclose(file);
	file = NULL;


	return 0;
}