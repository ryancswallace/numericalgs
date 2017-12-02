/**
 * 
 * strassen.c
 * 
 * Walter Martin
 * Ryan Wallace
 * Computer Science 124
 * Harvard University
 * Programming Assingment 2
 * 
 */

 /** TODO
 *
 * written report
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#define INT_LEN 12

// generates random integer inputfile
void genRand (int dimension) {
	FILE *fp;
	fp = fopen("randominput", "w");
	if (fp == NULL) {
		printf("error creating randominput file\n");
	}

	// random numbers are on [min, max]
	int min = 0;
	int max = 1;
	int lines = 2 * dimension * dimension;
	srand((int)time(NULL));

	for (int i = 0; i < lines; i++) {
		fprintf(fp, "%d\n", rand() % (max - min + 1 ) + min);
	}

	fclose(fp);
}

// returns allocated memory for a rows x columns matrix
int** getMatrix (int rows, int columns) {
    // allocate an array of row pointers to pointers to doubles
    int** matrixMemory = calloc (rows, sizeof(int*));
    if (matrixMemory == NULL) {
        return NULL;
    }
    
    // for each pointer in the array
    for (int i = 0; i < rows; i++) {
        // allocate a pointer to a double
        matrixMemory[i] = calloc (columns, sizeof(int*));
        if (matrixMemory[i] == NULL) {
            return NULL;
        }
    }
    
    // return the allocated for the matrix
    return (matrixMemory);
}

// multiplies the submatrix A from rows AR1 to AR2 and columns AC1 to AC2 with
// submatirx B from rows BR1 to BR2 and columns BC1 to BC2 
int** matrixMultiply (int** A, int** B, int AR1, int AR2, int AC1, int AC2, int BR1, int BR2, int BC1, int BC2) {
    // allocate memory for the product
    int** product = getMatrix (AR2 - AR1, BC2 - BC1);
    if (product == NULL) {
        return NULL;
    }

    // compute the product
    for (int i = 0; AR1 + i < AR2; i++) {
        for (int j = 0; BC1 + j < BC2; j++) {
            for (int c = 0; AR1 + c < AR2; c++) {
                product[i][j] += (A[AR1 + i][AC1 + c] * B[BR1 + c][BC1 + j]);
            }
        }
    }
    
    // return the the product
    return (product);
}

int** matrixAddMem (int** A, int** B, int** C, int AR1, int AR2, int AC1, int AC2, int BR1, int BR2, int BC1, int BC2) {
	// compute the sum
    for (int i = 0; AR1 + i < AR2; i++) {
        for (int j = 0; BC1 + j < BC2; j++) {
            C[i][j] = A[AR1 + i][AC1 + j] + B[BR1 + i][BC1 + j];
        }
    }
    
    // return the the sum
    return (C);
}

int** matrixSubMem (int** A, int** B, int** C, int AR1, int AR2, int AC1, int AC2, int BR1, int BR2, int BC1, int BC2) {
    // compute the difference
    for (int i = 0; AR1 + i < AR2; i++) {
        for (int j = 0; BC1 + j < BC2; j++) {
            C[i][j] = A[AR1 + i][AC1 + j] - B[BR1 + i][BC1 + j];
        }
    }

    // return the the difference
    return (C);
}

int** strassen_mix (int** A, int** B, int dim, int threshold, int AR1, int AR2, int AC1, int AC2, int BR1, int BR2, int BC1, int BC2) {
	// error checking
	if (threshold < 1) {
		printf("threshold must be >0\n");
		return NULL;
	}
	if (dim <= threshold) {
		int** result = getMatrix(dim, dim);
		result = matrixMultiply (A, B, AR1, AR2, AC1, AC2, BR1, BR2, BC1, BC2);
		return result;
	}

	int** result = getMatrix(AR2 - AR1, AC2 - AC1);

	dim = AR2 - AR1;
	int halfdim = dim / 2;

	// buffers for MatrixAddMem and MatrixSubMem
	int** C = getMatrix (AR2 - AR1, AC2 - AC1);
    if (C == NULL) {
        return NULL;
    }

    int** D = getMatrix (AR2 - AR1, AC2 - AC1);
    if (D == NULL) {
        return NULL;
    }

	// F - H
	int** FmH = matrixSubMem (B, B, C, BR1, halfdim + BR1, halfdim + BC1, dim + BC1, halfdim + BR1, dim + BR1, halfdim + BC1, dim + BC1);

	// P1
	int** P1 = strassen_mix (A, FmH, halfdim, threshold, AR1, halfdim + AR1, AC1, halfdim + AC1, 0, halfdim, 0, halfdim);

	// A + B
	int** ApB = matrixAddMem (A, A, C, AR1, halfdim + AR1, AC1, halfdim + AC1, AR1, halfdim + AR1, halfdim + AC1, dim + AC1);

	// P2
	int** P2 = strassen_mix (ApB, B, halfdim, threshold, 0, halfdim, 0, halfdim, halfdim + BR1, dim + BR1, halfdim + BC1, dim + BC1);

	// C + D
	int** CpD = matrixAddMem (A, A, C, halfdim + AR1, dim + AR1, AC1, halfdim + AC1, halfdim + AR1, dim + AR1, halfdim + AC1, dim + AC1);

	// P3
	int** P3 = strassen_mix (CpD, B, halfdim, threshold, 0, halfdim, 0, halfdim, BR1, halfdim + BR1, BC1, halfdim + BC1);

	// G - E
	int** GmE = matrixSubMem (B, B, C, halfdim + BR1, dim + BR1, BC1, halfdim + BC1, BR1, halfdim + BR1, BC1, halfdim + BC1);

	// P4
	int** P4 = strassen_mix (A, GmE, halfdim, threshold, halfdim + AR1, dim + AR1, halfdim + AC1, dim + AC1, 0, halfdim, 0, halfdim);

	// A + D
	int** ApD = matrixAddMem (A, A, C, AR1, halfdim + AR1, AC1, halfdim + AC1, halfdim + AR1, dim + AR1, halfdim + AC1, dim + AC1);

	// E + H
	int** EpH = matrixAddMem (B, B, D, BR1, halfdim + BR1, BC1, halfdim + BC1, halfdim + BR1, dim + BR1, halfdim + BC1, dim + BC1);

	// P5
	int** P5 = strassen_mix (ApD, EpH, halfdim, threshold, 0, halfdim, 0, halfdim, 0, halfdim, 0, halfdim);

	// B - D
	int** BmD = matrixSubMem (A, A, C, AR1, halfdim + AR1, halfdim + AC1, dim + AC1, halfdim + AR1, dim + AR1, halfdim + AC1, dim + AC1);

	// G + H
	int** GpH = matrixAddMem (B, B, D, halfdim + BR1, dim + BR1, BC1, halfdim + BC1, halfdim + BR1, dim + BR1, halfdim + BC1, dim + BC1);

	// P6
	int** P6 = strassen_mix (BmD, GpH, halfdim, threshold, 0, halfdim, 0, halfdim, 0, halfdim, 0, halfdim);

	// A - C
	int** AmC = matrixSubMem (A, A, C, AR1, halfdim + AR1, AC1, halfdim + AC1, halfdim + AR1, dim + AR1, AC1, halfdim + AC1);

	// E + F
	int** EpF = matrixAddMem (B, B, D, BR1, halfdim + BR1, BC1, halfdim + BC1, BR1, halfdim + BR1, halfdim + BC1, dim + BC1);

	// P7
	int** P7 = strassen_mix (AmC, EpF, halfdim, threshold, 0, halfdim, 0, halfdim, 0, halfdim, 0, halfdim);

	// copy Pi values into result matrix
	for(int i = 0; i < dim; i++) {
		for(int j = 0; j < dim; j++) {
			if(i < halfdim) {
				if(j < halfdim) {
					result[i][j] = P5[i][j] + P4[i][j] - P2[i][j] + P6[i][j];
				}
				else {
					result[i][j] = P1[i][j - halfdim] + P2[i][j - halfdim];
				}
			}
			else {
				if(j < halfdim) {
					result[i][j] = P3[i - halfdim][j] + P4[i - halfdim][j];
				}
				else {
					result[i][j] = P5[i - halfdim][j - halfdim] + P1[i - halfdim][j - halfdim] - P3[i - halfdim][j - halfdim] - P7[i - halfdim][j - halfdim];
				}
			}
		}
	}
	return result;
}

	// argument validation
	if (argc != 4) {
		printf("Use: ./strassen random dimension inputfile\n");
		return 1;
	}

	int dimension = atoi(argv[2]);
	if (dimension < 1) {
		printf("Use: ./strassen random dimension inputfile\n");
		return 1;
	}

	int random = atoi(argv[1]);
	if (random == 1) {
		if (strcmp(argv[3], "randominput") != 0) {
			printf("inputfile must be \"randominput\" when random is 1\n");
			return 1;
		}
		genRand (dimension);
	}
	else if (random != 0) {
		printf("Use: ./strassen random dimension inputfile\n");
		return 1;
	}

	// open the input file
	FILE *fp;
	fp = fopen(argv[3], "r");
	if (fp == NULL) {
		printf("error opening inputfile\n");
	}
	
	// cutoff is slightly dependent on matrix size
	int cutoff = 32;
	if(dimension > 256) {
		cutoff = 64;
	}
	
	// determine padding, enough so that every division above cutoff is even
	int power = floor(log2(dimension));
	if(dimension % ((int) pow(2, power)) != 0) {
		power++;
	}
	int pad_dim = (int) pow(2, power);
	int temp_power = power;
	for(int i = 3; i <= cutoff; i++) {
		while(i * ((int) pow(2, temp_power)) > dimension) {
			temp_power--;
		}
		temp_power++;
	
		if(i * ((int) pow(2, temp_power)) < pad_dim) {
			pad_dim = i * ((int) pow(2, temp_power));
		}
	}
	
	int padding = pad_dim - dimension;

	// allocate memory for two input matrices, A and B, and product C
	int** A = getMatrix (pad_dim, pad_dim);
	int** B = getMatrix (pad_dim, pad_dim);
	int** C = getMatrix (pad_dim, pad_dim);

	// read input into A and B matrices
	int j = 0;
	for (int i = 0; i < dimension; i++) {
		for (j = 0; j < dimension; j++) {
			char buffer[12];
			if (fgets (buffer, 12, fp) != NULL) {
				A[i][j] = atoi(buffer);
		    }
		    else {
			    printf("dimension and file size mismatch\n");
				return 1;
		    }
		}
	}

	for (int i = 0; i < dimension; i++) {
		for (j = 0; j < dimension; j++) {
			char buffer[INT_LEN];
			if (fgets (buffer, INT_LEN, fp) != NULL) {
		    	B[i][j] = atoi(buffer);
		    }
		    else {
		    	printf("dimension and file size mismatch\n");
				return 1;
		    }
		}
	}

	// establish timing
	clock_t start, end;
	double time_elapsed;

	// multiplication
	start = clock();
	// C = matrixMultiply(A, B, 0, pad_dim, 0, pad_dim, 0, pad_dim, 0, pad_dim);
	C = strassen_mix(A, B, pad_dim, cutoff, 0, pad_dim, 0, pad_dim, 0, pad_dim, 0, pad_dim);
	if (C == NULL) {
		printf("error in strassen calculation\n");
		return 1;
	}
	end = clock();
	time_elapsed = ((double) (end - start)) / CLOCKS_PER_SEC;

	// close input file
	fclose(fp);
	
	// prints the diagonal entries of a matrix as specified in assignment
	for (int i = 0; i < dimension; i++) {
		 printf("%d\n", C[i][i]);
	}
	printf("\n");

	// free allocated memory
	for (int i = 0; i < pad_dim; i++) {
    	free (A[i]);
    	free (B[i]);
    	free (C[i]);
    }
    free (A);
    free (B);
    free (C);
	
	return 0;
}