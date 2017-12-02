/**
*
* kk.c
* 
* Ryan Wallace	
* Computer Science 124
* Harvard University
* Programming Assingment 3
* April 23, 2016
*
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <stdbool.h>
#include <limits.h>

#define LONGLONG_LEN 21

int parent (int i) {
	return (i/2);
}

int left (int i) {
	return 2 * i;
}

int right (int i) {
	return 2 * i + 1;
}

void heapify (unsigned long long* H, int N, int* size) {
	unsigned long long l = left (N);
	unsigned long long r = right (N);
	int largest;

	if (l <= *size && H[l] > H[N]) {
		largest = l;
	}
	else {
		largest = N;
	}

	if (r <= *size && H[r] > H[largest]) {
		largest = r;
	}

	if (largest != N) {
		unsigned long long buffer = H[largest];
		H[largest] = H[N];
		H[N] = buffer;
		heapify (H, largest, size);
	}

	return;
}

void insert (unsigned long long* H, unsigned long long v, int* size) {
	*size += 1;
	H[*size] = v;
	int N = *size;
	while (N != 1 && H[parent (N)] < H[N]) {
		unsigned long long buffer = H[parent (N)];
		H[parent (N)] = H[N];
		H[N] = buffer;
		N = parent (N);
	}
	return;
}

unsigned long long extract (unsigned long long* H, int* size) {
	unsigned long long max = H[1];
	H[1] = H[*size];
	(*size)--;
	heapify (H, 1, size);
	return max;
}

unsigned long long KK (unsigned long long* A, int lines) {
	// construct the Max-Heap
	unsigned long long* H = calloc (lines + 1, sizeof (unsigned long long));
	int* size = malloc (sizeof (int));
	*size = 0;
	for (int i = 1; i <= lines; i++) {
		insert (H, A[i], size);
	}

	// perform the Karmarkar-Karp algorithm
	while (*size > 1) {
		unsigned long long max1 = extract (H, size);
		unsigned long long max2 = extract (H, size);
		insert (H, max1 - max2, size);
	}

	return (extract (H, size));
}

int main (int argc, char** argv) {
	if (argc != 2) {
		printf("use: ./kk inputfile\n");
		return 1;
	}

	// reads in input file to A
	int lines = 100;
	char filename[512];
	strcpy (filename, argv[1]);

   	FILE *fp = fopen(filename, "r");
	if (fp == NULL) {
		printf("error opening inputfile\n");
		return 1;
	}
	
	unsigned long long A[lines + 1];
	for (int line = 1; line <= lines; line++) {
		char buffer[LONGLONG_LEN];
		if (fgets (buffer, sizeof(buffer), fp) != NULL) {
			A[line] = atoll(buffer);
		}
		else {
		 	printf ("error reading inputfile\n");
		 	return 1;
		}
	}

	// Karmarkar-Karp
	unsigned long long R = KK (A, lines);
	printf("Residue: %llu\n", R);

	return 0;
}

