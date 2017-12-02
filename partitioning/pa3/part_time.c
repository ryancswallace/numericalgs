/**
*
* part_time.c
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
#include "mt64.h"

#define LONGLONG_LEN 21
#define E 2.718281828459

#define TRIAL 1
#define MAX_ITER 25000
#define MAX_TIME 1

// establish timing
clock_t start, end;
double kk_time = 0, rrb_time = 0, hcb_time = 0, sab_time = 0, rrp_time = 0, hcp_time = 0, sap_time = 0;

// size of problem
static int lines = 100, files = 50;

// generates random integer inputfile
void genRand (int lines, int files) {
	init_genrand64((unsigned long long)time(NULL));
	unsigned long long min = 1;
	unsigned long long max = pow (10, 12);
	
	for (int file = 0; file < files; file++) {
		FILE *fp;
		char filename[sizeof ("100.txt")];
		sprintf(filename, "%d.txt", file);
    	fp = fopen(filename, "w");
		
		if (fp == NULL) {
			printf("error creating randominput file\n");
		}

		for (int i = 0; i < lines; i++) {
			unsigned long long r = genrand64_int64 ();
			unsigned long long mod = r % max;
			fprintf(fp, "%llu\n", mod + min);
		}

		fclose(fp);
	}
}

int parent (int i) {
	return (i/2);
}

int left (int i) {
	return 2 * i;
}

int right (int i) {
	return 2 * i + 1;
}

void heapify (long long* H, int N, int* size) {
	long long l = left (N);
	long long r = right (N);
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
		long long buffer = H[largest];
		H[largest] = H[N];
		H[N] = buffer;
		heapify (H, largest, size);
	}

	return;
}

void insert (long long* H, long long v, int* size) {
	*size += 1;
	H[*size] = v;
	int N = *size;
	while (N != 1 && H[parent (N)] < H[N]) {
		long long buffer = H[parent (N)];
		H[parent (N)] = H[N];
		H[N] = buffer;
		N = parent (N);
	}
	return;
}

long long extract (long long* H, int* size) {
	long long max = H[1];
	H[1] = H[*size];
	(*size)--;
	heapify (H, 1, size);
	return max;
}

short* rand_bits () {
	short* bits = calloc (lines, sizeof(short) + 1);
	for (int i = 1; i <= lines; i++) {
		if (rand() % 2) {
			bits[i] = 1;
		}
		else {
			bits[i] = -1;
		}
	}

	return bits;
}

short* neighbor_bits (short* B) {
	// get two random distinct indices
	short i = (rand() % lines) + 1;
	short j = (rand() % lines) + 1;
	while (i == j) {
		j = (rand() % lines) + 1;
	}

	short* Bp = calloc (lines, sizeof(short) + 1);
	for (int c = 1; c <= lines; c++) {
		Bp[c] = B[c];
	}

	Bp[i] = -Bp[i];

	if (rand() % 2) {
		B[j] = -B[j];
	}

	return Bp;
}

short* rand_part () {
	short* part = calloc (lines, sizeof(short) + 1);
	for (int i = 1; i <= lines; i++) {
		short part_v = (rand() % lines) + 1;
			part[i] = part_v;
		}

	return part;
}

short* neighbor_part (short* P) {
	short i = (rand() % lines) + 1;
	short j = (rand() % lines) + 1;
	while (i == j) {
		j = (rand() % lines) + 1;
	}

	short* Pp = calloc (lines, sizeof(short) + 1);
	for (int c = 1; c <= lines; c++) {
		Pp[c] = P[c];
	}

	Pp[i] = j;

	return Pp;
}

long long* part_to_arr (long long* A, short* part) {
	long long* arr = calloc (lines, sizeof (long long) + 1);
	for (int i = 1; i <= lines; i++) {
		arr[i] = A[i];
	}

	for (short p = 1; p <= lines; p++) {
		long long sum = 0;
		bool mem[lines + 1];
		for (short i = 1; i <= lines; i++) {
			if (p == part[i]) {
				sum += A[i];
				mem[i] = true;
			}
			else {
				mem[i] = false;
			}
		}

		bool first = true;
		for (int i = 1; i <= lines; i++) {
			if (mem[i] == true) {
				if (first) {
					arr[i] = sum;
					first = false;
				} 
				else {
					arr[i] = 0;
				}
			}
		}
	}

	return arr;
}

long long kk (long long* A) {
	start = clock();

	// construct the Max-Heap
	long long* H = calloc (lines + 1, sizeof (long long));
	int* size = malloc (sizeof (int));
	*size = 0;
	for (int i = 1; i <= lines; i++) {
		insert (H, A[i], size);
	}

	// perform the Karmarkar-Karp algorithm
	while (*size > 1) {
		long long max1 = extract (H, size);
		long long max2 = extract (H, size);
		insert (H, max1 - max2, size);
	}

	long long result = extract (H, size);
	long long res = result;

	end = clock();
	kk_time += ((double) (end - start)) / CLOCKS_PER_SEC;

	return (res);
}

long long kk_p (long long* A) {
	// construct the Max-Heap
	long long* H = calloc (lines + 1, sizeof (long long));
	int* size = malloc (sizeof (int));
	*size = 0;
	for (int i = 1; i <= lines; i++) {
		insert (H, A[i], size);
	}

	// perform the Karmarkar-Karp algorithm
	while (*size > 1) {
		long long max1 = extract (H, size);
		long long max2 = extract (H, size);
		insert (H, max1 - max2, size);
	}

	long long res = extract (H, size);
	
	return (res);
}

long long rrb (long long* A) {
	rrb_time = 0;
	start = clock();

	short* B = rand_bits ();
	long long S = 0;
	for (int i = 1; i <= lines; i++) {
		long long Si = A[i];
		if (B[i] < 0) {
			Si = -Si;
		}
		S += Si;
	}
	S = llabs(S);

	end = clock();
	rrb_time += ((double) (end - start)) / CLOCKS_PER_SEC;

	while (rrb_time < MAX_TIME) {
		start = clock();

		short* B = rand_bits ();
		long long Sp = 0;
		for (int i = 1; i <= lines; i++) {
			long long Si = A[i];
			if (B[i] < 0) {
				Si = -Si;
			}
			Sp += Si;
		}
		Sp = llabs(Sp);

		if (Sp < S) {
			S = Sp;
		}

		end = clock();
		rrb_time += ((double) (end - start)) / CLOCKS_PER_SEC;
	}
	rrb_time = 0;
	return (S);
}

long long hcb (long long* A) {
	hcb_time = 0;
	start = clock();

	short* B = rand_bits ();
	long long S = 0;
	for (int i = 1; i <= lines; i++) {
		long long Si = A[i];
		if (B[i] < 0) {
			Si = -Si;
		}
		S += Si;
	}
	S = llabs(S);
	
	end = clock();
	hcb_time += ((double) (end - start)) / CLOCKS_PER_SEC;

	while (hcb_time < MAX_TIME) {
		start = clock();

		short* Bp = neighbor_bits (B);
		long long Sp = 0;
		for (int i = 1; i <= lines; i++) {
			long long Si = A[i];
			if (B[i] < 0) {
				Si = -Si;
			}
			Sp += Si;
		}
		Sp = llabs(Sp);

		if (Sp < S) {
			S = Sp;
			B = Bp;
		}

		end = clock();
		hcb_time += ((double) (end - start)) / CLOCKS_PER_SEC;
	}
	hcb_time = 0;
	return (S);
}

long long sab (long long* A) {
	sab_time = 0;
	start = clock();

	short* B = rand_bits ();
	long long S = 0;
	for (int i = 1; i <= lines; i++) {
		long long Si = A[i];
		if (B[i] < 0) {
			Si = -Si;
		}
		S += Si;
	}
	S = llabs(S);
	long long Spp = S;

	end = clock();
	sab_time += ((double) (end - start)) / CLOCKS_PER_SEC;
	
	long long iteration = 0;
	while (sab_time < MAX_TIME) {
		start = clock ();

		short* Bp = neighbor_bits (B);
		long long Sp = 0;
		for (int i = 1; i <= lines; i++) {
			long long Si = A[i];
			if (B[i] < 0) {
				Si = -Si;
			}
			Sp += Si;
		}
		Sp = llabs(Sp);

		if (Sp < S) {
			S = Sp;
			B = Bp;
		}
		else {
			double T = pow (10, 10) * pow (0.8, iteration / 300);
			double p = pow (E, -(Sp - S)/T);
			if ((rand() / (double) RAND_MAX) < p) {
				S = Sp;
				B = Bp;
			}
		}

		if (S < Spp) {
			Spp = S;
		}

		iteration++;
		end = clock();
		sab_time += ((double) (end - start)) / CLOCKS_PER_SEC;
	}

	sab_time = 0;
	return (Spp);
}

long long rrp (long long* A) {
	rrp_time = 0;
	start = clock();

	short* P = rand_part ();
	long long* M = part_to_arr (A, P);
	long long S = kk_p (M);
	S = llabs(S);

	end = clock();
	rrp_time += ((double) (end - start)) / CLOCKS_PER_SEC;

	while (rrp_time < MAX_TIME) {
		// printf("IN RRP WHILE\n");
		start = clock();

		short* Pp = rand_part ();
		long long* Mp = part_to_arr(A, Pp);
		long long Sp = kk_p (Mp);
		Sp = llabs(Sp);

		if (Sp < S) {
			S = Sp;
		}

		end = clock();
		rrp_time += ((double) (end - start)) / CLOCKS_PER_SEC;
	}

	rrp_time = 0;
	return S;
}

long long hcp (long long* A) {
	hcp_time = 0;
	start = clock();

	short* P = rand_part ();
	long long* M = part_to_arr (A, P);
	long long S = kk_p (M);
	S = llabs(S);

	end = clock();
	hcp_time += ((double) (end - start)) / CLOCKS_PER_SEC;

	while (hcp_time < MAX_TIME) {
		start = clock();
		short* Pp = neighbor_part (P);
		long long* Mp = part_to_arr(A, Pp);
		long long Sp = kk_p (Mp);
		Sp = llabs(Sp);

		if (Sp < S) {
			S = Sp;
			P = Pp;
		}

		end = clock();
		hcp_time += ((double) (end - start)) / CLOCKS_PER_SEC;
	}

	hcp_time = 0;
	return S;
}

long long sap (long long* A) {
	sap_time = 0;
	start = clock();
	
	short* P = rand_part ();
	long long* M = part_to_arr (A, P);
	long long S = kk_p (M);
	S = llabs(S);
	long long Spp = S;

	end = clock();
	sap_time += ((double) (end - start)) / CLOCKS_PER_SEC;

	long long iteration = 0;
	while (sap_time < MAX_TIME) {
		start = clock();

		short* Pp = neighbor_part (P);
		long long* Mp = part_to_arr(A, Pp);
		long long Sp = kk_p (Mp);
		Sp = llabs(Sp);

		if (Sp < S) {
			S = Sp;
			P = Pp;
		}

		else {
			double T = pow (10, 10) * pow (0.8, iteration / 300);
			double p = pow (E, -(Sp - S)/T);
			if ((rand() / (double) RAND_MAX) < p) {
				S = Sp;
				P = Pp;
			}
		}

		if (S < Spp) {
			Spp = S;
		}

		iteration++;
		end = clock();
		sap_time += ((double) (end - start)) / CLOCKS_PER_SEC;
	}

	sap_time = 0;
	return Spp;
}

int main (int argc, char** argv) {
	// arrays to store residues found by each method for each file
	long long kk_a[files];
	long long rrb_a[files];
	long long hcb_a[files];
	long long sab_a[files];
	long long rrp_a[files];
	long long hcp_a[files];
	long long sap_a[files];

	// generates files lines-line files
	genRand (lines, files);

	for (int file = 0; file < files; file++) {
		
		printf("FILE %d\n", file);
		
		// reads in input file to A
		char filename[512];
		sprintf(filename, "%d.txt", file);
    	FILE *fp = fopen(filename, "r");
		if (fp == NULL) {
			printf("error opening inputfile\n");
		}

		long long A[lines + 1];
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
		for (int iter = 0; iter < TRIAL; iter++) {
			kk_a[file] = kk (A);
		}
		printf("KK: ");
		printf("%lld\n", kk_a[file]);

		//*** Binary Values ***//
		// Repeated Random
		long long rrb_sum = 0;
		for (int iter = 0; iter < TRIAL; iter++) {
			long long R = rrb (A);
			rrb_sum += R;
		}
		rrb_a[file] = rrb_sum / TRIAL;
		printf("rrb avg: ");
		printf("%lld\n", rrb_a[file]);

		// Hill Climbing
		long long hcb_sum = 0;
		for (int iter = 0; iter < TRIAL; iter++) {
			long long R = hcb (A);
			hcb_sum += R;
		}
		hcb_a[file] = hcb_sum / TRIAL;
		printf("hcb avg: ");
		printf("%lld\n", hcb_a[file]);

		// Simulated Annealing
		long long sab_sum = 0;
		for (int iter = 0; iter < TRIAL; iter++) {
			long long R = sab (A);
			sab_sum += R;
		}
		sab_a[file] = sab_sum / TRIAL;
		printf("sab avg: ");
		printf("%lld\n", sab_a[file]);

		//*** Prepartitioning ***//
		// Repeated Random
		long long rrp_sum = 0;
		for (int iter = 0; iter < TRIAL; iter++) {
			long long R = rrp (A);
			rrp_sum += R;
		}
		rrp_a[file] = rrp_sum / TRIAL;
		printf("rrp avg: ");
		printf("%lld\n", rrp_a[file]);

		// Hill Climbing
		long long hcp_sum = 0;
		for (int iter = 0; iter < TRIAL; iter++) {
			long long R = hcp (A);
			hcp_sum += R;
		}
		hcp_a[file] = hcp_sum / TRIAL;
		printf("hcp avg: ");
		printf("%lld\n", hcp_a[file]);

		// Simulated Annealing
		long long sap_sum = 0;
		for (int iter = 0; iter < TRIAL; iter++) {
			long long R = sap (A);
			sap_sum += R;
		}
		sap_a[file] = sap_sum / TRIAL;
		printf("sap avg: ");
		printf("%lld\n", sap_a[file]);
	}

	// calculate averages over all files
	long long kk_tot = 0;
	long long rrb_tot = 0;
	long long hcb_tot = 0;
	long long sab_tot = 0;
	long long rrp_tot = 0;
	long long hcp_tot = 0;
	long long sap_tot = 0;

	for (int f = 0; f < files; f++) {
		kk_tot += kk_a[f];
		rrb_tot += rrb_a[f];
		hcb_tot += hcb_a[f];
		sab_tot += sab_a[f];
		rrp_tot += rrp_a[f];
		hcp_tot += hcp_a[f];
		sap_tot += sap_a[f];
	}

	long long kk_final = kk_tot / files;
	long long rrb_final = rrb_tot / files;
	long long hcb_final = hcb_tot / files;
	long long sab_final = sab_tot / files;
	long long rrp_final = rrp_tot / files;
	long long hcp_final = hcp_tot / files;
	long long sap_final = sap_tot / files;

	printf("\nFINAL AVERAGES:\n");
	printf("kk total avg: ");
	printf("%lli\n", kk_final);
	printf("rrb total avg: ");
	printf("%lli\n", rrb_final);
	printf("hcb total avg: ");
	printf("%lli\n", hcb_final);
	printf("sab total avg: ");
	printf("%lli\n", sab_final);
	printf("rrp total avg: ");
	printf("%lli\n", rrp_final);
	printf("hcp total avg: ");
	printf("%lli\n", hcp_final);
	printf("sap total avg: ");
	printf("%lli\n", sap_final);

	return 0;
}
