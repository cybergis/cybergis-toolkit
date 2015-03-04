/******************************************************************************
 * data.c: read problem instance data; define macros to access matrix         *
 * Author: Yan Y. Liu <yanliu@illinois.edu>                                   *
 * Date: 2014/08/17                                                           *
 * Copyright and license of this source file are specified in LICENSE.TXT     *
 * under the root directory of this software package.                         *
 ******************************************************************************/
#ifndef DATA_C
#define DATA_C

/** Read standard datasets
 * data is read into *v (value matrix), *w (weight matrix)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "data.h"

// global value and weight matrix
VTYPE *v = NULL; // value matrix
VTYPE *vt = NULL; // value matrix - transposed for fast memory access
WTYPE *w = NULL; // weight matrix
VTYPE *wt = NULL; // weight matrix - transposed for fast memory access
WTYPE *b = NULL; // capacity vector
int n = 0; // number of items
int m = 0; // number of bins

void readData(char * fileName)
{
	FILE *f;
	if ((f=fopen(fileName, "r")) == NULL) {
		fprintf(stderr, "Error read data from file %s\n", fileName);
		exit(0);
	}
	fscanf(f, "%d", &m);
	fscanf(f, "%d", &n);
	if (n<=0 || m<=0) {
		fprintf(stderr, "Invalid size of bins (%d) or items (%d)\n", m, n);
		exit(0);
	}
	v = (VTYPE *)malloc(m * n * sizeof(VTYPE));
	vt = (VTYPE *)malloc(n * m * sizeof(VTYPE));
	w = (WTYPE *)malloc(m * n * sizeof(WTYPE));
	wt = (WTYPE *)malloc(n * m * sizeof(WTYPE));
	b = (WTYPE *)malloc(m * sizeof(WTYPE));
	memset(v, 0, m * n * sizeof(VTYPE));
	memset(vt, 0, n * m * sizeof(VTYPE));
	memset(w, 0, m * n * sizeof(WTYPE));
	memset(wt, 0, n * m * sizeof(WTYPE));
	memset(b, 0, m * sizeof(WTYPE));
	int i, j, val;
	for (i=0; i<m; i++) {
		for (j=0; j<n; j++) {
			fscanf(f, "%d", &val);
			*(v + n * i + j) = val;
			*(vt + m * j + i) = val;
		} 
	}
	for (i=0; i<m; i++) {
		for (j=0; j<n; j++) {
			fscanf(f, "%d", &val);
			*(w + n * i + j) = val;
			*(wt + m * j + i) = val;
		} 
	}

	for (i=0; i<m; i++)
		fscanf(f, "%d", b + i);
	/* debug 
	fprintf(stdout, "GAP of %d items and %d bins\n", n, m);
	fprintf(stdout, "============value matrix===========\n");
	for (i=0; i<m; i++) {
		for (j=0; j<n; j++)
			fprintf(stdout, "%d ", M(v, i, j));
		fprintf(stdout, "\n");
	}
	fprintf(stdout, "============weight matrix===========\n");
	for (i=0; i<m; i++) {
		for (j=0; j<n; j++)
			fprintf(stdout, "%d ", M(w, i, j));
		fprintf(stdout, "\n");
	}
	fprintf(stdout, "============capacity vector===========\n");
	for (i=0; i<m; i++) {
		fprintf(stdout, "%d ", b[i]);
	}
	fprintf(stdout, "\n");
	*/
}

#endif
