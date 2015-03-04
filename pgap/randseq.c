/******************************************************************************
 * randseq.h: generating unique random sequence                               *
 * Author: Yan Y. Liu <yanliu@illinois.edu>                                   *
 * Date: 2014/08/17                                                           *
 * Copyright and license of this source file are specified in LICENSE.TXT     *
 * under the root directory of this software package.                         *
 ******************************************************************************/
#ifndef RANDSEQ_C
#define RANDSEQ_C
/* generate a random sequence of size n, the problem size */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "myrng.h"
#include "randseq.h"

// shuffle the random sequence memory using Fisher Yates algorithm
// http://en.wikipedia.org/wiki/Fisher%E2%80%93Yates_shuffle
int randseq_shuffle(int *randseq, int size)
{
	int i, j, tmp;
	for (i=0; i<size-1; i++) {
		j = MYRANDI(size-i);
		tmp = randseq[i];
		randseq[i] = randseq[j];
		randseq[j] = tmp;
	}
	return 1;
}
// init the random sequence memory
int randseq_init(int **randseq_holder, int size) {
	if (size <= 0) return 0;
	int i;
	int *randseq; 
	if (*randseq_holder == NULL) {
		randseq = (int *)malloc(sizeof(int) * size);
		if (randseq == NULL) {
			fprintf(stderr, "randseq_init(): ERROR: out of memory in getting %d integers.\n", size);
			exit(1);
		}
		*randseq_holder = randseq;
		for (i=0; i<size; i++) randseq[i] = i;
		randseq_shuffle(randseq, size);
	}
	return 1;
}
// free memory. NOTE: randseq will NOT be set to NULL after this call. Set it by yourself!!!
int randseq_finalize(int *randseq, int size)
{
	if (randseq != NULL) free(randseq);
	randseq = NULL;
	return 1;
}
int randseq_verify(int *randseq, int size)
{
	int i; char ihash[size];
	memset(ihash, 0, sizeof(char) * size);
	for (i=0; i<size; i++) {
		if (randseq[i]<0 || randseq[i]>=size || ihash[randseq[i]] > 0) {
			fprintf(stderr, "randseq[%d/%d] = %d. either dup or out of bound\n", i, size, randseq[i]);
			return 0;
		}
		else ihash[randseq[i]] ++;
	}
	return 1;
}

#endif
