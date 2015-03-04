/******************************************************************************
 * mysprng.c: SPRNG wrapper                                                   *
 * Author: Yan Y. Liu <yanliu@illinois.edu>                                   *
 * Date: 2014/08/17                                                           *
 * Copyright and license of this source file are specified in LICENSE.TXT     *
 * under the root directory of this software package.                         *
 ******************************************************************************/
#ifndef MYSPRNG_C
#define MYSPRNG_C

#include <stdio.h>
#include <stdlib.h>
#include "mysprng.h"

#ifdef PGAMODE
#include "mpi.h"
#endif

// global sprng stream
int * mysprng_stream;
// wrappers for sprng()
double mysprng() {
    if (mysprng_stream == NULL) {
        fprintf(stderr, "SPRNG stream is not initialized\n");
        exit(0); // fatal error
    } else {
        return sprng(mysprng_stream);
    }    
}
// this func must be called after MPI_Init()
void mysprng_init() {
    int gseed = make_sprng_seed();
#ifdef PGAMODE
    int myrank, net_size;
    MPI_Comm_rank( MPI_COMM_WORLD, &myrank );
    MPI_Comm_size( MPI_COMM_WORLD, &net_size);
    mysprng_stream = init_sprng(DEFAULT_RNG_TYPE, myrank, net_size, gseed, SPRNG_DEFAULT);
#else
    mysprng_stream = init_sprng(DEFAULT_RNG_TYPE, 0, 1, gseed, SPRNG_DEFAULT);
#endif
    return;
}
void mysprng_free() {
    if (mysprng_stream ==  NULL) return;
    free_sprng(mysprng_stream);
}

#endif
