/******************************************************************************
 * gsl_sprng.c: rewrite of gsl-sprng.h to use sprng streams                   *
 * Author: Yan Y. Liu <yanliu@illinois.edu>                                   *
 * Date: 2014/08/17                                                           *
 * Copyright and license of this source file are specified in LICENSE.TXT     *
 * under the root directory of this software package.                         *
 ******************************************************************************/
#ifndef GSL_SPRNG_C
#define GSL_SPRNG_C

#include <sys/time.h>
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include "gsl-sprng.h"
/***************************/
/* gsl definition of sprng */
/***************************/
void sprng_set(void * vstate, unsigned long int seed)
{
    gsl_sprng_state_t *mysprng = (gsl_sprng_state_t*) vstate;

    if (mysprng->stream != NULL) {
        free_sprng(mysprng->stream);
    }
    mysprng->stream = init_sprng(
                DEFAULT_RNG_TYPE,
                mysprng->streamnum,
                mysprng->nstreams,
                seed,
                SPRNG_DEFAULT);
}

unsigned long sprng_get(void * vstate)
{  
    return( (long) isprng( ((gsl_sprng_state_t*) vstate)->stream ) );
}

double sprng_get_double(void * vstate)
{
    return( (double) sprng( ((gsl_sprng_state_t*) vstate)->stream ) );
}

/* global sprng 2.0 type for gsl initialization */
const gsl_rng_type gsl_rng_sprng20 = {
    "sprng20",        /* name */
    0x7fffffffUL,     /* RAND_MAX */
    0,                /* RAND_MIN */
    0,                /* size of state - not sure about this */
    &sprng_set,          /* initialisation */
    &sprng_get,          /* get integer RN */
    &sprng_get_double    /* get double RN */
};

/***************************/
/*      gsl_sprng API      */
/***************************/

// if stream=NULL, create a new stream
// otherwise, just use the stream
gsl_rng *gsl_sprng_alloc (int streamnum, int nstreams, int seed, int * stream)
{
    gsl_sprng_state_t *  gsl_sprng_state = 
   	    (gsl_sprng_state_t*)  malloc( sizeof( gsl_sprng_state_t ) );

    if (gsl_sprng_state == NULL) {
        GSL_ERROR_VAL (
	   	  "failed to allocate space for gsl_sprng_state_t struct",
                     GSL_ENOMEM, 0);
        return NULL; // not reachable
    }
    gsl_sprng_state->streamnum = streamnum;
    gsl_sprng_state->nstreams = nstreams;
    // set seed
    if (seed == -1) {
        struct timeval tv; struct timezone tz;
        gettimeofday(&tv, &tz);
        seed = (int)(tv.tv_sec * streamnum);
    } else if (seed == 0) {
        seed = gsl_rng_default_seed; // zero, too, in fact
    }
    gsl_sprng_state->seed = seed;
    // sprng init
    if (stream == NULL) {
        stream = init_sprng(
                DEFAULT_RNG_TYPE,
                streamnum,
                nstreams,
                seed,
                SPRNG_DEFAULT);
    }
    gsl_sprng_state->stream = stream;
    // gsl_rng init
    gsl_rng * r = (gsl_rng *) malloc (sizeof (gsl_rng));
    if (r == NULL) {
        GSL_ERROR_VAL ("gsl_sprng: failed to allocate space for gsl_rng struct",
                        GSL_ENOMEM, 0);
        free(gsl_sprng_state); // not reachable
        return NULL;
    }
    r->type = &gsl_rng_sprng20;
    r->state = gsl_sprng_state;
    
    return r;
}

void gsl_sprng_free (gsl_rng * r)
{
    if (r==NULL) return;
    gsl_sprng_state_t *  gsl_sprng_state = (gsl_sprng_state_t *) r->state;
    if (gsl_sprng_state != NULL) {
        int numAvail = free_sprng(gsl_sprng_state->stream);
        fprintf(stdout, "gsl_sprng: free random stream. %d streams left\n", numAvail);
        free(r->state);
    }
    free(r);
}

void gsl_sprng_print (gsl_rng *r)
{
    if (r==NULL) return;
    gsl_sprng_state_t *  gsl_sprng_state = (gsl_sprng_state_t *) r->state;
    if (gsl_sprng_state == NULL) return;
    printf("gsl_sprng: streamnum=%d, nstreams=%d, seed=%d\n", gsl_sprng_state->streamnum, gsl_sprng_state->nstreams, gsl_sprng_state->seed);
    int * stream = gsl_sprng_state->stream;
    if (stream == NULL) return;
    print_sprng(stream);
}

int gsl_sprng_fread (char * fpath, gsl_rng * r)
{
    // delete old stream
    if (r!=NULL) {
        gsl_sprng_free(r);
        r = NULL;
    }
    // allocate memory to load data
    char buf[MAX_PACKED_LENGTH];
    // read stream from file to array
    FILE * f;
    int size;
    if ((f=fopen(fpath, "r")) == NULL) {
        fprintf(stderr, "gsl_sprng_fread(): failed to open %s for reading\n", fpath);
        return 0;
    }
    fread(&size,sizeof(int),1,f);
    fprintf(stdout, "sprng_fread[%s]: size=%d\n", fpath, size);
    int numRead = fread(buf, sizeof(char), size, f);
    if (numRead != size) {
        fprintf(stderr, "gsl_sprng_fread(): failed to read stream of size %d\n", size);
        fclose(f);
        return 0;
    }
    fclose(f);
    // load stream from array
    int * stream = unpack_sprng(buf);
    if (stream == NULL) {
        fprintf(stderr, "gsl_sprng_fread(): failed to read stream of size %d\n", size);
        return 0;
    }
    print_sprng(stream);
    // create new gsl_rng
    r = gsl_sprng_alloc(0, 0, 0, stream);
    if (r==NULL) return 0;
    return 1;
}

int gsl_sprng_fwrite (char * fpath, const gsl_rng * r)
{
    if (r==NULL) return 0;
    gsl_sprng_state_t *  gsl_sprng_state = (gsl_sprng_state_t*) r->state;
    if (gsl_sprng_state == NULL) return 0;
    int * stream = gsl_sprng_state->stream;
    if (stream==NULL) return 0;
    // write state to array
    char *buf;
    int size = pack_sprng(stream, &buf);
    // write array to file
    FILE * f;
    if ((f=fopen(fpath, "w")) == NULL) {
        fprintf(stderr, "gsl_sprng_fwrite(): failed to open %s for writing\n", fpath);
        return 0;
    }
    // our format: filesize, filecontent
    fwrite(&size,sizeof(int),1,f);
    int numWrote = fwrite(buf, sizeof(char), size, f);
    if (numWrote != size) {
        fprintf(stderr, "gsl_sprng_fwrite(): failed to write stream of size %d\n", size);
    }
    fclose(f);
    // sprng uses malloc() to allocate buf, we must free it properly
    free(buf);
    
    return 1;
}

#endif
