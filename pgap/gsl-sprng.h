/******************************************************************************
 * gsl_sprng.h: rewrite of gsl-sprng.h to use sprng streams                   *
 * Author: Yan Y. Liu <yanliu@illinois.edu>                                   *
 * Date: 2014/08/17                                                           *
 * Copyright and license of this source file are specified in LICENSE.TXT     *
 * under the root directory of this software package.                         *
 ******************************************************************************/

/* 
 * gsl-sprng.h: rewrite of gsl-sprng.h to use sprng streams
 * 
 * This header file is based on  Darren Wilkinson's version for simple interface,
 * which does not use streams. We also refer to IceCube's C++ extension based on
 * sprng 4.0.
 * 
 * Darren Wilkinson  http://www.staff.ncl.ac.uk/d.j.wilkinson/
 * IceCube project: http://icecube.wisc.edu/
 * 
 * To use, just add the line:
 * #include "gsl-sprng.h"
 * immediately after the line:
 * #include <gsl/gsl_rng.h>
 * near the start of your code. 
 * Make sure you alloc the rng on each processor. If you wish to
 * set a seed, you should set it to be the same on each processor.
 */

#ifndef GSL_SPRNG_H
#define GSL_SPRNG_H

#define USE_MPI
#include "sprng.h"

/***************************/
/* gsl definition of sprng */
/***************************/

typedef struct
{
    int streamnum;
    int nstreams;
    int seed;
    int *stream;
} gsl_sprng_state_t;

/* seed is gsl_rng_default_seed
 * either 0 by default or GSL_RNG_SEED environment variable
 */
void sprng_set(void * vstate, unsigned long int seed);
/* the get() function to get integer random number */
unsigned long sprng_get(void * vstate);
/* the get() function to get double random number */
double sprng_get_double(void * vstate);

/***************************/
/*      gsl_sprng API      */
/***************************/
/* replace all gsl_rng_* func with corresponding one below */

/* initialize gsl_rng */
gsl_rng *gsl_sprng_alloc (int streamnum, int nstreams, int seed, int * stream);
/* finalize gsl_rng */
void gsl_sprng_free (gsl_rng * r);
/* print out gsl_rng stream info for debug purpose */
void gsl_sprng_print (gsl_rng * r);
/* load gsl_rng from a file */
int gsl_sprng_fread (char * fpath, gsl_rng * r);
/* store gsl_rng to a file */
int gsl_sprng_fwrite (char * fpath, const gsl_rng * r);
#endif

