/******************************************************************************
 * mysprng.h: header for SPRNG wrapper                                        *
 * Author: Yan Y. Liu <yanliu@illinois.edu>                                   *
 * Date: 2014/08/17                                                           *
 * Copyright and license of this source file are specified in LICENSE.TXT     *
 * under the root directory of this software package.                         *
 ******************************************************************************/
#ifndef MYSPRNG_H
#define MYSPRNG_H

#include "sprng.h"

#ifdef PGAMODE
#define USE_MPI
#endif

double mysprng();
void mysprng_init();
void mysprng_free();

#endif
