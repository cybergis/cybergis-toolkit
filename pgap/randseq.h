/******************************************************************************
 * randseq.h: header for generating unique random sequence                    *
 * Author: Yan Y. Liu <yanliu@illinois.edu>                                   *
 * Date: 2014/08/17                                                           *
 * Copyright and license of this source file are specified in LICENSE.TXT     *
 * under the root directory of this software package.                         *
 ******************************************************************************/
#ifndef RANDSEQ_H
#define RANDSEQ_H

/* generate a random sequence of size n, the problem size */
int randseq_init(int **randseq_holder, int size);
int randseq_shuffle(int *randseq, int size);
int randseq_finalize(int *randseq, int size);
int randseq_verify(int *randseq, int size);
#endif
