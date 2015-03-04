/******************************************************************************
 * addr.h: id management                                                      *
 * Author: Yan Y. Liu <yanliu@illinois.edu>                                   *
 * Date: 2014/08/17                                                           *
 * Copyright and license of this source file are specified in LICENSE.TXT     *
 * under the root directory of this software package.                         *
 ******************************************************************************/
#ifndef ADDR_C
#define ADDR_C

#include "addr.h"

IDTYPE globalId = 0; // every process has a global id
IDTYPE localCount = 0; // locally the count indicates thread id

#endif
