/******************************************************************************
 * data.h: read problem instance data; define macros to access matrix         *
 * Author: Yan Y. Liu <yanliu@illinois.edu>                                   *
 * Date: 2014/08/17                                                           *
 * Copyright and license of this source file are specified in LICENSE.TXT     *
 * under the root directory of this software package.                         *
 ******************************************************************************/
#ifndef DATA_H
#define DATA_H

/* read data into memory
 * data format conforms OR-LIB format
 */

typedef int VTYPE; // type of value
typedef int WTYPE; // type of weight

// matrix access: size of a row is implicitly decided as n
#define M(matrix, i, j) (matrix[n * i + j])
// transposed matrix access
#define MT(matrix, i, j) (matrix[m * i + j])
#define Mr(matrix, i) (matrix[n * i])

void readData(char *fileName);

#endif
