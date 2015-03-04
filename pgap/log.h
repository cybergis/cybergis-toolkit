/******************************************************************************
 * log.h: header for logging support                                          *
 * Author: Yan Y. Liu <yanliu@illinois.edu>                                   *
 * Date: 2014/08/17                                                           *
 * Copyright and license of this source file are specified in LICENSE.TXT     *
 * under the root directory of this software package.                         *
 ******************************************************************************/
#ifndef LOG_H
#define LOG_H

#define LOG_BUF_SIZE 1024*64
#define LOG_ENTRY_BUF_SIZE 1024*16
#define LOG_INTERVAL 20

#include "addr.h"
void glog_init(char *d, char *p, IDTYPE i);
void glog(char *format, ...);


#endif
