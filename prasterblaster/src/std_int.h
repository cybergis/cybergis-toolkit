/*!
 * Copyright 0000 <Nobody>
 * @file
 * @author David Matthew Mattli <dmattli@usgs.gov>
 *
 * @section LICENSE
 *
 * This software is in the public domain, furnished "as is", without
 * technical support, and with no warranty, express or implied, as to
 * its usefulness for any purpose.
 *
 * @section DESCRIPTION
 *
 *
 *
 */

#ifndef SRC_STD_INT_H_
#define SRC_STD_INT_H_

#ifndef HAVE_CONFIG_H
//#include "../config.h"
#endif

#if HAVE_STDINT_H == 0
#include <cstdint>
#else
#include <stdint.h>
#endif

#endif  // SRC_STD_INT_H_
