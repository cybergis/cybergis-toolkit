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
 * Header for high-level API
 *
 */
#ifndef SRC_DEMOS_PRASTERBLASTER_PIO_H_
#define SRC_DEMOS_PRASTERBLASTER_PIO_H_

#include <string>

#include "src/configuration.h"
#include "src/demos/sptw.h"
#include "src/utils.h"

namespace librasterblaster {
/**
 * @brief write_rasterchunk writes chunk to the PTIFF ptiff using the sptw library
 *
 * @param ptiff PTIFF target for output
 * @param chunk RasterChunk to be written to file
 *
 */
PRB_ERROR write_rasterchunk(sptw::PTIFF *ptiff,
                            RasterChunk *chunk);

/**
 * @brief prasterblasterpio performs a complete, potentially parallel, raster
 * reprojection job using the parameters specified by conf.
 *
 * @param conf Configuration class that describes the reprojection task
 *
 *
 */
PRB_ERROR prasterblasterpio(librasterblaster::Configuration conf);
}

#endif  // SRC_DEMOS_PRASTERBLASTER_PIO_H_

