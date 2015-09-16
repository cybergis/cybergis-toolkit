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

#include <string>

/**
 * @brief rastercompare takes two raster files and compares the pixels,
 * verifying they are within a certain delta
 *
 * @param control_filename filename of a known good raster
 *
 * @param test_filename filename of a raster to be compared to control_filename
 *
 */
int rastercompare(std::string golden_filename, std::string test_filename);

