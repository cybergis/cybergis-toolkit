/**
 * Copyright 0000 <Nobody>
 * @file
 * @author David Matthew Mattli <dmattli@usgs.gov>
 *
 * @section LICENSE
 *
 * This software is in the public domain, furnished "as is", without
 * echnical support, and with no warranty, express or implied, as to
 * its usefulness for any purpose.
 *
 * @section DESCRIPTION
 *
 * The Configuration class represents the configuration of a reprojection task.
 *
 */


#ifndef SRC_CONFIGURATION_H_
#define SRC_CONFIGURATION_H_

#include <string>

#include "src/resampler.h"
#include "src/reprojection_tools.h"

using std::string;

namespace librasterblaster {
/// Configuration class
/**
 * This class implements parsing of command-line arguments
 */
struct Configuration {
  /**
   * @brief Default constructor
   */
  Configuration();
  /**
   * @brief Constructor, takes two arguments
   *
   * @param argc Number of arguments in array
   * @param argv Pointer to array of arguments
   *
   * Constructs a Configuration object by parsing the provided command-line
   * arguments and filling in the values of the class. If a value isn't provided
   * in the cli arguments, the variable is set to the default.
   *
   */
  Configuration(int argc, char *argv[]);

  /**
   * @brief This is set to the filesystem path to the input file. The default value is
   * "".
   */
  string input_filename;
  /**
   * @brief This is set to the spatial reference system string. The default value is "".
   */
  string input_srs;
  /**
   * @brief This is set to the filesystem path to the output file. The default value is
   * "".
   */
  string output_filename;
  /**
   * @brief This is set to the spatial reference system of the output file. The default
   * value is "".
   */
  string output_srs;
  /**
   * @brief This is set to the spatial reference system of the input file. The default
   * value is "".
   */
  string source_srs;
  /**
   * @brief This is set to the resampler method specified by the user. The default
   * value is NEAREST.
   */
  RESAMPLER resampler;
  /**
   * @brief This is set to the fillvalue specified by the user.
   */
  string fillvalue;
  /**
   * @brief Value to set in nodata areas of the raster. The interpretation of this
   * string will depend on the pixel type.
   */
  string nodata_value;
  /**
   * @brief The maximum size of output raster partitions. The default value is 0;
   */
  int partition_size;
  /**
   * @brief Tile size
   */
  int tile_size;
  /**
   * @brief Name of timing information file
   */
  string timing_filename;
};
}

#endif  // SRC_CONFIGURATION_H_
