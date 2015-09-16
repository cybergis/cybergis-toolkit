//
// Copyright 0000 <Nobody>
// @file
// @author David Matthew Mattli <dmattli@usgs.gov>
//
// @section LICENSE
//
// This software is in the public domain, furnished "as is", without
// technical support, and with no warranty, express or implied, as to
// its usefulness for any purpose.
//
// @section DESCRIPTION
//
// The Configuration class represents the configuration of a reprojection task.
//
//

#include <getopt.h>

#include "configuration.h"

namespace librasterblaster {

/** \cond DOXYHIDE **/
struct option longopts[] = {
  {"t_srs", required_argument, NULL, 'p'},
  {"s_srs", required_argument, NULL, 's'},
  {"partition-size", required_argument, NULL, 'n'},
  {"resampler", required_argument, NULL, 'r'},
  {"dstnodata", required_argument, NULL, 'f'},
  {"tile-size", required_argument, NULL, 'x'},
  {"timing-file", required_argument, NULL, 'c'},
  {"output-ratio", required_argument, NULL, 'o'},
  {0, 0, 0, 0}
};
/** \endcode **/

Configuration::Configuration() {
  resampler = NEAREST;
  partition_size = -1;
  tile_size = 1024;
  timing_filename = "";
  cell_dimension_ratio = 1.0;
}

Configuration::Configuration(int argc, char *argv[]) {
  signed char c = 0;
  partition_size = 1;
  std::string arg = "";
  resampler = BILINEAR;
  tile_size = 1024;
  timing_filename = "";
  cell_dimension_ratio = 1.0;

  while ((c = getopt_long(argc,
                          argv,
                          "p:r:f:n:x:c",
                          longopts, NULL)) != -1) {
    switch (c) {
      case 0:
        // getopt_long() set a variable, just keep going
      case 'p':
        output_srs = optarg;
        break;
      case 's':
        source_srs = optarg;
        break;
      case 'n':
        partition_size = std::stoi(optarg);
        break;
      case 'r':
        arg = optarg;
        if (arg == "mean") {
          resampler = MEAN;
        } else if (arg == "nearest") {
          resampler = NEAREST;
        } else if (arg == "bilinear") {
          resampler = BILINEAR;
        } else if (arg == "bicubic") {
          resampler = BICUBIC;
        } else if (arg == "lanczos") {
          resampler = LANCZOS;
        }
        break;
      case 'f':
        fill_value = optarg;
        break;
      case 'x':
        tile_size = std::stoi(optarg);
        break;
      case 'c':
        timing_filename = optarg;
        break;
      case 'o':
        cell_dimension_ratio = std::stof(optarg);
        break;
      default:
        fprintf(stderr, "%s: option '-%c' is invalid: ignored\n",
                argv[0], optopt);
        break;
    }
    if (argc > 2) {
      output_filename = argv[argc-1];
      input_filename = argv[argc-2];
    }
  }
  return;
}
}
