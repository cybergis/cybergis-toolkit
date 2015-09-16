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
// Implements main() function for prasterblasterpio command-line tool.
//
//

#include <mpi.h>

#include "../configuration.h"
#include "../demos/prasterblaster-pio.h"

using librasterblaster::Configuration;
using librasterblaster::prasterblasterpio;
using librasterblaster::RasterChunk;
/*! \page prasterblasterpio

\htmlonly
USAGE:
\endhtmlonly

\verbatim
prasterblaster [--t_srs target_srs] [--s_srs source_srs]
               [-r resampling_method] [-n partition_size]
               [--dstnodata no_data_value]
               source_file destination_file

\endverbatim

\section prasterblasterpio_description DESCRIPTION

<p>

The prasterblasterpio demo program implements parallel raster reprojection and
demonstrates the use of librasterblaster. The implementation can be found in
prasterblaster-pio.cc.

</p>
 */
int main(int argc, char *argv[]) {
  // Give MPI_Init first run at the command-line arguments
  MPI_Init(&argc, &argv);

  // Initialize Configuration object
  Configuration conf(argc, argv);

  int ret = prasterblasterpio(conf);
  MPI_Finalize();

  return ret;
}
