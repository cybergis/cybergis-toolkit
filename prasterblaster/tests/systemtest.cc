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
 * System level tests
 *
 */

#include <vector>

#include <mpi.h>
#include <unistd.h>

#include "gtest/gtest.h"
#include "../src/configuration.h"
#include "../src/reprojection_tools.h"
#include "../src/resampler.h"
#include "../src/demos/prasterblaster-pio.h"

#include "../src/demos/sptw.h"

#include "rastercompare.h"

using librasterblaster::Configuration;
using librasterblaster::PRB_NOERROR;

#define STR_EXPAND(tok) #tok
#define STR(tok) STR_EXPAND(tok)

namespace {
TEST(SystemTest, GLOBALVEG) {
  const std::string golden_rasters[] = { "aea", "cea", "eck4", "eck6", "gall",
                                         "gnom", "laea", "merc", "mill",
                                         "moll", "sinu", "vandg"};
  const int gold_count = std::extent<decltype(golden_rasters)>::value;
  Configuration conf;

  GDALAllRegister();

  for (unsigned int i = 0; i < gold_count; ++i) {
    const std::string gold_name = STR(__PRB_SRC_DIR__) "/tests/testdata/veg_" + golden_rasters[i] + ".tif";
    const std::string test_name = STR(__PRB_SRC_DIR__) "/tests/testdata/veg_test_" + golden_rasters[i] + ".tif";

    // Open golden raster to extract metadata
    GDALDataset *golden = static_cast<GDALDataset*>(GDALOpen(gold_name.c_str(),
                                                             GA_ReadOnly));

    if (golden == NULL) {
      ADD_FAILURE() << "Failed to open golden raster: " + gold_name;
    }
    ASSERT_TRUE(golden != NULL);

    conf.input_filename = STR(__PRB_SRC_DIR__) "/tests/testdata/veg.tif";
    conf.output_filename = test_name;
    conf.resampler = librasterblaster::NEAREST;
    conf.output_srs = golden->GetProjectionRef();
    conf.tile_size = 16;
    conf.partition_size = 1;

    GDALClose(golden);

    int ret = prasterblasterpio(conf);
    ASSERT_EQ(PRB_NOERROR, ret);

    int raster_compare_ret = rastercompare(gold_name, test_name);

    ASSERT_EQ(0, raster_compare_ret);
    unlink(test_name.c_str());
  }
  SUCCEED();
}
}  // namespace

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);
  ::testing::InitGoogleTest(&argc, argv);
  int ret = RUN_ALL_TESTS();
  MPI_Finalize();

  return ret;
}
