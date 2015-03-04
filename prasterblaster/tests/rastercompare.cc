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
//
//

#include <string>

#include <gdal_priv.h>
#include <gdal.h>

#include "tests/rastercompare.h"

using std::string;

int rastercompare(string control_filename, string test_filename) {
  const double delta = 0.001;
  GDALAllRegister();

  GDALDataset *control =
      static_cast<GDALDataset*>(GDALOpen(control_filename.c_str(),
                                         GA_ReadOnly));
  GDALDataset *test =
      static_cast<GDALDataset*>(GDALOpen(test_filename.c_str(),
                                         GA_ReadOnly));

  if (control == NULL) {
    fprintf(stderr, "Error opening control raster!\n");
    return 1;
  }

  if (test == NULL) {
    fprintf(stderr, "Error opening test raster!\n");
    return 1;
  }

  // Basic checks
  if (control->GetRasterYSize() != test->GetRasterYSize()) {
    printf("Control raster has %d rows and test has %d\n",
           control->GetRasterYSize(),
           test->GetRasterYSize());
    return 1;
  }

  if (control->GetRasterXSize() != test->GetRasterXSize()) {
    printf("Control raster has %d columns and test has %d\n",
           control->GetRasterXSize(),
           test->GetRasterXSize());
    return 1;
  }

  // Check pixel values
  const int band_count = control->GetRasterCount();
  GDALRasterBand *band = control->GetRasterBand(1);
  int block_x_size, block_y_size;
  band->GetBlockSize(&block_x_size, &block_y_size);

  double *control_pixels, *test_pixels;

  control_pixels = new double[band_count
                              * control->GetRasterXSize()
                              * sizeof(*control_pixels)];
  test_pixels = new double[band_count
                           * control->GetRasterXSize()
                           * sizeof(*test_pixels)];

  // Loop over rows of blocks
  //   Loop over columns of blocks
  //     Loop over rows in block
  //      Loop over columns in block
  int bad_pixels = 0;
  const int y_size = control->GetRasterYSize();
  const int x_size = control->GetRasterXSize();
  for (int y = 0; y < y_size; ++y) {
    // Read a row
    control->RasterIO(GF_Read,
                      0,
                      y,
                      control->GetRasterXSize(),
                      1,
                      static_cast<void*>(control_pixels),
                      control->GetRasterXSize(),
                      1,
                      GDT_Float64,
                      band_count,
                      NULL,
                      0,
                      0,
                      0);
    test->RasterIO(GF_Read,
                   0,
                   y,
                   test->GetRasterXSize(),
                   1,
                   static_cast<void*>(test_pixels),
                   test->GetRasterXSize(),
                   1,
                   GDT_Float64,
                   band_count,
                   NULL,
                   0,
                   0,
                   0);
    if (memcmp(control_pixels,
               test_pixels,
               x_size * sizeof(*control_pixels)) != 0) {
      for (int x = 0; x < x_size; ++x) {
        if (fabs(control_pixels[x] - test_pixels[x]) > delta) {
          printf("Values at (%d, %d) are too different! %f vs %f\n",
                 x, y, control_pixels[x], test_pixels[x]);
          bad_pixels++;
        }
      }
    }
  }

  // Cleanup
  delete control_pixels;
  delete test_pixels;
  GDALClose(control);
  GDALClose(test);
  if (bad_pixels == 0) {
    return 0;
  } else {
    if (bad_pixels == 1) {
      printf("One pixel was different\n");
    } else {
      printf("%d pixels were different\n", bad_pixels);
    }
    return 1;
  }
}


