///
/// Copyright 0000 <Nobody>
/// @file
/// @author David Matthew Mattli <dmattli@usgs.gov>
///
/// @section LICENSE
///
/// This software is in the public domain, furnished "as is", without
/// technical support, and with no warranty, express or implied, as to
/// its usefulness for any purpose.
///
/// @section DESCRIPTION
///
/// This file demonstrates how to use librasterblaster to implement serial
/// raster reprojection. This implementation uses GDAL through librasterblaster
/// to write to a tiff file in a single process.
///
///

#include <gdal.h>

#include <vector>

#include "src/configuration.h"
#include "src/reprojection_tools.h"

#include "src/utils.h"

#define STR_EXPAND(tok) #tok
#define STR(tok) STR_EXPAND(tok)

using librasterblaster::PRB_ERROR;
using librasterblaster::PRB_NOERROR;
using librasterblaster::Area;
using librasterblaster::RasterChunk;

int main(void) {
  const string input_filename = STR(__PRB_SRC_DIR__)
      "/tests/testdata/veg.tif";
  const string output_filename = STR(__PRB_SRC_DIR__)
      "/tests/testdata/veg_test_moll.tif";
  const int tile_size = 256;
  const std::string output_srs = "+proj=moll";
  const int partition_size = 4;
  const librasterblaster::RESAMPLER resampler = librasterblaster::MEAN;
  const string fillvalue = "4";

  // Very important!
  GDALAllRegister();

  // Open the input raster. Read the GDAL API tutorial to understand these GDAL
  // calls.
  GDALDataset *input_raster =
      static_cast<GDALDataset*>(GDALOpen(input_filename.c_str(),
                                         GA_ReadOnly));

  if (input_raster == NULL) {
    fprintf(stderr, "Error opening input file!\n");
    return 1;
  }

  // Create the output raster file. This function performs the minbox operation
  // to determine the size of the output raster file
  PRB_ERROR err = librasterblaster::CreateOutputRaster(input_raster,
                                                       output_filename,
                                                       output_srs,
                                                       tile_size);

  if (err != PRB_NOERROR) {
    fprintf(stderr, "Error creating output file!\n");
    return 1;
  }

  // Open new output raster file
  GDALDataset *output_raster =
      static_cast<GDALDataset*>(GDALOpen(output_filename.c_str(),
                                         GA_Update));

  std::vector<Area> partitions;
  partitions = librasterblaster::BlockPartition(0,
                                                1,
                                                output_raster->GetRasterYSize(),
                                                output_raster->GetRasterXSize(),
                                                tile_size,
                                                partition_size);

  for (std::vector<Area>::iterator i = partitions.begin();
       i != partitions.end();
       ++i) {
    // Each partition represents an area of the output raster. The input minbox
    // has to be calculated for each partition. This call to CreateRasterChunk
    // performs the minbox operation and initializes the RasterChunk.
    RasterChunk *in_chunk = RasterChunk::CreateRasterChunk(input_raster,
                                                                output_raster,
                                                                *i);

    // The RasterChunk is initialized but we still have to read the pixel values
    // from the file. This is done by calling ReadRasterChunk.
    PRB_ERROR chunk_err = RasterChunk::ReadRasterChunk(input_raster, in_chunk);

    if (chunk_err != PRB_NOERROR) {
      fprintf(stderr, "Error reading from input!\n");
      return 1;
    }

    RasterChunk *out_chunk = RasterChunk::CreateRasterChunk(output_raster,
                                                            *i);

    if (out_chunk == NULL) {
      fprintf(stderr, "Error initializing output chunk!\n");
      return 1;
    }

    bool ret = librasterblaster::ReprojectChunk(in_chunk,
                                                out_chunk,
                                                fillvalue,
                                                resampler);
    if (ret == false) {
      fprintf(stderr, "Error performing reprojection!\n");
      return 1;
    }

    // Finally write the output chunk to the output raster file
    chunk_err = RasterChunk::WriteRasterChunk(output_raster, out_chunk);

    delete in_chunk;
    delete out_chunk;
  }

  printf("Reprojection successful! Open the file %s to view the output\n",
         output_filename.c_str());
  return 0;
}
