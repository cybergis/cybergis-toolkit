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
// The RasterChunk class represents an in-memory, georeferenced section of
// raster.
//
//

#ifndef SRC_RASTERCHUNK_H_
#define SRC_RASTERCHUNK_H_

#include <string>
#include <gdal_priv.h>
#include "utils.h"

namespace librasterblaster {
/// A class representing an in-memory part of a raster.
class RasterChunk {
 public:

  RasterChunk() {
    pixels = NULL;
  }
  /**
   * @brief
   * This function creates a RasterChunk.
   *
   * @param ds Dataset to create chunk from
   * @param chunk_area The inclusive area that the chunk should represent.
   *
   */
  RasterChunk(GDALDataset *ds, Area chunk_area);

  /**
   * @brief
   * Copy constructor
   *
   * @param s Source RasterChunk to be copied
   */
  RasterChunk(const RasterChunk &s);

  /// RasterChunk destructor
  /**
   * This destructor frees the memory, if any, allocated at pixels_.
   */
  ~RasterChunk() {
    if (this->pixels != NULL) {
      free(this->pixels);
    }
  }

  /**
   * @brief The comparison operators compare the upper-left corner of the raster
   * locations. This is mostly meaningful when a raster is decomposed into
   * non-overlapping RasterChunks.
   *
   *
   */
  bool operator==(const RasterChunk &s);
  bool operator!=(const RasterChunk &s);
  bool operator<(const RasterChunk &s);
  bool operator>(const RasterChunk &s);
  bool operator<=(const RasterChunk &s);
  bool operator>=(const RasterChunk &s);

  /**
   * @brief 
   * This function reads the pixel values from the GDALDataset into the
   * RasterChunk.
   *
   * @param ds The raster to read from.
   *
   */
  PRB_ERROR Read(GDALDataset *ds);
  
  /**
   * @brief
   * This function writes the pixel values from the RasterChunk into the
   * file behind the GDALDataset.
   *
   * @param ds The GDALDataset to write to.
   * @param chunk The RasterChunk to write to the file.
   *
   */
  PRB_ERROR Write(GDALDataset *ds);

  Coordinate ChunkToRaster(Coordinate chunk_coordinate);
  Coordinate RasterToChunk(Coordinate raster_coordinate);

  std::string projection;
  /// Location of the chunk, in raster coordinates
  /** 
   * This variable represents the upper-left location of the raster chunk, in
   * the raster coordinates.
   */
  Coordinate raster_location;
  /// Upper-left corner, in projected coordinates
  /**
   * This variable represents the upper-left location of the raster chunk in
   * projected coordinates.
   */
  Coordinate ul_projected_corner;
  /// Size of pixel, in meters
  /**
   * This variable represents the size of the pixel in meters.
   */
  double pixel_size;  // in meters
  /// Number of rows
  int row_count;
  /// Number of columns
  int column_count;
  /// Datatype of pixel values
  GDALDataType pixel_type;
  /// Number of bands
  int band_count;
  /// GDAL geotransform
  double geotransform[6];
  /// Pointer to pixel values
  void *pixels;
};
}


#endif  // SRC_RASTERCHUNK_H_
