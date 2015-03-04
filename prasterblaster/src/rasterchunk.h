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

#include <gdal_priv.h>

#include <string>

#include "src/utils.h"

namespace librasterblaster {
/// A class representing an in-memory part of a raster.
class RasterChunk {
 public:
  /**
   * @brief
   * This function creates a RasterChunk.
   *
   * @param ds Dataset to create chunk from
   * @param chunk_area The inclusive area that the chunk should represent.
   *
   */
  static RasterChunk* CreateRasterChunk(GDALDataset *ds, Area chunk_area);

  /**
   * @brief
   *
   * This function creates a RasterChunk from the "destination" dataset by
   *   calculating the area in "destination" that corresponds to 
   *   "source_area" in "source" .
   *
   * @param destination Dataset to create the chunk from.
   * @param source Dataset that is used to find the area
   * @param source_area Area that is used to calculate area from 
   *                    destination
   *
   *
   */
  static RasterChunk* CreateRasterChunk(GDALDataset *destination,
                                        GDALDataset *source,
                                        Area source_area);

  /**
   * @brief
   * Copy constructor
   *
   * @param s Source RasterChunk to be copied
   */
  RasterChunk(const RasterChunk &s);

  /**
   * @brief 
   * This function reads the pixel values from the GDALDataset into the
   * RasterChunk.
   *
   * @param ds The raster to read from.
   * @param chunk The chunk to read to.
   *
   */

  RasterChunk& operator=(const RasterChunk &s);

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

  static PRB_ERROR ReadRasterChunk(GDALDataset *ds, RasterChunk *chunk);
  /**
   * @brief
   * This function writes the pixel values from the RasterChunk into the
   * file behind the GDALDataset.
   *
   * @param ds The GDALDataset to write to.
   * @param chunk The RasterChunk to write to the file.
   *
   */
  static PRB_ERROR WriteRasterChunk(GDALDataset *ds, RasterChunk *chunk);
  /// RasterChunk constructor
  RasterChunk() {
    this->pixels_ = NULL;
  }
  /// RasterChunk destructor
  /**
   * This destructor frees the memory, if any,  allocated at pixels_.
   */
  ~RasterChunk() {
    if (this->pixels_ != NULL) {
      free(this->pixels_);
    }
  }
  Coordinate ChunkToRaster(Coordinate chunk_coordinate);
  Coordinate RasterToChunk(Coordinate raster_coordinate);

  std::string projection_;
  /// Location of the chunk, in raster coordinates
  /** 
   * This variable represents the upper-left location of the raster chunk, in
   * the raster coordinates.
   */
  Coordinate raster_location_;
  /// Upper-left corner, in projected coordinates
  /**
   * This variable represents the upper-left location of the raster chunk in
   * projected coordinates.
   */
  Coordinate ul_projected_corner_;
  /// Size of pixel, in meters
  /**
   * This variable represents the size of the pixel in meters.
   */
  double pixel_size_;  // in meters
  /// Number of rows
  int row_count_;
  /// Number of columns
  int column_count_;
  /// Datatype of pixel values
  GDALDataType pixel_type_;
  /// Number of bands
  int band_count_;
  /// GDAL geotransform
  double geotransform_[6];
  /// Pointer to pixel values
  void *pixels_;
};
}


#endif  // SRC_RASTERCHUNK_H_
