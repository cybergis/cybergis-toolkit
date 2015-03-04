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

#include <gdal.h>
#include <string.h>

#include "src/reprojection_tools.h"
#include "src/rasterchunk.h"
#include "src/utils.h"

namespace librasterblaster {
RasterChunk* RasterChunk::CreateRasterChunk(GDALDataset *ds, Area chunk_area) {
  RasterChunk *temp = new RasterChunk;
  double gt[6];

  if (chunk_area.ul.x == -1.0) {  // Create a chunk with a single value for
                                  // resampling no data area.
          chunk_area.ul.x
              = chunk_area.ul.y
              = chunk_area.lr.x
              = chunk_area.lr.y
              = 0.0;
  }

  ds->GetGeoTransform(gt);
  ds->GetGeoTransform(temp->geotransform_);

  temp->projection_ = ds->GetProjectionRef();
  temp->raster_location_ = chunk_area.ul;
  temp->ul_projected_corner_ = Coordinate(gt[0]+(chunk_area.ul.x*gt[1]),
                                          gt[3]-(chunk_area.ul.y*gt[1]),
                                          UNDEF);
  temp->pixel_size_ = gt[1];
  temp->row_count_ = chunk_area.lr.y - chunk_area.ul.y + 1;
  temp->column_count_ = chunk_area.lr.x - chunk_area.ul.x + 1;
  temp->pixel_type_ = ds->GetRasterBand(1)->GetRasterDataType();
  temp->band_count_ = ds->GetRasterCount();
  temp->pixels_ = NULL;

  size_t buffer_size = (temp->row_count_ * temp->column_count_);
  temp->pixels_ = static_cast<unsigned char*>
      (calloc(buffer_size, GDALGetDataTypeSize(temp->pixel_type_)/8));

  if (temp->pixels_ == NULL) {
    fprintf(stderr, "Allocation error!\n");
    delete temp;
    return NULL;
  }

  return temp;
}

RasterChunk* RasterChunk::CreateRasterChunk(GDALDataset *input_raster,
                                            GDALDataset *output_raster,
                                            Area output_area) {
  // The RasterMinbox function calculates what part of the input raster
  // matches the given output partition.
  Area in_area = librasterblaster::RasterMinbox(output_raster,
                              input_raster,
                              output_area);
  return CreateRasterChunk(input_raster, in_area);
}

RasterChunk::RasterChunk(const RasterChunk &s) {
  projection_ = s.projection_;
  raster_location_ = s.raster_location_;
  ul_projected_corner_ = s.ul_projected_corner_;
  pixel_size_ = s.pixel_size_;
  row_count_ = s.row_count_;
  column_count_ = s.column_count_;
  pixel_type_ = s.pixel_type_;
  band_count_ = s.band_count_;
  memcpy(geotransform_, s.geotransform_, 6*sizeof(double));

  size_t pixel_buffer_size = static_cast<size_t>(row_count_)
      * static_cast<size_t>(column_count_)
      * static_cast<size_t>(GDALGetDataTypeSize(pixel_type_)/8)
      * static_cast<size_t>(band_count_);

  memcpy(pixels_, s.pixels_, pixel_buffer_size);
}

RasterChunk& RasterChunk::operator=(const RasterChunk &s) {
  if (this == &s) {
    return *this;
  }
  projection_ = s.projection_;
  raster_location_ = s.raster_location_;
  ul_projected_corner_ = s.ul_projected_corner_;
  pixel_size_ = s.pixel_size_;
  row_count_ = s.row_count_;
  column_count_ = s.column_count_;
  pixel_type_ = s.pixel_type_;
  band_count_ = s.band_count_;
  memcpy(geotransform_, s.geotransform_, 6*sizeof(double));

  size_t pixel_buffer_size = static_cast<size_t>(row_count_)
      * static_cast<size_t>(column_count_)
      * static_cast<size_t>(GDALGetDataTypeSize(pixel_type_)/8)
      * static_cast<size_t>(band_count_);

  memcpy(pixels_, s.pixels_, pixel_buffer_size);

  return *this;
}

bool RasterChunk::operator==(const RasterChunk &s) {

  // Maybe do some sort of normalization with the projections before
  // comparing. For now do a character-by-character comparison.
  if (projection_ != s.projection_
      || raster_location_ != s.raster_location_
      || ul_projected_corner_ != s.ul_projected_corner_
      || pixel_size_ != s.pixel_size_
      || row_count_ != s.row_count_
      || column_count_ != s.column_count_
      || pixel_type_ != s.pixel_type_
      || band_count_ != s.band_count_) {
    return false;
  }

  for (int i = 0; i < 6; ++i) {
    if (geotransform_[i] != s.geotransform_[i]) {
      return false;
    }
  }
  // Ignore pixel values for now
  return true;
}

bool RasterChunk::operator!=(const RasterChunk &s) {
  if (*this == s) {
    return false;
  }

  return true;
}

bool RasterChunk::operator<(const RasterChunk &s) {
  if (raster_location_.x < s.raster_location_.x
      || raster_location_.y < s.raster_location_.y) {
    return true;
  }

  return false;
}

bool RasterChunk::operator>(const RasterChunk &s) {
  if (raster_location_.x > s.raster_location_.x
      || raster_location_.y > s.raster_location_.y) {
    return true;
  }

  return false;
}

bool RasterChunk::operator<=(const RasterChunk &s) {
  if (*this == s) {
    return true;
  }

  if (*this < s) {
    return true;
  }

  return false;
}

bool RasterChunk::operator>=(const RasterChunk &s) {
  if (*this == s) {
    return true;
  }

  if (*this > s) {
    return true;
  }

  return false;
}

PRB_ERROR RasterChunk::ReadRasterChunk(GDALDataset *ds, RasterChunk *chunk) {
  // Read area of raster

  if (ds == NULL) {
    return PRB_BADARG;
  }

  if (ds->RasterIO(GF_Read,
                   chunk->raster_location_.x,
                   chunk->raster_location_.y,
                   chunk->column_count_,
                   chunk->row_count_,
                   chunk->pixels_,
                   chunk->column_count_,
                   chunk->row_count_,
                   chunk->pixel_type_,
                   chunk->band_count_,
                   NULL,
                   0, 0, 0) != CE_None) {
    return PRB_IOERROR;
  }

  return PRB_NOERROR;
}

PRB_ERROR RasterChunk::WriteRasterChunk(GDALDataset *ds, RasterChunk *chunk) {
  if (ds->RasterIO(GF_Write,
                   chunk->raster_location_.x,
                   chunk->raster_location_.y,
                   chunk->column_count_,
                   chunk->row_count_,
                   chunk->pixels_,
                   chunk->column_count_,
                   chunk->row_count_,
                   chunk->pixel_type_,
                   chunk->band_count_,
                   NULL, 0, 0, 0) != CE_None) {
    // Error!
    fprintf(stderr, "Error writing RasterChunk %p\n", chunk->pixels_);
    return PRB_IOERROR;
  }

  ds->FlushCache();

  return PRB_NOERROR;
}

Coordinate RasterChunk::ChunkToRaster(Coordinate chunk_coordinate) {
  return Coordinate(chunk_coordinate.x + raster_location_.x,
                    chunk_coordinate.y + raster_location_.y,
                    chunk_coordinate.units);
}

Coordinate RasterChunk::RasterToChunk(Coordinate raster_coordinate) {
  return Coordinate(raster_coordinate.x - raster_location_.x,
                    raster_coordinate.y - raster_location_.y,
                    raster_coordinate.units);
}
}
