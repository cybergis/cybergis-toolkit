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

#include <cstring>
#include <gdal.h>

#include "reprojection_tools.h"
#include "rasterchunk.h"
#include "utils.h"

namespace librasterblaster {
RasterChunk::RasterChunk(GDALDataset *ds, Area chunk_area) {
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
  ds->GetGeoTransform(geotransform);

  projection = ds->GetProjectionRef();
  raster_location = chunk_area.ul;
  ul_projected_corner = Coordinate(gt[0]+(chunk_area.ul.x*gt[1]),
                                          gt[3]-(chunk_area.ul.y*gt[1]));

  pixel_size = gt[1];
  row_count = chunk_area.lr.y - chunk_area.ul.y + 1;
  column_count = chunk_area.lr.x - chunk_area.ul.x + 1;
  pixel_type = ds->GetRasterBand(1)->GetRasterDataType();
  band_count = ds->GetRasterCount();

  size_t buffer_size = row_count * column_count;
  pixels = static_cast<uint8_t*>
      (calloc(buffer_size, GDALGetDataTypeSize(pixel_type)/8));

  if (pixels == NULL) {
    fprintf(stderr, "Allocation error!\n");
  }
}

bool RasterChunk::operator==(const RasterChunk &s) {

  // Maybe do some sort of normalization with the projections before
  // comparing. For now do a character-by-character comparison.
  if (projection != s.projection
      || raster_location != s.raster_location
      || ul_projected_corner != s.ul_projected_corner
      || pixel_size != s.pixel_size
      || row_count != s.row_count
      || column_count != s.column_count
      || pixel_type != s.pixel_type
      || band_count != s.band_count) {
    return false;
  }

  for (int i = 0; i < 6; ++i) {
    if (geotransform[i] != s.geotransform[i]) {
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
  if (raster_location.x < s.raster_location.x
      || raster_location.y < s.raster_location.y) {
    return true;
  }

  return false;
}

bool RasterChunk::operator>(const RasterChunk &s) {
  if (raster_location.x > s.raster_location.x
      || raster_location.y > s.raster_location.y) {
    return true;
  }

  return false;
}

bool RasterChunk::operator<=(const RasterChunk &s) {
  if (*this == s) {
      return true;
  }

  return *this < s;
}

bool RasterChunk::operator>=(const RasterChunk &s) {
  if (*this == s) {
    return true;
  }

  return *this > s;
}

PRB_ERROR RasterChunk::Read(GDALDataset *ds) {
  // Read area of raster

  if (ds == NULL) {
    return PRB_BADARG;
  }

  if (ds->RasterIO(GF_Read,
                   raster_location.x,
                   raster_location.y,
                   column_count,
                   row_count,
                   pixels,
                   column_count,
                   row_count,
                   pixel_type,
                   band_count,
                   NULL,
                   0, 0, 0) != CE_None) {
    return PRB_IOERROR;
  }

  return PRB_NOERROR;
}

PRB_ERROR RasterChunk::Write(GDALDataset *ds) {
  if (ds->RasterIO(GF_Write,
                   raster_location.x,
                   raster_location.y,
                   column_count,
                   row_count,
                   pixels,
                   column_count,
                   row_count,
                   pixel_type,
                   band_count,
                   NULL, 0, 0, 0) != CE_None) {
    fprintf(stderr, "Error while writing RasterChunk %p\n", pixels);
    return PRB_IOERROR;
  }

  ds->FlushCache();

  return PRB_NOERROR;
}

Coordinate RasterChunk::ChunkToRaster(Coordinate chunk_coordinate) {
  return Coordinate(chunk_coordinate.x + raster_location.x,
                    chunk_coordinate.y + raster_location.y);
}

Coordinate RasterChunk::RasterToChunk(Coordinate raster_coordinate) {
  return Coordinate(raster_coordinate.x - raster_location.x,
                    raster_coordinate.y - raster_location.y);
}
}
