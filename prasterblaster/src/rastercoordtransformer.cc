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
//
/* ! \mainpage RasterCoordTransformer
 *
 *
 */

#include <ogr_api.h>
#include <ogr_spatialref.h>
#include <gdal.h>
#include <gdal_priv.h>

#include <cmath>

#include "src/reprojection_tools.h"
#include "src/resampler.h"

namespace librasterblaster {
RasterCoordTransformer::
RasterCoordTransformer(string source_projection,
                       Coordinate source_ul,
                       double source_pixel_size,
                       int source_row_count,
                       int source_column_count,
                       string destination_projection,
                       Coordinate destination_ul,
                       double destination_pixel_size) {
  init(source_projection,
       source_ul,
       source_pixel_size,
       source_row_count,
       source_column_count,
       destination_projection,
       destination_ul,
       destination_pixel_size);
  return;
}

RasterCoordTransformer::~RasterCoordTransformer() {
  return;
}

void RasterCoordTransformer::init(string source_projection,
                                  Coordinate source_ul,
                                  double source_pixel_size,
                                  int source_row_count,
                                  int source_column_count,
                                  string destination_projection,
                                  Coordinate destination_ul,
                                  double destination_pixel_size) {
  source_ul_ = source_ul;
  source_pixel_size_ = source_pixel_size;
  destination_ul_ = destination_ul;
  destination_pixel_size_ = destination_pixel_size;

  OGRSpatialReference source_sr, dest_sr, *geo_sr;
  char *source_wkt = strdup(source_projection.c_str());
  char *dest_wkt = strdup(destination_projection.c_str());

  source_sr.SetFromUserInput(source_projection.c_str());
  dest_sr.SetFromUserInput(destination_projection.c_str());
  geo_sr = source_sr.CloneGeogCS();

  OGRCoordinateTransformation *t =
      OGRCreateCoordinateTransformation(&source_sr,
                                        &dest_sr);
  src_to_geo = OGRCreateCoordinateTransformation(&source_sr,
                                                 geo_sr);
  geo_to_src = OGRCreateCoordinateTransformation(geo_sr,
                                                 &source_sr);
  free(source_wkt);
  free(dest_wkt);

  if (t != NULL) {
    ctrans = t;
  } else {
    printf("BAD!\n\n");
    return;
  }

  Coordinate ul, lr;
  ul = source_ul_;
  src_to_geo->Transform(1, &ul.x, &ul.y);
  lr.x = source_ul_.x + (source_pixel_size_ * source_column_count);
  lr.y = source_ul_.y - (source_pixel_size_ * source_row_count);
  src_to_geo->Transform(1, &lr.x, &lr.y);
  ul.x = -179.99;
  lr.x = 179.99;
  maximum_geographic_area_.ul = ul;
  maximum_geographic_area_.lr = lr;
  return;
}

Area RasterCoordTransformer::
Transform(Coordinate source, bool area_check) {
  Area value;
  Coordinate temp1, temp2;

  // PROJ.4 will malloc a temporary Z value if one is
  // not provided. By passing in a local variable
  // we prevent these unnecessary allocations.
  double unused;

  temp1.x = temp1.y = 0.0;

  value.ul = temp1;
  value.lr = temp1;

  temp1.x = (static_cast<double>(source.x) * source_pixel_size_) + source_ul_.x;
  temp1.y = source_ul_.y - (static_cast<double>(source.y) * source_pixel_size_);
  temp2.x = temp1.x;
  temp2.y = temp1.y;

  src_to_geo->TransformEx(1, &temp2.x, &temp2.y, &unused);
  geo_to_src->TransformEx(1, &temp2.x, &temp2.y, &unused);

  if ((area_check && (fabs(temp1.y - temp2.y) > 0.01))
      || fabs(temp1.x - temp2.x) > 0.01) {
    // Point is outside defined projection area, return no-value
    value.ul.x = -1.0;
    value.lr.x = -1.0;
    return value;
  }

  temp1.x = (static_cast<double>(source.x) * source_pixel_size_) + source_ul_.x;
  temp1.y = source_ul_.y - (static_cast<double>(source.y) * source_pixel_size_);
  temp2 = temp1;

  // Now we are going to assign temp1 as the UL of our pixel and
  // temp2 as LR
  temp2.x += sqrt(2 * source_pixel_size_ * source_pixel_size_);
  temp2.y -= sqrt(2 * source_pixel_size_ * source_pixel_size_);

  ctrans->TransformEx(1, &temp1.x, &temp1.y, &unused);
  ctrans->TransformEx(1, &temp2.x, &temp2.y, &unused);

  // temp1/temp2 now contain coords to input projection
  // Now convert to points in the raster coordinate space.
  temp1.x -= destination_ul_.x;
  temp1.y = destination_ul_.y - temp1.y;
  temp1.x /= destination_pixel_size_;
  temp1.y /= destination_pixel_size_;
  temp2.x -= destination_ul_.x;
  temp2.y = destination_ul_.y - temp2.y;
  temp2.x /= destination_pixel_size_;
  temp2.y /= destination_pixel_size_;

  value.ul = temp1;
  value.lr = temp2;

  // Check that entries are valid
  if (value.ul.x < 0.0
      || value.lr.x < 0.0
      || value.ul.y < 0.0
      || value.lr.y < 0.0) {
    value.ul.x = -1.0;
    value.lr.x = -1.0;
    return value;
  }

  // Now validate and round pixel values
  // Truncate values
  value.ul.x = floor(fabs(value.ul.x));
  value.ul.y = floor(fabs(value.ul.y));
  value.lr.x = floor(fabs(value.lr.x));
  value.lr.y = floor(fabs(value.lr.y));

  if (value.ul.x > value.lr.x) {
    value.lr.x = value.ul.x;
  }

  if (value.ul.y > value.lr.y) {
    value.lr.y = value.ul.y;
  }

  return value;
}

bool RasterCoordTransformer::ready() {
  return true;
}
}
