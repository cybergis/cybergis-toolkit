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

/*
 *
 *
 */
#ifndef SRC_RASTERCOORDTRANSFORMER_H_
#define SRC_RASTERCOORDTRANSFORMER_H_

#include <ogr_spatialref.h>

#include <string>

#include "src/utils.h"

using std::string;


namespace librasterblaster {
/// Raster Coordinate transformation class
/*
 * This class implements the transformation of raster coordinates between two raster spaces with different projections and scales.
 */
class RasterCoordTransformer {
 public:
  // ! A constructor
  /* ! 

    This constructor takes six parameters. The first three
    describe the source raster: a projection, upper-left
    coordinate in projected coordinates, and a pixel size in meters.

    The second set of three parameters are the same as the first
    three but for the destination raster.

  */
  RasterCoordTransformer(string source_projection,
                         Coordinate source_ul,
                         double source_pixel_size,
                         int source_row_count,
                         int souce_column_count,
                         string destination_projection,
                         Coordinate destination_ul,
                         double destination_pixel_size);


  ~RasterCoordTransformer();

  // ! A normal member taking a single argument and returning an Area struct.
  /*

    This function takes a coordinate in the source raster space
    and maps is to an area in the destination raster space. The
    returned area consists of two points: an upper-left
    coordinate, and a lower-right. The coordinates are
    _inclusive_. So if the area consists of a single point, the
    upper-left == lower-right.

    \param source a Coordinate struct that specifies the point in the source raster space to map to the destination raster space.
  */
  Area Transform(Coordinate source, bool area_check = true);

  // ! A normal member function taking no arguments
  /*
    This function returns a boolean value indicating whether the
    RasterCoordTransformer constructed corrrectly and is ready to use.
   */
  bool ready();

 private:
  void init(string source_projection,
            Coordinate source_ul,
            double source_pixel_size,
            int source_row_count,
            int source_column_count,
            string destination_projection,
            Coordinate destination_ul,
            double destination_pixel_size);

  OGRCoordinateTransformation *ctrans, *src_to_geo, *geo_to_src;
  Area maximum_geographic_area_;
  Coordinate source_ul_;
  double source_pixel_size_;
  Coordinate destination_ul_;
  double destination_pixel_size_;
};
}

#endif  // SRC_RASTERCOORDTRANSFORMER_H_
