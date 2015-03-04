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
 * Helper functions to create and manipulate projections and projected rasters.
 *
 */

#ifndef SRC_REPROJECTION_TOOLS_H_
#define SRC_REPROJECTION_TOOLS_H_

#include <string>
#include <vector>

#include "src/rastercoordtransformer.h"
#include "src/resampler.h"
#include "src/std_int.h"
#include "src/utils.h"

/// Container namespace for librasterblaster project
namespace librasterblaster {
/** \cond DOXYHIDE **/
int simplerandom(int i);
/** \cond DOXYHIDE **/

#define GEN_RESAMPLER_CASES(C_PIXEL_TYPE)   \
        switch (resampler) {               \
        case MIN: \
          return ReprojectChunkType<C_PIXEL_TYPE>(source, \
                                                   destination, \
                                                   static_cast<C_PIXEL_TYPE>(fvalue), \
                                                   &(Min<C_PIXEL_TYPE>)); \
          break; \
        case MAX: \
          return ReprojectChunkType<C_PIXEL_TYPE>(source, \
                                                   destination, \
                                                   static_cast<C_PIXEL_TYPE>(fvalue), \
                                                   &(Max<C_PIXEL_TYPE>)); \
          break; \
    case NEAREST: \
    default: \
          return ReprojectChunkType<C_PIXEL_TYPE>(source, \
                                                   destination, \
                                                   static_cast<C_PIXEL_TYPE>(fvalue), \
                                                   NULL); \
      } \
      break; \


/**
 * @brief Creates an output raster based on a input and a new projection
 * 
 * This function creates a new raster file at the path
 * output_filename, with projection specified by output_srs. The
 * minbox of in is calculated with the new projection.
 *
 * @param in The GDALDataset that represents the input file.
 * @param output_filename The path the new raster will be created at.
 * @param output_tile_size Size in pixels of one dimension of the tiles of the output raster
 * @param output_srs String with a projection specification (WKT, proj4, EPSG) suitable for
 *        OGRSpatialReference->SetFromUserInput()
 * @param output_pixel_size Side of square in meters of projected coordinates that each pixel represents
 *
 */
PRB_ERROR CreateOutputRaster(GDALDataset *in,
                             string output_filename,
                             string output_srs,
                             int output_tile_size);
/**
 * @brief Creates an output raster based on an input raster, a new projection,
 * and a maximum pixel dimension. This is to be used when the dimensions of the
 * output raster are important e.g. when generating a sample output.
 *
 * @param in The GDALDataset that represents the input file.
 * @param output_filename The path the new raster will be created at.
 * @param output_srs String with a projection specification (WKT, proj4, EPSG) suitable for
 *        OGRSpatialReference->SetFromUserInput()
 * @param output_tile_size Size in pixels of one dimension of the tiles of the output raster.
 *        The output raster will be at most (output_tile_size X output_tile_size) size in pixels.
 */
PRB_ERROR CreateOutputRaster(GDALDataset *in,
                             string output_filename,
                             string output_srs,
                             int output_tile_size,
                             int output_max_dimension);

PRB_ERROR CreateOutputRasterFile(GDALDataset *in,
                                 string output_filename,
                                 string output_srs,
                                 int64_t output_columns,
                                 int64_t output_rows,
                                 double output_pixel_size,
                                 Area output_projected_area,
                                 int output_tile_size);
/** \cond DOXYHIDE **/
bool partition_compare(Area a, Area b);
/** \endcond **/
std::vector<Area> BlockPartition(int rank,
                                 int process_count,
                                 int row_count,
                                 int column_count,
                                 int tile_size,
                                 int partition_size);
/** \cond DOXYHIDE **/
void SearchAndUpdate(Area input_area,
                     string input_srs,
                     string output_srs,
                     double input_ulx,
                     double input_uly,
                     double input_pixel_size,
                     Area *output_area);
/** \endcond **/
Area ProjectedMinbox(Coordinate input_ul_corner,
                     string input_srs,
                     double input_pixel_size,
                     int input_row_count,
                     int input_column_count,
                     string output_srs);
/**
 * @brief RasterMinbox finds the equivalent minbox in the source raster of the
 *        given area in the destination raster
 *
 * @param source Dataset in which you want a minbox
 * @param destination Dataset which you are providing an area for
 * @param destination_raster_area Area in destination that you want mapped 
 *        to a minbox in source
 *
 */
Area RasterMinbox(GDALDataset *source,
                  GDALDataset *destination,
                  Area destination_raster_area);

Area RasterMinbox2(string source_projection,
                  Coordinate source_ul,
                  double source_pixel_size,
                  int source_row_count,
                  int source_column_count,
                  string destination_projection,
                  Coordinate destination_ul,
                  double destination_pixel_size,
                  int destination_row_count,
                  int destination_column_count,
                  Area destination_raster_area);
/**
 * \brief This function takes two RasterChunk pointers and performs
 *        reprojection and resampling
 * \param source Pointer to the RasterChunk to reproject from
 * \param destination Pointer to the RasterChunk to reproject to
 * \param fillvalue std::string that will be interpreted to be the fill value
 * \param resampler The resampler that should be used
 *
 * @return Returns a bool indicating success or failure.
 */

bool ReprojectChunk(RasterChunk *source,
                    RasterChunk *destination,
                    string fillvalue,
                    RESAMPLER resampler);
/** @cond DOXYHIDE **/
template <class pixelType>
bool ReprojectChunkType(RasterChunk *source,
                        RasterChunk *destination,
                        pixelType fillvalue,
                        pixelType (*resampler)(RasterChunk*,
                                               Area)) {
  Coordinate temp1, temp2;
  std::vector<char> inraster, outraster;


  Area pixelArea;

  RasterCoordTransformer rt(destination->projection_,
                            destination->ul_projected_corner_,
                            destination->pixel_size_,
                            destination->row_count_,
                            destination->column_count_,
                            source->projection_,
                            source->ul_projected_corner_,
                            source->pixel_size_);


  for (int chunk_y = 0; chunk_y < destination->row_count_; ++chunk_y)  {
    for (int chunk_x = 0; chunk_x < destination->column_count_; ++chunk_x) {
      temp1.x = chunk_x;
      temp1.y = chunk_y;

      pixelArea = rt.Transform(temp1);

      if (pixelArea.ul.x == -1.0 || (pixelArea.ul.x > source->column_count_ - 1)
          || (pixelArea.lr.y > source->row_count_ - 1)) {
        // The pixel is outside of the projected area
        reinterpret_cast<pixelType*>(destination->pixels_)
            [chunk_x + chunk_y * destination->column_count_] = fillvalue;
        continue;
      }

      temp1 = pixelArea.ul;
      temp2 = pixelArea.lr;

      int64_t ul_x = static_cast<int64_t>(temp1.x);
      int64_t ul_y = static_cast<int64_t>(temp1.y);
      int64_t lr_x = static_cast<int64_t>(temp2.x);
      int64_t lr_y = static_cast<int64_t>(temp2.y);

      if (ul_x < 0) {
        ul_x = 0;
      }

      if (ul_y < 0) {
        ul_y = 0;
      }

      if (lr_x > (source->column_count_ - 1)) {
        lr_x = source->column_count_ - 1;
      }

      if (ul_y > (source->row_count_ - 1)) {
        ul_y = source->row_count_ - 1;
      }

      // Perform resampling...
      int64_t dest_offset = chunk_x + chunk_y * destination->column_count_;
      int64_t src_offset = ul_x + ul_y * source->column_count_;

      if ((resampler == NULL) || ((ul_x <= lr_x) || (lr_y <= ul_x))) {
        // ul/lr do not enclose an area, use NN
        if (ul_x + 1 > source->column_count_ || ul_y + 1 > source->row_count_) {
          // TODO(dmattli) FIX THIS
        }
            reinterpret_cast<pixelType*>(destination->pixels_)[dest_offset] =
                reinterpret_cast<pixelType*>(source->pixels_)[src_offset];
        continue;
      }

      Area ia = Area(ul_x, ul_y, lr_x, lr_y);
      reinterpret_cast<pixelType*>(destination->pixels_)[dest_offset] =
          resampler(source, ia);
    }
  }

  return true;
}
/** @endcond **/
}


#endif  // SRC_REPROJECTION_TOOLS_H_
