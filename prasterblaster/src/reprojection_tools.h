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

#include <cmath>
#include <cstdint>
#include <functional>
#include <string>
#include <vector>

#include "rastercoordtransformer.h"
#include "resampler.h"
#include "utils.h"

/// Container namespace for librasterblaster project
namespace librasterblaster {
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
                             int output_tile_size,
                             double output_ratio = 1.0f,
                             double output_no_data_value = NAN);
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
                             double output_ratio,
                             int output_max_dimension);

PRB_ERROR CreateOutputRasterFile(GDALDataset *in,
                                 string output_filename,
                                 string output_srs,
                                 int64_t output_columns,
                                 int64_t output_rows,
                                 double output_pixel_size,
                                 Area output_projected_area,
                                 int output_tile_size,
                                 double output_no_data_value);

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

bool ReprojectChunk(RasterChunk& source,
                    RasterChunk& destination,
                    string fill_value,
                    RESAMPLER resampler);

/** @cond DOXYHIDE **/

template<typename T>
std::function<T(RasterChunk&, Area, float)> GetResampler(RESAMPLER type);

template <typename T>
bool ReprojectChunkType(RasterChunk& source,
                        RasterChunk& destination,
                        T fill_value,
                        std::function<T(RasterChunk&, Area, float)> resampler,
                        int filter_support);
/** @endcond **/

}


#endif  // SRC_REPROJECTION_TOOLS_H_
