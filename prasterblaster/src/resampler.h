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
// Several resampling implementations for use with the ReprojectChunk
//
//

#ifndef SRC_RESAMPLER_H_
#define SRC_RESAMPLER_H_

#include <gdal.h>

#include "src/rasterchunk.h"
#include "src/utils.h"

namespace librasterblaster {
/**
 * @brief Supported resampling algorithms
 *
 */
enum RESAMPLER {
  NEAREST, /** @brief Nearest-neighbor */
  MIN,     /** @brief Minimum value */
  MAX,     /** @brief Maximum value */
  MEAN,    /** @brief Arithmetic mean */
};

/** @cond DOXYHIDE */
template <typename T>
T Max(RasterChunk *input,
      Area pixel_area) {
  T *pixels = static_cast<T*>(input->pixels_);
  T temp = pixels[static_cast<int>(pixel_area.ul.y)
                  * input->column_count_
                  + static_cast<int>(pixel_area.ul.x)];
  T temp2 = 0;
  for (int x = pixel_area.ul.x; x <= pixel_area.lr.x; ++x) {
    for (int y = pixel_area.ul.y; y <= pixel_area.lr.y; ++y) {
      temp2 = pixels[y * input->column_count_  + x];

      if (temp2 > temp) {
        temp = temp2;
      }
    }
  }

  return temp;
}

template <typename T>
T Min(RasterChunk *input,
      Area pixel_area) {
  T *pixels = static_cast<T*>(input->pixels_);
  T temp = pixels[static_cast<int>(pixel_area.ul.y)
                  * input->column_count_
                  + static_cast<int>(pixel_area.ul.x)];
  T temp2 = 0;
  for (int x = pixel_area.ul.x; x <= pixel_area.lr.x; ++x) {
    for (int y = pixel_area.ul.y; y <= pixel_area.lr.y; ++y) {
      temp2 = pixels[y * input->column_count_  + x];

      if (temp2 < temp) {
        temp = temp2;
      }
    }
  }

  return temp;
}

template <typename T>
T Mean(RasterChunk *input,
       Area pixel_area) {
  T temp = 0;
  T *pixels = static_cast<T*>(input->pixels_);
  for (int x = pixel_area.ul.x; x <= pixel_area.lr.x; ++x) {
    for (int y = pixel_area.ul.y; y <= pixel_area.lr.y; ++y) {
      temp += pixels[y * input->column_count_ + x];
    }
  }

  temp /= (pixel_area.lr.x - pixel_area.ul.x)
      * (pixel_area.lr.y - pixel_area.ul.y);

  return temp;
}

template <typename T>
T Median(Coordinate input_ul,
         Coordinate input_lr,
         int input_column_count,
         T* input_pixels) {
}

template <typename T>
T Mode(Coordinate input_ul,
       Coordinate input_lr,
       int input_column_count,
       T* input_pixels) {
}

template <typename T>
T Sum(Coordinate input_ul,
      Coordinate input_lr,
      int input_column_count,
      T* input_pixels) {
}


template <typename T>
T Bilinear(Coordinate input_ul,
           Coordinate input_lr,
           int input_column_count,
           T* input_pixels) {
}
}
/** @cond DOXYHIDE */

#endif  // SRC_RESAMPLER_H_

