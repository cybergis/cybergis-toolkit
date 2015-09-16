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

#include <cmath>
#include <functional>

#include <gdal.h>

#include "rasterchunk.h"
#include "utils.h"

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
  BILINEAR,/** @brief Bilinear */
  BICUBIC, /** @brief Bicubic (Catmull-Rom spline) */
  LANCZOS, /** @brief Lanczos */
};

inline float bilinear_filter(float x) {
  if (x < 0.0)
    x = -x;

  if (x < 1.0)
    return 1.0 - x;

  return 0.0;
}

inline float bicubic_filter(float x) {
  const float a = -0.5;

  if (x < 0.0)
    x = -x;

  if (x < 1.0)
    return ((a + 2.0) * x - (a + 3.0)) * x*x + 1;

  if (x < 2.0)
    return (((x - 5.0) * x + 8) * x - 4.0) * a;

  return 0.0;
}

inline float sinc_filter(float x) {
  if (x == 0.0)
    return 1.0;

  x = x * M_PI;
  return sin(x) / x;
}

inline float lanczos_filter(float x) {
  if (-3.0 <= x && x <= 3.0) {
    return sinc_filter(x) * sinc_filter(x / 3);
  }

  return 0.0;
}

/** @cond DOXYHIDE */
template <typename T>
T Max(RasterChunk& input, Area pixel_area, float) {
  T *pixels = static_cast<T*>(input.pixels);
  T temp = pixels[static_cast<int>(pixel_area.ul.y)
                  * input.column_count
                  + static_cast<int>(pixel_area.ul.x)];
  T temp2 = 0;
  for (int y = pixel_area.ul.y; y <= pixel_area.lr.y; ++y) {
    for (int x = pixel_area.ul.x; x <= pixel_area.lr.x; ++x) {
      temp2 = pixels[y * input.column_count + x];

      if (temp2 > temp) {
        temp = temp2;
      }
    }
  }

  return temp;
}

template <typename T>
T Min(RasterChunk& input, Area pixel_area, float) {
  T *pixels = static_cast<T*>(input.pixels);
  T temp = pixels[static_cast<int>(pixel_area.ul.y)
    * input.column_count
    + static_cast<int>(pixel_area.ul.x)];

  T temp2 = 0;

  for (int y = pixel_area.ul.y; y <= pixel_area.lr.y; ++y) {
    for (int x = pixel_area.ul.x; x <= pixel_area.lr.x; ++x) {
      temp2 = pixels[y * input.column_count + x];

      if (temp2 < temp) {
        temp = temp2;
      }
    }
  }

  return temp;
}

template <typename T>
T Mean(RasterChunk& input, Area pixel_area, float) {
  T temp = 0;
  T *pixels = static_cast<T*>(input.pixels);

  int cells = 0;

  for (int y = pixel_area.ul.y; y <= pixel_area.lr.y; ++y) {
    for (int x = pixel_area.ul.x; x <= pixel_area.lr.x; ++x) {
      temp += pixels[(int64_t) y * input.column_count + (int64_t) x];
      cells++;
    }
  }

  temp /= cells;

  return temp;
}

template <typename T, typename F>
T Filter(RasterChunk& input, Area pixel_area, F filter, int support, float scale_factor) {
  T temp = 0;
  T *pixels = static_cast<T*>(input.pixels);

  float ss = support / scale_factor;

  // Assume center of source cell is in the center of the input pixel area
  float x_center = pixel_area.ul.x + (pixel_area.lr.x - pixel_area.ul.x) / 2.0;
  float y_center = pixel_area.ul.y + (pixel_area.lr.y - pixel_area.ul.y) / 2.0;

  float total_weight = 0.0;

  for (int y = pixel_area.ul.y; y <= pixel_area.lr.y; ++y) {
    float y_filter_point = (y + 0.5 - y_center) * ss;
    float y_weight = filter(y_filter_point) * ss;

    for (int x = pixel_area.ul.x; x <= pixel_area.lr.x; ++x) {
      float x_filter_point = (x + 0.5 - x_center) * ss;
      float x_weight = filter(x_filter_point) * ss;

      float weight = x_weight * y_weight;

      temp += pixels[(int64_t) y * input.column_count + (int64_t) x] * weight;
      total_weight += weight;
    }
  }

  return temp / total_weight;
}

template<typename T>
T Bilinear(RasterChunk& input, Area pixel_area, float scale_factor) {
  return Filter<T>(input, pixel_area, bilinear_filter, 1, scale_factor);
}

template<typename T>
T Bicubic(RasterChunk& input, Area pixel_area, float scale_factor) {
  return Filter<T>(input, pixel_area, bicubic_filter, 2, scale_factor);
}

template<typename T>
T Lanczos(RasterChunk& input, Area pixel_area, float scale_factor) {
  return Filter<T>(input, pixel_area, lanczos_filter, 3, scale_factor);
}
}

/** @cond DOXYHIDE */

#endif  // SRC_RESAMPLER_H_

