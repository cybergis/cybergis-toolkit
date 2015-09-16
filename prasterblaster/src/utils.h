/**
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
 * Helper utilities and definitions
 *
 */

#ifndef SRC_UTILS_H_
#define SRC_UTILS_H_

#include <cmath>

namespace librasterblaster {
/** An enum of possible return values */
enum PRB_ERROR {
  PRB_NOERROR, /*!< No error occurred */
  PRB_IOERROR, /*!< Error communicating or performing file I/O */
  PRB_BADARG,  /*!< Bad argument provided */
  PRB_PROJERROR, /*!< Error with projection specification */
};

/*! Coordinate struct
   This class provides a more readable way of storing and passing coordinate
   parameters for the Transformer class. It stores x and y as doubles and units
   corresponds to the GCTP enumeration as defined in constants.h
*/
struct Coordinate {
  /*! Default Constructor
    0's all attributes.
  */
  Coordinate():x(0.0), y(0.0) {}
  /*! Full Constructor
    Sets all attributes in the Coordinate according to the parameters.
  */
  Coordinate(double xx, double yy)
      : x(xx), y(yy) {}

  bool operator==(const Coordinate &c) {
    return (x == c.x && y == c.y);
  }

  bool operator!=(const Coordinate &c) { return !(*this == c); }

  double x;
  double y;
};

/** 
 * @class 
 */
struct Area {
  /// A constructor
  /** 
   * @brief This constructor initializes points to zero.
   */
  Area():ul(), lr() {}
  /// A constructor
  /**
   * @brief This constructor allows both coordinates to be given initial values.
   * @param ulx Value for the upper-left x 
   * @param uly Value for the upper-left y
   * @param ulx Value for the lower-right x 
   * @param ulx Value for the lower-right y
   */
Area(double ulx,
     double uly,
     double lrx,
     double lry) : ul(ulx, uly), lr(lrx, lry){}
  /// Upper-left coordinate of area
  Coordinate ul;
  /// Lower-right coordinate of area
  Coordinate lr;
};
}

#endif  // SRC_UTILS_H_
