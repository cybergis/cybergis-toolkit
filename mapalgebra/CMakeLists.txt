cmake_minimum_required(VERSION 2.8)

project("CyberGIS Map Algebra")

# Check for MPI
find_package(MPI REQUIRED)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})

# Check for GDAL
find_package(GDAL)
if (NOT GDAL_FOUND)
  find_program(GDAL_CONFIG gdal-config)
  if(GDAL_CONFIG)
    message("-- Setting with gdal-config")
    exec_program(${GDAL_CONFIG} ARGS --prefix OUTPUT_VARIABLE GDAL_PREFIX)
    set(GDAL_INCLUDE_DIR "${GDAL_PREFIX}/include")
    set(GDAL_LIBRARY "${GDAL_PREFIX}/lib/libgdal.so")
  else()
    message(FATAL_ERROR "GDAL NOT FOUND")
  endif()
endif()
include_directories(SYSTEM ${GDAL_INCLUDE_DIR})

add_subdirectory(src)

###### Installing
install (PROGRAMS "${CMAKE_CURRENT_BINARY_DIR}/src/mapalg"
                  "${CMAKE_CURRENT_BINARY_DIR}/src/mapalg-kernel" DESTINATION bin)

