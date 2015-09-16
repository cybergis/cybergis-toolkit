cmake_minimum_required(VERSION 2.8)

set(CMAKE_DISABLE_SOURCE_CHANGES ON)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)

project(pRasterBlaster)

option(BundleGDAL "Bundle private GDAL and PROJ" OFF)

# Set release build as default
if (NOT CMAKE_BUILD_TYPE)
		message(STATUS "No build type selected, default to Release")
		set(CMAKE_BUILD_TYPE "Release")
endif()

# Set path to additional CMake modules
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

# Set path to bundled GDAL and proj
set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} ${CMAKE_SOURCE_DIR}/src/gdal)

include_directories(${CMAKE_SOURCE_DIR}/src/gtest/include)

include(FindCXX11)

# Set CXXFLAGS
set(CMAKE_CXX_FLAGS "-D__PRB_SRC_DIR__=${CMAKE_SOURCE_DIR} ${CMAKE_CXX_FLAGS} -W -Wall -Wextra -Wcast-align -Wpointer-arith -Wsign-compare -Wformat=2 -Wno-format-y2k  -Wmissing-braces -Wparentheses -Wtrigraphs -Wstrict-aliasing=2")

# MPI
find_package(MPI REQUIRED)
set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})
include_directories(${MPI_INCLUDE_PATH})

# 4.0 is required for BigTIFF support
# It will be used as a backup if GDAL does not have internal libtiff
find_package(TIFF 4.0.0)

find_package(GDAL REQUIRED)
find_package(Proj REQUIRED)

if(BundleGDAL)
    set(CMAKE_INSTALL_RPATH "\$ORIGIN/../lib/prasterblaster:\$ORIGIN/../lib")
endif()

# include_directories(${GDAL_INCLUDE_DIR} ${PROJ_INCLUDE_DIR} ${TIFF_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/src/gdal/gdal-1.11.0/frmts/gtiff/libtiff/)
include_directories(${GDAL_INCLUDE_DIR} ${TIFF_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/src/gdal/gdal-1.11.0/frmts/gtiff/libtiff/)

add_library(sptw SHARED src/demos/sptw.cc)
add_library(rasterblaster SHARED src/configuration.cc src/rastercoordtransformer.cc 
  src/reprojection_tools.cc src/rasterchunk.cc)
add_library(prasterblaster SHARED src/demos/prasterblaster-pio.cc)

target_link_libraries(sptw ${GDAL_LIBRARY} ${MPI_LIBRARIES} ${TIFF_LIBRARY})
target_link_libraries(rasterblaster ${GDAL_LIBRARY} ${PROJ_LIBRARY})
target_link_libraries(prasterblaster rasterblaster sptw ${MPI_LIBRARIES})

add_executable(prasterblasterpio src/demos/prasterblaster-main.cc)
target_link_libraries(prasterblasterpio sptw rasterblaster prasterblaster ${MPI_CXX_LIBRARIES})

add_executable(prasterblaster-simple src/demos/prasterblaster-simple.cc)
target_link_libraries(prasterblaster-simple rasterblaster prasterblaster)

add_subdirectory(src/gtest)
add_executable(tests tests/systemtest.cc tests/check_reprojection_tools.cc tests/rastercompare.cc)
target_link_libraries(tests gtest rasterblaster sptw prasterblaster)

# Add a target to generate API documentation with Doxygen
find_package(Doxygen 1.8)
if(DOXYGEN_FOUND)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
  add_custom_target(doc
    ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating API documentation with Doxygen" VERBATIM
    )
endif(DOXYGEN_FOUND)

# Add a target to run the 'tests' binary
add_custom_target(check ./tests
                  DEPENDS tests)

# Installation
install(TARGETS prasterblasterpio sptw rasterblaster prasterblaster
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib)

if(BundleGDAL)
    # Figure out soname for gdal and proj
    execute_process(COMMAND objdump -p ${GDAL_LIBRARY} COMMAND grep SONAME COMMAND awk "{print \$2}" OUTPUT_VARIABLE GDAL_LIBRARY_FN OUTPUT_STRIP_TRAILING_WHITESPACE)
    execute_process(COMMAND objdump -p ${PROJ_LIBRARY} COMMAND grep SONAME COMMAND awk "{print \$2}" OUTPUT_VARIABLE PROJ_LIBRARY_FN OUTPUT_STRIP_TRAILING_WHITESPACE)

    # Resolve symlinks
    get_filename_component(GDAL_LIBRARY ${GDAL_LIBRARY} REALPATH)
    get_filename_component(PROJ_LIBRARY ${PROJ_LIBRARY} REALPATH)

    # Install the libraries using their soname filenames to make sure loader finds them
    install(FILES ${GDAL_LIBRARY} DESTINATION lib/prasterblaster RENAME ${GDAL_LIBRARY_FN})
    install(FILES ${PROJ_LIBRARY} DESTINATION lib/prasterblaster RENAME ${PROJ_LIBRARY_FN})
endif()
