include_directories(${CMAKE_CURRENT_SOURCE_DIR} 
                    "${CMAKE_CURRENT_SOURCE_DIR}/utils")

add_subdirectory(utils)

add_executable(mapalg mapalg.cc)
target_link_libraries(mapalg utils)
target_link_libraries(mapalg ${MPI_CXX_LIBRARIES})
target_link_libraries(mapalg ${GDAL_LIBRARIES})

add_executable(mapalg-kernel mapalg-kconvolve.cc)
target_link_libraries(mapalg-kernel utils)
target_link_libraries(mapalg-kernel ${MPI_CXX_LIBRARIES})
target_link_libraries(mapalg-kernel ${GDAL_LIBRARIES})

