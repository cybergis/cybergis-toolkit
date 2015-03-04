cmake_minimum_required(VERSION 2.6)
 
# Get user input of IPM library path (e.g. cmake -DIPM_DIR=YOUR/IPM/FILEPATH)
set(IPM_DIR "" CACHE FILEPATH "Path to IPM directory")
 
 
# Project name
project(PABM)
 
# Add source files
add_executable(PABM main.c sim.c mpi.c env.c agent.c sugar.c lb.c migrate.c)
 
# PABM uses two packages, MPI is required and IPM is optional
# MPI is bundled with CMake, so direclty use find_package(). If it is not found, the error is handled by CMake
find_package(MPI REQUIRED)
if(MPI_FOUND)
        set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
        set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})
endif()
 
 
# Include directories and link against MPI
include_directories(${MPI_INCLUDE_PATH})
target_link_libraries(PABM ${MPI_LIBRARIES})

install(TARGETS PABM DESTINATION bin)
 
# Add a sample test (4 nodes mpi test):
enable_testing()
add_test(test1 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${CMAKE_CURRENT_BINARY_DIR}/PABM)
 
 
# Optional IPM is not bundled with CMake, so we need to find it through user
# input directory
 
# Clear variables cache
unset(IPM_LIBRARY CACHE)
unset(IPM_INCLUDE_DIR CACHE)
 
# find_library() and find_path() are key functions to find IPM library
find_library(IPM_LIBRARY ipm PATHS ${IPM_DIR}/lib)
find_path(IPM_INCLUDE_DIR ipm.h PATHS ${IPM_DIR}/include)
 
if(IPM_LIBRARY AND IPM_INCLUDE_DIR)
        # Link IPM flags if it is found
        set(CMAKE_EXE_LINKER_FLAGS "-L${IPM_DIR}/lib -lipm")
else()
        # Warining message if IPM not found
        message("WARNING: IPM is not found, please specify the directory in the command line, or you will not see the IPM outputs.")
endif()
