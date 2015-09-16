# Map-Algebra

## Summary
It is a parallel geospatial programming models training package. The
package contains parallel map algebra code using CUDA, MPI, and Parallel
I/O.
	

Authors: 
	Yan Y. Liu <yanliu@illinois.edu>
	Garrett Nickel <gmnicke2@illinois.edu>
	Yiming Huang <yhuan125@illinois.edu>
	Mingze Gao <mgao16@illinois.edu>
	Sunwoo Kim <kim392@illinois.edu>
    Jeff Terstriep <jefft@illinois.edu>
	
## Requirements
    MPI
    CMake 2.8
    gdal 1.11.2

## Content
- src/
	(parallel map algebra using MPI)
	CMakeList.txt: script to build Makefile using cmake command
	mapalg.cc: main program
	mapalg-kconvolve.cc: map algebra - kernel convolution
	util.{h,cc}: utility functions
	data.{h,cc}: GDAL-based raster I/O functions with data 
		decomposition strategy (block-wise)

	jenkins_build.sh: Provided example Jenkins build script to SSH and use

- test/
	decomp.cc: illustration of how to do block-wise data decomposition 
	rastercp.cc: illustration of using GDAL raster I/O for raster copy
