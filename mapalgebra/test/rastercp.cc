/** rastercp.cc: Illustration of GDAL C API-based raster copy
 * Author: Yan Y. Liu <yanliu@illinois.edu>
 * Date: 06/18/2014
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "util.h"
#include "data.h"
#include "gdal.h"
#include "cpl_conv.h"
#include "cpl_string.h"

// test read and write a raster file in data blocks - raster copy
// g++ -DDEBUG -I. -I../code/ -I$GDAL_HOME/include -o rastercp rastercp.cc ../code/util.cc ../code/data.cc -L$GDAL_HOME/lib -lgdal -lm
int main(int argc, char **argv) {
	char * fn = argv[1];
	char * ofn = argv[2];
	int np = 8; // np, also num of data blocks
	int i, j;
	
	double georef[6]; // georef data structure for a raster
	char prj[2048]; // store projection wkt
	double nodata;
	int x, y; // size of raster on x and y dim
	GDALDatasetH rin;

	// get input raster info
	rin = raster_open(fn, georef, prj, &nodata, &x, &y);
	raster_info(rin, georef, prj, nodata, x, y);
	// determine block sizes
	int *offsetx, *offsety, *sizex, *sizey;
	offsetx = (int *) malloc(sizeof(int) * np);
	memset(offsetx, 0, sizeof(int) * np);
	offsety = (int *) malloc(sizeof(int) * np);
	memset(offsety, 0, sizeof(int) * np);
	sizex = (int *) malloc(sizeof(int) * np);
	memset(sizex, 0, sizeof(int) * np);
	sizey = (int *) malloc(sizeof(int) * np);
	memset(sizey, 0, sizeof(int) * np);
	for (i=0; i<np; i++) {
		get_block(i, np, x, y, &offsetx[i], &offsety[i], &sizex[i], &sizey[i]);
	}
	// find max sizex and sizey
	int maxx=0, maxy=0;
	for (i=0; i<np; i++) {
		if (maxx < sizex[i]) maxx = sizex[i];
		if (maxy < sizey[i]) maxy = sizey[i];
	}
	// allocate data block memory
	float *raster = (float *)malloc(sizeof(float) * maxx * maxy * np);
	memset(raster, 0, sizeof(float) * maxx * maxy * np);
	// read blocks from input raster
	float *block;
	for (i=0; i<np; i++) {
		block = raster + i * (maxx * maxy);
		raster_read(rin, block, offsetx[i], offsety[i], sizex[i], sizey[i]);
	}
	// close raster
	raster_close(rin);
	// processing ...

	// create output raster
	GDALDatasetH rout;
	rout = raster_create(ofn, x, y, georef, prj, nodata);
	for (i=0; i<np; i++) {
		block = raster + i * (maxx * maxy);
		raster_write(rout, block, offsetx[i], offsety[i], sizex[i], sizey[i]);
	}
	raster_close(rout);
	
	// free memory
	free(offsetx);
	free(offsety);
	free(sizex);
	free(sizey);
	free(raster);
}
