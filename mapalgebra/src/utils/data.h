#ifndef DATA_H
#define DATA_H
/** data.h: data I/O functions
 * Author: Yan Y. Liu <yanliu@illinois.edu>
 * Date: 06/18/2014
 */
#include "gdal.h"
#include "cpl_conv.h"
#include "cpl_string.h"

GDALDatasetH raster_open(char * pszFilename, double *georef, char *prj, double *nodata, int *x, int *y, int band=1);
void raster_info(GDALDatasetH hDataset, double *georef, char *prj, double nodata, int x, int y);
int raster_read(GDALDatasetH hDataset, float *buf, int offsetx, int offsety, int sizex, int sizey);
void raster_close(GDALDatasetH hDataset);
GDALDatasetH raster_create(char *fn, int x, int y, double *georef, char *prj, double nodata);
int raster_write(GDALDatasetH hDstDS, float *buf, int offsetx, int offsety, int sizex, int sizey);

#endif
