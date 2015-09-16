#ifndef DATA_CC
#define DATA_CC
/** data.cc: data I/O functions using GDAL C API 
 * Author: Yan Y. Liu <yanliu@illinois.edu>
 * Date: 06/18/2014
 */

#include "data.h"
#include "util.h"

// open a raster file
GDALDatasetH raster_open(char * pszFilename, double *georef, char *prj,
                         double *nodata, int *x, int *y, int band)
{
	GDALDatasetH  hDataset;

	GDALAllRegister();

	// open the dataset
	hDataset = GDALOpen( pszFilename, GA_ReadOnly );
	if( hDataset == NULL ) {
		fprintf(stderr, "Error open raster %s\n", pszFilename);
		exit(1);
	}

	// get projection info
	const char * prjInfo = GDALGetProjectionRef( hDataset );
	strcpy(prj, prjInfo);

	// get geo info
	if (GDALGetGeoTransform( hDataset, georef) == CE_Failure) {
		fprintf(stderr, "Info: fetch geo transformation info on raster %s\n",
                pszFilename);
	}

	*x = GDALGetRasterXSize( hDataset );
	*y = GDALGetRasterYSize( hDataset );
    int maxBand = GDALGetRasterCount(hDataset);
    if (band>maxBand || band<=0) {
        fprintf(stderr, "Error: invalid band index, number of bands is %d\n",
                maxBand);
        exit(-1);
    }

	GDALRasterBandH hBand;
	hBand = GDALGetRasterBand( hDataset, band );
	int status;
	*nodata = GDALGetRasterNoDataValue(hBand, &status);

	return hDataset;
}


// print raster info
void raster_info(GDALDatasetH hDataset, double *georef, char *prj, 
                 double nodata, int x, int y)
{
	GDALDriverH   hDriver;

	hDriver = GDALGetDatasetDriver( hDataset );
	fprintf(stdout, "Data driver: %s%s\n", GDALGetDriverShortName( hDriver ),
            GDALGetDriverLongName( hDriver ) );
	fprintf(stdout, "Size: %dx%dx%d\n", x, y, GDALGetRasterCount( hDataset ) );
	fprintf(stdout, "Projection: \n%s\n", prj);
	fprintf(stdout, "Origin = (%.6f,%.6f)\n", georef[0], georef[3] );
	fprintf(stdout, "Pixel Size = (%.6f,%.6f)\n", georef[1], georef[5] );
	fprintf(stdout, "Band 1 NoData = %.6f\n", nodata );
}


// read data block into memory
int raster_read(GDALDatasetH hDataset, float *buf, int offsetx, int offsety,
                int sizex, int sizey)
{
	GDALRasterBandH hBand;
	hBand = GDALGetRasterBand( hDataset, 1 );
	float *pafScanline;
	int nSize = sizex * sizey;
	pafScanline = (float *) CPLMalloc(sizeof(float)*nSize);
	GDALRasterIO(hBand, GF_Read, offsetx, offsety, sizex, sizey, 
                 pafScanline, sizex, sizey, GDT_Float32, 0, 0 );
	int i;
	for (i=0; i<nSize; i++) {
		buf[i] = pafScanline[i];
	}
	CPLFree(pafScanline);
	return 1;
}

// close a raster dataset
void raster_close(GDALDatasetH hDataset)
{
	GDALClose(hDataset);
}
// create output raster
GDALDatasetH raster_create(char *fn, int x, int y, double *georef, 
                           char *prj, double nodata)
{
	GDALDriverH gtiff_driver = NULL;
	GDALDatasetH ds = NULL;
	char **options = NULL;
	const char *format = "GTiff";

	GDALAllRegister();
	// get data driver
	gtiff_driver = GDALGetDriverByName( format );
	if (gtiff_driver == NULL) {
		fprintf(stderr, "Error get raster data driver %s\n", format);
		exit(1);
	}

	options = CSLSetNameValue(options, "BIGTIFF", "YES");
	options = CSLSetNameValue(options, "INTERLEAVE", "PIXEL");
	options = CSLSetNameValue(options, "COMPRESS", "NONE");
	// create raster
	ds = GDALCreate( gtiff_driver, fn, x, y, 1, GDT_Float32, options);

	// Clean up options
	CSLDestroy(options);

	// set geo info
	GDALSetGeoTransform( ds, georef );
	// set projection
	if (GDALSetProjection( ds, prj) == CE_Failure) {
		fprintf(stderr, "Error set output raster projection %s\n", fn);
		exit(1);
	}
	GDALRasterBandH hBand;
	hBand = GDALGetRasterBand( ds, 1 );
	GDALSetRasterNoDataValue(hBand, nodata);
	return ds;
}

// write data block from memory
int raster_write(GDALDatasetH hDstDS, float *buf, int offsetx, int offsety,
                 int sizex, int sizey)
{
	GDALRasterBandH hBand;
	hBand = GDALGetRasterBand( hDstDS, 1 );
	GDALRasterIO(hBand, GF_Write, offsetx, offsety, sizex, sizey, buf,
                 sizex, sizey, GDT_Float32, 0, 0 );
	return 1;
}
#endif
