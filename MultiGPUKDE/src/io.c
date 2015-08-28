/**
 * io.c: methods used to read input points and output GeoTiff
 * Authors: Yizhao Gao <ygao29@illinois.edu>
 * Date: {08/19/2015}
 */

#include <stdio.h>
#include <stdlib.h>
#include <gdal.h>
#include <ogr_srs_api.h>
#include <ogr_api.h>
#include <cpl_conv.h>

int getNumPoints(FILE * file)
{
	float x, y;
	int count = 0;

	rewind(file);

	while(fscanf(file, "%f,%f\n", &x, &y) != EOF)
	{
		count ++;
	}

	return count;
}

//Read points, while reading points, find the blockID it lies in and find the number of points in each block
int * readFile(FILE * file, float * xCor, float * yCor, int * blockID, int nPoints, float xMin, float yMax, float blockSizeDis, int gridX, int gridY)
{
	int * nPointsInB;	
	if(NULL == (nPointsInB = (int *)malloc(sizeof(int) * gridX * gridY)))
	{
		printf("ERROR: Out of Memory at line %d at file %s\n", __LINE__, __FILE__);
		exit(1);
	}
	for(int i = 0; i < gridX * gridY; i++)
	{
		nPointsInB[i] = 0;
	}

	int blockIDY, blockIDX;

	rewind(file);
	for(int i = 0; fscanf(file, "%f,%f\n", xCor + i, yCor + i) != EOF; i++)
	{
		blockIDX = (int)((xCor[i] - xMin) / blockSizeDis);
		blockIDY = (int)((yMax - yCor[i]) / blockSizeDis);

		if(blockIDX < 0 || blockIDX >= gridX || blockIDY < 0 || blockIDY >= gridY)
		{
			blockID[i] = -1;
		}
		else
		{
			blockID[i] = blockIDY * gridX + blockIDX;	
			nPointsInB[blockID[i]] ++;
		}
	}

	return nPointsInB;
}

void writeGeoTiff(char * fileName, float * grid, int nX, int nY, float xMin, float yMax, float cellSize)
{
	GDALAllRegister();
	OGRRegisterAll();

	GDALDatasetH hDstDS;
	GDALDriverH hDriver;
	GDALRasterBandH hBand;
	double adfGeoTransform[6];

	char *papszOptions[] = {"COMPRESS=LZW",NULL};
	const char *pszFormat="GTiff";

	if(NULL == (hDriver = GDALGetDriverByName(pszFormat)))
	{
		printf("ERROR: hDriver is null cannot output using GDAL\n");
		exit(1);
	}

	hDstDS = GDALCreate(hDriver, fileName, nX, nY, 1, GDT_Float32, papszOptions);

	adfGeoTransform[0] = xMin;
	adfGeoTransform[1] = cellSize;
	adfGeoTransform[2] = 0;
	adfGeoTransform[3] = yMax;
	adfGeoTransform[4] = 0;
	adfGeoTransform[5] = -cellSize;

	GDALSetGeoTransform(hDstDS,adfGeoTransform);

	hBand=GDALGetRasterBand(hDstDS,1);
	GDALSetRasterNoDataValue(hBand,-1);
	GDALRasterIO(hBand, GF_Write, 0, 0, nX, nY, grid, nX, nY, GDT_Float32, 0, 0 );
	
	GDALClose(hDstDS);
}

void writeGeoTiff(char * fileName, int * grid, int nX, int nY, float xMin, float yMax, float cellSize)
{
	GDALAllRegister();
	OGRRegisterAll();

	GDALDatasetH hDstDS;
	GDALDriverH hDriver;
	GDALRasterBandH hBand;
	double adfGeoTransform[6];

	char *papszOptions[] = {"COMPRESS=LZW",NULL};
	const char *pszFormat="GTiff";

	if(NULL == (hDriver = GDALGetDriverByName(pszFormat)))
	{
		printf("ERROR: hDriver is null cannot output using GDAL\n");
		exit(1);
	}

	hDstDS = GDALCreate(hDriver, fileName, nX, nY, 1, GDT_Int32, papszOptions);

	adfGeoTransform[0] = xMin;
	adfGeoTransform[1] = cellSize;
	adfGeoTransform[2] = 0;
	adfGeoTransform[3] = yMax;
	adfGeoTransform[4] = 0;
	adfGeoTransform[5] = -cellSize;

	GDALSetGeoTransform(hDstDS,adfGeoTransform);

	hBand=GDALGetRasterBand(hDstDS,1);
	GDALSetRasterNoDataValue(hBand,-1);
	GDALRasterIO(hBand, GF_Write, 0, 0, nX, nY, grid, nX, nY, GDT_Int32, 0, 0 );
	
	GDALClose(hDstDS);
}
