/**
 * io.h: header file for io.c
 * Authors: Yizhao Gao <ygao29@illinois.edu>
 * Date: {08/19/2015}
 */

#ifndef IOH
#define IOH


/**
 * NAME:	getNumPoints
 * DESCRIPTION:	get the number of points in an input file
 * PARAMETERS:
 * 	FILE * file:	A file pointer to the input file
 * RETURN:
 * 	TYPE:	int
 * 	VALUE:	the number of points in the input file
 */
int getNumPoints(FILE * file);


/**
 * NAME:	readFile
 * DESCRIPTION:	read all the x, y cooridinates of each point, and find the block it belongs to
 * PARAMETERS:
 * 	FILE * file:		A file pointer to the input file
 * 	float * xCor:		An array of X coordinates of all the points
 * 	float * yCor:		An array of Y coordinates of all the points
 *	int * blockID:		An array of the block ID each point belongs to
 *	float * xMin:		The western bound of data
 *	float * yMax:		The northern bound of data
 *	float blockSizeDis:	Block size in terms of distance
 * 	int gridX:		Number of blocks along X dimension
 * 	int gridY:		Number of blocks along Y dimension
 * RETURN:
 * 	TYPE:	int *
 * 	VALUE:	An array storing the number of points in each block 
 */
int * readFile(FILE * file, float * xCor, float * yCor, int * blockID, int nPoints, float xMin, float yMax, float blockSizeDis, int gridX, int gridY);


/**
 * NAME:	writeGeoTiff
 * DESCRIPTION:	write a raster to a GeoTiff
 * PARAMETERS:
 * 	char * fileName:	Output file Name
 *	float * grid:		An array of raster values at each cell
 *	int nX:			Number of cells in x dimension
 *	int nY:			Number of cells in y dimension
 *	float xMin:		Western bound of raster
 *	float yMax:		Northern bound of raster
 *	flaot cellSize:		Cell size 
 * RETURN:	void
 */
void writeGeoTiff(char * fileName, float * grid, int nX, int nY, float xMin, float yMax, float cellSize);


/**
 * NAME:	writeGeoTiff
 * DESCRIPTION:	write a raster to a GeoTiff
 * PARAMETERS:
 * 	char * fileName:	Output file Name
 *	int * grid:		An array of raster values at each cell
 *	int nX:			Number of cells in x dimension
 *	int nY:			Number of cells in y dimension
 *	float xMin:		Western bound of raster
 *	float yMax:		Northern bound of raster
 *	flaot cellSize:		Cell size 
 * RETURN:	void
 */
void writeGeoTiff(char * fileName, int * grid, int nX, int nY, float xMin, float yMax, float cellSize);


#endif

