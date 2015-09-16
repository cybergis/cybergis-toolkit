/**
 * kde.h: header file for kde.u
 * Authors: Yizhao Gao <ygao29@illinois.edu>
 * Date: {08/19/2015}
 */

#ifndef KDEH
#define KDEH


/**
 * NAME:	performKDE
 * DESCRIPTION:	calculating KDE
 * PARAMETERS:
 * 	float * xCor:		An array of X coordinates of all the points
 * 	float * yCor:		An array of Y coordinates of all the points
 * 	int nPoints:		Total number of points
 * 	int * pointIndex:	An array storing the ending index of points in each block
 * 	int gridX:		Number of blocks along X dimension
 * 	int gridY:		Number of blocks along Y dimension
 * 	int blockBandwidth:	Bandwidth in terms of the number of blocks
 *	int cellX:		Number of cells in x dimension
 *	int cellY:		Number of cells in y dimension
 *	float xMin:		Western bound of raster
 *	float yMax:		Northern bound of raster
 *	float cellSize:		Cell size
 *	float bandwidth:	Bandwidth for KDE calculation 
 * RETURN:
 * 	TYPE:	float *
 * 	VALUE:	An array of density value at each cell
 */
float * performKDE(float * xCor, float * yCor, int nPoints, int * pointIndex, int gridX, int gridY, int blockBandwidth, int cellX, int cellY, float xMin, float yMax, float cellSize, float bandwidth);


#endif
