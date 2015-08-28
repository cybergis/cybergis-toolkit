/**
 * decompose.h: header file of decompose.c
 * Authors: Yizhao Gao <ygao29@illinois.edu>
 * Date: {08/19/2015}
 */

#ifndef DCH
#define DCH

/**NAME:	spatialDomainDecom
 * DESCRIPTION:	Estimate the computational intensity surface, and determine the extent of each sub-region to be processed by teach node
 * PARAMETERS:
 * 	int nXRegion:		Number of splits along X dimension
 * 	int nYRegion:		Number of splits along Y dimension
 * 	int gridX:		Number of blocks along X dimension
 * 	int gridY:		Number of blocks along Y dimension
 * 	int blockBandwidth:	Bandwidth in terms of the number of blocks
 * 	int * nPointsInB:	An array storing the number of points in each block
 * RETURN:
 * 	TYPE:	int *
 * 	VALUE:	The extents of all the decomposed sub-regions. The first four elements are the xMin, xMax, yMin and yMax (in term of blocks) for the first sub-region. It is followed by those of the second, third ... (if any) region.
 */
int * spatialDomainDecom(int nXRegion, int nYRegion, int gridX, int gridY, int blockBandwidth, int * nPointsInB);


/**NAME:	indexPoints
 * DESCRIPTION:	Index and reorder points based on the block they fall into
 * PARAMETERS:
 * 	float * &allX:		An array of X coordinates of all the points. The pointer will be pointing to reordered points after execution.
 * 	float * &allY:		An array of Y coordinates of all the points. The pointer will be pointing to reordered points after execution.
 * 	int * allBlockID: 	An array of the blocks each point falls into.
 * 	int &allPoints:		The total number of points, will be changed to the number of points that are indexed after execution.
 * 	int dataGridX:		Number of data blocks along X dimension
 * 	int dataGridY:		Number of data blocks along Y dimension
 * 	int * allPointsInB:	An array storing the number of points in each block
 * RETURN:
 * 	TYPE:	int *
 * 	VALUE:	An array storing the ending index of points in each block
 */
int * indexPoints(float * &allX, float * &allY, int * allBlockID, int &allPoints, int dataGridX, int dataGridY, int * allPointsInB);


/**
 * NAME:	getNPoints
 * DESCRIPTION:	get the number of points (including neighboring areas) within a spatial region
 * PARAMETERS:
 * 	int * allPointsInB:	An array storing the number of points in each block
 * 	int dataGridX:		Number of data blocks along X dimension
 * 	int dataGridY:		Number of data blocks along Y dimension
 * 	int blockBandwidth:	Bandwidth in term of the number of blocks
 * 	int * spatialRegion:	The xMin, xMax, yMin and yMax (in terms of block) for the spatial region
 * RETURN:
 * 	TYPE:	int
 * 	VALUE:	The number of points within a spatial region itself and its neighboring areas
 */
int getNPoints(int * allPointsInB, int dataGridX, int dataGridY, int blockBandwidth, int * spatialRegion);


/**
 * NAME:	getDataInRegion
 * DESCRIPTION:	get the points within a spatial region and index these points
 * PARAMETERS:
 * 	int * allPointsInB:	An array storing the number of points in each block
 * 	int * allPointsIndex:	An array storing the ending index of points in each block
 * 	float * allX:		An array of X coordinates of all the points
 * 	float * allY:		An array of Y coordinates of all the points
 * 	int dataGridX:		Number of data blocks along X dimension
 * 	int dataGridY:		Number of data blocks along Y dimension
 * 	int blockBandwidth:	Bandwidth in term of the number of blocks
 * 	int * spatialRegion:	The xMin, xMax, yMin and yMax (in terms of block) for the spatial region
 * 	float * &regionX:	An array used to store and send back X coordinates of points in the region
 * 	float * &regionY:	An array used to store and send back Y coordinates of points in the region
 * 	float * &regionIndex:	An array used to store and send back the spatial index of points in the region
 * RETURN:
 * 	TYPE:	int
 * 	VALUE:	The number of points within a spatial region itself and its neighboring areas
 */
int getDataInRegion(int * allPointsInB, int * allPointsIndex, float * allX, float * allY, int dataGridX, int dataGridY, int blockBandwidth, int * spatialRegion, float * &regionX, float * &regionY, int * &regionIndex);


/**
 * NAME:	mergeInto
 * DESCRIPTION:	merge the density output in one region into the final output
 * PARAMETERS:
 * 	float * fullFile:	An array of the density values at each map cell for the whole study area
 * 	int xFull:		Number of colums in the output raster
 * 	float * partFile:	An array of the density values at each map cell for the region to be merged into the whole study area
 * 	int partXmin:		The minium cell ID in x dimension of the region
 * 	int partYmin:		The minium cell ID in y dimension of the region
 *	int xPart:		The number of cells in x dimension of the region
 *	int yPart:		The number of cells in y dimension of the region
 * RETURN:	void
 */
void mergeInto(float * fullFile, int xFull, float * partFile, int partXMin, int partYMin, int xPart, int yPart);

#endif
