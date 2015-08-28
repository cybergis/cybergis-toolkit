/**
 * decompose.c: methods related to spatial domain decomposition
 * Authors: Yizhao Gao <ygao29@illinois.edu>
 * Date: {08/19/2015}
 */

#include <stdio.h>
#include <stdlib.h>
#include "io.h"

int * spatialDomainDecom(int nXRegion, int nYRegion, int gridX, int gridY, int blockBandwidth, int * nPointsInB)
{
	int * allRegions; //xMin, xMax, yMin, yMax
	if(NULL == (allRegions = (int *)malloc(sizeof(int) * nXRegion * nYRegion * 4)))
	{
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}

	int * compInten;
	if(NULL == (compInten = (int *)malloc(sizeof(int) * gridX * gridY)))
	{
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}

	long intensity;
	long totalInten = 0;
	long * totalIntenXRegion;
	if(NULL == (totalIntenXRegion = (long *)malloc(sizeof(long) * nXRegion)))
	{
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}

	int dataGridX = gridX + 2 * blockBandwidth;
	for(int blockIDY = 0; blockIDY < gridY; blockIDY ++)
	{
		for(int blockIDX = 0; blockIDX < gridX; blockIDX ++)
		{
			intensity = 0;
			for(int i = 0; i < blockBandwidth * 2 + 1; i ++)
			{
				for(int j = 0; j < blockBandwidth * 2 + 1; j ++)
				{
					intensity += nPointsInB[(blockIDY + i) * dataGridX + blockIDX + j];
				}
			}
			compInten[blockIDY * gridX + blockIDX] = intensity;
			totalInten += intensity;
		}
	}

	
//	writeGeoTiff("/gpfs_scratch/tsccsj/KDE/compInten.tif", compInten, gridX, gridY, -2380000, 1420000, 32000);

	int * splitX;
	int * splitY;
	int splitIndex = 1;
	if(NULL == (splitX = (int *)malloc(sizeof(int) * (nXRegion + 1))))
	{
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}
	splitX[0] = 0;
	splitX[nXRegion] = gridX;
	if(NULL == (splitY = (int *)malloc(sizeof(int) * (nYRegion + 1))))
	{
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}
	splitY[0] = 0;
	splitY[nYRegion] = gridY;

	if(nXRegion > 1)
	{
		intensity = 0;
		int colIntensity;
		for(int blockIDX = 0; blockIDX < gridX; blockIDX ++)
		{
			colIntensity = 0;
			for(int blockIDY = 0; blockIDY < gridY; blockIDY ++)
			{
				colIntensity += compInten[blockIDY * gridX + blockIDX];
			}

			if((intensity + colIntensity) * nXRegion > totalInten * splitIndex)
			{
				splitX[splitIndex] = blockIDX;
				totalIntenXRegion[splitIndex - 1] = intensity;
				splitIndex ++;
				if(splitIndex >= nXRegion)
				{
					break;
				}
			}

			intensity += colIntensity;
		}

		if(splitIndex != nXRegion)
		{
			printf("ERROR: Can't find correct splitting points\nsplitIndex=%d\tnXRegion=%d\n", splitIndex, nXRegion);
			exit(1);
		}

		totalIntenXRegion[nXRegion - 1] = totalInten;
		for(int i = nXRegion - 1; i > 0; i--)
		{
			totalIntenXRegion[i] = totalIntenXRegion[i] - totalIntenXRegion[i - 1];
		}
	}
	else
	{
		totalIntenXRegion[0] = totalInten;
	}

	int regionXMin;
	int regionXMax;
	for(int regionX = 0; regionX < nXRegion; regionX++)
	{
		regionXMin = splitX[regionX];
		regionXMax = splitX[regionX + 1];

//		printf("%ld\t%ld\n", totalIntenXRegion[regionX], totalInten);

		if(nYRegion > 1)
		{
			splitIndex = 1;	
			intensity = 0;
			for(int blockIDY = 0; blockIDY < gridY; blockIDY ++)
			{
				for(int blockIDX = regionXMin; blockIDX < regionXMax; blockIDX ++)
				{
					intensity += compInten[blockIDY * gridX + blockIDX];
				}
				if(intensity * nYRegion > totalIntenXRegion[regionX] * splitIndex)
				{
					splitY[splitIndex] = blockIDY;
					splitIndex ++;
					if(splitIndex >= nYRegion)
					{
						break;
					}
				}
			}

			if(splitIndex != nYRegion)
			{
				printf("ERROR: Can't find correct splitting points\nsplitIndex=%d\tnYRegion=%d\n", splitIndex, nYRegion);
				exit(1);
			}
		}

		for(int regionY = 0; regionY < nYRegion; regionY ++)
		{
			allRegions[(regionY * nXRegion + regionX) * 4] = regionXMin;
			allRegions[(regionY * nXRegion + regionX) * 4 + 1] = regionXMax;
			allRegions[(regionY * nXRegion + regionX) * 4 + 2] = splitY[regionY];
			allRegions[(regionY * nXRegion + regionX) * 4 + 3] = splitY[regionY + 1];
		}

	}
	free(splitX);
	free(splitY);
	free(compInten);
	free(totalIntenXRegion);

	return allRegions;
}

int * indexPoints(float * &allX, float * &allY, int * allBlockID, int &allPoints, int dataGridX, int dataGridY, int * allPointsInB)
{
	int totalPoints = 0;
	int tempInt;

	int * indexTable;
	if(NULL == (indexTable = (int *)malloc(sizeof(int) * dataGridX * dataGridY)))
	{
		printf("ERROR: Out of memory at line %d at file %s\n", __LINE__, __FILE__);
		exit(1);
	}


	for(int i = 0; i < dataGridX * dataGridY; i++)
	{
		indexTable[i] = totalPoints;
		totalPoints += allPointsInB[i];
	}

	float * newX;
	float * newY;

	if(NULL == (newX = (float *)malloc(sizeof(float) * totalPoints)))
	{
		printf("ERROR: Out of memory at line %d at file %s\n", __LINE__, __FILE__);
		exit(1);
	}

	if(NULL == (newY = (float *)malloc(sizeof(float) * totalPoints)))
	{
		printf("ERROR: Out of memory at line %d at file %s\n", __LINE__, __FILE__);
		exit(1);
	}

	for(int i = 0; i < allPoints; i ++)
	{
		tempInt = allBlockID[i];
		if(tempInt > -1)
		{
			newX[indexTable[tempInt]] = allX[i];
			newY[indexTable[tempInt]] = allY[i];
			indexTable[tempInt] ++;
		}
	}

	free(allX);
	free(allY);

	allX = newX;
	allY = newY;
	allPoints = totalPoints;


	return indexTable;
}

int getNPoints(int * allPointsInB, int dataGridX, int dataGridY, int blockBandwidth, int * spatialRegion)
{
	int nPoints = 0;
	for(int blockIDY = spatialRegion[2]; blockIDY < spatialRegion[3] + 2 * blockBandwidth; blockIDY ++)
	{
		for(int blockIDX = spatialRegion[0]; blockIDX < spatialRegion[1] + 2 * blockBandwidth; blockIDX ++)
		{
			nPoints += allPointsInB[blockIDY * dataGridX + blockIDX];
		}
	}

	return nPoints;
}

int getDataInRegion(int * allPointsInB, int * allPointsIndex, float * allX, float * allY, int dataGridX, int dataGridY, int blockBandwidth, int * spatialRegion, float * &regionX, float * &regionY, int * &regionIndex)
{
	int regionGridXMin = spatialRegion[0];
	int regionGridXMax = spatialRegion[1] + 2 * blockBandwidth;
	int regionGridYMin = spatialRegion[2];
	int regionGridYMax = spatialRegion[3] + 2 * blockBandwidth;

	int regionDataGridX = regionGridXMax - regionGridXMin;
	int regionDataGridY = regionGridYMax - regionGridYMin;
	
	if(NULL == (regionIndex = (int *)malloc(sizeof(int) * regionDataGridX * regionDataGridY)))
	{
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}

	int regionNPoints = 0;
	for(int blockIDY = 0; blockIDY < regionDataGridY; blockIDY ++)
	{
		for(int blockIDX = 0; blockIDX < regionDataGridX; blockIDX ++)
		{
			regionNPoints += allPointsInB[(blockIDY + regionGridYMin) * dataGridX + (blockIDX + regionGridXMin)]; 
			regionIndex[blockIDY * regionDataGridX + blockIDX] = regionNPoints;
		}
	}

	
	if(NULL == (regionX = (float *)malloc(sizeof(float) * regionNPoints)))
	{
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}
	if(NULL == (regionY = (float *)malloc(sizeof(float) * regionNPoints)))
	{
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}

	int idAll, idAllEnd;
	int idReg = 0;

	for(int blockIDY = regionGridYMin; blockIDY < regionGridYMax; blockIDY ++)
	{
		if(blockIDY == 0 && regionGridXMin == 0)
		{
			idAll = 0;
		}
		else
		{
			idAll = allPointsIndex[blockIDY * dataGridX + regionGridXMin - 1];
		}

		idAllEnd = allPointsIndex[blockIDY * dataGridX + regionGridXMax - 1];
		for(; idAll < idAllEnd; idAll ++)
		{
			regionX[idReg] = allX[idAll];
			regionY[idReg] = allY[idAll];
			idReg ++;
		}
	}

	return regionNPoints;
}

void mergeInto(float * fullFile, int xFull, float * partFile, int partXMin, int partYMin, int xPart, int yPart)
{
	for(int cellIDY = 0; cellIDY < yPart; cellIDY ++)
	{
		for(int cellIDX = 0; cellIDX < xPart; cellIDX ++)
		{
			fullFile[(cellIDY + partYMin) * xFull + cellIDX + partXMin] = partFile[cellIDY * xPart + cellIDX];
		}
	}
}
