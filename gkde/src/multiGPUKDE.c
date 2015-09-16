/**
 * multiGPUKDE.c: main method for GKDE
 * Authors: Yizhao Gao <ygao29@illinois.edu>
 * Date: {08/19/2015}
 */

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include <sys/time.h>
#include "io.h"
#include "decompose.h"
#include "kde.h"


#define BLOCKSIZE 16
int main(int argc, char ** argv)
{
	int rank, size;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	if(argc != 11)
	{
		if(rank == 0)
		{
			printf("ERROR: Incorrect input parameters, to run the code:\nMultiGPUKDE InputFile OutputFile xMin xMax yMin yMax cellSize bandwidth numberOfRegionsInX numberOfRegionsInY\n");
		}
		MPI_Finalize();
		exit(1);
	}

	float xMin = atof(argv[3]); 
	float xMax = atof(argv[4]);
	float yMin = atof(argv[5]); 
	float yMax = atof(argv[6]); 
	float cellSize = atof(argv[7]); 
	float bandwidth = atof(argv[8]);
	int nXRegion = atoi(argv[9]);
	int nYRegion = atoi(argv[10]);

//	printf("Size: %d\n", size);

	if(nXRegion * nYRegion > size)
	{
		if(rank == 0)
		{
			printf("ERROR: The number of regions %d(x) * %d(y) should be no bigger than the number of processors %d\n", nXRegion, nYRegion, size);
		}
		MPI_Finalize();
		exit(1);
	}
		
	
	FILE * file;

	int allPoints = 0;
	float * allX;
	float * allY;
	int * allBlockID;

	int * allPointsInB;
	int * allPointsIndex;
	
	int * allRegions;


	int nXCell = ceil((xMax - xMin)/cellSize);
	int nYCell = ceil((yMax - yMin)/cellSize);
	xMax = xMin + nXCell * cellSize;
	yMax = yMin + nYCell * cellSize;

	int gridX = ceil((float) nXCell / BLOCKSIZE);
	int gridY = ceil((float) nYCell / BLOCKSIZE);

	float blockSizeDis = BLOCKSIZE * cellSize;
	int blockBandwidth = ceil(bandwidth / blockSizeDis);

	int dataGridX = gridX + 2 * blockBandwidth;
	int dataGridY = gridY + 2 * blockBandwidth;


	int myPoints;
	float * myX;
	float * myY;
	int * myIndexTable;
	float * myDensity;
	int myRegion[4];
	int myCellXMin, myCellXMax, myCellYMin, myCellYMax;
	float myXMin, myXMax, myYMin, myYMax;
	int myXCell, myYCell, myGridX, myGridY, myDataGridX, myDataGridY;


	
	struct timeval tBegin;
	struct timeval tEnd;

	if(rank == 0)
	{

		printf("**********************************************\nMulti-GPU Kernel Density Estimation\nCyberInfrastructure and Geospatial Information\n(CIGI) Labortory(http://www.cigi.illinois.edu)\nDeveloped by Yizhao Gao (ygao29@illinois.edu)\n**********************************************\n");
		printf("Running on %d (%d * %d) nodes\n", nXRegion * nYRegion, nXRegion, nYRegion);
		printf("Study area:\n\txMin: %f, xMax: %f, yMin: %f, yMax %f\n", xMin, xMax, yMin, yMax);
		printf("\tOutput raster size: %d(row) * %d(column)\n", nYCell, nXCell);


		
		gettimeofday(&tBegin, NULL);
		if(NULL == (file = fopen(argv[1], "r")))
		{
			printf("ERROR: Can't open the input file.\n");
			//MPI_Abort(MPI_COMM_WORLD, MPI_ERR_FILE);
			exit(1);
		}

		allPoints = getNumPoints(file);

		printf("Number of input points: %d\n", allPoints);
		
		if(NULL == (allX = (float *)malloc(sizeof(float) * allPoints)))
		{
			printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
			exit(1);
		}
		if(NULL == (allY = (float *)malloc(sizeof(float) * allPoints)))
		{
			printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
			exit(1);
		}
		if(NULL == (allBlockID = (int *)malloc(sizeof(int) * allPoints)))
		{
			printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
			exit(1);
		}

		allPointsInB = readFile(file, allX, allY, allBlockID, allPoints, xMin - blockSizeDis, yMax + blockSizeDis, blockSizeDis, dataGridX, dataGridY);

		fclose(file);

		gettimeofday(&tEnd, NULL);
		printf("Time - Input:\t%lfms\n", ((&tEnd)->tv_sec - (&tBegin)->tv_sec) * 1000 + (double)((&tEnd)->tv_usec - (&tBegin)->tv_usec) / 1000);

		tBegin = tEnd;
	
		allPointsIndex = indexPoints(allX, allY, allBlockID, allPoints, dataGridX, dataGridY, allPointsInB);
		free(allBlockID);

		gettimeofday(&tEnd, NULL);
		printf("Time - Spatial index:\t%lfms\n", ((&tEnd)->tv_sec - (&tBegin)->tv_sec) * 1000 + (double)((&tEnd)->tv_usec - (&tBegin)->tv_usec) / 1000);
		
		tBegin = tEnd;

		allRegions = spatialDomainDecom(nXRegion, nYRegion, gridX, gridY, blockBandwidth, allPointsInB);
/*
		for(int i = 0; i < nXRegion * nYRegion; i++)
		{
			printf("Region %d:\t%d\t%d\t%d\t%d\n", i, allRegions[4 * i],  allRegions[4 * i + 1], allRegions[4 * i + 2], allRegions[4 * i + 3]);
		}
*/

		for(int i = 1; i < nXRegion * nYRegion; i++)
		{
//			myPoints = getNPoints(allPointsInB, dataGridX, dataGridY, blockBandwidth, allRegions + 4 * i);

			myDataGridX = allRegions[4 * i + 1] - allRegions[4 * i] + 2 * blockBandwidth;
			myDataGridY = allRegions[4 * i + 3] - allRegions[4 * i + 2] + 2 * blockBandwidth;

			myPoints = getDataInRegion(allPointsInB, allPointsIndex, allX, allY, dataGridX, dataGridY, blockBandwidth, allRegions + 4 * i, myX, myY, myIndexTable);	

			MPI_Send(&myPoints, 1, MPI_INT, i, 0, MPI_COMM_WORLD);

			MPI_Send(allRegions + i * 4, 4, MPI_INT, i, 1, MPI_COMM_WORLD);

			MPI_Send(myX, myPoints, MPI_FLOAT, i, 2, MPI_COMM_WORLD);
			MPI_Send(myY, myPoints, MPI_FLOAT, i, 3, MPI_COMM_WORLD);
			MPI_Send(myIndexTable, myDataGridX * myDataGridY, MPI_INT, i, 4, MPI_COMM_WORLD);

			free(myX);
			free(myY);
			free(myIndexTable);
		}

		myRegion[0] = allRegions[0];	
		myRegion[1] = allRegions[1];	
		myRegion[2] = allRegions[2];	
		myRegion[3] = allRegions[3];
		myDataGridX = myRegion[1] - myRegion[0] + 2 * blockBandwidth;
		myDataGridY = myRegion[3] - myRegion[2] + 2 * blockBandwidth;

		myPoints = getDataInRegion(allPointsInB, allPointsIndex, allX, allY, dataGridX, dataGridY, blockBandwidth, allRegions, myX, myY, myIndexTable);	
		
		//writeGeoTiff("/gpfs_scratch/tsccsj/KDE/allInB.tif", allPointsInB, dataGridX, dataGridY, xMin - blockSizeDis, yMax + blockSizeDis, blockSizeDis);
		gettimeofday(&tEnd, NULL);
		printf("Time - spliting data and computation to nodes:\t%lfms\n", ((&tEnd)->tv_sec - (&tBegin)->tv_sec) * 1000 + (double)((&tEnd)->tv_usec - (&tBegin)->tv_usec) / 1000);

		free(allPointsInB);
		free(allPointsIndex);
		free(allX);
		free(allY);
	}

	else if(rank < nXRegion * nYRegion)
	{
		MPI_Recv(&myPoints, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
/*
		if(myPoints < 0)
		{
			MPI_Finalize();
			exit(1);
		}
*/
		MPI_Recv(&myRegion, 4, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

		myDataGridX = myRegion[1] - myRegion[0] + 2 * blockBandwidth;
		myDataGridY = myRegion[3] - myRegion[2] + 2 * blockBandwidth;

		if(NULL == (myIndexTable = (int *)malloc(sizeof(int) * myDataGridX * myDataGridY)))
		{
			printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
			exit(1);
		}

		if(NULL == (myX = (float *)malloc(sizeof(float) * myPoints)))
		{
			printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
			exit(1);
		}
		if(NULL == (myY = (float *)malloc(sizeof(float) * myPoints)))
		{
			printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
			exit(1);
		}

		MPI_Recv(myX, myPoints, MPI_FLOAT, 0, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(myY, myPoints, MPI_FLOAT, 0, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(myIndexTable, myDataGridX * myDataGridY, MPI_INT, 0, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}

	if(rank < nXRegion * nYRegion)
	{

		myCellXMin = myRegion[0] * BLOCKSIZE;
		myCellXMax = myRegion[1] * BLOCKSIZE;
		if(myCellXMax > nXCell)
		{
			myCellXMax = nXCell;
		}
		myCellYMin = myRegion[2] * BLOCKSIZE;
		myCellYMax = myRegion[3] * BLOCKSIZE;
		if(myCellYMax > nYCell)
		{
			myCellYMax = nYCell;
		}

		myXMin = xMin + myCellXMin * cellSize; 
		myXMax = xMin + myCellXMax * cellSize; 
		myYMin = yMax - myCellYMax * cellSize; 
		myYMax = yMax - myCellYMin * cellSize;

		myXCell = myCellXMax - myCellXMin;
		myYCell = myCellYMax - myCellYMin; 

		myGridX = myRegion[1] - myRegion[0];
		myGridY = myRegion[3] - myRegion[2];


/*
		if(rank == 0)
		{
			printf("Rank %d: x: %f to %f\ty: %f to %f\n", rank, myXMin, myXMax, myYMin, myYMax); 
			for(int i = 0; i < 1000; i ++)
			{
				printf("%f\t%f\n", myX[i], myY[i]);
			}
		}
*/
//		printf("Region %d\tThe region I get: x:\t%d - \t%d\t y\t%d - \t%d\n", rank, myCellXMin, myCellXMax, myCellYMin, myCellYMax);
//		printf("Region %d\tThe region I get: x:\t%f - \t%f\t y\t%f - \t%f\n", rank, myXMin, myXMax, myYMin, myYMax);

//		printf("Rank %d: x: %f to %f\ty: %f to %f\n", rank, myXMin, myXMax, myYMin, myYMax); 

		gettimeofday(&tBegin, NULL);

		myDensity = performKDE(myX, myY, myPoints, myIndexTable, myGridX, myGridY, blockBandwidth, myXCell, myYCell, myXMin, myYMax, cellSize, bandwidth);

		gettimeofday(&tEnd, NULL);
		printf("****************************************\nNode %d:\n\tStudy area:\txMin: %f, xMax: %f, yMin: %f, yMax %f\n\tOutput raster size: %d(row) * %d(column)\n\tNumber of input points: %d\n  Time - KDE on this node:\t%lfms\n", rank, myXMin, myXMax, myYMin, myYMax, myYCell, myXCell, myPoints, ((&tEnd)->tv_sec - (&tBegin)->tv_sec) * 1000 + (double)((&tEnd)->tv_usec - (&tBegin)->tv_usec) / 1000);


		free(myIndexTable);
		free(myX);
		free(myY);

	}

	if(rank == 0)
	{
		float * fullDensity;
		float * densityBuffer;
		if(NULL == (fullDensity = (float *)malloc(sizeof(float) * nXCell * nYCell)))
		{
			printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
			exit(1);
		}
		if(NULL == (densityBuffer = (float *)malloc(sizeof(float) * nXCell * nYCell)))
		{
			printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
			exit(1);
		}

		gettimeofday(&tBegin, NULL);

		mergeInto(fullDensity, nXCell, myDensity, myCellXMin, myCellYMin, myXCell, myYCell);
		free(myDensity);

		for(int i = 1; i < nXRegion * nYRegion; i++)
		{
			
			myCellXMin = allRegions[4 * i] * BLOCKSIZE;
			myCellXMax = allRegions[4 * i + 1] * BLOCKSIZE;
			if(myCellXMax > nXCell)
			{
				myCellXMax = nXCell;
			}
			myCellYMin = allRegions[4 * i + 2] * BLOCKSIZE;
			myCellYMax = allRegions[4 * i + 3] * BLOCKSIZE;
			if(myCellYMax > nYCell)
			{
				myCellYMax = nYCell;
			}

			myXCell = myCellXMax - myCellXMin;
			myYCell = myCellYMax - myCellYMin; 

//			printf("Rank 0 is about to receive %d from Rank %d\n", myXCell * myYCell, i);
			MPI_Recv(densityBuffer, myXCell * myYCell, MPI_FLOAT, i, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

			mergeInto(fullDensity, nXCell, densityBuffer, myCellXMin, myCellYMin, myXCell, myYCell);
		}
		free(allRegions);
		free(densityBuffer);
	
		gettimeofday(&tEnd, NULL);
		printf("****************************************\nTime - Merge result:\t%lfms\n", ((&tEnd)->tv_sec - (&tBegin)->tv_sec) * 1000 + (double)((&tEnd)->tv_usec - (&tBegin)->tv_usec) / 1000);


		tBegin = tEnd;

/*
		for(int i = 0; i < nYCell; i++)
		{
			for(int j = 0; j < nXCell; j++)
			{
				printf("\t%f", fullDensity[i * nXCell + j]);
			}
			printf("\n");
		}
*/
//Output
		writeGeoTiff(argv[2], fullDensity, nXCell, nYCell, xMin, yMax, cellSize);

		gettimeofday(&tEnd, NULL);
		printf("Time - Output GeoTiff:\t%lfms\n", ((&tEnd)->tv_sec - (&tBegin)->tv_sec) * 1000 + (double)((&tEnd)->tv_usec - (&tBegin)->tv_usec) / 1000);

		free(fullDensity);

		printf("Finished !\n");
	}
	else if(rank < nXRegion * nYRegion)
	{	
//		printf("Rank %d is about to send %d to Rank 0\n", rank, myXCell * myYCell);	
		MPI_Send(myDensity, myXCell * myYCell, MPI_FLOAT, 0, 5, MPI_COMM_WORLD);
		free(myDensity);
	}

	MPI_Finalize();

	return 0;
}

