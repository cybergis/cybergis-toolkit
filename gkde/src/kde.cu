/**
 * kde.cu: kde calculation methods, including the CUDA kernel function
 * Authors: Yizhao Gao <ygao29@illinois.edu>
 * Date: {08/19/2015}
 */

#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>

#define BLOCKSIZE 16

__global__ void kdeKernel(float * dDensity, float * dX, float * dY, int * dIndex, int nXCell, int nYCell, float xMin, float yMax, float cellSize, float bandwidth2, int blockBandwidth)
{
	__shared__ float sX[BLOCKSIZE * BLOCKSIZE];
	__shared__ float sY[BLOCKSIZE * BLOCKSIZE];

	int i = blockIdx.y * blockDim.y + threadIdx.y;
	int j = blockIdx.x * blockDim.x + threadIdx.x;

	int idInThread = threadIdx.y * blockDim.x + threadIdx.x;

	float cellX = xMin + cellSize * (j + 0.5); 
	float cellY = yMax - cellSize * (i + 0.5);

	float density = 0.0f;
	float dist2;

	int pointProcessed;
	int pointToProcess;
	int endPoint;

	for(int k = 0; k < 1 + 2 * blockBandwidth; k ++)
	{
		int dataBID = (blockIdx.y + k) * (gridDim.x + 2 * blockBandwidth) + blockIdx.x;
		if(dataBID < 1)
		{
			pointProcessed = 0;
		}	
		else
		{
			pointProcessed = dIndex[dataBID - 1];
		}
		endPoint = dIndex[dataBID + 2 * blockBandwidth];

		pointToProcess = BLOCKSIZE * BLOCKSIZE;

		for(; pointProcessed < endPoint; pointProcessed += BLOCKSIZE * BLOCKSIZE)
		{
			if(pointProcessed + pointToProcess > endPoint)
			{	
				pointToProcess = endPoint - pointProcessed;
			}

			if(idInThread < pointToProcess)
			{
				sX[idInThread] = dX[pointProcessed + idInThread];
				sY[idInThread] = dY[pointProcessed + idInThread];
			}
			__syncthreads();

			for(int m = 0; m < pointToProcess; m++)
			{
				dist2 = (cellX - sX[m]) * (cellX - sX[m]) + (cellY - sY[m]) * (cellY - sY[m]);
				if(dist2 < bandwidth2)
				{
					density += (1 - dist2/bandwidth2);
				}
			}
			
			__syncthreads();
		}
	}

	
	if(i < nYCell && j < nXCell && i > -1 && j > -1)
	{
		dDensity[i * nXCell + j] = density * 2 / (M_PI * bandwidth2) * cellSize * cellSize;
		//dDensity[i * nXCell + j] = density;
	}
}

float * performKDE(float * xCor, float * yCor, int nPoints, int * pointIndex, int gridX, int gridY, int blockBandwidth, int cellX, int cellY, float xMin, float yMax, float cellSize, float bandwidth)
{
	cudaError_t err;

	dim3 dimBlock (BLOCKSIZE, BLOCKSIZE);
	dim3 dimGrid (gridX, gridY);

	float * myDensity;
	
	if(NULL == (myDensity = (float *) malloc(sizeof(float) * cellX * cellY)))
	{
		printf("ERROR: Out of memory at %d in file %s!\n", __LINE__, __FILE__);
		exit(1);
	}

	int dataGridX = gridX + 2 * blockBandwidth;
	int dataGridY = gridY + 2 * blockBandwidth;

	float * dX;
	float * dY;
	int * dIndex;
	float * dDensity;


	err = cudaMalloc((void **)&dX, sizeof(float) * nPoints);
	if(err != cudaSuccess)
	{
		printf("%s in %s at line %d\n", cudaGetErrorString(err), __FILE__, __LINE__);
		exit(1);
	}
	err = cudaMalloc((void **)&dY, sizeof(float) * nPoints);
	if(err != cudaSuccess)
	{
		printf("%s in %s at line %d\n", cudaGetErrorString(err), __FILE__, __LINE__);
		exit(1);
	}

	err = cudaMalloc((void **)&dIndex, sizeof(int) * dataGridX * dataGridY);
	if(err != cudaSuccess)
	{
		printf("%s in %s at line %d\n", cudaGetErrorString(err), __FILE__, __LINE__);
		exit(1);
	}
	err = cudaMalloc((void **)&dDensity, sizeof(float) * cellX * cellY);
	if(err != cudaSuccess)
	{
		printf("%s in %s at line %d\n", cudaGetErrorString(err), __FILE__, __LINE__);
		exit(1);
	}

	err = cudaMemcpy(dX, xCor, sizeof(float) * nPoints, cudaMemcpyHostToDevice);
	if(err != cudaSuccess)
	{
		printf("%s in %s at line %d\n", cudaGetErrorString(err), __FILE__, __LINE__);
		exit(1);
	}
	err = cudaMemcpy(dY, yCor, sizeof(float) * nPoints, cudaMemcpyHostToDevice);
	if(err != cudaSuccess)
	{
		printf("%s in %s at line %d\n", cudaGetErrorString(err), __FILE__, __LINE__);
		exit(1);
	}
	err = cudaMemcpy(dIndex, pointIndex, sizeof(int) * dataGridX * dataGridY, cudaMemcpyHostToDevice);
	if(err != cudaSuccess)
	{
		printf("%s in %s at line %d\n", cudaGetErrorString(err), __FILE__, __LINE__);
		exit(1);
	}


	//KDE kernel
	kdeKernel<<<dimGrid,dimBlock>>>(dDensity, dX, dY, dIndex, cellX, cellY, xMin, yMax, cellSize, bandwidth * bandwidth, blockBandwidth);

	err = cudaMemcpy(myDensity, dDensity, sizeof(float) * cellX * cellY, cudaMemcpyDeviceToHost);
	if(err != cudaSuccess)
	{
		printf("%s in %s at line %d\n", cudaGetErrorString(err), __FILE__, __LINE__);
		exit(1);
	}
	
	cudaFree(dIndex);
	cudaFree(dDensity);

	cudaFree(dX);
	cudaFree(dY);

	return myDensity;
}
