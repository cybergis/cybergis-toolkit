/** decomp.cc: Illustration of block-wise raster decomposition
 * Author: Yan Y. Liu <yanliu@illinois.edu>
 * Date: 06/18/2014
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "util.h"
// test data decomposition functions
// g++ -DDEBUG -I. -I../code/ -o decomp decomp.cc ../code/util.cc
int main(int argc, char **argv) {
	int np=atoi(argv[1]);
	int x=atoi(argv[2]);
	int y=atoi(argv[3]);
	int i;
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
	int numCells = 0;
	for (i=0; i<np; i++) {
		numCells += sizex[i] * sizey[i];
	}
	printf("rank\toffsetx\toffsety\tsizex\sizey\n");
	for (i=0; i<np; i++) {
		printf("%d\t%d\t%d\t%d\t%d\n", i, offsetx[i], offsety[i], sizex[i], sizey[i]);
	}
	free(offsetx);
	free(offsety);
	free(sizex);
	free(sizey);
	assert( numCells == x*y);
}
