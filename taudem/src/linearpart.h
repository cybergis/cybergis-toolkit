/*  Taudem parallel linear partition classes

  David Tarboton, Kim Schreuders, Dan Watson, Ahmet A. Yildirim
  Utah State University  
  May 23, 2010
  
 */

/*  Copyright (C) 2010  David Tarboton, Utah State University

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License 
version 2, 1991 as published by the Free Software Foundation.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

A copy of the full GNU General Public License is included in file 
gpl.html. This is also available at:
http://www.gnu.org/copyleft/gpl.html
or from:
The Free Software Foundation, Inc., 59 Temple Place - Suite 330, 
Boston, MA  02111-1307, USA.

If you wish to use or incorporate this program (or parts of it) into 
other software that does not meet the GNU General Public License 
conditions contact the author to request permission.
David G. Tarboton  
Utah State University 
8200 Old Main Hill 
Logan, UT 84322-8200 
USA 
http://www.engineering.usu.edu/dtarb/ 
email:  dtarb@usu.edu 
 */

//  This software is distributed from http://hydrology.usu.edu/taudem/
#ifndef LINEARPART_H
#define LINEARPART_H

#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <stdint.h>
#include <cstring>
#include <cmath>

#include <queue>
#include <limits>
#include <vector>
#include <memory>
#include <iostream>
#include <exception>

#include "mpi.h"
#include "partition.h"

template<class datatype>
class linearpart : public tdpartition {
private:
    enum BorderType {
        LEFT_BORDER,
        RIGHT_BORDER,
        TOP_BORDER,
        BOTTOM_BORDER,
        TOPLEFT_BORDER,
        TOPRIGHT_BORDER,
        BOTTOMLEFT_BORDER,
        BOTTOMRIGHT_BORDER
    };

protected:
    // Member data inherited from partition
    //long totalx, totaly;
    //long nx, ny;
    //double dx, dy;
    int rank, size;
    MPI_Comm myComm;
    int coordX, coordY;
    int sizeX, sizeY;
    int leftRank, rightRank;
    int topRank, bottomRank;
    int topleftRank, toprightRank;
    int bottomleftRank, bottomrightRank;
    int numOfNeighbours;

    int rawOffsetX = 0;
    int rawOffsetY = 0;
    int rawWidth;
    int rawHeight;

    MPI_Datatype MPI_type;
    datatype noData;

    std::unique_ptr<datatype[]> rawData;
    datatype *gridData;

    BorderType *borderTypes;
    MPI_Datatype leftrighttype;
    MPI_Datatype topbottomtype;
    datatype **gridDataBorderPointers;
    datatype **borderPointers;
    datatype **tmpBorderPointers;
    int **tmpIntBorderPointers;
    long *borderLengths;
    MPI_Datatype *dataTypes;
    int *neighbourRanks;

public:
    linearpart() : tdpartition() {}

    linearpart(long totalx, long totaly, double dx_in, double dy_in, MPI_Datatype dt, datatype nd) {
        init(totalx, totaly, dx_in, dy_in, dt, nd);
    }

    ~linearpart();

    void init(long totalx, long totaly, double dx_in, double dy_in, MPI_Datatype MPIt, datatype nd);

    bool isInPartition(int x, int y) const;
    bool hasAccess(int x, int y) const;

    void share();

    void passBorders();
    void addBorders();
    void clearBorders();

    int ringTerm(int isFinished);

    bool globalToLocal(int globalX, int globalY, int &localX, int &localY);
    void localToGlobal(int localX, int localY, int &globalX, int &globalY);

    void transferPack(int** neighbourBuffers, int**neighbourCountArr);

    void *getGridPointer() const {
        return gridData;
    }

    int getGridPointerStride() const {
        return rawWidth * sizeof(datatype);
    }

    bool isNodata(int x, int y) const;
    void setToNodata(int x, int y);

    datatype getData(int x, int y) const;
    datatype getData(int x, int y, datatype &val) const;
    void setData(int x, int y, datatype val);
    void addToData(int x, int y, datatype val);

    void savedxdyc(tiffIO &obj);
    void getdxdyc(int y, double &val_dxc, double &val_dyc);

    bool hasLeftNeighbour();
    bool hasRightNeighbour();
    bool hasTopNeighbour();
    bool hasBottomNeighbour();
    bool hasTopLeftNeighbour();
    bool hasTopRightNeighbour();
    bool hasBottomLeftNeighbour();
    bool hasBottomRightNeighbour();

    int getNeighbourCount();

private:
    void findClosestFactors(int number, int &firstFactor, int &secondFactor);

};

//Destructor.  Just frees up memory.

template<class datatype>
linearpart<datatype>::~linearpart() {
    if (gridDataBorderPointers)
        delete[] gridDataBorderPointers;

    if (borderPointers)
        delete[] borderPointers;

    if (dataTypes)
        delete[] dataTypes;

    if (neighbourRanks)
        delete[] neighbourRanks;

    if (borderTypes)
        delete[] borderTypes;

    if (tmpBorderPointers) {
        for (unsigned int i = 0; i < numOfNeighbours; ++i) {
            delete[] tmpBorderPointers[i];
        }

        delete[] tmpBorderPointers;
    }

    if (tmpIntBorderPointers) {
        for (unsigned int i = 0; i < numOfNeighbours; ++i) {
            delete[] tmpIntBorderPointers[i];
        }

        delete[] tmpIntBorderPointers;
    }

    if (borderLengths) {
        delete[] borderLengths;
    }
}

template<class datatype>
bool linearpart<datatype>::hasLeftNeighbour() {
    return leftRank != -1;
}

template<class datatype>
bool linearpart<datatype>::hasRightNeighbour() {
    return rightRank != -1;
}

template<class datatype>
bool linearpart<datatype>::hasTopNeighbour() {
    return topRank != -1;
}

template<class datatype>
bool linearpart<datatype>::hasBottomNeighbour() {
    return bottomRank != -1;
}

template<class datatype>
bool linearpart<datatype>::hasTopLeftNeighbour() {
    return topleftRank != -1;
}

template<class datatype>
bool linearpart<datatype>::hasTopRightNeighbour() {
    return toprightRank != -1;
}

template<class datatype>
bool linearpart<datatype>::hasBottomLeftNeighbour() {
    return bottomleftRank != -1;
}

template<class datatype>
bool linearpart<datatype>::hasBottomRightNeighbour() {
    return bottomrightRank != -1;
}

template<class datatype>
int linearpart<datatype>::getNeighbourCount() {
    return numOfNeighbours;
}

template<class datatype>
void linearpart<datatype>::findClosestFactors(int number, int &firstFactor, int &secondFactor) {
    int sroot = (int) sqrt((double) number);
    int i = 1;
    firstFactor = 1;
    int mindiff = INT_MAX;
    while (i <= number) {
        if (number % i == 0) {
            int diff = abs(i - sroot);
            if (diff < mindiff) {
                mindiff = diff;
                firstFactor = i;
            }
        }
        i++;
    }

    secondFactor = number / firstFactor;
}

//Init routine.  Takes the total number of rows and columns in the ENTIRE grid to be partitioned,
//dx and dy for the grid, MPI datatype (should match the template declaration), and noData value.

template<class datatype>
void linearpart<datatype>::init(long totalx, long totaly, double dx_in, double dy_in, MPI_Datatype MPIt, datatype nd) {
    MPI_Comm_rank(MCW, &rank);
    MPI_Comm_size(MCW, &size);

    leftRank = -1;
    rightRank = -1;
    topRank = -1;
    bottomRank = -1;
    topleftRank = -1;
    toprightRank = -1;
    bottomleftRank = -1;
    bottomrightRank = -1;
    numOfNeighbours = 0;
    gridDataBorderPointers = NULL;
    borderPointers = NULL;
    dataTypes = NULL;
    neighbourRanks = NULL;
    tmpBorderPointers = NULL;
    borderLengths = NULL;
    borderTypes = NULL;
    tmpIntBorderPointers = NULL;

    this->totalx = totalx;
    this->totaly = totaly;

    // ahmet: variables of dxA and dyA are never used throughout this file
    // keeping them for possible later use
    dxA = dx_in;
    dyA = dy_in;

    MPI_type = MPIt;
    noData = nd;

    if (decompType == DECOMP_BLOCK) {
        int firstFactor, secondFactor;
        // find the closest pair of factors of processor size
        findClosestFactors(size, firstFactor, secondFactor);

        // largest factor goes to largest dimension
        if (totaly > totalx) {
            sizeY = std::max(firstFactor, secondFactor);
            sizeX = std::min(firstFactor, secondFactor);
        } else {
            sizeY = std::min(firstFactor, secondFactor);
            sizeX = std::max(firstFactor, secondFactor);
        }
    } else if (decompType == DECOMP_ROW) {
        sizeX = 1;
        sizeY = size;
    } else {
        sizeX = size;
        sizeY = 1;
    }

    int dim = 2;
    int dimsize[2] = {sizeY, sizeX};
    // not circular neighbourship for both dimensions
    int periods[2] = {0, 0};
    // no reorder is required. otherwise need to pass myComm
    // to all mpi calls
    int reorder = 0;

    // create cartesian topology
    MPI_Cart_create(MPI_COMM_WORLD, dim, dimsize, periods, reorder, &myComm);

    int coords[2];
    // get the coordinate
    MPI_Cart_coords(myComm, rank, dim, coords);

    coordX = coords[1];
    coordY = coords[0];

    // Find the neighbours

    // find the rank of left neighbour
    if (coordX != 0) {
        int nrank;
        int ncoords[2] = {coordY, coordX - 1};
        MPI_Cart_rank(myComm, ncoords, &nrank);
        leftRank = nrank;
        ++numOfNeighbours;
    }

    // find the rank of right neighbour
    if (coordX != sizeX - 1) {
        int nrank;
        int ncoords[2] = {coordY, coordX + 1};
        MPI_Cart_rank(myComm, ncoords, &nrank);
        rightRank = nrank;
        ++numOfNeighbours;
    }

    // find the rank of top neighbour
    if (coordY != 0) {
        int nrank;
        int ncoords[2] = {coordY - 1, coordX};
        MPI_Cart_rank(myComm, ncoords, &nrank);
        topRank = nrank;
        ++numOfNeighbours;
    }

    // find the rank of bottom neighbour
    if (coordY != sizeY - 1) {
        int nrank;
        int ncoords[2] = {coordY + 1, coordX};
        MPI_Cart_rank(myComm, ncoords, &nrank);
        bottomRank = nrank;
        ++numOfNeighbours;
    }

    // find the rank of left-top neighbour
    if (coordX != 0 && coordY != 0) {
        int nrank;
        int ncoords[2] = {coordY - 1, coordX - 1};
        MPI_Cart_rank(myComm, ncoords, &nrank);
        topleftRank = nrank;
        ++numOfNeighbours;
    }

    // find the rank of right-top neighbour
    if (coordX != sizeX - 1 && coordY != 0) {
        int nrank;
        int ncoords[2] = {coordY - 1, coordX + 1};
        MPI_Cart_rank(myComm, ncoords, &nrank);
        toprightRank = nrank;
        ++numOfNeighbours;
    }

    // find the rank of left-bottom neighbour
    if (coordX != 0 && coordY != sizeY - 1) {
        int nrank;
        int ncoords[2] = {coordY + 1, coordX - 1};
        MPI_Cart_rank(myComm, ncoords, &nrank);
        bottomleftRank = nrank;
        ++numOfNeighbours;
    }

    // find the rank of right-bottom neighbour
    if (coordX != sizeX - 1 && coordY != sizeY - 1) {
        int nrank;
        int ncoords[2] = {coordY + 1, coordX + 1};
        MPI_Cart_rank(myComm, ncoords, &nrank);
        bottomrightRank = nrank;
        ++numOfNeighbours;
    }

    // calculate global coordinates on DEM surface
    cx = (totalx / sizeX) * coordX;
    cy = (totaly / sizeY) * coordY;
    nx = coordX != (sizeX - 1) ? totalx / sizeX : totalx - cx;
    ny = coordY != (sizeY - 1) ? totaly / sizeY : totaly - cy;

    // Figure out the bounds needed to include borders.
    rawWidth = nx;
    rawHeight = ny;

    // Left border
    if (coordX > 0) {
        rawOffsetX = 1;
        rawWidth++;
    }

    // Right border
    if (coordX < sizeX - 1) {
        rawWidth++;
    }

    // Top border
    if (coordY > 0) {
        rawOffsetY = 1;
        rawHeight++;
    }

    // Bottom border
    if (coordY < sizeY - 1) {
        rawHeight++;
    }

    //printf("rank %d @ %d %d: total %d %d n %d %d c %d %d raw %d %d %d %d\n", rank, coordX, coordY, totalx, totaly, nx, ny, cx, cy, rawOffsetX, rawOffsetY, rawWidth, rawHeight);

    // create vector type for left and right borders if necessary
    if (leftRank != -1 || rightRank != -1) {
        MPI_Type_vector(ny, 1, rawWidth, MPI_type, &leftrighttype);
        MPI_Type_commit(&leftrighttype);
    }

    // maybe not necessary but lets create
    // vector type for top and bottom border as well
    if (topRank != -1 || bottomRank != -1) {
        MPI_Type_contiguous(nx, MPI_type, &topbottomtype);
        MPI_Type_commit(&topbottomtype);
    }
    
    // Allocate memory for data and fill with noData value.  Catch exceptionfs
    uint64_t numPixels = (uint64_t) rawWidth * rawHeight;

    try {
        rawData = std::unique_ptr<datatype[]>(new datatype[numPixels]);
        gridData = rawData.get() + rawOffsetY * rawWidth + rawOffsetX;
    } catch (std::bad_alloc &) {
        //  DGT added clause below to try trap for insufficient memory in the computer.
        fprintf(stdout, "Memory allocation error during partition initialization in process %d.\n", rank);
        fprintf(stdout, "NCols: %ld, NRows: %ld, NCells: %ld\n", nx, ny, numPixels);
        fflush(stdout);
        MPI_Abort(MCW, -999);
    }

    std::fill(rawData.get(), rawData.get() + numPixels, noData);

    // store all pointers, mpi types and neighbour ranks
    // to be used in border-exchange operations
    if (numOfNeighbours > 0) {
        gridDataBorderPointers = new datatype *[numOfNeighbours];
        borderPointers = new datatype *[numOfNeighbours];
        dataTypes = new MPI_Datatype[numOfNeighbours];
        neighbourRanks = new int[numOfNeighbours];
        tmpBorderPointers = new datatype *[numOfNeighbours];
        borderLengths = new long[numOfNeighbours];
        tmpIntBorderPointers = new int *[numOfNeighbours];
        borderTypes = new BorderType[numOfNeighbours];

        int neighbourIndex = 0;
        if (topRank != -1) {
            gridDataBorderPointers[neighbourIndex] = gridData;
            borderPointers[neighbourIndex] = gridData - rawWidth;
            dataTypes[neighbourIndex] = topbottomtype;
            neighbourRanks[neighbourIndex] = topRank;
            tmpBorderPointers[neighbourIndex] = new datatype[nx];
            tmpIntBorderPointers[neighbourIndex] = new int[nx];
            borderLengths[neighbourIndex] = nx;
            borderTypes[neighbourIndex] = TOP_BORDER;

            ++neighbourIndex;
        }

        if (bottomRank != -1) {
            gridDataBorderPointers[neighbourIndex] = gridData + ((ny - 1) * rawWidth);
            borderPointers[neighbourIndex] = gridData + ny * rawWidth;
            dataTypes[neighbourIndex] = topbottomtype;
            neighbourRanks[neighbourIndex] = bottomRank;
            tmpBorderPointers[neighbourIndex] = new datatype[nx];
            tmpIntBorderPointers[neighbourIndex] = new int[nx];
            borderLengths[neighbourIndex] = nx;
            borderTypes[neighbourIndex] = BOTTOM_BORDER;

            ++neighbourIndex;
        }

        if (leftRank != -1) {
            gridDataBorderPointers[neighbourIndex] = gridData;
            borderPointers[neighbourIndex] = gridData - rawOffsetX;
            dataTypes[neighbourIndex] = leftrighttype;
            neighbourRanks[neighbourIndex] = leftRank;
            tmpBorderPointers[neighbourIndex] = new datatype[ny];
            tmpIntBorderPointers[neighbourIndex] = new int[ny];
            borderLengths[neighbourIndex] = ny;
            borderTypes[neighbourIndex] = LEFT_BORDER;

            ++neighbourIndex;
        }

        if (rightRank != -1) {
            gridDataBorderPointers[neighbourIndex] = gridData + nx - 1;
            borderPointers[neighbourIndex] = gridData + nx;
            dataTypes[neighbourIndex] = leftrighttype;
            neighbourRanks[neighbourIndex] = rightRank;
            tmpBorderPointers[neighbourIndex] = new datatype[ny];
            tmpIntBorderPointers[neighbourIndex] = new int[ny];
            borderLengths[neighbourIndex] = ny;
            borderTypes[neighbourIndex] = RIGHT_BORDER;

            ++neighbourIndex;
        }

        if (topleftRank != -1) {
            gridDataBorderPointers[neighbourIndex] = gridData;
            borderPointers[neighbourIndex] = gridData - rawWidth - 1;
            dataTypes[neighbourIndex] = MPI_type;
            neighbourRanks[neighbourIndex] = topleftRank;
            tmpBorderPointers[neighbourIndex] = new datatype[1];
            tmpIntBorderPointers[neighbourIndex] = new int[1];
            borderLengths[neighbourIndex] = 1;
            borderTypes[neighbourIndex] = TOPLEFT_BORDER;

            ++neighbourIndex;
        }

        if (toprightRank != -1) {
            gridDataBorderPointers[neighbourIndex] = gridData + nx - 1;
            borderPointers[neighbourIndex] = gridData + nx - rawWidth;
            dataTypes[neighbourIndex] = MPI_type;
            neighbourRanks[neighbourIndex] = toprightRank;
            tmpBorderPointers[neighbourIndex] = new datatype[1];
            tmpIntBorderPointers[neighbourIndex] = new int[1];
            borderLengths[neighbourIndex] = 1;
            borderTypes[neighbourIndex] = TOPRIGHT_BORDER;

            ++neighbourIndex;
        }

        if (bottomleftRank != -1) {
            gridDataBorderPointers[neighbourIndex] = gridData + ((ny - 1) * rawWidth);
            borderPointers[neighbourIndex] = gridData + ny * rawWidth - 1;
            dataTypes[neighbourIndex] = MPI_type;
            neighbourRanks[neighbourIndex] = bottomleftRank;
            tmpBorderPointers[neighbourIndex] = new datatype[1];
            tmpIntBorderPointers[neighbourIndex] = new int[1];
            borderLengths[neighbourIndex] = 1;
            borderTypes[neighbourIndex] = BOTTOMLEFT_BORDER;

            ++neighbourIndex;
        }

        if (bottomrightRank != -1) {
            gridDataBorderPointers[neighbourIndex] = gridData + (ny * rawWidth) - 1;
            borderPointers[neighbourIndex] = gridData + ny * rawWidth;
            dataTypes[neighbourIndex] = MPI_type;
            neighbourRanks[neighbourIndex] = bottomrightRank;
            tmpBorderPointers[neighbourIndex] = new datatype[1];
            tmpIntBorderPointers[neighbourIndex] = new int[1];
            borderLengths[neighbourIndex] = 1;
            borderTypes[neighbourIndex] = BOTTOMRIGHT_BORDER;
        }
    }
}

//Returns true if (x,y) is in partition
template<class datatype>
bool linearpart<datatype>::isInPartition(int x, int y) const {
    return x >= 0 && x < nx && y >= 0 && y < ny;
}

//Returns true if (x,y) is in or on borders of partition
template<class datatype>
bool linearpart<datatype>::hasAccess(int x, int y) const {
    x += rawOffsetX;
    y += rawOffsetY;

    return x >= 0 && x < rawWidth && y >= 0 && y < rawHeight;
}

//Shares border information between adjacent processes
template<class datatype>
void linearpart<datatype>::share() {
    MPI_Status status;
    if (size <= 1) return; //if there is only one process, we're all done sharing

    MPI_Request *requests = (MPI_Request *) malloc(sizeof (MPI_Request) * numOfNeighbours * 2);

    for (unsigned int i = 0; i < numOfNeighbours; ++i) {
        MPI_Irecv(borderPointers[i], 1, dataTypes[i], neighbourRanks[i], 0, MCW, &requests[i]);
    }

    for (unsigned int i = 0; i < numOfNeighbours; ++i) {
        MPI_Isend(gridDataBorderPointers[i], 1, dataTypes[i], neighbourRanks[i], 0, MCW,
                &requests[i + numOfNeighbours]);
    }

    MPI_Status *statuses = (MPI_Status *) malloc(sizeof (MPI_Status) * numOfNeighbours * 2);

    MPI_Waitall(numOfNeighbours * 2, requests, statuses);

    delete requests;
    delete statuses;
}

//Swaps border information between adjacent processes.  In this way, no data is
//overwritten.  If this function is called a second time, the original state is
//restored.

template<class datatype>
void linearpart<datatype>::passBorders() {
    if (size <= 1) return; //if there is only one process, we're all done sharing

    // swaps the borders and store the data in the temp-data buffer of each border
    MPI_Request *requests = (MPI_Request *) malloc(sizeof (MPI_Request) * numOfNeighbours * 2);

    for (unsigned int i = 0; i < numOfNeighbours; ++i) {
        MPI_Irecv(tmpBorderPointers[i], borderLengths[i], MPI_type, neighbourRanks[i], 0, MCW, &requests[i]);
    }

    for (unsigned int i = 0; i < numOfNeighbours; ++i) {
        MPI_Isend(borderPointers[i], 1, dataTypes[i], neighbourRanks[i], 0, MCW,
                &requests[i + numOfNeighbours]);
    }

    MPI_Status *statuses = (MPI_Status *) malloc(sizeof (MPI_Status) * numOfNeighbours * 2);

    MPI_Waitall(numOfNeighbours * 2, requests, statuses);

    delete requests;
    delete statuses;

    // now transfer the tmp-border data to the corresponding border
    for (unsigned int i = 0; i < numOfNeighbours; ++i) {
        if (dataTypes[i] == leftrighttype) {
            for (int x = 0; x < borderLengths[i]; x++) {
                borderPointers[i][x*rawWidth] = tmpBorderPointers[i][x]; 
            }
        } else {
            memmove(borderPointers[i], tmpBorderPointers[i], borderLengths[i] * sizeof(datatype));
        }
    }
}

//Swaps border information between adjacent processes,
//then adds the values from received borders to the local copies.

template<class datatype>
void linearpart<datatype>::addBorders() {
    //Start by calling passBorders to get information.
    passBorders();

    uint64_t i;
    for (i = 0; i < nx; i++) {
        //Add the values passed in from other process
        if (hasAccess(i, -1)) {
            if (isNodata(i, -1) || isNodata(i, 0))
                setData(i, 0, noData);
            else
                addToData(i, 0, getData(i, -1));
        }

        if (hasAccess(i, ny)) {
            if (isNodata(i, ny) || isNodata(i, ny - 1))
                setData(i, ny - 1, noData);
            else
                addToData(i, ny - 1, getData(i, ny));
        }
    }

    for (i = 0; i < ny; i++) {
        //Add the values passed in from other process
        if (hasAccess(0, i)) {
            if (isNodata(-1, i) || isNodata(0, i))
                setData(0, i, noData);
            else
                addToData(0, i, getData(-1, i));
        }

        if (hasAccess(nx, i)) {
            if (isNodata(nx, i) || isNodata(nx - 1, i))
                setData(nx - 1, i, noData);
            else
                addToData(nx - 1, i, getData(nx, i));
        }
    }

    if (hasAccess(-1, -1)) {
        if (isNodata(-1, -1) || isNodata(0, 0))
            setData(0, 0, noData);
        else
            addToData(0, 0, getData(-1, -1));
    }

    if (hasAccess(nx, -1)) {
        if (isNodata(nx, -1) || isNodata(nx - 1, 0))
            setData(nx - 1, 0, noData);
        else
            addToData(nx - 1, 0, getData(nx, -1));
    }

    if (hasAccess(-1, ny)) {
        if (isNodata(-1, ny) || isNodata(0, ny - 1))
            setData(0, ny - 1, noData);
        else
            addToData(0, ny - 1, getData(-1, ny));
    }
    
    if (hasAccess(nx, ny)) {
        if (isNodata(nx, ny) || isNodata(nx - 1, ny - 1))
            setData(nx - 1, ny - 1, noData);
        else
            addToData(nx - 1, ny - 1, getData(nx, ny));
    }
}

//Clears borders (sets them to zero).

template<class datatype>
void linearpart<datatype>::clearBorders() {
    for (int x = -1; x <= nx; x++) {

        if (hasAccess(x, -1)) setData(x, -1, 0);
        if (hasAccess(x, ny)) setData(x, ny, 0);
    }

    for (int y = -1; y <= ny; y++) {
        if (hasAccess(-1, y)) setData(-1, y, 0);
        if (hasAccess(nx, y)) setData(nx, y, 0);
    }
}

template<class datatype>
int linearpart<datatype>::ringTerm(int isFinished) {
    if (size == 1)
        return isFinished;

    int ringBool;
    MPI_Allreduce(&isFinished, &ringBool, 1, MPI_INT, MPI_MIN, MCW);

    return ringBool;
}

//Converts global coordinates (for the whole grid) to local coordinates (for this
//partition).  Function returns TRUE only if the coordinates are contained
//in this partition.

template<class datatype>
bool linearpart<datatype>::globalToLocal(int globalX, int globalY, int &localX, int &localY) {
    double minWidthPerProc = totalx / sizeX;
    double minHeightPerProc = totaly / sizeY;

    localX = globalX - coordX * minWidthPerProc;
    localY = globalY - coordY * minHeightPerProc;

    return isInPartition(localX, localY);
}

//Converts local coordinates (for this partition) to the whole grid.

template<class datatype>
void linearpart<datatype>::localToGlobal(int localX, int localY, int &globalX, int &globalY) {
    double minWidthPerProc = totalx / sizeX;
    double minHeightPerProc = totaly / sizeY;

    globalX = coordX * minWidthPerProc + localX;
    globalY = coordY * minHeightPerProc + localY;
}

template<class datatype>
void linearpart<datatype>::transferPack(int** neighbourBuffers, int**neighbourCountArr) {
    if (size == 1) return;
    
    MPI_Request *requests = (MPI_Request *) malloc(sizeof (MPI_Request) * numOfNeighbours * 2);

    for (unsigned int i = 0; i < numOfNeighbours; ++i) {
        MPI_Isend(neighbourBuffers[i], *neighbourCountArr[i], MPI_INT, neighbourRanks[i], 0, MCW,
                &requests[i]);
    }

    MPI_Status *statuses = (MPI_Status *) malloc(sizeof (MPI_Status) * numOfNeighbours * 2);

    MPI_Waitall(numOfNeighbours, requests, statuses);

    for (unsigned int i = 0; i < numOfNeighbours; ++i) {
        MPI_Probe(neighbourRanks[i], 0, MCW, &statuses[i]); // Blocking function this only returns when there is a message to receive
        MPI_Get_count(&statuses[i], MPI_INT, neighbourCountArr[i]);

        MPI_Irecv(tmpIntBorderPointers[i], *neighbourCountArr[i], MPI_INT, neighbourRanks[i], 0, MCW, &requests[i + numOfNeighbours]);
    }

    MPI_Waitall(numOfNeighbours, requests + numOfNeighbours, statuses + numOfNeighbours);


    for (unsigned int i = 0; i < numOfNeighbours; ++i) {
        int bufsize = *neighbourCountArr[i];
        if (bufsize > 0) {
            memcpy(neighbourBuffers[i], tmpIntBorderPointers[i], sizeof (int) * bufsize);
        }
    }

    delete requests;
    delete statuses;
}

//Returns true if grid element (x,y) is equal to noData.
template<class datatype>
bool linearpart<datatype>::isNodata(int x, int y) const {
    return getData(x, y) == noData;
}

//Sets the element in the grid to noData.
template<class datatype>
void linearpart<datatype>::setToNodata(int x, int y) {
    setData(x, y, noData);
}

//Returns the element in the grid with coordinate (x,y).
template<class datatype>
datatype linearpart<datatype>::getData(int x, int y) const {
    assert(hasAccess(x, y));
    return gridData[x + y * rawWidth];
}

// FIXME: get rid of this ugly function
template<class datatype>
datatype linearpart<datatype>::getData(int x, int y, datatype &val) const {
    val = getData(x, y);
    return val;
}

template<class datatype>
void linearpart<datatype>::savedxdyc(tiffIO &obj) {
    dxc = new double[ny];
    dyc = new double[ny];
    double minHeightPerProc = totaly / sizeY;
    int globalY = coordY * minHeightPerProc;

    for (int i = 0; i < ny; i++) {
        dxc[i] = obj.getdxc(globalY);
        dyc[i] = obj.getdyc(globalY);
        ++globalY;
    }
}

template<class datatype>
void linearpart<datatype>::getdxdyc(int y, double &val_dxc, double &val_dyc) {
    if (y >= 0 && y < ny) {
        val_dxc = dxc[y];
        val_dyc = dyc[y];
    }
}

//Sets the element in the grid to the specified value.
template<class datatype>
void linearpart<datatype>::setData(int x, int y, datatype val) {
    assert(isInPartition(x, y));
    gridData[x + y * rawWidth] = val;
}

//Increments the element in the grid by the specified value.
template<class datatype>
void linearpart<datatype>::addToData(int x, int y, datatype val) {
    assert(isInPartition(x, y));
    gridData[x + y * rawWidth] += val;
}

#endif
