#Binning
import pysal
import numpy as np

# delta to get buckets right
DELTA = 0.000001

QUEEN = 1
ROOK = 2

# constants for bucket sizes
BUCK_SM = 8
BUCK_LG = 80
SHP_SMALL = 1000


def bbcommon(bb, bbother):
    """
    Checks for overlaps of bounding boxes. First, east-west, then north-south.
    Element 0 is west, element 2 is east, element 1 is north?, element 3 is
    south?
    All four checks must be false for chflag to be true, meaning the two
    bounding boxes do not overlap.
    """
    chflag = 0
    if not ((bbother[2] < bb[0]) or (bbother[0] > bb[2])):
        if not ((bbother[3] < bb[1]) or (bbother[1] > bb[3])):
            chflag = 1
    return chflag


def bin_shapefile(shpFile, wtype='rook', n_cols=10, n_rows=10, buff=1.0001,
        shapeDensity=4):


    shpFileObject = pysal.open(shpFile)

    if shpFileObject.type != pysal.cg.Polygon: #This is really slow...
        return False

    shapebox = shpFileObject.bbox      # bounding box

    numPoly = len(shpFileObject)
    shapes = [[]] * numPoly

    bucketmin = np.int(np.ceil(np.sqrt((numPoly * 1. /  shapeDensity))))
    #print bucketmin

    # bucket size
    #if (numPoly < SHP_SMALL):
    #    bucketmin = numPoly / BUCK_SM + 2
    #else:
    #    bucketmin = numPoly / BUCK_LG + 2
    # bucket length
    lengthx = ((shapebox[2] + DELTA) - shapebox[0]) / bucketmin
    lengthy = ((shapebox[3] + DELTA) - shapebox[1]) / bucketmin

    # initialize buckets
    columns = [set() for i in range(bucketmin)]
    rows = [set() for i in range(bucketmin)]

    minbox = shapebox[:2] * \
        2                                  # minx,miny,minx,miny
    binWidth = [lengthx, lengthy] * \
        2                              # lenx,leny,lenx,leny
    bbcache = {}
    poly2Column = [set() for i in range(numPoly)]
    poly2Row = [set() for i in range(numPoly)]
    for i in range(numPoly):
        shpObj = shpFileObject.get(i)
        bbcache[i] = shpObj.bounding_box[:]
        shapes[i] = shpObj
        projBBox = [int((shpObj.bounding_box[:][j] -
                         minbox[j]) / binWidth[j]) for j in xrange(4)]
        
        for j in range(projBBox[0], projBBox[2] + 1):
            columns[j].add(i)
            poly2Column[i].add(j)
        for j in range(projBBox[1], projBBox[3] + 1):
            rows[j].add(i)
            poly2Row[i].add(j)
    # loop over polygons rather than bins
    #The union calls are what is really slow here.
    w = {}

    for polyId in xrange(numPoly):
        idRows = poly2Row[polyId]
        idCols = poly2Column[polyId]
        rowPotentialNeighbors = set()
        colPotentialNeighbors = set()
        for row in idRows:
            rowPotentialNeighbors = rowPotentialNeighbors.union(rows[row])
        for col in idCols:
            colPotentialNeighbors = colPotentialNeighbors.union(
                columns[col])
        potentialNeighbors = rowPotentialNeighbors.intersection(
            colPotentialNeighbors)
        if polyId not in w:
            w[polyId] = set()
        for j in potentialNeighbors:
            if polyId < j:
                if bbcommon(bbcache[polyId], bbcache[j]):
                    w[polyId].add(j)

    results = {}
    results['n_polygons'] = numPoly
    results['potential_neighbors'] = w
    results['shapes'] = shapes
    shpFileObject.close()
    return results


if __name__ == '__main__':

    import pysal as ps

    import time 
    results = {}
 
    sf =  "1000x10.shp"
    #sf =  ps.examples.get_path("nat.shp")
    w = bin_shapefile(sf,shapeDensity=4) # default 
      
    for density in range(1,101):
        t1 = time.time()
        w4 = bin_shapefile(sf,shapeDensity = density)
        if w4['potential_neighbors'] != w['potential_neighbors']:
            print 'no', density
        t2 = time.time()
        #print 'density %d: %f'%(density,t2 - t1)
        results[density] = t2-t1



