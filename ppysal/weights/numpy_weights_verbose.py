#Contiguity using apply_async

import pysal as ps
from collections import defaultdict
import multiprocessing as mp
import time
import sys
import ctypes

import numpy as np
from numpy.random import randint


def check_contiguity(checks,lock,weight_type='ROOK'):
    cid = mp.current_process()._name
    geoms = np.frombuffer(sgeoms)
    geoms.shape = (2,geoms.shape[0] / 2)
    offsets = np.frombuffer(soffsets)  #This is float, but should be int...
    contmatrix = np.frombuffer(scontmatrix)
    contmatrix.shape = (len(offsets), len(offsets))
    if weight_type == 'ROOK':
        for polys in checks:
            potential_neigh = polys.tolist()
            vertices = {}
            for poly in polys:
                vstart = 0
                vend = offsets[poly]
                if poly - 1 > 0:
                    vstart = offsets[int(poly) - 1]
                vertices[poly] = geoms[:,vstart:vend]
            for k, v in vertices.iteritems():
                potential_neigh.remove(k)
                root_geom = v
                for neigh in potential_neigh:
                    test_geom = vertices[neigh]
                    #If the geoms share a common vertex, we need to test for a common edge.
                    xintersects = np.intersect1d(root_geom[0], test_geom[0])
                    if len(xintersects) > 1:
                        yintersects = np.intersect1d(root_geom[1], test_geom[1])
                        if len(yintersects) > 1:
                            #We have two shared points - are they adjacent in the poly geom, i.e. an edge?
                            x1root = np.where(root_geom[0] == xintersects[0])[0]
                            x2root = np.where(root_geom[0] == xintersects[1])[0]
                            if np.absolute(x1root - x2root).any() == 1:
                                x1test = np.where(test_geom[0] == xintersects[0])[0]
                                x2test = np.where(test_geom[0] == xintersects[1])[0]
                                if np.absolute(x1test - x2test).any() == 1:
                                    with lock:
                                        contmatrix[k, neigh] += 1
                                        contmatrix[neigh, k] += 1

def global_pointers(_cgeoms, _coffsets, _contmatrix):
    global sgeoms
    global soffsets
    global scontmatrix
    sgeoms = _cgeoms
    soffsets = _coffsets
    scontmatrix = _contmatrix

if __name__ == "__main__":
    if len(sys.argv) > 1:
        cores = int(sys.argv[1])
    else:
        cores = mp.cpu_count()

    #print "This version uses apply_async with a callback function and {0} cores.".format(cores)

    #fnames = ['1024_lattice.shp', '10000_lattice.shp', '50176_lattice.shp', '100489_lattice.shp', '1000_poly.shp', '10000_poly.shp', '50000_poly.shp', '100000_poly.shp']
    fnames = ['2500_poly.shp']
    for fname in fnames:
        ta = time.time()  #Global time keeper

        t1 = time.time()
        #Phase 1: Bin the shapefile
        shpFileObject = ps.open(fname)
        t2 = time.time()
        print "Reading the shapefile took {} seconds".format(t2-t1)

        t1 = time.time()
        if shpFileObject.type != ps.cg.Polygon:
            break
        t2 = time.time()
        print "Checking the geometry took {} seconds".format(t2-t1)

        t1 = time.time()
        shapebox = shpFileObject.bbox      # bounding box
        numPoly = len(shpFileObject)
        t2 = time.time()
        print "Getting the BBox and length took {} seconds".format(t2-t1)

        t1 = time.time()
        t3 = time.time()
        ranseq = sorted([randint(0,numPoly) for r in xrange(5)])
        geomx = []
        geomy = []
        bboxes = np.empty((numPoly, 4))
        pieces = 0
        total_perim = 0
        lens = np.empty(numPoly)
        t4 = time.time()
        for g in xrange(numPoly):
            shpobj = shpFileObject.get(g)
            x, y = zip(*shpobj.vertices)
            geomx += x
            geomy += y
            lens[g] = shpobj.len
            bboxes[g][:] = shpobj.bounding_box[:]  #Add 0.3 seconds for 5625 super inefficient!
            if g in ranseq:
                pieces += lens[g] - 1
                total_perim += shpobj.perimeter
        cellsize = total_perim / pieces * 1.
        cellsize *= 2  #This needs to be tests, is a cell size of l better ot l*c
        geoms = np.empty((2, len(geomx)))
        geoms[0] = geomx
        geoms[1] = geomy
        del geomx, geomy
        t2 = time.time()
        print "***THIS IS ALL READ TIME***"
        print "Flattening vertices and cellsize computation required {} seconds".format(t2 - t1)
        print "     Within this {} seconds were used for allocation".format(t4-t3)
        print "***DONE READING***"
        print "Processing with a cell size of {} units".format(cellsize)

        t1 = time.time()
        xdimension = abs(int((shapebox[2] - shapebox[0]) / cellsize))
        ydimension = abs(int((shapebox[3] - shapebox[1]) / cellsize))
        #Partition the space into a regular grid
        xmesh = np.linspace(shapebox[0], shapebox[2], xdimension)
        ymesh = np.linspace(shapebox[1], shapebox[3], ydimension)
        xv, yv = np.meshgrid(xmesh,ymesh)
        memship = np.empty((numPoly, 5), dtype=np.int)
        #Intersect the BBoxes with the meshgrid
        memship[:,2] = np.searchsorted(yv[:,0], bboxes[:,1], side='left')
        memship[:,3] = np.searchsorted(yv[:,0], bboxes[:,3], side='left')
        memship[:,0] = np.searchsorted(xv[0], bboxes[:,0], side='left')
        memship[:,1] = np.searchsorted(xv[0], bboxes[:,2], side='left')
        #Fix floating point inaccuracies, i.e. all the 0s and all the max + 1 values
        ystart = memship[:,2]
        ystart[ystart == 0] = 1
        xstart = memship[:,0]
        xstart[xstart == 0] = 1
        ystop = memship[:,3]
        ystop[ystop == len(yv[:,0] + 1)] = len(yv[:,0])
        xstop = memship[:,1]
        xstop[xstop == len(xv[0]) + 1] = len(xv[0])
        #Add the keys
        memship[:,4] = indices = np.arange(len(bboxes))
        #Lexicographical sort on xstart, ystart, xend, yend
        ind = np.lexsort((memship[:,0], memship[:,2], memship[:,1], memship[:,3]))
        sortmem = memship[ind]
        t2 = time.time()
        print "Getting buckets and generating data structure took {} seconds.".format(t2-t1)

        t1 = time.time()
        potential_neighbors = {}
        #Can this be vectorized or use itertools?
        for i in xrange(1, len(xv[0])):
            stepback = {}  #A list of x and y crossers that we need to deincrement x for
            crosseridx = np.where((sortmem[:,0]==i) & (sortmem[:,1]!=sortmem[:,0]))
            crosseridy = np.where((sortmem[:,0]==i)\
                                  & (sortmem[:,2]!=sortmem[:,3])\
                                  & (sortmem[:,1]!=sortmem[:,0]))
            yrollback = sortmem[crosseridy, 2]
            for j in xrange(1, len(yv[:,0])):
                #Step over all y cells in the x column
                yidx = np.logical_and(sortmem[:,0] == i, sortmem[:,2] == j)
                if len(sortmem[yidx, -1]) > 0:
                    potential_neighbors[(i,j)] = sortmem[yidx, -1]
                #Same idea as below, but with all j in this i - using bitwise operators
                # should be safe as arrays are all boolean checks.
                idx = np.where((sortmem[:,2]==j) & (sortmem[:,2]!=sortmem[:,3]) & (sortmem[:,0]==i))
                sortmem[idx,2] = (j + 1)

            #We know that all the values are sorted, so if start != end, increment
            # start until it start == end.  Then the poly is added to all
            # row / column pairs between start and end.
            sortmem[crosseridx, 0] = (i + 1)
            #Rollback the y crossers for the new x.
            sortmem[crosseridy,2] = yrollback

        t2 = time.time()
        print "Extracting vectors to polygon membership lists too {} seconds".format(t2-t1)

        t1 = time.time()
        #Can I get a vertex count from a shapefile header?
        # If so no need for lists to arrays, just allocate and pack.
        cgeoms = mp.RawArray(ctypes.c_double, geoms.size)
        npgeoms = np.frombuffer(cgeoms)
        npgeoms.shape = (2, geoms.shape[1])
        npgeoms[:] = geoms
        coffsets = mp.RawArray(ctypes.c_int, lens.size * 2)
        npoffsets = np.frombuffer(coffsets)
        npoffsets[:] = np.cumsum(lens)
        contmatrix = mp.RawArray(ctypes.c_int, (lens.size * lens.size * 2))
        npcontmatrix = np.frombuffer(contmatrix)
        npcontmatrix.shape = (lens.size, lens.size)
        npcontmatrix[:] = 0
        global_pointers(cgeoms, coffsets, contmatrix)

        t2 = time.time()
        print "Creating ctype shared memory vertices took {} seconds".format(t2-t1)


        '''
        t1 = time.time()
        cores = mp.cpu_count()
        pool = mp.Pool(cores)
        t2 = time.time()
        print "Initializing the pool of workers took {} seconds".format(t2 - t1)
        '''

        t1 = time.time()
        #We don't care what 'cell' polys are in, only that they
        # might be neighbors.
        neighbor_checks = [v for v in potential_neighbors.itervalues()]
        starts = range(0,len(neighbor_checks), len(neighbor_checks) / cores)
        stops = starts[1:]
        if len(stops) == 1:
            stops.append(len(neighbor_checks))
        offsets = [ range(z[0],z[1]) for z in zip(starts, stops)]
        t2 = time.time()
        print "Computing decomposition took {} seconds".format(t2-t1)

        t1 = time.time()
        jobs = []
        lock = mp.Lock()
        for offset in offsets:
            checks = [neighbor_checks[j] for j in offset]
            job = mp.Process(target=check_contiguity, args=(checks,lock, 'ROOK'))
            jobs.append(job)

        for job in jobs:
            job.start()
        for job in jobs:
            job.join()
        t2 = time.time()
        print "Multicore contiguity check took {} seconds".format(t2-t1)

        t1 = time.time()
        w = {}
        nonzero =np.transpose(np.nonzero(npcontmatrix))
        for i in range(numPoly):
            neigh = nonzero[nonzero[:,0] == i]
            w[i] = neigh[:,1].tolist()
        t2 = time.time()
        print "Generating a W from a sparse matrix took {} seconds".format(t2-t1)

        tb = time.time()
        print "Total processing time was {} seconds".format(tb-ta)

