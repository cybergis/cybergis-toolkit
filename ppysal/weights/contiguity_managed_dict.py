#Contiguity using map_async

import pysal
from binning import bin_shapefile, bbcommon
from collections import defaultdict
import multiprocessing as mp
import time
import sys

def pcheck_joins(mdict,x,step, weight_type='ROOK',polygon_ids = []):

    #PolyId Setup
    polygon_ids = range(x, x+step)
    
    if x+step > len(shapes):
        polygon_ids = range(x, len(shapes))
        
    weight_type = weight_type.upper()
    
    if not polygon_ids:
        polygon_ids = xrange(len(shapes))
    
    if weight_type == 'ROOK':
        # check for a shared edge
        edgeCache = {}
        for polyId in polygon_ids:
            if polyId not in edgeCache:
                iEdges ={}
                iVerts = shapes[polyId].vertices
                nv = len(iVerts)
                ne = nv - 1
                for i in range(ne):
                    l = iVerts[i]
                    r = iVerts[i+1]
                    iEdges[(l,r)] = []
                    iEdges[(r,l)] = []
                edgeCache[polyId] = iEdges
            nbrs = potential_neighbors[polyId]
            if polyId not in mdict:
                mdict[polyId] = set()           
            for j in nbrs:
                join = False
                if j not in edgeCache:
                    jVerts = shapes[j].vertices
                    jEdges = {}
                    nv = len(jVerts)
                    ne = nv - 1
                    for e in range(ne):
                        l = jVerts[e]
                        r = jVerts[e+1]
                        jEdges[(l,r)] = []
                        jEdges[(r,l)] = []
                    edgeCache[j] = jEdges
                for edge in edgeCache[j]:
                    if edge in edgeCache[polyId]:
                        join = True
                        d = mdict[polyId]
                        d.add(j)
                        mdict[polyId] = d
                        if j not in mdict:
                            mdict[j] = set()                        
                        k = mdict[j]
                        k.add(polyId)
                        mdict[j] = k
                        break
    else:
        print 'unsupported weight type'
        return None 

def managed_dict(res,cores):   

    t1 = time.time()
    #cores = mp.cpu_count()
    pool = mp.Pool(cores)    
    
    step = len(res['shapes']) / cores
    manager = mp.Manager()
    mdict = manager.dict() #The w
    
    for x in range(len(res['shapes'])):
        mdict[x] = set()
    jobs = [pool.Process(target=pcheck_joins, args=(mdict,x,step)) for x in range(0,len(res['shapes']), len(res['shapes'])/cores)]
    
    for job in jobs:
        job.start()
    for job in jobs:
        job.join()
    t2 = time.time()

    return (t2 - t1)

if __name__ == "__main__":

    cores = int(sys.argv[1])
    print "This version uses a managed dictionary using {0} cores.".format(cores)
    
    fnames = ['1024_lattice.shp', '10000_lattice.shp', '50176_lattice.shp', '100489_lattice.shp', '1000_poly.shp', '10000_poly.shp', '50000_poly.shp', '100000_poly.shp']
    
    for fname in fnames:
        res = bin_shapefile('TestData/'+fname)
        global shapes
        global potential_neighbors 
        shapes = res['shapes']
        potential_neighbors = res['potential_neighbors']   
        
        t = managed_dict(res,cores)
        
        print "{0} required {1} seconds.".format(fname, t)
