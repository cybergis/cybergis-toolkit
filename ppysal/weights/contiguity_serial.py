#Contiguity using map_async

import pysal
from binning import bin_shapefile, bbcommon
from collections import defaultdict
import multiprocessing as mp
import time

def check_joins(weight_type='ROOK', polygon_ids = []):
    t1 = time.time()
    
    w = {}
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
            if polyId not in w:
                w[polyId] = set()
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
                        w[polyId].add(j)
                        if j not in w:
                            w[j] = set()
                        w[j].add(polyId)
                        break
        
    else:
        print 'unsupported weight type'
        return None
    
    t2 = time.time()
    return (t2 - t1)

if __name__ == "__main__":

    print "This version runs in serial."
    
    fnames = ['1024_lattice.shp', '10000_lattice.shp', '50176_lattice.shp', '100489_lattice.shp', '1000_poly.shp', '10000_poly.shp', '50000_poly.shp', '100000_poly.shp']
    
    for fname in fnames:
        res = bin_shapefile('TestData/'+fname)
        global shapes
        global potential_neighbors 
        shapes = res['shapes']
        potential_neighbors = res['potential_neighbors']   
        
        t = check_joins()
        
        print "{0} required {1} seconds.".format(fname, t)
