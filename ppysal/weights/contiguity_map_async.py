#Contiguity using map_async

import pysal
from binning import bin_shapefile, bbcommon
from collections import defaultdict
import multiprocessing as mp
import time
import sys

def check_joinb(iterable, weight_type='ROOK'):
    
    polygon_ids = iterable
    
    weight_type = weight_type.upper()
    w = {}

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
                        d = w[polyId]
                        d.add(j)
                        w[polyId] = d
                        if j not in w:
                            w[j] = set()                       
                        k = w[j]
                        k.add(polyId)
                        w[j] = k
                        break
        return w
    else:
        print 'unsupported weight type'
        return None 

def pool_map(res,cores):
    t9 = time.time()
    #cores = mp.cpu_count()
    ddict = defaultdict(set)
    
    def callback_dict_map(w):
        for l in w:
            for key, value in l.items():
                for v in value:
                    ddict[key].add(v)     
 
    #Get offsets
    n = len(res['shapes'])
    starts = range(0,n,n/cores)
    ends = starts[1:]
    ends.append(n)  
    offsets = [ range(z[0],z[1]) for z in zip(starts, ends) ]     
    
    cores = mp.cpu_count()
    pool = mp.Pool()
       
    r = pool.map_async(check_joinb, offsets,callback=callback_dict_map)
    ta = time.time()
    r.wait() #This is killing speed...
    
    #Cleanup
    pool.close()
    pool.join()
    
    t10 = time.time()
    
    return (t10 - t9)

if __name__ == "__main__":

    cores = int(sys.argv[1])
    print "This version uses map_async with a callback function and {0} cores.".format(cores)

    fnames = ['1024_lattice.shp', '10000_lattice.shp', '50176_lattice.shp', '100489_lattice.shp', '1000_poly.shp', '10000_poly.shp', '50000_poly.shp', '100000_poly.shp']
    
    for fname in fnames:
        res = bin_shapefile('TestData/'+fname)
        global shapes
        global potential_neighbors 
        shapes = res['shapes']
        potential_neighbors = res['potential_neighbors']   
        
        t = pool_map(res,cores)
        
        print "{0} required {1} seconds.".format(fname, t)
