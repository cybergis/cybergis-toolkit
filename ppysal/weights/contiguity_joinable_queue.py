import pysal
from binning import bin_shapefile, bbcommon
from collections import defaultdict
import multiprocessing as mp
import time
import sys

def pcheck_joins2(q,resultq, weight_type='ROOK'):
    while True:
        work = q.get()
        if work == None:
            #print "Got the pill."
            q.task_done()
            break
        
        #Unpack the args from q
        polygon_ids = work
        mdict = {}
        weight_type = weight_type.upper()
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
                    mdict[polyId] = []
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
                            d.append(j)
                            mdict[polyId] = d
                            if j not in mdict:
                                mdict[j] = []
                            k = mdict[j]
                            k.append(polyId)
                            mdict[j] = k
                            break
        #Put the resultant dict back into the queue and alert that the work is done.           
        resultq.put(mdict)
        q.task_done()
        return

def joinable_queue(res,cores):
    
    t1 = time.time()
    #cores = mp.cpu_count()

    #Create a joinable queue from which to draw cells and a solution queue to get results
    q = mp.JoinableQueue()
    resultq = mp.Queue()

    #Start up a number of child workers equal to the number of cores
    #This is a great way to manage a web service.
    jobs = [mp.Process(target=pcheck_joins2, args=(q, resultq)) for x in range(cores)]
    for job in jobs:
        job.start()
    
    n = len(res['shapes'])
    starts = range(0,n,n/cores)
    ends = starts[1:]
    ends.append(n)        
    offsets = [ range(z[0],z[1]) for z in zip(starts, ends) ]
    
    #As the jobs are loaded, they start, so we avoid some of the packing overhead.
    for i in offsets:
        q.put_nowait(i)

    #Load a poison pill into the queue to kill the children when work is done
    for i in range(cores):
        q.put_nowait(None)    

    results = []
    for i in range(len(offsets)):
        results.append(resultq.get())
    
    ddict = defaultdict(set)
    for d in (results):
        for key, value in d.items():
            for v in value:
                ddict[key].add(v)

    for job in jobs:
        job.join()
    for job in jobs:
        job.join()
    t2 = time.time()
    return (t2 - t1)

if __name__ == "__main__":

    cores = int(sys.argv[1])
    print "This version uses a joinable processing queue and {0} cores.".format(cores)
    
    fnames = ['1024_lattice.shp', '10000_lattice.shp', '50176_lattice.shp', '100489_lattice.shp', '1000_poly.shp', '10000_poly.shp', '50000_poly.shp', '100000_poly.shp']
    
    for fname in fnames:
        res = bin_shapefile('TestData/'+fname)
        global shapes
        global potential_neighbors
        shapes = res['shapes']
        potential_neighbors = res['potential_neighbors']
        t = joinable_queue(res,cores)
        
        print "{0} required {1} seconds.".format(fname, t)
