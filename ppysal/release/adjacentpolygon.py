import collections
import itertools
import sys
import numpy as np
import pysal as ps
from mpi4py import MPI

import time

from globalsort import globalsort



if __name__ == '__main__':
    t1 = time.time()
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    #Phase I: Compute Hi bounds and sort the points - this get the points local to the cores as well
    if rank == 0:
        fname = sys.argv[1]
        print "Using {} cores for {} polygons".format(comm.size, fname.split('.')[0])
        t2 = time.time()
        shpfileobj = ps.open(fname)
        geoms = []
        x = []
        y = []
        for i, poly in enumerate(shpfileobj):
            for j in poly.vertices[:-1]:
                geoms.append(i)
                x.append(j[0])
                y.append(j[1])
        nvertices = len(x)
        pts = np.empty((nvertices, 3))
        pts[:,0] = x
        pts[:,1] = y
        pts[:,2] = geoms
        t3 = time.time()
        print "File I/O required {} seconds".format(t3 - t2)
    else:
        nvertices = None
   
    npivots = int(sys.argv[2])

    nvertices = comm.bcast(nvertices)
    shape = (nvertices, 3)
    comm.Barrier()
    
    if rank == 0:
        local_hi = globalsort(comm, rank, shape, pts=pts,
                              axis='y', samplesize = npivots)
    else:
        local_hi = globalsort(comm, rank, shape, pts=None,
                              axis='y', samplesize = npivots)

    ''' 
    for i in range(comm.size):
        if rank == i:
            print i, local_hi[local_hi[:,0].argsort()].shape
            sys.exit()
    '''
    comm.Barrier()

    if rank == 0:
        t4 = time.time()
        print "Global sort took {} seconds".format(t4 - t3)
        
    local_hi = local_hi[np.lexsort((local_hi[:,0], local_hi[:,1]))]
    
    '''
    for i in range(comm.size):
        if i == rank:
            print len(local_hi)
    '''

    #if rank == 0:
        #a = local_hi
        #ua, uind = np.unique(np.ascontiguousarray(a[:,:2]).view(np.dtype((np.void,a[:,:2].dtype.itemsize * a[:,:2].shape[1]))),return_inverse=True)
        #for i in range(np.max(uind) + 1):
        #print local_hi
    coincident = []
    seed = local_hi[0][:2]
    clist = set([])
    for i in local_hi:
        if np.array_equal(i[:2], seed):
            clist.add(i[2])
        else:
            coincident.append(clist)
            clist = set([i[2]])
            seed = i[:2]

    coincident.append(clist)  #Have to get the final iteration

    neighbors = collections.defaultdict(set)
    for n in coincident:
        for c in n:
            neighbors[c] = neighbors[c].union(n)
    comm.Barrier()

    '''
    for i in range(comm.size):
        if i == rank:
            print i, neighbors
            comm.Barrier()  # Just to get prints to be pretty
    '''
    if rank == 0:
        t5 = time.time()
        print "Computing local coincident points took {} seconds".format(t5 - t4)
    
    neighbors_list = comm.gather(neighbors, root=0)
    if rank == 0:
        neighbors = neighbors_list[0]
        for n in neighbors_list[1:]:
            for k, v in n.iteritems():
                try:
                    neighbors[k] = neighbors[k].union(v)
                except KeyError:
                    neighbors[k] = v
        
        for k, v in neighbors.iteritems():
            v.remove(k)

        t6 = time.time()
        print "Collecting and parsing neighbors took {} seconds".format(t6 - t5)
        t6a = time.time()
        w_mpi = ps.W(neighbors)
        t7 = time.time()
        print "Generating the PySAL W Object took {} seconds".format(t7-t6a)
        print "Total runtime was {} seconds".format(t7 - t1)
        
        
        '''
        t8 = time.time()
        w_ps = ps.queen_from_shapefile(sys.argv[1])
        t9 = time.time()
        print "Serial W generation took {} seconds".format(t9 - t8)
        for k, v in w_mpi.neighbors.iteritems():
            #print k, sorted(v), sorted(w_ps.neighbors[k])
            assert(sorted(v) == sorted(w_ps.neighbors[k]))
        t10 = time.time()
        print "Assertion that PySAL matches pPySAL took {} seconds".format(t10 - t9)
        '''
