import collections
import itertools
import sys
import numpy as np
import pysal as ps
from mpi4py import MPI

from globalsort import globalsort



if __name__ == '__main__':
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    #Phase I: Compute Hi bounds and sort the points - this get the points local to the cores as well
    if rank == 0:
        fname = sys.argv[1]
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
    else:
        nvertices = None

    nvertices = comm.bcast(nvertices)
    shape = (nvertices, 3)
    comm.Barrier()

    if rank == 0:
        local_hi = globalsort(comm, rank, shape, pts=pts, axis='y')
    else:
        local_hi = globalsort(comm, rank, shape, pts=None, axis='y')

    '''
    for i in range(comm.size):
        if rank == i:
            print i, local_hi[local_hi[:,0].argsort()]
            sys.exit()
    '''
    comm.Barrier()

    #Compute coincident x geometries
    local_hi = local_hi[np.lexsort((local_hi[:,0], local_hi[:,1]))]
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

    #sys.exit()

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
        w_mpi = ps.W(neighbors)

        w_ps = ps.queen_from_shapefile(sys.argv[1])

        for k, v in w_mpi.neighbors.iteritems():
            #print k, sorted(v), sorted(w_ps.neighbors[k])
            assert(sorted(v) == sorted(w_ps.neighbors[k]))



