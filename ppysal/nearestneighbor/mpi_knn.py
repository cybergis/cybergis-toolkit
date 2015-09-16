#Standard Lib. Imports
from copy_reg import pickle
import itertools
import multiprocessing as mp
import sys
import time
from types import MethodType
import cPickle
from mpi4py import MPI

import numpy as np
import scipy.spatial
from scipy.spatial import kdtree

#Monkey patch the kdtree to make it pickable
kdtree.node = kdtree.KDTree.node
kdtree.leafnode = kdtree.KDTree.leafnode
kdtree.innernode = kdtree.KDTree.innernode


#MPI Boilerplate
comm = MPI.COMM_WORLD
rank = comm.Get_rank()


if rank == 0:
    shape = (int(sys.argv[1]),2)
    npoints = shape[0]
    rstate = np.random.RandomState(123456)
    print "Using {} cores for {} points".format(comm.size, npoints)
    t1 = time.time()
    points = rstate.rand(shape[0], shape[1])
    t2 = time.time()
    print "Generating points array took {} seconds.".format(t2 - t1)

    t3 = time.time()
    kdt_class = kdtree.KDTree(points)
    kdt = MPI.pickle.dumps(kdt_class)
    t4 = time.time()
    print "KDTree generation required {} seconds".format(t4 - t3)

else:
    kdt = None
    npoints = None
    points = None

kdt_pickle = comm.bcast(kdt, root=0)
kdt = MPI.pickle.loads(kdt_pickle)
npoints = comm.bcast(npoints, root=0)

'''
#Sanity check - are the memory addresses all different? - YES
for r in range(comm.size):
    if rank == r:
        print kdt
'''
comm.Barrier()
if rank == 0:
    t5 = time.time()
    print "Communication of the KDTRee took {} seconds".format(t5 - t4)

#Compute the offsets in the kdtree.data structure
quotient, remainder = divmod(npoints, comm.size)
scattersize = list([(quotient + remainder)  ]) +\
              [quotient for i in range(comm.size - 1)]
scatteroffsets = [0] + (np.cumsum(scattersize)[:-1].tolist())
comm.Barrier()

'''
for r in range(comm.size):
    if r == rank:
        print scatteroffsets
'''

start = scatteroffsets[rank]
if rank == comm.size - 1:
    stop = None
else:
    stop = scatteroffsets[rank + 1]

local_pts = kdt.data[start:stop]

'''
for r in range(comm.size):
    if rank == r:
        print rank, local_pts.shape
'''
ni, di = kdt.query(local_pts, k=2)
nn_result = np.column_stack((ni, di[:,1]))

'''
for r in range(comm.size):
    if rank == r:
        print nn_result
'''
if rank == 0:
    t6 = time.time()
    print "KDQuery took {} seconds".format(t6 - t5)
    print "Total runtime was {} seconds".format(t6 - t1)
