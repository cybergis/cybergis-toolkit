
#Standard Lib. Imports
import collections
from copy_reg import pickle
import itertools
import multiprocessing as mp
import sys
import time
from types import MethodType
from mpi4py import MPI

import numpy as np

import pysal as ps
from pysal.cg import RTree, Rect
from pysal.cg.standalone import get_shared_segments

#MPI Boilerplate
comm = MPI.COMM_WORLD
rank = comm.Get_rank()

if rank == 0:
    t1 = time.time()
    shp = ps.open(sys.argv[1])
    t2 = time.time()
    print "File I/O took {} seconds".format(t2 - t1)

    tree_class = RTree()
    polys = []
    for i, poly in enumerate(shp):
        vertices = poly.vertices
        b = poly.bounding_box
        bbox = Rect(b.left, b.lower, b.right, b.upper)
        tree_class.insert(i,bbox)
        polys.append((b, set(vertices)))

    quotient, remainder = divmod(len(polys), comm.size)
    scattersize = list([(quotient + remainder)  ]) +\
              [quotient for i in range(comm.size - 1)]
    scatteroffsets = [0] + (np.cumsum(scattersize)[:-1].tolist()) + [None]
    #chunks = [polys[scatteroffsets[i]:scatteroffsets[i+1]] for i in range(len(scatteroffsets)-1)]
    t3 = time.time()
    print "Generating tree took {} seconds.".format(t3 - t2)
    tree = MPI._p_pickle.dumps(tree_class)
else:
    tree = None
    polys = None
    scatteroffsets = None

tree_pickle = comm.bcast(tree, root=0)
scatteroffsets = comm.bcast(scatteroffsets, root=0)
polys = comm.bcast(polys, root=0)
tree = MPI._p_pickle.loads(tree_pickle)

comm.Barrier()
if rank == 0:
    t4 = time.time()
    print "Communication of the rTRee took {} seconds".format(t4 - t3)

'''
for r in range(comm.size):
    if rank == r:
        print dir(tree), len(chunks)
'''

localw = collections.defaultdict(set)
polystart = scatteroffsets[rank]
cpid = polystart
for bbox, poly1 in polys[scatteroffsets[rank]:scatteroffsets[rank+1]]:
    potential_neighbors = tree.intersection(bbox)
    for pid in potential_neighbors:
        if pid == cpid:  #Self neighbor
            continue
        #Only queen case for now
        poly2 = polys[pid][1]
        common = poly1.intersection(poly2)
        if len(common) > 0:
            localw[pid].add(cpid)
            localw[cpid].add(pid)
    cpid += 1

if rank == 0:
    t5 = time.time()
    print "Locals RTree query took {} seconds".format(t5 - t4)

'''
for r in range(comm.size):
    if r == rank:
        print w
'''

gatherw = comm.gather(localw, root=0)
comm.Barrier()

if rank == 0:
    t6 = time.time()
    print "Gathering locals Ws took {} seconds".format(t6 - t5)

if rank == 0:
    w = gatherw[0]
    for localw in gatherw[1:]:
        for k, v in localw.iteritems():
            w[k] = w[k].union(v)
    t7 = time.time()
    print "Merge took {} seconds".format(t7 - t6)
    W = ps.W(w)
    t8 = time.time()
    print "PS W creation took {} seconds".format(t8 - t7)
    print "Total runtime was {} seconds".format(t8 - t1)
