#!/usr/bin/env python

import sys
import pysal as ps
import collections
import time

from mpi4py import MPI

#MPI Boilerplate
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
ncores = comm.Get_size()
status = MPI.Status()


def genvertices(shp, pcount):
    """
    A generator that yields chunks of vertices

    Assuming that the data is spatially autocorrelated within the shapefile
    this should generate a fair number of neighbors and avoid sparse dicts.

    This would be something interesting to test: what function does SA have
    on the performance of the algorithm?

    """
    stepper = int(500)  #Maybe too big for MPI?
    numPoly = len(shp)

    vcount = 0

    vertices = collections.defaultdict(set)  # 2.7%
    while vcount < stepper and pcount < numPoly:
        newvertices = shp.get(pcount).vertices # 33%
        vcount += len(newvertices)
        for v in newvertices[:-1]: # 8%
            vertices[v].add(pcount) # 19.2%
        pcount += 1
    yield pcount,vertices
    vcount = 0
    vertices = collections.defaultdict(set)  # 2.7%

if rank == 0:
    print "Using {} cores".format(ncores)
    t1 = time.time()
    shp = ps.open(sys.argv[1])
    numPoly = len(shp)
    nworkers = ncores - 1
    closed_workers = 0
    pcount = 0
    while closed_workers < nworkers:
        data = comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)
        source = status.Get_source()
        tag = status.Get_tag()
        if tag == 1:
            pcount,vertices = next(genvertices(shp, pcount))
            if vertices:
                comm.send(vertices, dest=source, tag=4)
            #All work is done.
            else:
                comm.send(None, dest=source, tag=3)
        elif tag == 3:
            closed_workers += 1
    t2 = time.time()
    print t2 - t1

elif rank == ncores - 1:
    nworkers = ncores - 2
    mergew = collections.defaultdict(set)
    while True:
        w = comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)
        source = status.Get_source()
        tag = status.Get_tag()
        #print "Joiner getting data from {}".format(source)
        if tag == 2:
            #Merge the dicts
            for k, v in w.iteritems():
                mergew[k] = mergew[k] | v
            #print "Rank {} merged {} entries yielding a w with {} keys".format(rank, len(w), len(masterw.keys()) )
        elif tag == 3:
            nworkers -= 1
            if nworkers == 0:
                masterw = ps.W(mergew)
                print masterw.n
                break

    comm.send(None, dest=0, tag=3)

else:
    while True:
        comm.send(None, dest=0, tag=1)
        vertices = comm.recv(source=0, tag=MPI.ANY_TAG, status=status)
        tag = status.Get_tag()
        if tag == 4:
            #Step over the vertices and generate a local W
            w = collections.defaultdict(set)
            for neighbors in vertices.itervalues():
                for neighbor in neighbors:
                    w[neighbor] = w[neighbor] | neighbors
            comm.send(w, dest=ncores - 1, tag=2)  # The highest rank is the joiner
        elif tag == 3:
            comm.send(None, dest=ncores - 1, tag=3)
            break
    comm.send(None, dest=0, tag=3)

