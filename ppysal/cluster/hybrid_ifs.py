import ctypes
from functools import partial
import multiprocessing as mp
import random
import sys

from mpi4py import MPI

import pysal as ps
from pysal.region.components import is_component

import numpy as np
from numpy.random import RandomState

from ifs import IFS, initshared_soln
from ls import LocalSearch, initshared_localsoln

def assignenclaves(rng, w, z):
    start, stop = rng
    for x in range(start,stop):
        soln_column = soln_space[x][1:]
        enclaves = np.where(soln_column == 0)[0].tolist()
        #print soln_column[enclaves], enclaves
        for e in enclaves:
            nid = set()
            potentialneighbors = w[e].keys()
            if len(potentialneighbors) == 0:
                print "ERROR NO NEIGHBORS!!!"
            for p in potentialneighbors:
                nid.add(soln_column[p])
            nid.discard(0)
            try:
                reg = random.sample(nid, 1)[0]
                soln_space[x,e + 1] = reg
            except:
                #Situation where an enclave's only neighbors are enclaves
                enclaves.append(e)
            #print e, reg
            #soln_column[e] = regi
           
            #Replace the region count with the Obj. function value
            wss = objective_func(soln_column, z)
            soln_space[x][0] = wss
    return

def objective_func(regions, z):
    """
    Computes the global objective function value
    """
    wss = 0
    for r in range(1, max(regions)):
        ids = np.where(regions == r)[0]
        m = z[ids]
        var = m.var()
        wss += np.sum(var * len(ids))
    return wss


def checkcontiguity(idx, w):
    """
    Check the contiguity of a solution in the shared memory space. Called by
    the test script to validate IFS generation.
    """
    valid = True
    start, stop = idx
    
    for idx in range(start,stop):
        soln_column = soln_space[idx]
        soln = soln_column[1:]
        nregions = int(soln_column[0])
        for i in xrange(1, nregions + 1):
            ids = np.where(soln == i)[0]
            if is_component(w, ids) != True:
                valid = False
    if valid == True:
        pass
    else:
        print "FAILURE:", soln_column.reshape(10,10)

def f(idx, w, rank):
    print "I am shared memory worker {}, managed by {}, and I see a W object with {} entries.".format(idx, rank, w.n)

initialo = np.zeros(100)
finalo = np.zeros(100)
preg = np.zeros(100)
for a in range(100):
    #MPI Boilerplate
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    nmanagers = comm.Get_size()
    status = MPI.Status()
    host = MPI.Get_processor_name()
    info = MPI.INFO_NULL

    nlocalcores = mp.cpu_count()  #One core is manager

    if rank == 0:
        """
        The rank 0 process is the master manager.  This process:

        1. Reads the data from the shapefile or DB
        2. Generates the W Object
        3. Sends the W object and attribute vector to all children
        """
        w = ps.lat2W(8,8)
        random_int = RandomState(123456789)
        attribute = random_int.random_sample((w.n, 1))
        numifs =  8

        data = {'w':w,
                'numifs':numifs}

        print "I have {} cores in a shared memory space".format(nlocalcores)
    else:
        data = None

    #Broadcast 2 sets of data, a list of Python objects and an array of attribute information
    data = comm.bcast(data, root=0) #Inefficient Python object, better to get full, pass and reform?
    if rank != 0:
        w = data['w']
        numifs = data['numifs']
        attribute = np.empty((w.n, 1), dtype=np.float)
    comm.Bcast([attribute, MPI.DOUBLE])
    """
    for r in range(nmanagers):
        if r == rank:
            print "I am manager {} with {}".format(rank, w)
    """

    solution_lock = mp.Lock()
    csoln_space = mp.Array(ctypes.c_float, numifs * (w.n + 1), lock=solution_lock)
    soln_space = np.frombuffer(csoln_space.get_obj(), dtype=np.float32)
    soln_space[:] = 0.0
    soln_space.shape = (-1, w.n + 1)
    initshared_soln(csoln_space)
    #Create a put/get memory window on each machine
    window = MPI.Win.Create(soln_space, soln_space.size, info, comm)

    jobs = []
    for i in xrange(nlocalcores):
        p = IFS(attribute, w, lock=solution_lock, pid=i, floor=7)
        jobs.append(p)
        p.start()
    for j in jobs:
        j.join()

    #Local CMAX
    localmax = np.max(soln_space[:,0])
    #This is a Python type gather, I should move to a np array gather, maybe?
    globalmax = comm.allgather(localmax)
    max_globalmax = nregions = max(globalmax)

    #Synchronize after IFS are generate without enclaves assignment
    comm.Barrier()

    group = window.Get_group()
    group.Free()

    if localmax <  max_globalmax:
        print "The local solutions on rank {} are inferior to the global best p. Updating...".format(rank)
        #Another node found a better maximum number of regions.
        idxchoices = [i for i,x in enumerate(globalmax) if x == max_globalmax]
        idx = random.choice(idxchoices)
        #Using one sided communication, get a better solution space
        window.Lock(idx)
        window.Get(soln_space, idx)
        window.Unlock(idx)

    for r in xrange(nmanagers):
        if rank == r:
            newlocalmax = max(soln_space[:,0])
            print "My soln space looks like: {}".format(soln_space[:,0])
            print "I am manager {} and I initially had {} regions.  I now have {} regions.".format(rank, localmax, newlocalmax)

    #Synchronize now that all solutions are standardized
    comm.Barrier()

    #Assign enclaves and replace region count with wss (at index 0)
    starts = range(0,numifs, 2)
    stops = starts[1:] + [numifs]
    offsets = zip(starts, stops)
    pool = mp.Pool(nlocalcores)
    for i in offsets:
        result = pool.apply(assignenclaves, args=(i, w, attribute))

    pool.close()
    pool.join()

    if len(np.where(soln_space == 0)[0]) > 0:
        print "Error in enclaves assignment! Some unit was unassigned."

    pool = mp.Pool(nlocalcores)
    for i in offsets:
        result = pool.apply(checkcontiguity, args=(i, w))

    pool.close()
    pool.join()

    comm.Barrier()

    initialbest = min(soln_space[:,0])

    #Local Search
    initshared_localsoln(csoln_space)
    jobs = []
    for i in xrange(nlocalcores):
        p = LocalSearch(attribute, w, nregions,
                lock=solution_lock,
                pid=i,
                floor=7,
                intensification=0.5,
                maxfailures = 50,
                maxiterations = 15)
        jobs.append(p)

    for j in jobs:
        j.start()
    for j in jobs:
        j.join()


    objfunc_values = soln_space[:,0]
    cbest = min(objfunc_values)
    cbest_idx = np.where(objfunc_values == cbest)[0]
    sol = soln_space[cbest_idx]
    sol = sol.ravel()[1:]

    print """Manager {}.
             Initial Best: {} 
             Local Search Best: {}
             My solution looks like:
    {}
    """.format(rank, initialbest, cbest, sol.reshape(-1, 8))
    initialo[a] = initialbest
    finalo[a] = cbest
    preg[a] = len(np.unique(sol))

print initialo
print finalo
print preg
