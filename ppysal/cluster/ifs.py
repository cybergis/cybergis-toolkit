import ctypes
import multiprocessing as mp
import random

import pysal as ps
from pysal.region.components import is_component

import numpy as np
from numpy.random import RandomState

class IFS(mp.Process):

    def __init__(self, attribute, w, iterations=100, floor=3, lock=None, pid=None, minp = 0):
        mp.Process.__init__(self)
        self.attribute = attribute
        self.w = w
        self.iterations = iterations
        self.floor = floor
        self.lock = lock
        self.identity = pid
        self.currentiteration = 0
        self.nregions = 0

    def __repr__(self):
        return """
        The current state of core {} is:
            totaliterations: {}
            iterations remaining: {}
            floor: {}
            current number of regions: {}
            current region membership: {}
            current enclaves: {}
    """.format(self.identity,
               self.iterations,
               self.iterations - self.currentiteration,
               self.floor,
               None,
               None,
               None)

    def run(self):
        solnspace = np.frombuffer(shared_solnspace.get_obj(), dtype=np.float32)
        solnspace.shape = (-1, self.w.n + 1)
        while self.currentiteration < self.iterations:
            self.currentiteration += 1

            regions = []
            enclaves = []
            #All regions start as candidates
            #This is random seed selection - can look at different selection methods
            candidates = set(self.w.id_order)
            while candidates:
                seed = random.choice(tuple(candidates))
                #Remove the seed from the candidates list and add to the region
                candidates.remove(seed)
                region = set([seed])
                #Grow the region
                while len(region) < self.floor:
                    potentialaddition = set([])
                    for i in region:
                        potentialaddition = potentialaddition.union(self.w.neighbors[i])
                    #Remove all potential regions are have been assigned
                    additions = potentialaddition.intersection(candidates)
                    if additions:
                        addedunit = random.choice(tuple(additions))
                        region.add(addedunit)
                        candidates.remove(addedunit)
                    else:
                        enclaves.append(seed)
                        break
                if len(region) >= self.floor:
                    regions.append(region)

            pooreridx = np.where(solnspace[:,0] < len(regions))[0]
            currentmax = np.max(solnspace[:,0])
            if len(regions) < currentmax:
                #This solution is not the current highest p
                self.currentiteration -= 1
            elif len(pooreridx) > 0:
                #A solution with a worse total number of regions exists, replace it.
                rmap = np.zeros(self.w.n)
                for i, reg in enumerate(regions):
                    for r in reg:
                        rmap[r] = i + 1

                with self.lock:
                    idx = pooreridx[0]
                    solnspace[idx][0] = len(regions)
                    solnspace[idx][1:] = rmap
            else:
                #The solution is worse than all others
                self.currentiteration += (self.iterations - self.currentiteration)
            """
            if self.currentiteration % 25 == 0:
                print self.__repr__()
            """

def checkcontiguity(soln_column, w):
    """
    Check the contiguity of a solution in the shared memory space. Called by
    the test script to validate IFS generation.
    """
    soln = soln_column[1:]
    nregions = soln_column[0]

    for i in xrange(1, nregions + 1):
        ids = np.where(soln == i)[0]
        try:
            assert is_component(w, ids) == True
        except:
            print "Failure in connected component check"

def initshared_soln(_solnspace):
    global shared_solnspace
    shared_solnspace = _solnspace


def test(cores=None):
    """
    """
    #Test data
    w = ps.lat2W(10, 10)
    random_int = RandomState(123456789)
    attribute = random_int.random_sample((w.n, 2))

    #mp Boilerplate
    if cores == None:
        cores = mp.cpu_count()
    numifs = 20

    #Locking solution space
    solution_lock = mp.Lock()
    csoln_space = mp.Array(ctypes.c_int32, numifs * (w.n + 1), lock=solution_lock)
    soln_space = np.frombuffer(csoln_space.get_obj(), dtype=np.int32)
    soln_space[:] = 0
    soln_space.shape = (-1, w.n + 1)
    initshared_soln(csoln_space)

    jobs = []
    for i in xrange(cores):
        p = IFS(attribute, w, lock=solution_lock, pid=i)
        jobs.append(p)
        p.start()
    for j in jobs:
        j.join()

    for i in range(numifs):
        checkcontiguity(soln_space[i], w)

    """
    for i in range(numifs):
        print soln_space[i][1:].reshape(-1,10)
        print
    """
    print "Generated solution space with {} regions per solution".format(soln_space[:,0])

if __name__ == '__main__':
    test()
