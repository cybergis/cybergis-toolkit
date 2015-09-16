import collections
import copy
import multiprocessing as mp
from operator import gt, lt
import random
from random import randint, uniform, shuffle

import pysal as ps
from pysal.region.components import check_contiguity

import numpy as np
from numpy.random import RandomState

class LocalSearch(mp.Process):
    """
    Class Attributes
    ----------------
    cbest                   float   the current best known solution

    Instance Attributes
    --------------------

    failures                int     the current number of failures for this iteration
    intensificationsize     int     the size of the soln space to propagate the best soln
    tabulist                deque   of tuples in the form (unit, move to region)
    
    """

    def __init__(self ,attribute, w, nregions, lock = None, pid=None, floor=3,
            maxfailures=100, maxiterations=15, intensification=0.5):
        mp.Process.__init__(self)
        self.index = pid
        self.lock = lock
        self.w = w
        self.z = attribute
        self.floor = floor
        self.maxiterations = maxiterations
        self.wss = 0
        self.cbest = float('inf')

        #Shared memory setup
        self.solnspace = np.frombuffer(shared_solnspace.get_obj(), dtype=np.float32)
        self.solnspace.shape = (-1, self.w.n + 1)
        self.solnspacesize = self.solnspace[:,0].size
        self.nregions = nregions
        #Work on a copy of the shared memory space.
        self.solncolumn = np.empty(self.w.n)
        self.unitchooser = np.arange(len(self.solncolumn))

        #Setup for intensification and diversification
        self.intensificationsize = int(self.solnspacesize * intensification)

        #Tabu parameters
        self.failures = 0
        self.maxfailures = maxfailures + int(maxfailures * uniform(-1.1, 1.2))
        self.maxiterations = 15
        self.tabulength = self.computetabulength()
        self.tabulist = collections.deque(maxlen=self.tabulength)

        #Seed the python random number generator
        self.randstate = RandomState(self.index)
        random.seed(self.index)
        
    def __repr__(self):
        return """
        The current state of index {} is:
            working on soln {} / {}
            current region membership: 
{}
            current tabu list length: {}
            current maximum number of failures: {}
            current obj. func value: {}
            current iteration: {}
        """.format(self.index, self.index, self.solnspacesize,
                self.solncolumn.reshape(8,8), self.tabulength, self.maxfailures, self.wss,
                self.maxiterations)

    def localsearch(self):
        '''
        '''
        swapping = True
        swap_iteration = 0
        total_moves = 0

        sln = self.solncolumn
        ids = np.arange(len(sln))

        k = int(self.nregions)
        changed_regions = [1] * k
        nr = range(k)
        while swapping:
            moves_made = 0
            regionIds = [r+1 for r in nr if changed_regions[r]]
            shuffle(regionIds)
            changed_regions = [0] * k
            
            for seed in regionIds:
                local_swapping = True
                local_attempts = 0
                while local_swapping:
                    local_moves = 0
                    members = ids[sln == seed].tolist()
                    neighbors = set()

                    for member in members:
                        neighbors |= set(self.w.neighbors[member])
                   
                    neighbors -= set(members)
                    candidates = []

                    for neighbor in neighbors:
                        rid = sln[neighbor]
                        block = ids[sln == rid].tolist()
                        if len(block) <= self.floor:
                            continue
                        if check_contiguity(self.w, block, neighbor):
                            candidates.append(neighbor)
                        
                    if not candidates:
                        local_swapping = False
                    else:
                        nc = len(candidates)
                        best = None
                        cv = 0.0  #TODO Needs to be a high positive number with an aspiration func
                        for area in candidates:
                            current_internal = members
                            rid = sln[area]
                            current_outter = ids[sln == rid].tolist()
                            currentwss = self.objective_func([current_internal, current_outter])
                            new_internal = copy.copy(current_internal)
                            new_outter = copy.copy(current_outter)
                            new_internal.append(area)
                            new_outter.remove(area)
                            newwss = self.objective_func([new_internal, new_outter])
                            change = newwss - currentwss
                            old_region = int(sln[area])
                            if (area, old_region) in self.tabulist:
                                continue
                            elif change < cv:
                                best = area
                                cv = change
                            else:
                                pass
                                #Aspiration function here

                        if best:
                            area = best
                            moves_made += 1
                            old_region = int(sln[area])
                            #changed_regions is 0 based
                            changed_regions[seed - 1] = 1
                            changed_regions[old_region - 1] = 1
                            #print "Moving area: {} from {} to {}".format(area, old_region, seed)
                            sln[area] = seed 
                            self.tabulist.appendleft((area, old_region))
                            self.failures = 0
                        else:
                            self.failures += 1
                            local_swapping = False
                            #swapping = False

                if self.failures >= self.maxfailures:
                    swapping = False
            if moves_made == 0:
                swapping = False

        new_obj = self.objective_func(sln=sln)
        diversify = False
        #print sln.reshape(8,8), globalobj
        with self.lock:
            current_obj_vector = self.solnspace[:,0]
            if new_obj < current_obj_vector.all():
                sortedidx = np.argsort(current_obj_vector)[::-1]
                idx = sortedidx[0:self.intensificationsize]
                for i in idx:
                    self.solnspace[i][1:] = sln
                    self.solnspace[i][0] = new_obj
            elif new_obj < current_obj_vector[self.index]:
                self.solnspace[self.index][1:] = sln
                self.solnspace[self.index][0] = new_obj
            else:
                diversify = True
        
        if diversify:
            pass

        #Manual 1 iteration breaker
        #self.maxiterations = 0
        return
    
    def run(self):
        #Populate the initial objective function value
        with self.lock:
            self.wss = self.solncolumn[0]
        while self.maxiterations > 0: 
            with self.lock:
                #Populate the local working space
                self.solncolumn[:] = self.solnspace[self.index][1:]
                #Compute the current objective function value
                self.wss = self.solnspace[self.index,0]
            self.failures = 0  #Reset the failure counter before each iteration
            self.localsearch()
            ''' 
            #This is a constant contiguity check that can be removed once validated.
            cont = True
            for i in range(1, int(self.nregions) + 1):
                region = np.where(self.solncolumn == i)[0].tolist()

                if test_region(self.w, region) == False:
                    cont = False
            if cont == False:
                print "ERROR: ", self.__repr__()
            '''
            # Uncomment to enable cycling around the soln space
            #Increment the index counter to step around the solution space
            self.index += 1
            if self.index >= self.solnspacesize:
                self.index = 0
            self.maxiterations -= 1
            
    def computetabulength(self):
        '''Talliard 1990'''
        smin = (self.nregions - 1) * 0.9
        smax = (self.nregions - 1) * 1.1
        tabu_length = 6 + (randint(0, int(smax - smin)))
        return tabu_length

    def objective_func(self, regions=None, sln=None):
        """
        Computes the  objective function value

        Parameters
        ----------
        regions     list of regionsids, if regions is none, computed for all regions

        Returns
        -------
        wss         float wss
        """
        wss = 0
        if regions == None:
            computespace = range(1, int(self.nregions) + 1)
            for r in computespace:
                ids = np.where(sln == r)[0]
                m = self.z[ids]
                var = m.var()
                wss += np.sum(var * len(ids))
        else:
            computespace = regions
            for r in computespace:
                m = self.z[r]
                var = m.var()
                wss += np.sum(var * len(m))
        return wss 


        

def initshared_localsoln(_solnspace):
    global shared_solnspace
    shared_solnspace = _solnspace

def test_region(w,neighbors):
    d={}
    g=Graph()
    for i in neighbors:
        d[i]=[j for j in w.neighbors[i] if (j in neighbors)]
    for i in d:
        for j in d[i]:
            g.add_edge(i,j,1.0)
    cc=g.connected_components(op=gt)
    if len(cc)==1:
        return True
    else:
        return False

class Graph(object):
    def __init__(self):
        self.nodes=set()
        self.edges={}
        self.cluster_lookup={}
        self.no_link={}

    def add_edge(self,n1,n2,w):
        self.nodes.add(n1)
        self.nodes.add(n2)
        self.edges.setdefault(n1,{}).update({n2:w})
        self.edges.setdefault(n2,{}).update({n1:w})

    def connected_components(self,threshold=0.9, op=lt):
        nodes = set(self.nodes)
        components,visited =[], set()
        while len(nodes) > 0:
            connected, visited = self.dfs(nodes.pop(), visited, threshold, op)
            connected = set(connected)
            for node in connected:
                if node in nodes:
                    nodes.remove(node)
            subgraph=Graph()
            subgraph.nodes = connected
            subgraph.no_link = self.no_link
            for s in subgraph.nodes:
                for k,v in self.edges.get(s,{}).iteritems():
                    if k in subgraph.nodes:
                        subgraph.edges.setdefault(s,{}).update({k:v})
                if s in self.cluster_lookup:
                    subgraph.cluster_lookup[s] = self.cluster_lookup[s]
            components.append(subgraph)
        return components
    
    def dfs(self, v, visited, threshold, op=lt, first=None):
        aux=[v]
        visited.add(v)
        if first is None:
            first = v
        for i in (n for n, w in self.edges.get(v,{}).iteritems() \
                  if op(w, threshold) and n not in visited):
            x,y=self.dfs(i,visited,threshold,op,first)
            aux.extend(x)
            visited=visited.union(y)
        return aux, visited
