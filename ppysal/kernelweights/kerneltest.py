#Standard Lib. Imports
from copy_reg import pickle
import itertools
import multiprocessing as mp
import time
from types import MethodType

#PySAL
import pysal as ps
from pysal.weights.user import kernelW_from_shapefile as kws
from util import get_points_array_from_shapefile
from pysal.weights import W
from pysal.common import KDTree

#SciPy/NumPy
import numpy as np
import scipy.spatial


#Hardcoded paths to my data for testing
shapefile = '/home/jlaura/Downloads/julia/us_tracts_data.shp'
subset = '/home/jlaura/Downloads/julia/subset.shp'
tiny = '/home/jlaura/Downloads/julia/tiny.shp'
a = '/home/jlaura/Downloads/julia/25.shp'

#Hardcoded parameters
k = 2
function = 'triangular'
fixed = False
diagonal = False

ta = time.time()
#Get the centroids from the polys
t1 = time.time()
points = get_points_array_from_shapefile(shapefile)
t2 = time.time()
print "Generating points array took {} seconds.".format(t2 - t1)

#Pickle and unpickle a class method - avoids instancemethod error in mp
#http://stackoverflow.com/questions/1816958/cant-pickle-type-instancemethod-when-using-pythons-multiprocessing-pool-ma/7309686#7309686
def _pickle_method(method):
    func_name = method.im_func.__name__
    obj = method.im_self
    cls = method.im_class
    return _unpickle_method, (func_name, obj, cls)

def _unpickle_method(func_name, obj, cls):
    for cls in cls.mro():
        try:
            func = cls.__dict__[func_name]
        except KeyError:
            pass
        else:
            break
    return func.__get__(obj, cls)

def bqwrapper(datai):
    """
    Wraps the kdtree ball query for concurrent tree search.
    """
    return kdtbq(datai, r=bw[0])

def qwrapper(args):
    """
    Wraps the kdtree query for concurrent tree search
    """
    i = args[0]
    nids = args[1]
    datai = args[2]
    di, ni = kdtq(datai, k=len(nids))
    zi = np.array([dict(zip(ni, di))[nid] for nid in nids]) / bw[i]
    return zi


def loadkd(_kdtbq, _kdtq, _bw):
    """
    Initialize the child process the necessary instance methods
    and data structures to support concurrent tree search

    Since the children are in their own namespace, the global is local
    to each child.

    _kdtbq      (instancemethod) KDTree instance method query_ball_point
    _kdtq       (instancemethod) KDTree instance method query
    _bw         (ndarray) Kernel attribute bandwidth
    """
    global kdtbq
    kdtbq = _kdtbq

    global kdtq
    kdtq = _kdtq

    global bw
    bw = _bw

class Kernel(W):
    def __init__(self, data, bandwidth=None, fixed=True, k=2,
                 function='triangular', eps=1.0000001, ids=None,
                 diagonal=False, ncores=1):
        if issubclass(type(data), scipy.spatial.KDTree):
            self.kdt = data
            self.data = self.kdt.data
            data = self.data
        else:
            self.data = data
            self.kdt = KDTree(self.data)
        self.k = k + 1
        self.function = function.lower()
        self.fixed = fixed
        self.eps = eps
        self.ncores = ncores

        if bandwidth:
            try:
                bandwidth = np.array(bandwidth)
                bandwidth.shape = (len(bandwidth), 1)
            except:
                bandwidth = np.ones((len(data), 1), 'float') * bandwidth
            self.bandwidth = bandwidth
        else:
            self._set_bw()

        self._eval_kernel()
        neighbors, weights = self._k_to_W(ids)
        if diagonal:
            for i in neighbors:
                weights[i][neighbors[i].index(i)] = 1.0
        W.__init__(self, neighbors, weights, ids)

    def _k_to_W(self, ids=None):
        allneighbors = {}
        weights = {}
        if ids:
            ids = np.array(ids)
        else:
            ids = np.arange(len(self.data))
        for i, neighbors in enumerate(self.kernel):
            if len(self.neigh[i]) == 0:
                allneighbors[ids[i]] = []
                weights[ids[i]] = []
            else:
                allneighbors[ids[i]] = list(ids[self.neigh[i]])
                weights[ids[i]] = self.kernel[i].tolist()
        return allneighbors, weights

    def _set_bw(self):
        dmat, neigh = self.kdt.query(self.data, k=self.k)
        if self.fixed:
            # use max knn distance as bandwidth
            bandwidth = dmat.max() * self.eps
            n = len(dmat)
            self.bandwidth = np.ones((n, 1), 'float') * bandwidth
        else:
            # use local max knn distance
            self.bandwidth = dmat.max(axis=1) * self.eps
            self.bandwidth.shape = (self.bandwidth.size, 1)
            # identify knn neighbors for each point
            nnq = self.kdt.query(self.data, k=self.k)
            self.neigh = nnq[1]

    def _eval_kernel(self):
        t1 = time.time()
        # get points within bandwidth distance of each point
        kdtbq = self.kdt.query_ball_point
        kdtq = self.kdt.query
        bw = self.bandwidth
        if self.ncores > 1:
            pool = mp.Pool(processes=self.ncores, initializer=loadkd, initargs=(kdtbq,kdtq,bw))
        if not hasattr(self, 'neigh'):
            if self.ncores > 1:
                neighbors = pool.map(bqwrapper,self.data, chunksize = len(self.bandwidth) / self.ncores)
            else:
                neighbors = [kdtbq(self.data[i], r=bwi[0]) for i,
                            bwi in enumerate(self.bandwidth)]
            self.neigh = neighbors
        t2 = time.time()
        print "Ball Point Query took {} seconds.".format(t2 - t1)
        # get distances for neighbors
        bw = self.bandwidth

        #kdtq = self.kdt.query
        z = []
        t1 = time.time()
        if self.ncores > 1:
            iterable = [(i,nids, self.data[i]) for i, nids in enumerate(self.neigh)]
            z = pool.map(qwrapper, iterable)
        else:
            for i, nids in enumerate(self.neigh):
                di, ni = kdtq(self.data[i], k=len(nids))
                zi = np.array([dict(zip(ni, di))[nid] for nid in nids]) / bw[i]
                z.append(zi)
        t2 = time.time()
        print "Local query took: {} seconds".format(t2 - t1)
        zs = z
        # functions follow Anselin and Rey (2010) table 5.4
        if self.function == 'triangular':
            self.kernel = [1 - zi for zi in zs]
        elif self.function == 'uniform':
            self.kernel = [np.ones(zi.shape) * 0.5 for zi in zs]
        elif self.function == 'quadratic':
            self.kernel = [(3. / 4) * (1 - zi ** 2) for zi in zs]
        elif self.function == 'quartic':
            self.kernel = [(15. / 16) * (1 - zi ** 2) ** 2 for zi in zs]
        elif self.function == 'gaussian':
            c = np.pi * 2
            c = c ** (-0.5)
            self.kernel = [c * np.exp(-(zi ** 2) / 2.) for zi in zs]
        else:
            print 'Unsupported kernel function', self.function


pickle(MethodType, _pickle_method, _unpickle_method)
t1 = time.time()
w = Kernel(points, k=k, function=function, fixed=fixed,
            diagonal=diagonal, ncores=6)
t2 = time.time()

print "In total, Kernel weights took {} seconds".format(t2 - t1)
print
tb = time.time()
print "Total runtime was {} seconds".format(tb - ta)
