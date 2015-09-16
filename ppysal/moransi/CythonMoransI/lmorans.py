import numpy as np
import pysal
import time
from pysal.weights.spatial_lag import lag_spatial as slag
import scipy.stats as stats

import ctypes
import multiprocessing as mp


#cython
#import pyximport
#pyximport.install(setup_args={"include_dirs":np.get_include()},
                  #reload_support=True)
import shuffle1d
from shuffle1d import mpfori

PERMUTATIONS = 999

class Moran_Local:
    """Local Moran Statistics


    Parameters
    ----------
    y : n*1 array

    w : weight instance assumed to be aligned with y

    transformation : string
                     weights transformation,  default is row-standardized "r".
                     Other options include
                     "B": binary,
                     "D": doubly-standardized,
                     "U": untransformed (general weights),
                     "V": variance-stabilizing.

    permutations   : number of random permutations for calculation of pseudo
                     p_values

    Attributes
    ----------

    y            : array
                   original variable
    w            : W
                   original w object
    permutations : int
                   number of random permutations for calculation of pseudo
                   p_values
    I            : float
                   value of Moran's I
    q            : array (if permutations>0)
                   values indicate quadrat location 1 HH,  2 LH,  3 LL,  4 HL
    sim          : array (if permutations>0)
                   vector of I values for permuted samples
    p_sim        : array (if permutations>0)
                   p-value based on permutations (one-sided)
                   null: spatial randomness
                   alternative: the observed Ii is further away or extreme
                   from the median of simulated values. It is either extremelyi
                   high or extremely low in the distribution of simulated Is.
    EI_sim       : float (if permutations>0)
                   average value of I from permutations
    VI_sim       : float (if permutations>0)
                   variance of I from permutations
    seI_sim      : float (if permutations>0)
                   standard deviation of I under permutations.
    z_sim        : float (if permutations>0)
                   standardized I based on permutations
    p_z_sim      : float (if permutations>0)
                   p-value based on standard normal approximation from
                   permutations (one-sided)
                   for two-sided tests, these values should be multiplied by 2

    Examples
    --------
    >>> import pysal
    >>> import numpy as np
    >>> np.random.seed(10)
    >>> w = pysal.open(pysal.examples.get_path("desmith.gal")).read()
    >>> f = pysal.open(pysal.examples.get_path("desmith.txt"))
    >>> y = np.array(f.by_col['z'])
    >>> lm = Moran_Local(y, w, transformation = "r", permutations = 99)
    >>> lm.q
    array([4, 4, 4, 2, 3, 3, 1, 4, 3, 3])
    >>> lm.p_z_sim[0]
    0.46756830387716064

    Note random components result is slightly different values across
    architectures so the results have been removed from doctests and will be
    moved into unittests that are conditional on architectures
    """
    #@profile
    def __init__(self, y, w, transformation="r", permutations=PERMUTATIONS, cores=None):
        self.y = y
        n = len(y)
        self.n = n
        self.n_1 = n - 1
        z = y - y.mean()
        # setting for floating point noise
        orig_settings = np.seterr()
        np.seterr(all="ignore")
        sy = y.std()
        z /= sy
        np.seterr(**orig_settings)
        self.z = z
        w.transform = transformation
        self.w = w
        self.permutations = permutations
        self.den = sum(z * z)
        self.Is = self.calc(self.w, self.z)
        self.__quads()
        if permutations:
            self.rlisas = mpcrand(self, conditional=False)
            #self.__crand()
            sim = np.transpose(self.rlisas)
            above = sim >= self.Is
            larger = np.sum(above, axis = 0)
            #larger = sum(above)
            low_extreme = (self.permutations - larger) < larger
            larger[low_extreme] = self.permutations - larger[low_extreme]
            self.p_sim = (larger + 1.0) / (permutations + 1.0)
            self.sim = sim
            self.EI_sim = sim.mean()
            self.seI_sim = sim.std()
            self.VI_sim = self.seI_sim * self.seI_sim
            self.z_sim = (self.Is - self.EI_sim) / self.seI_sim
            self.p_z_sim = 1 - stats.norm.cdf(np.abs(self.z_sim))
    #@profile
    def calc(self, w, z):
        zl = slag(w, z)
        return self.n_1 * self.z * zl / self.den

    #@profile
    def __crand(self):
        """
        conditional randomization

        for observation i with ni neighbors,  the candidate set cannot include
        i (we don't want i being a neighbor of i). we have to sample without
        replacement from a set of ids that doesn't include i. numpy doesn't
        directly support sampling wo replacement and it is expensive to
        implement this. instead we omit i from the original ids,  permutate the
        ids and take the first ni elements of the permuted ids as the
        neighbors to i in each randomization.

        """
        z = np.array(self.z)
        lisas = np.zeros((self.n, self.permutations))
        n_1 = self.n - 1
        k = self.w.max_neighbors + 1
        nn = self.n - 1
        rids = np.tile(np.arange(nn), (self.permutations, 1))
        shuffle1d.shuffle2(rids, int(time.time() * 10e6))

        ids = np.arange(self.w.n)
        ido = self.w.id_order
        w = [self.w.weights[ido[i]] for i in ids]
        wc = np.array([self.w.cardinalities[ido[i]] for i in ids])

        lisas = shuffle1d.for_i_in_wn(z, ids, rids[:,-k:], wc, lisas, w)
        self.rlisas = (n_1 / self.den) * lisas

    def __quads(self):
        zl = slag(self.w, self.z)
        zp = self.z > 0
        lp = zl > 0
        pp = zp * lp
        np = (1 - zp) * lp
        nn = (1 - zp) * (1 - lp)
        pn = zp * (1 - lp)
        self.q = 1 * pp + 2 * np + 3 * nn + 4 * pn

def main(perms):
    np.random.seed(10)
    w = pysal.queen_from_shapefile(pysal.examples.get_path("NAT.shp"))
    f = pysal.open(pysal.examples.get_path('NAT.dbf'))
    y = np.array(f.by_col('POL90'))
    y = np.reshape(y, (3085,))
    t1 = time.time()
    lm = Moran_Local(y, w, transformation = "r", permutations = perms, cores=None)
    t2 = time.time()
    print "Local Morans (only) time was {}".format(t2-t1)
    print "I value: {}".format(lm.Is)
    print "p sum: {}".format(lm.p_z_sim)

def simple_test():
    #10 x 5 array
    lisas_c = mp.RawArray(ctypes.c_double, 12 * 5)
    lisas = np.frombuffer(lisas_c, dtype=np.float)
    lisas[:] = np.arange(60)
    lisas = lisas.reshape(12,5)
    print lisas
    ncores = mp.cpu_count()
    pool = mp.Pool(ncores)
    step = lisas.shape[0] / ncores
    starts = range(0, lisas.shape[0], step)
    stops = starts[1:]
    stops.append(lisas.shape[0])
    offsets =  zip(starts, stops)
    jobs =[mp.Process(target=shuffle1d.simple, args=(lisas,offsets[i][0], offsets[i][1])) for i in range(ncores)]
    for j in jobs:
        j.start()
    for j in jobs:
        j.join()
    pool.close()
    pool.join()
    print lisas

#@profile
def mpcrand(lm, conditional=True, cores=None):
    z = np.array(lm.z)
    #Has to be allocated as a C style array, not a numpy wrapped C style array.
    lisas_c = mp.RawArray(ctypes.c_double, lm.n * lm.permutations)
    lisas = np.frombuffer(lisas_c, dtype = np.float).reshape(lm.n, lm.permutations)
    n_1 = lm.n - 1
    k = lm.w.max_neighbors + 1
    nn = lm.n - 1
    rids = np.tile(np.arange(nn), (lm.permutations, 1))
    shuffle1d.shuffle2(rids, int(time.time() * 10e6))

    ids = np.arange(lm.w.n)
    ido = lm.w.id_order
    wc = np.array([lm.w.cardinalities[ido[i]] for i in ids])

    #Segment the jobs over available cores
    if cores is None:
       ncores = mp.cpu_count()
    else:
        ncores = cores
    pool = mp.Pool(ncores)
    
    #Compute the offsets - this is a decomposition by areal unit.
    step = lm.n / ncores
    starts = range(0, lm.n, step)
    stops = starts[1:]
    stops[-1] = stops[-1] + 1
    offsets = zip(starts[:-1], stops)
    #Split the weights into ncores lists 
    w = [lm.w.weights[ido[i]] for i in ids]

    if conditional == True:
        jobs =[mp.Process(target=shuffle1d.mpfori, args=(lisas, z, ids, rids, wc, w, offsets[i][0], offsets[i][1], time.time() * 10e6)) for i in range(ncores)]
    else:
        jobs = [mp.Process(target=shuffle1d.mpfori_nonconditional, args=(lisas, z, ids, rids, wc, w, offsets[i][0], offsets[i][1], time.time() * 10e6)) for i in range(ncores)]
    for j in jobs:
        j.start()
    for j in jobs: 
        j.join()
    return (n_1 / lm.den) * lisas



if __name__ == '__main__':
     main(999)
