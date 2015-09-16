import pysal as ps
import multiprocessing as mp
import numpy as np
from pysal.weights.spatial_lag import lag_spatial as slag
from pysal.common import *
import ctypes

PERMUTATIONS = 999

def localmoran_mp(y,w,transformation='r',permutations=PERMUTATIONS):
    n = len(y)
    n_1 = n - 1
    z = y - y.mean()
    orig_settings = np.seterr()
    np.seterr(all="ignore")
    sy = y.std() #This is super inefficient with numpy.  It doubles the memory space.
    z /= sy
    np.seterr(**orig_settings)
    w.transform = transformation
    den = sum(z * z)
    Is = calc(w, z, n_1, den)
    q = quads(w,z)


    #Multiprocessing Permutations
    global c_lisas
    c_lisas = mp.RawArray(ctypes.c_double,n * permutations)
    cores = mp.cpu_count()
    pool = mp.Pool()
    step = n / cores
    for start in range(0,n,step):
        pool.apply_async(crand, args=(z, w, permutations, n,start,step))
    pool.close()
    pool.join()
    lisas = np.frombuffer(c_lisas)
    lisas.shape = (n,permutations)
    lisas *= n_1 / den
    #Permutations
    #rlisas = crand(z, w, permutations, n, den, None, None) #We shouldn't return anything, write rlisas to shmem

    #Returned Info
    sim = np.transpose(lisas)
    above = sim >= Is
    larger = sum(above)
    low_extreme = (permutations - larger) < larger
    larger[low_extreme] = permutations - larger[low_extreme]
    p_sim = (larger + 1.) / (permutations + 1.)
    EI_sim = sim.mean()
    seI_sim = sim.std()
    VI_sim = seI_sim * seI_sim
    z_sim = (Is - EI_sim) / seI_sim
    p_z_sim = 1 - stats.norm.cdf(np.abs(z_sim))
    return  p_sim

def calc(w, z ,n_1, den):
    zl = slag(w,z)
    return n_1 * z * zl / den

#Send the permutations over the line to the child, not a split of the ids over all the cores.
def crand(z, weights, permutations, n,start,step):
    #pid = mp.current_process()._identity[0]
    lisas = np.frombuffer(c_lisas)
    lisas.shape = (n,permutations)
    lisas[start:start+step] = 0
    #lisas = np.zeros((n,permutations))
    rid = range(n - 1)
    prange = range(permutations) #999
    k = weights.max_neighbors + 1 # 10
    nn = n - 1 # 76
    rids = np.array([np.random.permutation(nn)[0:k] for i in prange]) #(999,10)
    ids = np.arange(weights.n) #(0 - 77)
    ido = weights.id_order #1-78
    w  = [weights.weights[ido[i]] for i in ids]
    wc = [weights.cardinalities[ido[i]] for i in ids]
    for i in range(start,start+step):
        idsi = ids[ids != i]
        np.random.shuffle(idsi)
        tmp = z[idsi[rids[:, 0:wc[i]]]]
        lisas[i] = z[i] * (w[i] * tmp).sum(1)

    #return (n_1 / den) * lisas

def quads(w,z):
    zl = slag(w,z)
    zp = z > 0
    lp = zl > 0
    pp = zp * lp
    np = (1 - zp) * lp
    nn = (1 - zp) * (1 - lp)
    pn = zp * (1 - lp)
    return 1 * pp + 2 * np + 3 * nn + 4 * pn

if __name__ == "__main__":
    w = ps.open(ps.examples.get_path("stl.gal")).read()
    f = ps.open(ps.examples.get_path("stl_hom.txt"))
    y = np.array(f.by_col['HR8893'])

    ###999###
    t1 = time.time()
    mi = ps.Moran_Local(y,  w, 'r', 999)
    t2 = time.time()
    print "PySAL serial (999): ", t2-t1
    t1 = time.time()
    mp_i = localmoran_mp(y, w, 'r', 999)
    t2 = time.time()
    print "PySAL MP (999): ", t2-t1
    #assert(mi.q.all() == mp_i.all())
    ###9999###
    t1 = time.time()
    mi = ps.Moran_Local(y,  w, 'r', 9999)
    t2 = time.time()
    print "PySAL serial (9999): ", t2-t1
    t1 = time.time()
    mp_i = localmoran_mp(y, w, 'r', 9999)
    t2 = time.time()
    print "PySAL MP (9999): ", t2-t1
    #assert(mi.I == mp_i)

    ###19999###
    t1 = time.time()
    mi = ps.Moran_Local(y,  w, 'r', 19999)
    t2 = time.time()
    print "PySAL serial (19999): ", t2-t1
    t1 = time.time()
    mp_i = localmoran_mp(y, w, 'r', 19999)
    t2 = time.time()
    print "PySAL MP (19999): ", t2-t1
    #assert(mi.I == mp_i)

    ###29999###
    t1 = time.time()
    mi = ps.Moran_Local(y,  w, 'r', 29999)
    t2 = time.time()
    print "PySAL serial (29999): ", t2-t1
    t1 = time.time()
    mp_i = localmoran_mp(y, w, 'r', 29999)
    t2 = time.time()
    print "PySAL MP (29999): ", t2-t1
    #assert(mi.I == mp_i)


