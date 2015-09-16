import pysal as ps
import multiprocessing as mp
import numpy as np
from pysal.weights.spatial_lag import lag_spatial as slag
from pysal.common import *
import ctypes

PERMUTATIONS = 999

def moran_mp(y,w,transformation='r',permutations=PERMUTATIONS):
    w.transform = transformation
    n,z2ss,z,EI,seI_norm,seI_rand   = get_moments(y,w)
    I = _calc(z,w,n,z2ss)
    z_norm = (I - EI) / seI_norm
    p_norm = 2.0 * (1 - stats.norm.cdf(np.abs(z_norm)))
    z_rand = (I - EI) / seI_rand
    p_rand = 2.0 * (1 - stats.norm.cdf(np.abs(z_rand)))

    global c_perm
    c_perm = mp.RawArray(ctypes.c_double, permutations)

    cores = mp.cpu_count()
    pool = mp.Pool()
    step = permutations / cores
    for start in range(0,permutations,step):
        pool.apply_async(calc, args=(z,w,n,z2ss,start,step))
    pool.close()
    pool.join()
    sim = np.frombuffer(c_perm)

    above = sim >= I
    larger = sum(above)
    if (permutations - larger) < larger:
        larger = permutations - larger
    p_sim = (larger + 1.) / (permutations + 1.)
    EI_sim = sum(sim) / permutations
    seI_sim = np.array(sim).std()
    VI_sim = seI_sim ** 2
    z_sim = (I - EI_sim) / seI_sim
    p_z_sim = 2.0 * (1 - stats.norm.cdf(np.abs(z_sim)))

    return I

def get_moments(y,w):
    n = len(y)
    z = y - y.mean()
    z2ss = sum(z*z)
    EI = -1 / (n-1)
    s1 = w.s1
    s0 = w.s0
    s2 = w.s2
    s02 = s0 * s0

    v_num = n * n * s1 - n * s2 + 3 * s0 * s0
    v_den = (n - 1) * (n + 1) * s0 * s0

    VI_norm = v_num / v_den - (1.0 / (n - 1)) ** 2
    seI_norm = VI_norm ** (1 / 2.)

    k = (1 / (sum(z ** 4)) * ((sum(z ** 2)) ** 2))
    vi = (1 / (((n - 1) ** 3) * s02)) * ((n * ((n * n - 3 * n + 3) * s1 - n * s2 + 3 * s02))
            - (k * ((n * n - n) * s1 - 2 * n * s2 + 6 * s02)))

    VI_rand = vi
    seI_rand = vi ** (1 / 2.)

    return n,z2ss,z,EI,seI_norm,seI_rand

def _calc(z,w,n,z2ss):
    zl = slag(w, z)
    inum = sum(z * zl)
    return n / w.s0 * inum / z2ss


def calc(z,w,n,z2ss,start,stop):
    pid=mp.current_process()._identity[0]
    shared_sim = np.frombuffer(c_perm)
    for i in range(start,start+stop):
        r_num = np.random.RandomState(pid + i)
        z = r_num.permutation(z)
        zl = slag(w, z)
        inum = sum(z * zl)
        shared_sim[i] =  n / w.s0 * inum / z2ss

if __name__ == "__main__":
    global w
    w = ps.open(ps.examples.get_path("stl.gal")).read()
    f = ps.open(ps.examples.get_path("stl_hom.txt"))
    y = np.array(f.by_col['HR8893'])

    ###999###
    t1 = time.time()
    mi = ps.Moran(y,  w, 'r', 999)
    t2 = time.time()
    print "PySAL serial (999): ", t2-t1
    t1 = time.time()
    mp_i = moran_mp(y, w, 'r', 999)
    t2 = time.time()
    print "PySAL MP (999): ", t2-t1
    assert(mi.I == mp_i)

    ###9999###
    t1 = time.time()
    mi = ps.Moran(y,  w, 'r', 9999)
    t2 = time.time()
    print "PySAL serial (9999): ", t2-t1
    t1 = time.time()
    mp_i = moran_mp(y, w, 'r', 9999)
    t2 = time.time()
    print "PySAL MP (9999): ", t2-t1
    assert(mi.I == mp_i)

    ###19999###
    t1 = time.time()
    mi = ps.Moran(y,  w, 'r', 19999)
    t2 = time.time()
    print "PySAL serial (19999): ", t2-t1
    t1 = time.time()
    mp_i = moran_mp(y, w, 'r', 19999)
    t2 = time.time()
    print "PySAL MP (19999): ", t2-t1
    assert(mi.I == mp_i)

    ###29999###
    t1 = time.time()
    mi = ps.Moran(y,  w, 'r', 29999)
    t2 = time.time()
    print "PySAL serial (29999): ", t2-t1
    t1 = time.time()
    mp_i = moran_mp(y, w, 'r', 29999)
    t2 = time.time()
    print "PySAL MP (29999): ", t2-t1
    assert(mi.I == mp_i)

