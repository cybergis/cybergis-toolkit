#!/usr/bin/env python

import sys
import os
import time
import numpy as np
import pysal as ps
from fj_vect import fisher_jenks as vFisher

from mpi4py import MPI

#Override sys.execpthook
_excepthook = sys.excepthook

def excepthook(t, v, tb):
    _excepthook(t, v, tb)
    if (not MPI.Is_finalized() and MPI.Is_initialized()):
        MPI.COMM_WORLD.Abort(1)

sys.excepthook = excepthook


def gen_structure_data(n, rho):
    eps = np.random.normal(size=(n*n))
    eps.shape = (n*n, 1)
    infile = 'Ident_{}_{}.npy'.format(n, rho)
    f = os.path.join('/home/jlaura/pPysal/geoda_cluster/fisher_jenks/identity_matrices',infile)
    inverse = np.load(f, mmap_mode='r')
    y = np.dot(inverse, eps)
    inverse = None
    eps = None
    variance = np.var(y)
    return np.array(y, dtype=np.float32), variance


def get_tss(k, y, pivots):
    """
    Total sum of squares around class means

    Returns sum of squares over all class means
    """
    tss = 0
    y = np.sort(y)
    classes = np.arange(k)
    starts = [0]
    for i in range(len(pivots)):
        starts.append(pivots[i])
    pivots.append(len(y))
    breaks = zip(starts, pivots)
    for b in breaks:
        ymean = y[b[0]:b[1]].mean()
        css = y[b[0]:b[1]] - ymean
        css *= css
        tss += np.sum(css)
    return tss


def fj_generate_sample(y, pct=0.10, random=True):
    n = y.size
    if random:
        ids = np.random.random_integers(0, n - 1, n * pct)
    else:
        ids = np.arange(int(n*pct))
    yr = y[ids]
    yr[-1] = max(y)  # make sure we have the upper bound
    yr[0] = min(y)  # make sure we have the min
    return yr

def fj(realizations, n, k, p):
    local_log = np.zeros((realizations, 28))
    #Realizations
    for r in range(realizations):
        #Fully computed
        orig_data, variance = gen_structure_data(n, p)
        data = np.copy(orig_data).ravel()
        t1 = time.time()
        pivots = vFisher(data, classes=k)
        t2 = time.time()
        fj_time = t2 - t1
        fj_tss = get_tss(k, data, pivots)
        del data
        #Quantiles
        data = np.copy(orig_data).ravel()
        t1 = time.time()
        quantiles = ps.esda.mapclassify.Quantiles(data, k=k) 
        t2 = time.time()
        q_time = t2-t1
        q_tss = quantiles.get_tss()
        ordered_time = np.empty(5)
        ordered_tss = np.empty(5)
        random_time = np.empty(5)
        random_tss = np.empty(5)
        #Since each sample size is a different method, compute those in a loop
        for i, s in enumerate([0.05, 0.1, 0.15, 0.2, 0.25]):
            del data
            data = np.copy(orig_data).ravel()
            ####Ordered Sample####
            t3 = time.time()
            sample = fj_generate_sample(data, pct=s, random=False)
            #Returns pivots in the sample, not pivots in the total - hence C-TSS is way off.
            fj_ordered = vFisher(sample, classes=k)
            t4 = time.time()
            data.sort()
            pivots = []
            for cl in range(0,k-1):
                    pivots.append(np.where(data==sample[fj_ordered[cl]])[0][-1])
            ordered_tss[i] = get_tss(k, data, pivots)
            ordered_time[i] = t4-t3

            #Cleanup
            del data
            data = np.copy(orig_data).ravel()

            ####Random Sample####
            t5 = time.time()
            sample = fj_generate_sample(data, pct=s)
            fj_random = vFisher(sample, classes=k)
            t6 = time.time()

            data.sort()
            pivots = []
            for cl in range(0,k-1):
                pivots.append(np.where(data==sample[fj_random[cl]])[0][-1])
            random_tss[i] = get_tss(k, data, pivots)
            random_time[i] = t6 - t5
        sampled = np.hstack((ordered_tss,ordered_time, random_tss, random_time))

        log = [n,k,p,variance, fj_tss, fj_time, q_tss, q_time]
        log = map(str, log)
        local_log[r] = np.hstack((np.array(log),sampled))

    return local_log

t0 = MPI.Wtime()
#Read in the command line args
REALIZATIONS = int(sys.argv[1])
n = int(sys.argv[2])
k = int(sys.argv[3])
rho = float(sys.argv[4])

#MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
ncores = comm.Get_size()

#Just like single machine - compute the load for each core assuming even compute time.
comps = REALIZATIONS / ncores

#Generate the parameters to pass to the children
logs = np.zeros((REALIZATIONS, 28))
if rank == 0:
    data = np.array([comps, n, k, rho])
else:
    data = None

#Broadcast the data array to all children
#with communication:
data = comm.bcast([data, MPI.DOUBLE], root=0)

#with computation:
_logs = fj(int(data[0][0]), int(data[0][1]), int(data[0][2]), data[0][3])

#Synchronize all the cores
comm.Barrier()

#with communication:
comm.Gather([_logs, MPI.DOUBLE],[logs, MPI.DOUBLE],  root=0)

#with fj_end: pass

#with writing:
if rank == 0:
    fname = "{}_{}_{}.txt".format(n,k,rho)
    with open(fname, 'w') as logfile:
        np.savetxt(logfile, logs)
    print "Completed {} realizations with n={}, k={}, and rho={} in {} seconds.".format(REALIZATIONS, n,k,rho, MPI.Wtime() - t0)
