import lmorans
import numpy as np
cimport numpy as np
cimport cython
from libc.stdlib cimport rand, RAND_MAX, srand
'''
The original implementation used j = randrange(i+1).
This was quite slow because the C code was having
to call a pyhton module (random).  This appears to have
a significant overhead.  Therefore, rand from the libc
standard library is used to keep the C function from
having to step back into Python.

One concern is that randomness of rand() if N is low.  The
algorithm works by:

1. randomly select an index in the array (n-1)
2. swap the first value in the array with the value from 1.
3. repeat with a random selection from n-2, then n-3, ... n - (n+1)

Therefore, random permutation has the potential to exhibit serial
correlation for for all these permutations as n decreases.
'''

np.import_array()  # Suppressed one warning that I get at compile time

DTYPE_INT = np.int
ctypedef np.int_t DTYPE_t
DTYPE_FL = np.float
ctypedef np.float_t DTYPE_f


cdef extern from "time.h":
    struct tm:
        int tm_mday
        int tm_mon
        int tm_year

    ctypedef long time_t
    time_t time(time_t *tloc)

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef inline void shuffle(long [:] a):
    '''
    From: http://en.wikipedia.org/wiki/Fisher–Yates_shuffle

    Parameters
    ----------
    a  1d ndarray of type int
    '''

    cdef int i, j, n
    n = a.shape[0]
    for i in range(n-1, 0, -1):
        j = rand()
        while (n < RAND_MAX and j >= RAND_MAX - (RAND_MAX % n)):
            j = rand()
        j %= i + 1
        #j = int(rand() / RAND_MAX * (i+1))
        a[i], a[j] = a[j], a[i]

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cpdef shuffle2(long [:,:] a, long time):
    srand(time)
    cdef unsigned int r,i,j
    cdef unsigned int ax0 = a.shape[0]
    cdef unsigned int ax1 = a.shape[1]
    cdef int ax1_1 = ax1 - 1
    for r in range(ax0):
            for i in range(ax1_1, 0, -1):
                j = rand()
                while (ax1 < RAND_MAX and j >= RAND_MAX - (RAND_MAX % ax0)):
                    j = rand()
                j %= i + 1
                a[r,i], a[r,j] = a[r,j], a[r,i]

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cpdef for_i_in_wn(np.ndarray[np.float_t, ndim=1] z,
                np.ndarray[np.int_t, ndim=1] ids,
                np.ndarray[np.int_t, ndim=2] rids,
                np.ndarray[np.int_t, ndim=1] wc,
                np.ndarray[np.float_t, ndim=2] lisas,
                list w):
    #Args
    #z vector of normalized observations
    #ids vector of ids
    #rids permutations x k array of shuffled ids
    #wc vector of cardinalities in ido order
    #lisas the lisa matrix to fill
    #w a list of lists of weights
    #n number of observations

    idsi = np.empty((ids.shape[0] - 1), dtype=DTYPE_INT)
    cdef long [:] idsi_view = idsi
    cdef np.ndarray[np.float_t, ndim=2] tmp
    cdef unsigned int k, c, i
    cdef unsigned int n = ids.shape[0]

    for i in range(n):
        k = wc[i]
        idsi[:i] = ids[:i]
        idsi[i:] = ids[i+1:]
        shuffle(idsi_view)
        tmp = z[idsi[rids[:,0:k]]]
        lisas[i] = z[i] * np.sum(w[i] * tmp, axis=1)
    return lisas

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cpdef mpfori(double [:,::1] lisas,
           np.ndarray[np.float_t, ndim=1] z,
           np.ndarray[np.int_t, ndim=1] ids,
           np.ndarray[np.int_t, ndim=2] rids,
           np.ndarray[np.int_t, ndim=1] wc,
           list w,
           int start,
           int stop,
           long time):
    #Seed the random number generator
    srand(time)
    cdef np.ndarray[np.int_t, ndim=1, mode='c'] idsi = np.zeros(ids.shape[0] - 1, dtype=DTYPE_INT)
    cdef double [::1]result = np.empty(lisas.shape[1], dtype=DTYPE_FL)
    cdef unsigned int k,j,i = 0
    
    for i in range(start,stop):
        k = wc[i]
        for j in range(0, i):
            idsi[j] = ids[j]
        for j in range(i+1, ids.shape[0]):
            idsi[j] = ids[j]
        shuffle(idsi) 
        result = z[i] * (np.sum(w[i] * z[idsi[rids[:,0:k]]], axis=1))
        with nogil:
            lisas[i] = result
    
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cpdef mpfori_nonconditional(double [:,::1] lisas,
           np.ndarray[np.float_t, ndim=1] z,
           np.ndarray[np.int_t, ndim=1] ids,
           np.ndarray[np.int_t, ndim=2] rids,
           np.ndarray[np.int_t, ndim=1] wc,
           list w,
           int start,
           int stop,
           long time):
    #Seed the random number generator
    srand(time)
    cdef double [::1]result = np.empty(lisas.shape[1], dtype=DTYPE_FL)
    cdef unsigned int k,j,i = 0
    
    for i in range(start,stop):
        k = wc[i]
        #Can we just jump around with the slice of rids?
        result = z[i] * (np.sum(w[i] * z[rids[:,0:k]], axis=1))
        with nogil:
            lisas[i] = result
 


cpdef arrtest(np.ndarray[np.int_t, ndim=1, mode='c'] a):
        cdef unsigned int i
        for i in range(a.shape[0]):
            a[i] += 10


cpdef simple(double [:,::1] a, unsigned start, unsigned stop):    
    cdef unsigned int i, j
    cdef unsigned int ay = a.shape[0]
    cdef unsigned int ax = a.shape[1]
    for i in range(start, stop):
        for j in range(ax):
            a[i,j] = a[i,j] + a[i,j]

def sample_without_replacement(int popsize, int sampsize, int perms):

    #From Knuth 3.4.2S - Seminumeric Algorithms

    cdef np.ndarray shuf_idsi = np.zeros((perms,sampsize), dtype=DTYPE_INT)

    cdef double u

    cdef int n = sampsize
    cdef int N = popsize
    cdef int t = 0
    cdef m = 0
    for p in range(perms):
        t = 0
        m = 0
        while m < n:
            u = (<double>rand()/<double>RAND_MAX)
            if (N-t)*u >= n-m:
                t += 1
            else:
                shuf_idsi[p,m] = t
                t += 1
                m += 1
    return shuf_idsi
'''
#Just for testing the distirbution of the algorithm
def fisher_yates(int time, int n=100000):
    cdef int i=0
    cdef int j=0
    cdef np.ndarray x
    #Ineffecient, but easy to work with
    d = {}
    srand(time)
    for i in range(n):
        x = np.arange(3, dtype=DTYPE_INT)
        a = shuffle(x)
        k = tuple(a.tolist())
        if k not in d.iterkeys():
            d[k] = 1
        else:
            d[k] += 1
    return d
'''


