from cython.parallel import prange

from libc.math cimport sqrt
cimport cython

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
@cython.nonecheck(False)
cpdef updatesubmatrix(double [:] x, double[:] y, double [:,:] a,
                   int n):

    cdef int j, i

    for j in range(n):
        for i in range(n):
            y [i] = y[i]  + x[i] * a[i,j]
    return

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
@cython.nonecheck(False)
cpdef naivechol(double [:,:] a ):
    cdef unsigned int i, j, k, n
    n = a.shape[0]
    for k in range(0,n):
        a[k,k] = sqrt(a[k,k])
        for i in range(k + 1, n):
            a[i,k] = a[i,k] / a[k,k]
        for j in range(k+1, n):
            _updateelement(a, j, k, n)
            #for i in range(j,n):
                #a[i,j] = a[i,j] - a[i,k] * a[j,k]
    return

@cython.boundscheck(False)
@cython.wraparound(False)
cdef _updateelement(double [:,:] a, int j, int k, int n):
    """
    Updates a single column in the matrix
    """
    for i in range(j,n):
        a[i,j] = a[i,j] - a[i,k] * a[j,k]

@cython.boundscheck(False)
@cython.wraparound(False)
#@cython.nonecheck(False)
cpdef updateelement(double [:,:] a, int j, int k, int n):
    """
    Updates a single column in the matrix

    """
    for i in range(j,n):
        a[i,j] = a[i,j] - a[i,k] * a[j,k]
    return


'''
This is not working...even with openmp

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cpdef pnaivechol(double [:,:] a ):
    cdef int i, j, k, n
    n = a.shape[0]
    for k in range(0,n):
        a[k,k] = sqrt(a[k,k])
        for i in range(k + 1, n):
            a[i,k] = a[i,k] / a[k,k]
        for j in range(k+1, n):
            for i in prange(j,n, nogil=True, schedule='static'):
                a[i,j] = a[i,j] - a[i,k] * a[j,k]
    return
'''
