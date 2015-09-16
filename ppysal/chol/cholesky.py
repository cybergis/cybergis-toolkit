from math import sqrt
import time
from cchol import naivechol as cnaive
from cchol import updateelement
import numpy as np
from numpy.random import randint
from numpy.linalg import cholesky

import matplotlib.pyplot as plt


'''
def blaschol(a):
    """
    Non-blocked BLAS style cholesky
    """
    n = a.shape[0]
    for j in range(n):
        # u_{j,j} = \sqrt{A_{j,j} - \sum_k=1^j-1 L_^{2}_{j,k}}
        if j >= 1:
            #Is np.dot() faster or slower than np.sum(a**2)?
            a[j,j] = sqrt(a[j,j] - np.dot(a[0:j:,0], a[0:j:,0]))
        else:
            a[j,j] = sqrt(a[j,j])

        if j < n:
            for i in range(
            #Compute elements of j+1:n of column j
            x = a[j:j+1,j+1:]
            y = a[j+1:j+2, j+1:]

            #A_{i,j} - \sum_{k=1}^{j-1} L_{i,k} * L_{j,k}
            # Element wise multiplication.
            #The summation is y=\alpha * A * x + \beta y, where x and y are vectors
            # A is a submatrix, \alpha = -1 and \beta = 1.0
            #This simplifies to y=\alpha * -1 * x + y
            #submatrix = a[j+1:,j+1:]
            #x =
            #y =
            #Multiply the submatrix by 1/a[j,j]
            a[j+1:] *= (1/a[j,j])

        print a
'''
def purepythoncholesky(A):
    """Performs a Cholesky decomposition of A, which must
    be a symmetric and positive definite matrix. The function
    returns the lower variant triangular matrix, L."""
    n = len(A)

    # Create zero matrix for L
    L = [[0.0] * n for i in xrange(n)]

    # Perform the Cholesky decomposition
    for i in xrange(n):
        for k in xrange(i+1):
            tmp_sum = sum(L[i][j] * L[k][j] for j in xrange(k))

            if (i == k): # Diagonal elements
                # LaTeX: l_{kk} = \sqrt{ a_{kk} - \sum^{k-1}_{j=1} l^2_{kj}}
                L[i][k] = sqrt(A[i][i] - tmp_sum)
            else:
                # LaTeX: l_{ik} = \frac{1}{l_{kk}} \left( a_{ik} - \sum^{k-1}_{j=1} l_{ij} l_{kj} \right)
                L[i][k] = (1.0 / L[k][k] * (A[i][k] - tmp_sum))
    return L

def naivechol(a):
    """
    A naive cholesky decomposition using an unoptimized
    Cholesky-Banachiewicz implementation.

    Reverse i, j to j, i for a Cholesky-Crout implementation.

    This implementation operates in-place and is O(n^3/6)
    """
    #Cast to float since the random SPD matrix comes in as int
    #Algorithm
    n = a.shape[0]
    for k in range(0, n):
        a[k,k] = sqrt(a[k, k])
        for i in range(k + 1, n):
            a[i, k] = a[i, k] / a[k, k]
        for j in range(k + 1, n):
            for i in range(j, n):
                a[i, j] = a[i, j] - a[i, k] * a[j, k]
    #Since modified in place, set the UT to 0, offest by 1
    # from the diagonal
    idx = np.triu_indices(n, 1)
    a[idx] = 0
    return a

def mixedcython(a):
    """
    This is a naive implementation with cython slotted in
    for the inner for loop.
    """
    n = a.shape[0]
    for k in range(0, n):
        a[k,k] = sqrt(a[k, k])
        for i in range(k + 1, n):
            a[i, k] = a[i, k] / a[k, k]
        for j in range(k + 1, n):
            updateelement(a, j, k, n)
    idx = np.triu_indices(n, 1)
    a[idx] = 0
    return a

#@profile
def vectorizedchol(a):
    """
    This is the naive chol, but vectorized by column.
    """
    n = a.shape[0]
    for k in range(0, n):
        #Have to recompute every time stepping 'down' diagonal
        a[k,k] = sqrt(a[k,k])
        #Update the column
        a[k+1:,k] = a[k+1:,k] / a[k,k]
        #Update the remainder of the lower triangle
        c = 0
        for j in range(k + 1, n):  # For each remaining column
            a[c+j:, j] = a[c+j:, j] - (a[j:, k] *a[j,k])
        c += 1
    #Since modified in place, set the UT to 0, offest by 1
    # from the diagonal
    idx = np.triu_indices(n, 1)
    a[idx] = 0
    return a

def cython_implementation(a):
    """
    Calls a naive Cython implementation.
    """
    cnaive(a)
    idx = np.triu_indices(a.shape[0], 1)
    a[idx] = 0
    return a

def main():
    li = [10, 25, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 1000,
              1500, 2000, 2500, 3000]
    x = np.array(li)
    y = np.zeros((6, x.shape[0]))

    for i, s in enumerate(li):
        print '\n'
        print "n = {}".format(s**2)
        print '---------------------------'
        a = np.asfortranarray(randint(0,10,size=(s, s)), dtype='d')
        #Conversion using copy from int64 to float64
        #Dot should ensure positive semi-definite
        data = np.dot(a, a.T)
        del a

        #Numpy Cholesky Decomposition (Dense Matrix)
        spd_matrix = data.copy()
        t1 = time.time()
        np_chol = cholesky(spd_matrix)
        t2 = time.time()
        print "Numpy Implementation took {} seconds".format(t2-t1)
        y[0][i] = t2-t1

        if s <= 200:
            #Naive
            spd_matrix = data.copy()
            t1 = time.time()
            out = naivechol(spd_matrix)
            t2 = time.time()
            np.testing.assert_allclose(out, np_chol, rtol=1e-07)
            print "Naive implementation took {} seconds".format(t2-t1)
            y[1][i] = t2-t1
        else:
            y[1][i] = np.nan

        #Cython
        spd_matrix = data.copy()
        t1 = time.time()
        out = cython_implementation(spd_matrix)
        t2 = time.time()
        #np.testing.assert_allclose(out, np_chol, rtol=1e-07)
        print "Cython implementation took {} seconds".format(t2-t1)
        y[2][i] = t2-t1

        if s <= 500:
            #Vectorized
            spd_matrix = data.copy()
            t1 = time.time()
            out = vectorizedchol(spd_matrix)
            t2 = time.time()
            #np.testing.assert_allclose(out, np_chol, rtol=1e-07)
            print "Vectorized implementation took {} seconds".format(t2-t1)
            y[3][i] = t2-t1
        else:
            y[3][i] = np.nan

        if s <= 500:
            #Mixed Naive / Cython
            spd_matrix = data.copy()
            t1 = time.time()
            out = mixedcython(spd_matrix)
            t2 = time.time()
            np.testing.assert_allclose(out, np_chol, rtol=1e-07)
            print "Mixed Cython/Python implementation took {} seconds".format(t2-t1)
            y[4][i] = t2 - t1
        else:
            y[4][i] = np.nan

        if s <=500:
            #new cholesky
            spd_matrix = data.copy()
            t1 = time.time()
            out = purepythoncholesky(spd_matrix)
            t2 = time.time()
            print "Pure Python Cholesky implementation took {} seconds".format(t2-t1)
            y[5][i] = t2-t1
        else:
            y[5][i] = np.nan

    #plot
    labels = ['NP', 'Naive', 'Cython', 'Vect', 'Mixed C/P']
    for i in range(y.shape[0]):
        plt.plot(x, y[i], label=labels[i])

    plt.show()
def prof():
    s = 500
    a = np.zeros(shape=(s,s), dtype=np.float64, order='F')
    a[:] = randint(0,10,size=(s, s))
    #Conversion using copy from int64 to float64
    a = a.astype(np.float64, copy=False, order='F')
    #Dot should ensure positive semi-definite
    data = np.dot(a, a.T)
    vectorizedchol(data)

def testvectorized():
    b = np.array([[9,9,6,-1],[9,30,5,22],[6,5,10,0],[-1,22,0,36.0]], order='F')
    print vectorizedchol(b)

def testblaschol():
    b = np.array([[9,9,6,-1],[9,30,5,22],[6,5,10,0],[-1,22,0,36.0]], order='F')
    blaschol(b)

if __name__ == '__main__':
    #prof()
    main()
    #testvectorized()
    #testblaschol()
