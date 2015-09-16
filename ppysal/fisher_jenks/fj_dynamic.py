'''
This is for the exploration of dynamically computing
the required portions of the diameter matrix for the
dynamic population of the error matrix.

The idea is that the error matrix is only (nxk), so
scalability in the memory domain will be better than
storing the entire (nxn) distance matrix.

Next steps:
    1. Intelligent balancing.
    2. Per core caching of preceeding values to avoid
        redundant computation.
'''

import multiprocessing as mp
import ctypes
import numpy as np

#np.set_printoptions(linewidth = 200)
def allocate(values, classes):
    numvalues = len(values)

    errctype = mp.RawArray(ctypes.c_double,numvalues * classes)
    errormatrix = np.frombuffer(errctype)
    errormatrix.shape = (classes, numvalues)

    arrRow = np.array([values])
    n = np.arange(numvalues) + 1
    errormatrix[0] =  ((np.cumsum(np.square(arrRow))) - \
                ((np.cumsum(arrRow)*np.cumsum(arrRow)) / (n)))

    pivotmatrix = np.ndarray((1, classes), dtype=np.float)
    val = np.asarray(values)


    initarr(errormatrix, val)
    return pivotmatrix

def err(row, y, stop, lenrow):
    if stop > lenrow:
        sop = lenrow
    while y < stop:
        err = sharederror[row-1][row-1:y+row]
        diam = np.zeros(err.shape)
        #diam is a column from the diameter matrix.
        for x in range(len(diam)):
            obs = sharedval.copy()
            #Diagonally symmetric, prepend zeros
            obs[:row + x] = 0
            n = len(obs) - row + x
            #Window of values to compute
            arrRow = obs[row+x:row+y+1]
            diam[x] = ((np.sum(np.square(arrRow))) - \
                ((np.sum(arrRow)*np.sum(arrRow)) / (len(arrRow))))
        sharedrow[y] = np.amin(diam + err)
        y += 1

def initarr(errormatrix_, val_):
    global sharederror
    sharederror = errormatrix_
    global sharedval
    sharedval = val_

def init_error_row(errorrow_):
    global sharedrow
    sharedrow = errorrow_

def main(values, classes, sort=True):
    if sort:
        values.sort()

    pivotmatrix = allocate(values, classes)
    cores = mp.cpu_count()

    #Naive partitioning
    row = 1
    jobs = []
    #Step through rows of the error matrix
    for x in sharederror[row:]:
        errorrow = x[row:]
        init_error_row(errorrow)

        #Computations is the cumsum of all computations per row
        computations = np.cumsum(np.array([i+1 for i in range(len(errorrow))]))
        step = computations[-1] / cores
        #Every step computations I need to farm out to another core.
        starts = [0]
        stops = []
        comp_counter = 0
        for j in [i+1 for i in range(len(errorrow))]:
            if comp_counter + j < step:
                comp_counter += j
            else:
                stops.append(j)
                starts.append(j)
                comp_counter = 0
        stops.append(len(errorrow))
        offsets = zip(starts,stops)
        #step = len(errorrow) // cores

        for y in offsets:
            start = y[0]
            stop = y[1]
            p = mp.Process(target=err,args=(row,start, stop, len(errorrow)))
            jobs.append(p)

        for j in jobs:
            j.start()
        for j in jobs:
            j.join()
        del jobs[:], p, j
        row += 1
    print sharederror
if __name__ == "__main__":
    values = [12,10.8, 11, 10.8, 10.8, 10.8, 10.6, 10.8, 10.3, 10.3, 10.3,10.4, 10.5, 10.2, 10.0, 9.9]
    main(values, 5)
