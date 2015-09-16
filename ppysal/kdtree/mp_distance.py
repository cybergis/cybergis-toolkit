import numpy as np
import pysal as ps
import time
import multiprocessing as mp

import kd_parallel


def g(tup):
    '''This function is called by a pool of
    workers to prep the neighbors data.'''
    lines = tup[0]
    idset = tup[1]
    offset = tup[1][0]
    neighbors = {}
    for i, row in enumerate(lines):
        row = row.tolist()
        row.remove(i + offset)
        neighbors[idset[i]] = list(row)
    return neighbors


def main():
    cores = mp.cpu_count()
    pool = mp.Pool(processes=cores)

    #sizes = [ 10**i for i in range(1,5) ]
    #data = np.random.random_integers(0,100,(sizes[-1],2))
    #for size in sizes:
        #print size
        #t1 = time.time()
        #wnn2 = ps.weights.Distance.knnW(data,2)
        #t2 = time.time()
        #print size, t2-t1

    ##sizes = [ 10**i for i in range(4,6) ]
    ##data = np.random.random_integers(0,100,(sizes[-1],2))
    ##for size in sizes:
        ##print size
        ##t1 = time.time()
        ##wnn2 = ps.weights.Distance.knnW(data,2)
        ##t2 = time.time()
        ##print size, t2-t1

    #print data.shape

    #kd = ps.common.KDTree(data)
    #nnq = kd.query(data,k=2+1, p=2)
    #info = nnq[1]
    ##Here we start the multiprocessing work
    #t3 = time.time()
    #slice_size = len(info) / cores
    #id_set = np.arange(len(info))

    ##Pack up an arg list
    #sections = []
    #for line in xrange(0,len(info),slice_size):
        #sections.append((info[line:line+slice_size], id_set[line:line+slice_size], line))

    ##print sections[0] #Just to show the tuple that is going passed

    #result_pool = pool.map(g, iterable=sections)

    ##Pack the results back into 1 dict
    #result = {}
    #map(result.update, result_pool)
    #t4 = time.time()
    #print size, t4-t3
    ##Pack the W
    #w = ps.weights.W(result)
    #print w

    #'''SERIAL'''
    #sizes = [ 10**i for i in range(2,3) ]
    #data = np.random.random_integers(0,100,(sizes[-1],2))
    #for size in sizes:
        #kd = ps.common.KDTree(data)
        #nnq = kd.query(data,k=2+1, p=2)
        #info = nnq[1] #This is the indices of the neighbors
        #neighbors = {}
        #idset = np.arange(len(info)) #Indices of the input point
        #for i, row in enumerate(info):
            #row = row.tolist()
            #row.remove(i)
            #neighbors[idset[i]] = list(row)
    #print neighbors

    '''SERIAL TO TEST PARALLEL KDQUERY'''
    #sizes = [ 10**i for i in range(3,4)]
    #data = np.random.random_integers(0,100,(sizes[-1],2))
    #for size in sizes:
        #kd = kd_parallel.KDTree(data)
        #nnq = kd.query(data,k=2+1, p=2)
        #info = nnq[1] #This is the indices of the neighbors
        #neighbors = {}
        #idset = np.arange(len(info)) #Indices of the input point
        #for i, row in enumerate(info):
            #row = row.tolist()
            #row.remove(i)
            #neighbors[idset[i]] = list(row)
    #print neighbors

    '''PARALLEL SINGLE STATEMENT w/ PARALLEL KDQUERY'''
    sizes = [10**i for i in range(1,5)]
    for size in sizes:
        data = np.random.random_integers(0, 10000, (size, 2))
        t1 = time.time()
        kd = kd_parallel.KDTree(data)
        t2 = time.time()
        nnq = kd.query(data, k=2+1, p=2)
        t3 = time.time()
        info = nnq[1] #This is the indices of the neighbors
        info = info[info[:, 0].argsort()]
        slice_size = len(info) // cores
        id_set = np.arange(len(info))
        sections = []
        for line in xrange(0,len(info),slice_size):
            sections.append((info[line:line+slice_size], id_set[line:line+slice_size]))
        result_pool = pool.map(g, iterable=sections)
        result = {}
        map(result.update, result_pool)
        w = ps.weights.W(result)
        t4 = time.time()
        print size, " KD Tree || time: {0}".format(t2-t1)," KD Query || time: {0}".format(t3-t2), " W prep time (multi-core): {0}".format(t4-t3), " Total time: {0}".format(t4-t1)

        t1 = time.time()
        kd = ps.common.KDTree(data)
        t2 = time.time()
        nnq = kd.query(data,k=2+1, p=2)
        t3 = time.time()
        info = nnq[1] #This is the indices of the neighbors
        slice_size = len(info) / cores
        id_set = np.arange(len(info))
        sections = []

        for line in xrange(0,len(info),slice_size):
            sections.append((info[line:line+slice_size], id_set[line:line+slice_size]))
        result_pool = pool.map(g, iterable=sections)
        result = {}
        map(result.update, result_pool)
        w = ps.weights.W(result)
        t4 = time.time()
        print size, " KD Tree Serial time: {0}".format(t2-t1)," KD Query Serial time: {0}".format(t3-t2), " W prep time (multi-core): {0}".format(t4-t3), " Total time: {0}".format(t4-t1)
        #print w.neighbors

    '''PARALLEL SINGLE STATEMENT'''
    #sizes = [ 10**i for i in range(2,3) ]
    #print sizes
    #data = np.random.random_integers(0,100,(sizes[-1],2))
    #print "Data generated: ", data
    #for size in sizes:
        #print size
        #t1 = time.time()
        #kd = ps.common.KDTree(data)
        #t2 = time.time()
        #nnq = kd.query(data,k=2+1, p=2)
        #t3 = time.time()
        #info = nnq[1] #This is the indices of the neighbors
        #print "Length of info: ", len(info)
        #slice_size = len(info) / cores
        #print "Slice Size: ", slice_size
        #id_set = np.arange(len(info))
        #sections = []
        #print info
        #for line in xrange(0,len(info),slice_size):
            #sections.append((info[line:line+slice_size], id_set[line:line+slice_size]))
        #print sections[1]
        #exit()
        #result_pool = pool.map(g, iterable=sections)
        #result = {}
        #map(result.update, result_pool)
        #w = ps.weights.W(result)
        #t4 = time.time()
        #print size, " KD Tree time: {}".format(t2-t1)," KD Query time: {}".format(t3-t2), " W prep time (multi-core): {}".format(t4-t3), " Total time: {}".format(t4-t1)
        #print w.neighbors[1]

if __name__ == "__main__":
    main()

