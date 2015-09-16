import sys
import numpy as np
from mpi4py import MPI


axis_lookup = {'x':0, 'y':1}



rstate = np.random.RandomState(123456)

#Phase I: Initialize
def globalsort(comm, rank, shape, pts=None, axis='y', pivots=[], samplesize=8):
    #axis = axis_lookup[axis]
    #local_pts = np.empty((shape[0] / comm.size, shape[1]))
    #Phase II: Scatter the data, local quicksort, return regular sample
    #comm.Scatter([pts, MPI.DOUBLE], [local_pts, MPI.DOUBLE])
    #local_pts = local_pts[local_pts[:, axis].argsort()]

    #'''
    axis = axis_lookup[axis]

    #Compute the data breaks, assign remainder to rank 0.  Could be smarter and balance more
    quotient, remainder = divmod(shape[0], comm.size)
    if rank == 0:
        local_pts = np.empty((quotient + remainder, shape[1]))
    else:
        local_pts = np.empty((quotient, shape[1]))

    #Compute the scatter sizes and scatter offsets
    scattersize = list([(quotient + remainder) * shape[1]]) +\
                  [quotient * shape[1] for i in range(comm.size-1)]
    scatteroffsets = [0] + (np.cumsum(scattersize)[:-1]).tolist()
    
    #Phase II: Scatter the data, local quicksort, return regular sample
    comm.Scatterv([pts,scattersize, scatteroffsets, MPI.DOUBLE],
                  [local_pts, MPI.DOUBLE])
    local_pts = local_pts[local_pts[:, axis].argsort()]
    #'''
    comm.Barrier()
    
    if len(pivots) == 0:
        #local_pts = np.sort(local_pts, axis=0)
        idx = np.linspace(0, len(local_pts) - 1,
                        samplesize + 1,
                        dtype=int)[1:]
        samples = local_pts[idx]
        comm.Barrier()

        #Phase III: Gather, merge, and broadcast samples
        if rank == 0:
            sampledpts = np.empty(samplesize * comm.size * shape[1])
        else:
            sampledpts = None
        comm.Gather([samples.ravel(), MPI.DOUBLE],
                    [sampledpts,MPI.DOUBLE],
                    root=0)

        if rank == 0:
            sampledpts = sampledpts.reshape(-1, shape[1])
            sampledpts = sampledpts[sampledpts[:, axis].argsort()]
            idx = np.linspace(0,
                            len(sampledpts) - 1,
                            comm.size + 1,
                            dtype=np.int)[1:-1]
            pivots = sampledpts[idx][:,axis]
        else:
            pivots = np.empty(comm.size - 1, dtype=np.float64)
        comm.Barrier()
        #Ndarrays need to be 'flat' to communicate
        comm.Bcast([pivots.ravel(), MPI.DOUBLE])

    comm.Barrier()
    #Phase IV: Partition local data into p classes using pivots
    indices = local_pts[:,axis].searchsorted(pivots)
    partitions = np.split(local_pts, indices, axis=0)
    partitionsizes = np.array([i.shape[0] * shape[1] for i in partitions], dtype='i')
    comm.Barrier()

    '''
    for i in range(comm.size):
        if rank == 0:
            print pivots
        if rank == i:
            print rank, partitionsizes
    '''

    #Phase V: Gather the ith partition onto this core
    sizes = np.zeros(comm.size, dtype='i')
    for i in range(comm.size):
        comm.Gather([partitionsizes[i], MPI.INT32_T],
                    [sizes, MPI.INT32_T], root=i)
    pointcount = np.sum(sizes)

    #Create the memory space to collect the ith lists
    localclasses = np.empty(pointcount, dtype=np.float64)

    for i in range(comm.size):
        comm.Gatherv([partitions[i].ravel(), MPI.DOUBLE],
                    [localclasses, (sizes,None), MPI.DOUBLE],
                    root=i)
    comm.Barrier()

    localclasses = localclasses.reshape(-1, shape[1])
    localclasses = localclasses[localclasses[:, axis].argsort()]
    for i in range(comm.size):
        if rank == i:
            return localclasses
    #for i in range(comm.size):
        #if rank == i:
            #print i, localclasses


if __name__ == '__main__':
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    shape = (40,3)

    if rank == 0:
       #Generate an array of points
        pts = rstate.rand(shape[0], shape[1])
        pts[:,2] = np.arange(shape[0])
        l = globalsort(comm, rank, shape=shape, pts=pts, axis='y')
    else:
        pts = None
        l = globalsort(comm, rank, shape=shape, axis='y')

    comm.Barrier()

    for i in range(comm.size):
        if rank == i:
            print i, l
