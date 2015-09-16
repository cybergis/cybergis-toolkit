import os
name = os.uname()[1]

try:
    from mpi4py import MPI

    t0 = MPI.Wtime()

    import pysal as ps
    import numpy as np
    import os
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    ncores = comm.Get_size()
   
    #Check that we have all cores available
    assert(ncores == 8)
    #Check numpy
    assert(np.__version__ == '1.7.1')
    #Check PySAL
    assert(ps.version == '1.6.0')

    t1 = MPI.Wtime()

    time = t1 - t0
    times = comm.gather(time, root=0)

    if rank == 0:
        time = max(times)
        with open('test1.log', 'a') as f:
            f.write("Total testing time for {} using {} cores was {} seconds.\n".format(name, ncores, t1-t0))
except:
    with open('test1.log', 'a') as f:
        f.write("{} failed.".format(name))

