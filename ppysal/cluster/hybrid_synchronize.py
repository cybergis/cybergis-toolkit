import multiprocessing as mp
import random

from mpi4py import MPI

#MPI Boilerplate
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
nmanagers = comm.Get_size()
status = MPI.Status()
host = MPI.Get_processor_name()

#Parallel function
def f(x):
    multiplier = random.randint(0,10)
    return x * multiplier

nlocalcores = mp.cpu_count() - 1  #One core is manager
localdata = [[1 for i in range(nlocalcores)] for i in range(nmanagers)]
data = [[1 for j in range(nlocalcores)] for i in range(nmanagers)]

vect = range(nlocalcores)
pool = mp.Pool()
localdata[rank] = pool.map(f, localdata[rank])
#update the collective 
data[rank] = localdata[rank]

#Collect the results on all managers
comm.Allgather(localdata, data)

for r in range(nmanagers):
    if r == rank:
        print "Rank {} gathered all data yielding: {}".format(rank, data)

