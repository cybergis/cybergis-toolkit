import multiprocessing as mp
from mpi4py import MPI

#MPI Boilerplate
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
ncores = comm.Get_size()
status = MPI.Status()
host = MPI.Get_processor_name()

nlocalcores = mp.cpu_count()

def f(pid,cid):
    print "I am worker {} managed by mpi process {}".format(pid, cid)

for r in range(ncores):
    if r == rank:
        print "I am a manager of rank {} with {} available shared memory cores.".format(rank, nlocalcores)

fakejobs = []
for i in range(nlocalcores):
    p = mp.Process(target=f, args=(i,rank))
    fakejobs.append(p)
    p.start()
