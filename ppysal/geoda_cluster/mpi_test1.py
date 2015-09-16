'''
Embarassingly parallel example to see how mpi4py works.
Optimal speedup is Amdahl's Law, with p ~= 1.

This uses 3 cores - 2 workers and 1 receiver.  Adding additional
workers will not improve runtimes.

mpirun -np 3 python mpi_test1.py
'''

from mpi4py import MPI
import timeit

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

N = 1000000

def parSum():
	if rank == 0: 
	  s = sum(range(N/2))
	  comm.send(s,dest=2,tag=11)
	elif rank == 1:
	  s = sum(range(N/2+1,N))
	  comm.send(s,dest=2,tag=11)
	elif rank == 2:
	  s1 = comm.recv(source=0, tag=11)
	  s2 = comm.recv(source=1, tag=11)
	  print s1+s2

def serSum():
    s = sum(range(N))

if rank == 0:
    print 'Parallel time:'
    tp = timeit.Timer("parSum()","from __main__ import parSum")
    print tp.timeit(number=10)

    print 'Serial time:'
    ts = timeit.Timer("serSum()","from __main__ import serSum")
    print ts.timeit(number=10)
