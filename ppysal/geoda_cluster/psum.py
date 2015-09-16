<<<<<<< HEAD
'''
iPython example from:
http://ipython.org/ipython-doc/stable/parallel/parallel_mpi.html#actually-using-mpi
'''
from mpi4py import MPI
import numpy as np
import time
def psum(a):
    time.sleep(10)
    locsum = np.sum(a)
    rcvBuf = np.array(0.0,'d')
    MPI.COMM_WORLD.Allreduce([locsum, MPI.DOUBLE],
        [rcvBuf, MPI.DOUBLE],
        op=MPI.SUM)
    return rcvBuf
=======
from mpi4py import MPI
import numpy as np

def psum(a):
	locsum = np.sum(a)
	rcvBuf = np.array(0.0, 'd')
	MPI.COMM_WORLD.Allreduce([locsum, MPI.DOUBLE]),
		[rcvBuf, MPI.DOUBLE], op=MPI.SUM)
	return rcvBuf
>>>>>>> 1605d6ec0cab655469ab2ec6f029b5ef0b08663a
