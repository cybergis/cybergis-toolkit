import pysal as ps
import numpy as np
import sys
import itertools
from mpi4py import MPI

#sys.excepthook = MPI.COMM_WORLD.Abort(1)

ns = [75]
rho = [0.9, 0.7, 0.5, 0.3, 0.1, 0, -0.1, -0.3, -0.5, -0.7, -0.9]
rho=[0.7, 0.9]
params = [x for x in itertools.product(ns, rho)]

def save_structured(n, rho):
    print "Generating..."
    outfile = 'Ident_{}_{}'.format(n, rho)
    ident = np.identity(n*n, dtype=np.float32)
    w = ps.weights.util.lat2W(nrows=n, ncols=n, rook=True, id_type='int')
    w.transform = 'r'
    inverse = np.linalg.inv((ident - rho * w.sparse))
    np.save(outfile, inverse)
    print 'Completed generating {}.'.format(outfile)

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
ncores = comm.Get_size()

if rank == 0:
    data = params
    chunks = [[] for i in range(ncores)]
    for i, chunk in enumerate(data):
        chunks[i%ncores].append(chunk)
else:
    data = None
    chunks = None

data = comm.scatter(chunks, root=0)

if rank != 0:
    for d in data:
        print d
        save_structured(d[0], d[1])

comm.Barrier()
