from mpi4py import MPI 
from array import array 

comm = MPI.COMM_WORLD 
size = comm.size

for r in range(size):
    if comm.rank == r:
        print "I am rank {}".format(r)

if comm.rank == 0: 
    # process zero exposes an array of 10 integers 
    memory = array('i', range(10)) 
    disp_unit = memory.itemsize 
else: 
    # other process do not expose memory 
    memory = None 
    disp_unit = 1 
    print('Creating window')
    #win = MPI.Win.Create(memory, disp_unit, comm=comm) 
    print('Window created') 
    # all processes get three integers from process zero 
    buf = array('i', [0]*3) 
    print('Created buffer')
    #win.Fence() 
    #win.Get(buf, 0) 
    #win.Fence() 
    #print("[%d] %r" % (comm.rank, buf)) 
