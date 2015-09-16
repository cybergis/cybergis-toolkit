import numpy as np
from mpi4py import MPI

def square(data):
    '''
    Squares each value in the array
    '''
    return data ** 2

#Create the communicator object
comm = MPI.COMM_WORLD
#Get the IDs of the cores
rank = comm.Get_rank()
#Get the number of cores
ncores = comm.Get_size()

#Create the sample data in the parent process
if rank == 0:
    data = np.arange(100)
    #Compute the data segmentation
    load_per_core = len(data) / ncores
    #Simple example, simple decomposition
    data = data.reshape(ncores, load_per_core)
else:
    data = None

#Local data that we will gather back to
local_data = np.zeros(100)
data = comm.scatter(data, root=0)

#Child processes run this function
_data = square(data)

#Wait for all the workers to be done
comm.Barrier()

local_data = comm.gather(_data)
if rank == 0:
    print local_data
