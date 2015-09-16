import sys
import time
import numpy as np
import pysal as ps

from pysal.esda.mapclassify import Fisher_Jenks, Fisher_Jenks_Sampled
from fj_refactored import fisher_jenks as pFisher
from fj_vect import fisher_jenks as vFisher

from mpi4py import MPI
#from mpi4py import ANY_SOURCE


def gen_structure_data(n, rho):
    # eps is not spatially correlated
    ident = np.identity(n*n) # identity matrix
    eps = np.random.random(n*n) #uniformly distributed error vector
    eps.shape = (n*n, 1)  #convert to vector
    w = ps.weights.util.lat2W(nrows = n, ncols = n, rook = True, id_type = 'int')  #Lattice
    w.transform = 'r' #Row standardize
    y = np.array((np.matrix(ident-rho*w.sparse).I) * np.matrix(eps))
    #Uncomment to compute Moran's I and check that the data is generated as expected.
    #morI = ps.Moran(np.array(y), w)
    #print "Sample Data: Rho: {} yielded Moran's I: {} with p-value {}".format(rho, morI.I, morI.p_sim)
    return np.array(y)

def get_tss(k, y, pivots):
    """
    Total sum of squares around class means

    Returns sum of squares over all class means
    """
    tss = 0
    y = np.sort(y)
    classes = np.arange(k)
    starts = [0]
    for i in range(len(pivots)):
        starts.append(pivots[i])
    pivots.append(len(y))
    breaks = zip(starts, pivots)
    for b in breaks:
        ymean = y[b[0]:b[1]].mean()
        css = y[b[0]:b[1]] - ymean
        css *= css
        tss += np.sum(css)
    return tss

def fj_generate_sample(y, pct=0.10, random=True):
    n = y.size
    if random:
        ids = np.random.random_integers(0, n - 1, n * pct)
    else:
        ids = np.arange(int(n*pct))
    yr = y[ids]
    yr[-1] = max(y)  # make sure we have the upper bound
    yr[0] = min(y)  # make sure we have the min
    return yr

def fj(realizations, n, p, k, s):
    local_log = np.zeros((realizations, 11))
    #Realizations
    for r in range(realizations):
	orig_data = gen_structure_data(n,p).flatten()
        variance = np.var(orig_data)
	data = orig_data.copy()
	#np.save('data.npy', data)

	t1 = time.time()
	pivots = vFisher(data, classes=k)
	t2 = time.time()
	fj_tss = get_tss(k, data, pivots)

	#Cleanup
	del data
	data = orig_data.copy()
	#data = np.load('data.npy')

	####Ordered Sample####
	t3 = time.time()
	sample = fj_generate_sample(data, pct=s, random=False)
	#Returns pivots in the sample, not pivots in the total - hence C-TSS is way off.
	fj_ordered = vFisher(sample, classes=k)

	t4 = time.time()
	data.sort()
	pivots = []
	for cl in range(0,k-1):
            pivots.append(np.where(data==sample[fj_ordered[cl]])[0][-1])
	fj_ordered_tss = get_tss(k, data, pivots)

	#Cleanup
	del data
	data = orig_data.copy()
	#data = np.load('data.npy')

	####Random Sample####
	t5 = time.time()
	sample = fj_generate_sample(data, pct=s)
	fj_random = vFisher(sample, classes=k)
	t6 = time.time()

	data.sort()
	pivots = []
	for cl in range(0,k-1):
	    pivots.append(np.where(data==sample[fj_random[cl]])[0][-1])
	fj_random_tss = get_tss(k, data, pivots)
	log = [n,k,p, s,variance, fj_tss, fj_random_tss, fj_ordered_tss, t2-t1, t4-t3, t6-t5]
	log = map(str, log)
	local_log[r] = np.array(log)
    return local_log

#MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
ncores = comm.Get_size()

HEADER = 'N,K,RHO,S,VARIANCE,FJ_TSS,RANDOM_TSS,ORDERED_TSS,FJ_TIME,RANDOM_TIME,ORDERED_TIME'

#Parameter space to explore
rho = [0.9, 0.7, 0.3, 0.1, 0, -0.1, -0.3, -0.7, -0.9]
ns = [25,50,75,100]#,125,150,175,200]#,225,250,275,300] #Large need to use HDF5 for full
ks = [5,7,9]
ss = [0.05, 0.10, 0.15, 0.20, 0.25]
REALIZATIONS = int(sys.argv[1])

#Compute the number of computations per core
comps = REALIZATIONS / ncores
#N
for i,n in enumerate(ns):
    #Rho
    for r,p in enumerate(rho):
	#K
	for j,k in enumerate(ks):
	    #Sample Size
	    for l,s in enumerate(ss):
		#Generate the data to pass to the children
		logs = np.zeros((REALIZATIONS, 11))
		if rank == 0:
		    data = np.array([comps, n, p, k, s])
		else:
		    data = None
		
		#Broadcast the data array to all children
		data = comm.bcast([data, MPI.DOUBLE], root=0)
		
		_logs = fj(int(data[0][0]), int(data[0][1]), data[0][2], int(data[0][3]), data[0][4])
		#comm.Barrier() #Synchronize all processes
	        comm.Gather([_logs, MPI.DOUBLE],[logs, MPI.DOUBLE],  root=0)
		if rank == 0:
		    with open('log.txt', 'a') as logfile:
			np.savetxt(logfile, logs)
		    print "Completed {} realizations with n={}, k={}, rho={}, and sample size={}".format(REALIZATIONS, n,k,p,s)	    
		    
'''
with open('log.txt', 'a') as log_file:
    log_file.write(HEADER + '\n')

    if comm.rank == 0:
	for l in logs:
           for e in l:
               log_file.write(e + '\n')
'''	
