import numpy
import time
import multiprocessing
import ctypes
import warnings

#Suppress the divide by zero errors
warnings.filterwarnings('ignore', category=RuntimeWarning)

def fisher_jenks(values, classes=5, cores=None, sort=True):
    '''Fisher-Jenks Optimal Partitioning of an ordered array into k classes

    Parameters
    ==========

    x: list of values to be partitioned
    k: number of classes to form
    
    Returns
    =======

    pivots: list (k-1) 
            since x order has to be respected, the first cluster begins 
            at element 0. The first value in pivots is the start of the second 
            class and the last value the start of the last class.

    Example
    =======
    >>> x = []
    >>> for y in range(500):
    ...     x.append(15)
    ...     x.append(11)
    ...     x.append(2)
    ...     x.append(26)
    ...     x.append(191)
    >>> p = fisher_jenks(x, 5)
    >>> p
    [500, 1000, 1500, 2000]
    >>> x = [12,10.8, 11, 10.8, 10.8, 10.8, 10.6, 10.8, 10.3, 10.3, 10.3,10.4, 10.5, 10.2, 10.0, 9.9]
    >>> p1 = fisher_jenks(x, 5)
    >>> p1
    [2, 7, 9, 15]
    '''

    def allocate(values, classes):
	'''This function allocates memory for the variance matrix, error matrix, 
	and pivot matrix.  It also moves the variance matrix and error matrix from
	numpy types to a ctypes, shared memory array.
	
	In this version the diameter and error matrices are allocated to size 
	n*k.  This way we are not storing the entire n^2 diameter matrix in memory.
	'''
	
	numClass = classes
	numVal = len(values)
	
	varCtypes = multiprocessing.RawArray(ctypes.c_double, classes*numVal)
	varMat = numpy.frombuffer(varCtypes)
	varMat.shape = (classes,numVal)
	
	print "Shape ", len(varMat[0])
	
	for x in range(0,len(values)):
	    varMat[x] = values[:]
	    varMat[x][0:x] = 0
    
	errCtypes = multiprocessing.RawArray(ctypes.c_double, classes*numVal)
	errorMat = numpy.frombuffer(errCtypes)
	errorMat.shape = (classes, numVal)
	
	pivotShape = (classes, numVal)
	pivotMat = numpy.ndarray(pivotShape, dtype=numpy.float)
	
	#Initialize the arrays as globals.
	initArr(varMat, errorMat)
	
	return pivotMat, numClass
    
    def initArr(varMat_, errorMat_):
	'''Initialize the ctypes arrays as global variables for multiprocessing'''
	global sharedVar
	sharedVar = varMat_
	
	global sharedErr
	sharedErr = errorMat_
    
    def initErrRow(errRow_):
	'''Initialize the sharedErrRow as a global for multiprocessing'''
	global sharedErrRow
	sharedErrRow = errRow_
    
    def fj(sharedVar,i, values, start):
	'''This function facilitates passing multiple rows to each process and
	then performing multiple vector calculations along individual rows.'''
	arr = sharedVar
	arr[i] = numpy.apply_along_axis(calcVar, 1, arr[i], len(values))
	arr[i][numpy.isnan(arr[i])] = 0
    
    def calcVar(arrRow, lenValues):
	'''This function calculates the diameter matrix.  It is called by fj.
	The return line performs the calculation.  All other lines prepare an
	order vector, prepended with zeros when necesary which stores n, the number 
	of elements summed for each index.'''
	
	lenN = (arrRow != 0).sum()
	n = numpy.arange(1, lenN+1)
	
	if lenN != lenValues:
	    n.resize(arrRow.shape[0]) 	
	    n[arrRow.shape[0]-lenN:] =  n[:lenN-arrRow.shape[0]] 
	    n[0:arrRow.shape[0]-lenN] = 0 
    
	return ((numpy.cumsum(numpy.square(arrRow))) - \
	        ((numpy.cumsum(arrRow)*numpy.cumsum(arrRow)) / (n)))
    
    def err(row,y,step, lenrow):
	'''This function computes the error on a segment of each error row, from the error matrix.  
	The function is provided with the row number, starting index, step size, and total row length.
	Since start + step could be greater than row length the first three lines check for this condition
	and set the stop variable to the appropriate value'''
    
	stop = (y+step)
	if stop+1 > lenrow:
	    stop = lenrow-1
	while y <= stop:
	    sharedErrRow[y] = numpy.amin(sharedErr[row-1][row-1:y+row] + sharedVar[:,y+row][row:y+row+1])
	    y+=1   
    
    if sort:
	values.sort()
    
    if cores == None:
	cores = multiprocessing.cpu_count()
	
    numVal = len(values)

    #Allocate all necessary memory
    pivotMat, k = allocate(values, classes)
    
    #Calculate the number of cores over which to multiprocess
    if cores > len(values):
	cores = len(values)
    step = numVal // cores

    t0 = time.time()
    jobs = []
    
    #Calculate the variance matrix
    for i in range(0,len(values),step):
	p = multiprocessing.Process(target=fj,args=(sharedVar,slice(i, i+step),values, i))
	jobs.append(p) 
    for job in jobs:
	job.start()
    for job in jobs:
	job.join()
    del jobs[:], p, job
    t1 = time.time()
    
    #Calculate the error matrix
    sharedErr[0] = sharedVar[0]
    
    #Error Matrix Calculation
    row = 1    
    for x in sharedErr[1:]:
	errRow = x[row:]
	#Initialize the errorRow as a global for multiprocessing writes
	initErrRow(errRow)
	
	#Calculate the step size for the errRow to be processed
	step = len(errRow) // cores
	
	for y in range(0,len(errRow), step+1):
	    p = multiprocessing.Process(target=err, args=(row,y, step, len(errRow)))
	    jobs.append(p)
	
	for job in jobs:
	    job.start()
	for job in jobs:
	    job.join()
	del jobs[:], p, job
	 
	row += 1
    
    t2 = time.time()
    
    #Calculate Pivots
    pivots = []
    j = k - 1
    col = numVal - 1
    
    while j > 0:
	ev = sharedErr[j, col]
	pivot_search = True
	right = col
	
	while pivot_search:
	    left_error = sharedErr[j-1, right-1]
	    right_error = sharedVar[right, col]
	    if left_error + right_error == ev:
		pivots.insert(0, right)
		col = right -1
		pivot_search = False
	    right -= 1
	j -=1
    t3 = time.time()
	

    return pivots

def _test():
    import doctest
    doctest.testmod(verbose=True)

#Not called by time_test.py
if __name__ == '__main__':
    x = []
    for y in range(500):
	x.append(15)
	x.append(11)
	x.append(2)
	x.append(26)
	x.append(191)
    p = fisher_jenks(x, 5)
    print p   
#_test()
    
    