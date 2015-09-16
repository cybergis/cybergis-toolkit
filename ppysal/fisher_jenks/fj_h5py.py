import numpy
import time
import multiprocessing
import ctypes
import warnings
import h5py

#Simple edit

#Suppress the divide by zero errors
warnings.filterwarnings('ignore', category=RuntimeWarning)

#We can not use nested functions with map.
def calcVar(arrRow):
    '''This function calculates the diameter matrix.  It is called by fj.
    The return line performs the calculation.  All other lines prepare an
    order vector, prepended with zeros when necesary which stores n, the number 
    of elements summed for each index.'''

    lenN = (arrRow != 0).sum()
    n = numpy.arange(1, lenN+1)
    
    if lenN != numVal:
	n.resize(arrRow.shape[0]) 	
	n[arrRow.shape[0]-lenN:] =  n[:lenN-arrRow.shape[0]] 
	n[0:arrRow.shape[0]-lenN] = 0 

    arrRow =  ((numpy.cumsum(numpy.square(arrRow))) - \
            ((numpy.cumsum(arrRow)*numpy.cumsum(arrRow)) / (n)))
    arrRow[numpy.isnan(arrRow)] = 0
    return arrRow

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
	numpy types to a ctypes, shared memory array.'''
	
	numClass = classes
	numVal = len(values)
	
	varMatFile = h5py.File('/Users/jay/github/pPysal/fisher_jenks/varMatFile.hdf5', 'w') #Overwrites existing
	varMat = varMatFile.create_dataset('VarMat', (numVal,numVal), 'float64')
	
	for x in range(0,len(values)):
	    varMat[x] = values[:]
	    #Must use ellipses with this code, not arr[row][index]
	    if x != 0:
		varMat[x,...,0:x] = 0

	errorMat = numpy.zeros(shape=(classes, numVal))
	
	pivotShape = (classes, numVal)
	pivotMat = numpy.ndarray(pivotShape, dtype=numpy.float)
	
	#Initialize the arrays as globals.
	#initArr(varMat, errorMat)
	
	return pivotMat, numClass, varMatFile, errorMat, varMat
    
    def initErrRow(errRow_):
	'''Initialize the sharedErrRow as a global for multiprocessing'''
	global sharedErrRow
	sharedErrRow = errRow_

    def err(row,y,step, lenrow):
	'''This function computes the error on a segment of each error row, from the error matrix.  
	The function is provided with the row number, starting index, step size, and total row length.
	Since start + step could be greater than row length the first three lines check for this condition
	and set the stop variable to the appropriate value'''

	stop = (y+step)
	if stop+1 > lenrow:
	    stop = lenrow-1
	while y <= stop:
	    sharedErrRow[y] = numpy.amin(sharedErr[row-1][row-1:y+row] + varMatSlice[:,y+row][row:y+row+1])
	    y+=1   
	print sharedErrRow
	
    if sort:
	values.sort()
    
    if cores == None:
	cores = multiprocessing.cpu_count()
    
    #Hacky, I think I can use a class to make this variable global,
    # without polluting the namespace.
    global numVal 
    numVal = len(values)

    #Allocate all necessary memory
    pivotMat, k, varMatFile, errorMat, varMat = allocate(values, classes)
    
    #Calculate the number of cores over which to multiprocess
    if cores > len(values):
	cores = len(values)
    pool = multiprocessing.Pool(cores)
    
    t0 = time.time()
    jobs = []

    for x in range(0,numVal,cores):
	iterable = [row for row in varMat[x:x+cores]]
	results = pool.map(calcVar, iterable)
	index = x
	for result in results:
	    varMat[index] = result
	    index += 1
	#Manual cleanup
	del results

    t1 = time.time()

    #Calculate the error matrix
    errorMat[0] = varMat[0]
    print varMat
    print "Computing Error Matrix"
    #Error Matrix Calculation
    row = 1
    for x in errorMat[1:]: #Get each row, save the first
	errRow = x[row:] #Work only on the indices which will hold valid values
	for y in range(0,len(errRow)): #Iterate through each index in errRow
	    errRow[y] = numpy.amin(errorMat[row-1][row-1:y+row] + varMat[:,y+row][row:y+row+1])
	row += 1 #Iterate the row counter to get valid slices
    
    t2 = time.time()
    
    #Calculate Pivots
    pivots = []
    j = k - 1
    col = numVal - 1
    
    while j > 0:
	ev = errorMat[j, col]
	pivot_search = True
	right = col
	
	while pivot_search:
	    left_error = errorMat[j-1, right-1]
	    right_error = varMat[right, col]
	    if left_error + right_error == ev:
		pivots.insert(0, right)
		col = right -1
		pivot_search = False
	    right -= 1
	j -=1
    t3 = time.time()
	
    print "Time: ", t3-t1
    
    #Close the file internally, so that it definitely gets closed.
    varMatFile.close()
    
    return pivots

def _test():
    import doctest
    doctest.testmod(verbose=True)

#Not called by time_test.py
if __name__ == '__main__':
    #_test()
    #x = [12,10.8, 11, 10.8, 10.8, 10.8, 10.6, 10.8, 10.3, 10.3, 10.3,10.4, 10.5, 10.2, 10.0, 9.9]
    #p1 = fisher_jenks(x, 5)
    #print p1  

    
    x = []
    for y in range(500):
	x.append(15)
	x.append(11)
	x.append(2)
	x.append(26)
	x.append(191)
    p = fisher_jenks(x, 5)    
    print p
