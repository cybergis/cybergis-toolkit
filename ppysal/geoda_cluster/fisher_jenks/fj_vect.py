import numpy
import time
import ctypes
import warnings

#Suppress the divide by zero errors
warnings.filterwarnings('ignore', category=RuntimeWarning)

def fisher_jenks(values, classes=5, sort=True):
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
	
	varMat = numpy.zeros(shape=(numVal, numVal))
	
	for x in range(0,len(values)):
	    varMat[x] = values[:]
	    varMat[x][0:x] = 0
	
	errorMat = numpy.zeros(shape=(classes, numVal))
	
	return varMat, errorMat, numClass

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
    
	row = ((numpy.cumsum(numpy.square(arrRow))) - \
	        ((numpy.cumsum(arrRow)*numpy.cumsum(arrRow)) / (n)))
	row[numpy.isnan(row)] = 0 #Set the nan to 0
	return row
    
    if sort:
	values.sort()
	
    numVal = len(values)

    #Allocate all arrays
    varMat, errorMat, k = allocate(values, classes)

    t0 = time.time()

    #Calculate the variance matrix
    varMat = numpy.apply_along_axis(calcVar, 1, varMat, len(values))
    t1 = time.time()
    
    #Calculate the error matrix
    errorMat[0] = varMat[0]
    
    #Error Matrix Calculation
    row = 1
    for x in errorMat[1:]: #Get each row, save the first
	errRow = x[row:] #Work only on the indices which will hold valid values
	for y in range(0,len(errRow)): #Iterate through each index in errRow
	    #print errorMat[row-1][row-1:y+row] 
	    #print sharedVar.asarray()[:,y+row][row:y+row+1]
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
	
    return pivots

def _test():
    import doctest
    doctest.testmod(verbose=True)

#Not called by time_test.py
if __name__ == '__main__':    
    _test()
    
    

