import pysal as ps
import sys
import time as t

t1 = t.time()
w = ps.rook_from_shapefile(sys.argv[1])
t2 = t.time()
print "Serial time was {} seconds.".format(t2-t1)
print w[0]
print w[2499]
