import sys
import time
import numpy
from fj_refactored import fisher_jenks

cores = [1,2,4,16,32]
classes = [5,6,7]
data_sizes = [500, 1000, 2500, 5000, 7500, 10000]

for c in cores:
    for d in data_sizes:
        for k in classes:
            data = numpy.random.ranf(size=d)
            try:
                t1 = time.time()
                #wrapped in try since we will blow out RAM at some point
                classification = fisher_jenks(data, k, c)
                t2 = time.time()
                print "Processed {0} data points in {1} classes using {2} cores. Total time: {3}".format(d, k, c, t2-t1)
            except KeyboardInterrupt:
                print "Aborting"
                sys.exit(1)
            except:
                print "FAILURE: {0} data points.".format(d)
