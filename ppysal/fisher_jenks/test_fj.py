import math
import sys
import time
import numpy as np
from fj_refactored import fisher_jenks, fj_generate_sample


def testfull():
    """
    Tests the fully enumerated Fisher-Jenks implementation
    """
    cores = [1,2,4,16,32]
    classes = [5,6,7]
    data_sizes = [500, 1000, 2500, 5000, 7500, 10000, 12500, 15000, 17500, 20000, 22500, 25000]

    for c in cores:
        for d in data_sizes:
            for k in classes:
                data = np.random.ranf(size=d)
                try:
                    t1 = time.time()
                    #wrapped in try since we will blow out RAM at some point
                    classification = fisher_jenks(data, k, c)
                    t2 = time.time()
                    print "Processed {0} data points in {1} classes using {2} cores. Total time: {3}".format(d, k, c, t2-t1)
                    data = None
                except KeyboardInterrupt:
                    print "Aborting"
                    sys.exit(1)
                except:
                    print "FAILURE: {0} data points.".format(d)


def testsample():
    """
    Tests the sampled Fisher-Jenks implementation
    """
    cores = [1,2,4,16,32]
    classes = [5,6,7]
    data_sizes = [10000, 20000, 40000, 80000, 160000, 320000, 640000, 1280000,
                  2560000, 5120000, 10240000, 20480000, 40960000, 81920000,
                  163840000, 327680000, 655360000]
    for c in cores:
        for d in data_sizes:
            for k in classes:
                #Generate the test data and save to disk
                data = np.random.ranf(size=d)
                nobs = len(data)
                np.save('testfile.npy', data)
                data = None

                #Compute the sample size as the sqrt of nobs
                sqrt = math.sqrt(nobs)
                if sqrt > 40000:
                    sqrt = 40000
                pct = sqrt / float(d)

                #Load the data back into memory as a mmapped file
                f = np.load('testfile.npy', mmap_mode='r+')
                t1 = time.time()
                data = fj_generate_sample(f, pct=pct)
                t2 = time.time()
                print "Randomly sampling {0} percent of {1} observations for a run size of {2} observations took {3} seconds.".format(pct, nobs, sqrt, t2 - t1)
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


if __name__ =='__main__':
    #Test the fully enumerated FJ
    testfull()

    #Test FJ using sampling
    testsample()
