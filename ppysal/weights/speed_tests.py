import sys, time
import pysal as ps
from pysal.weights._contW_binning import ContiguityWeights_binning
from pysal.weights._contW_binning import ContiguityWeightsLists
from pysal.weights._contW_rtree import ContiguityWeights_rtree

td = ['NAT.shp',
      'tract_shapes_pop_hsu_all.shp',
      'Merged__blockgroups_contUS_2010Census.shp']

methods = [ContiguityWeights_binning,
           ContiguityWeights_rtree,
           ContiguityWeightsLists]

for f in td:
    polys = ps.open(f)
    print 'Processing file {} which has {} polys.'.format(f, len(polys))
    print "Queen"
    for m in methods:
        t1 = time.time()
        w = m(polys,1)
        t2 = time.time()
        print "Method {} required {} seconds.".format(m, t2-t1)
    print "Rook"
    for m in methods:
        t1 = time.time()
        w = m(polys,2)
        t2 = time.time()
        print "Method {} required {} seconds.".format(m, t2-t1)

    exit()
