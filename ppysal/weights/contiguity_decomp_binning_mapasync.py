import pysal as ps
import numpy as np
import multiprocessing as mp
import time
from decomp_binning import binning, binning_map

sf = ps.open(ps.examples.get_path("nat.shp"))
#sf = ps.open(ps.examples.get_path("columbus.shp"))
#sf = ps.open(ps.examples.get_path("sids2.shp"))
bb = sf.bbox

ncpus = mp.cpu_count()

w = (bb[0] - bb[2]) / (ncpus)

w = np.abs(w)

bounds = np.arange(bb[0], bb[2] + w, w)[1:]

bounds[-1] += 0.0001 * w

east = bounds
west = [bb[0]]
west.extend(east[0:-1])
bins = {}
ids = {}
BBS = {}
for b in range(len(bounds)):
    bins[b] = []
    ids[b] = []
    BBS[b] = [west[b], bb[1], east[b], bb[3]]

shps = []

for i, shp in enumerate(sf):
    shps.append(shp)

    bbi = shp.bounding_box
    left = np.nonzero((bounds > bbi.left))[0][0]
    right = np.nonzero((bounds > bbi.right))[0][0]
    bids = range(left, right + 1)
    for bid in bids:
        bins[bid].append(shp)
        ids[bid].append(i)

sf.close()

t1 = time.time()
res = binning(shps, bb)
t2 = time.time()

print 'Sequential: ', t2 - t1

############## Map_Async #####################
t1 = time.time()

"""
I think this is a potential point of optimization.
We prebin the data into left and right but then pass full extents. Using
the code below, we can custom pass extents and potential realize better
performance.
"""

# note that trying to pass in unique bbs/extents in the map raises errors. for
# now we use the same extent for each mapping but this is inefficient

#Start to setup to pass all these args as a single argument iterable for map
#Get the total number of geometries
n = len(shps)

#Compute the offsets
starts = range(0, n, n / ncpus)

ends = starts[1:]
ends.append(n)
offsets = [range(z[0], z[1]) for z in zip(starts, ends)]

#Package the args list
args = [(shps[offset[0]:offset[-1] + 1], bb, offset) for offset in offsets]

#Multiprocessing call
pool = mp.Pool(ncpus)
results = pool.map_async(binning_map, args)
results.wait()

# combine results
# this has to be fixed for correct mapping of IDS
neighbors = {}
for res_list in results.get():
    for key, result in res_list.iteritems():
        if key not in neighbors:
            neighbors[key] = result
        else:
            neighbors[key] = neighbors[key].union(result)

t2 = time.time()
print res[1539], neighbors[1539]
print res[1540], neighbors[1540]
print res[1541], neighbors[1541]
print res[1542], neighbors[1542]

print res[1000], neighbors[1000]
print len(res)
print len(neighbors)
print 'Map Async: ', t2 - t1
############## Apply_Async #####################
