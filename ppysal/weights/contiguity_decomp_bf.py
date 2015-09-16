"""
Contiguity builder using brute force check for queen

start this from terminal with:
    ipcluster start -n 4

This is to explore parallelization prior to a binning approach

The idea here is to first subdivide the extent into ncp-1 regions, (where ncp
is number of compute units), find the shapes who have bounding boxes that are
in or overlap each region and then farm out a check for contiguity between the
pairs of shapes within each region (i.e., each core gets a region). then we
recombine after returning from the mapping.

Example run-times using 3 cores:

Sequential:  183.793218851
importing combinations from itertools on engine(s)
Parallel:  54.3720309734

"""
_author_ = "Serge Rey <sjsrey@gmail.com>"

from IPython.parallel import Client
from itertools import combinations

client = Client()

print client.ids

ncpus = len(client.ids)

import pysal as ps
import numpy as np

sf = ps.open(ps.examples.get_path("nat.shp"))
#sf = ps.open(ps.examples.get_path("columbus.shp"))
#sf = ps.open(ps.examples.get_path("sids2.shp"))
bb = sf.bbox


w = (bb[0] - bb[2]) / (ncpus-1)

w = np.abs(w)

bounds = np.arange(bb[0], bb[2]+w, w)[1:]
bounds[-1] += 0.0001*w
bins = {}
ids = {}
for b in range(len(bounds)):
    bins[b] = []
    ids[b] = []

shps = []

for i,shp in enumerate(sf):
    shps.append(shp)

    bbi = shp.bounding_box
    left = np.nonzero((bounds > bbi.left))[0][0]
    right = np.nonzero((bounds > bbi.right))[0][0]
    bids = range(left, right+1)
    for bid in bids:
        bins[bid].append(shp)
        ids[bid].append(i)



sf.close()

def bf_queen(shps, ids = []):
    """
    Brute force queen contiguity builder

    Arguments
    ---------

    shps: a list of pysal cg.Polygons

    ids: a list of integer ids for the shps


    Returns

    w: a dictionary with id as key and list of neighboring ids as values

    coords: a dictionary with (x,y) as the key and the value is a list of ids
    for shapes that have that (x,y) one or more times on their boundary

    """
    n = len(shps)
    w = {}
    coords = {}
    if not ids:
        ids = xrange(n)

    for i in range(n-1):
        si = shps[i]
        vertsi = si.vertices
        for vi in vertsi:
            if vi not in coords:
                coords[vi] = set([ids[i]])
            else:
                coords[vi] = coords[vi].union(set([ids[i]]))
        for j in range(i+1,n):
            sj = shps[j]
            vertsj = sj.vertices
            for vj in vertsj:
                if vj not in coords:
                    coords[vj] = set([ids[j]])
                else:
                    coords[vj] = coords[vj].union(set([ids[j]]))
    for coord in coords:
        if len(coords[coord]) > 1:
            pairs = combinations(coords[coord],2)
            for pair in pairs:
                l,r = pair
                if l not in w:
                    w[l] = set([r])
                else:
                    w[l] = w[l].union(set([r]))
                if r not in w:
                    w[r] = set([l])
                else:
                    w[r] = w[r].union(set([l]))
    return w, coords

def grid(bbox, n, delta=0.000001, dimension=10):
    bucket_min = dimension
    width = (( bbox[2] + delta) - bbox[0] ) / bucket_min
    height = (( bbox[3] + delta) - bbox[1] ) / bucket_min
    return  [width, height] * 2



if __name__ == '__main__':

    import time

    #t1 = time.time()
    #res, c = bf_queen(shps)
    #t2 = time.time()
    #
    #print 'Sequential: ',t2-t1
    #
    ## parallel
    #t1 = time.time()
    #view = client[0:-1]
    #with client[:].sync_imports():
    #    from itertools import combinations
    #results = view.map(bf_queen, bins.values(), ids.values())
    #t2 = time.time()
    #
    ## combine results
    #neighbors = {}
    #for i,result in enumerate(results.result):
    #    neigh,c = result
    #    for key in neigh:
    #        if key not in neighbors:
    #            neighbors[key] = neigh[key]
    #        else:
    #            neighbors[key] = neighbors[key].union(neigh[key])
    #
    #t3 = time.time()
    #print 'Parallel: ', t3-t1

    # vectorized parallel
    GRID_DIM = 80
    t1 = time.time()
    dims = grid(bb, GRID_DIM)
    bboxes = np.array([shp.bounding_box[:] for shp in shps])
    origin = np.array([bb[0],bb[1], bb[0], bb[1]])
    grid_cells = ((bboxes-origin) / dims).astype(int) # [r0,c0, r1,c1]
    grid_cells[:,[2,3]] += 1
    max_g = grid_cells.max(axis=1)
    n = len(shps)
    #global rows_2_polys
    #global cols_2_polys

    rows_2_polys = np.zeros((max_g[2],n), 'int')
    cols_2_polys = np.zeros((max_g[3],n), 'int')
    r_slices = grid_cells[:,[0,2]]
    c_slices = grid_cells[:,[1,3]]
    for i, row in enumerate(grid_cells):
        cols_2_polys[r_slices[i,0]:r_slices[i,1],i] = 1
        rows_2_polys[c_slices[i,0]:c_slices[i,1],i] = 1
    vertices = np.empty(len(bboxes), dtype=object)
    for i,shp in enumerate(shps):
        vertices[i] = np.array(shp.vertices)
    print type(vertices)
    #vertices = np.array([np.array(shp.vertices) for shp in shps])
    #vertices = [shp.vertices for shp in shps]
    import numpy

    def bb_check(i):
        r = numpy.dot(numpy.diag(rows_2_polys[:,i]), rows_2_polys).sum(axis=0)
        c = numpy.dot(numpy.diag(cols_2_polys[:,i]), cols_2_polys).sum(axis=0)
        potential_neighbors = numpy.nonzero((r>0) * (c>0))[0]
        # need to add: for neighbor in neighbors do explict check for shared vertex in shapes
        # potential neighbors have bounding boxes that share a commmon grid
        # row and a common grid column, but bounding boxes may not overlap
        bbox_i = bboxes[i]
        vertices_i = vertices[i]
        potential_neighbors  = potential_neighbors[potential_neighbors != i]
        neighbors = []
        for j in potential_neighbors:
            vertices_j = vertices[j]
            bbox_j = bboxes[j]
            if not ((bbox_j[2] < bbox_i[0]) or (bbox_j[0] > bbox_i[2])):
                if not ((bbox_j[3] < bbox_i[1]) or (bbox_j[1] > bbox_i[3])):
                    # bounding boxes overlap
                    com_v = [ vj for vj in vertices_j if vj in vertices_i ]
                    if len(com_v) > 0:
                        neighbors.append(j)
        return  neighbors

    view = client[0:-1]
    # put rows_2_poly and cols_2_polys in the namespaces of all the views
    # these are "small as in GRID_DIM x n  where GRID_DIM is the number of
    # rows in the grid == number of columns
    view['rows_2_polys'] = rows_2_polys
    view['cols_2_polys'] = cols_2_polys
    view['bboxes'] = bboxes
    view['vertices'] = vertices
    #view['shps'] = shps

    with client[:].sync_imports():
        import numpy
    neighbors = view.map_sync(bb_check, range(n))
    t2 = time.time()
    print 'Parallel vectorized: ', t2-t1

    #vectorized sequential

    t1 = time.time()
    dims = grid(bb, GRID_DIM)
    bboxes = np.array([shp.bounding_box[:] for shp in shps])
    origin = np.array([bb[0],bb[1], bb[0], bb[1]])
    grid_cells = ((bboxes-origin) / dims).astype(int) # [r0,c0, r1,c1]
    grid_cells[:,[2,3]] += 1
    max_g = grid_cells.max(axis=1)
    n = len(shps)
    rows_2_polys = np.zeros((max_g[2],n), 'int')
    cols_2_polys = np.zeros((max_g[3],n), 'int')
    r_slices = grid_cells[:,[0,2]]
    c_slices = grid_cells[:,[1,3]]
    for i, row in enumerate(grid_cells):
        cols_2_polys[r_slices[i,0]:r_slices[i,1],i] = 1
        rows_2_polys[c_slices[i,0]:c_slices[i,1],i] = 1


    neighbors = map(bb_check, range(n))
    t2 = time.time()
    print 'Serial vectorized: ', t2-t1


