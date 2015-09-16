
"""
Brute force contiguity builders used for testing results of parallel
algorithms
"""
_author_ = "Serge Rey <sjsrey@gmail.com>"


import pysal as ps
import numpy as np
from itertools import combinations

#sf = ps.open(ps.examples.get_path("nat.shp"))
#sf = ps.open(ps.examples.get_path("columbus.shp"))
#sf = ps.open(ps.examples.get_path("sids2.shp"))
#bb = sf.bbox


def bf_contiguity(shps, wttype = "QUEEN"):
    """
    Brute force contiguity builder

    Arguments
    ---------

    shps: list of pysal.cg.Polygon shapes

    wttype: string
            contiguity type

    Returns
    -------
    neighbors: dict
               key is id, value is list of neighboring ids

    """
    neighbors = {}
    if wttype.upper() == "QUEEN":
        vertices = {}
        for i, shp in enumerate(shps):
            si = set([i])
            for vertex in shp.vertices:
                if vertex not in vertices:
                    vertices[vertex] = set()
                vertices[vertex] = vertices[vertex].union(si)
        for vertex in vertices:
            pairs = combinations(vertices[vertex], 2)
            for l,r in pairs:
                if l not in neighbors:
                    neighbors[l] = set()
                if r not in neighbors:
                    neighbors[r] = set()
                neighbors[l] = neighbors[l].union([r])
                neighbors[r] = neighbors[r].union([l])
        return neighbors
    elif wttype.upper() == 'ROOK':
        edges = {}
        neighbors = {}
        for i, shp in enumerate(shps):
            neighbors[i] = set()
            nv = len(shp.vertices)
            for o in range(nv-1):
                d = o + 1
                edge = [shp.vertices[o], shp.vertices[d]]
                edge.sort()
                edge = tuple(edge)
                if edge not in edges:
                    edges[edge] = set()
                edges[edge] = edges[edge].union([i])
        checked = {}
        for edge in edges:
            pairs = combinations(edges[edge], 2)
            for pair in pairs:
                l,r = pair
                if pair not in checked and (r,l) not in checked:
                    neighbors[l] = neighbors[l].union([r])
                    neighbors[r] = neighbors[r].union([l])
                    checked[pair] = pair
                    checked[(r,l)] = (r,l)

        return neighbors
    else:
        print "Weight type not supported: ", wttype


def qf_shapefile(sf):
    shps = []
    f = ps.open(sf)
    for shp in f:
        shps.append(shp)
    f.close()
    return  bf_contiguity(shps, wttype = 'QUEEN')

def rf_shapefile(sf):
    shps = []
    f = ps.open(sf)
    for shp in f:
        shps.append(shp)
    f.close()
    return  bf_contiguity(shps, wttype = 'ROOK')


if __name__ == '__main__':

    sf = ps.examples.get_path("columbus.shp")
    queen_col = qf_shapefile(sf)
    rook_col = rf_shapefile(sf)
    wrc = ps.W(rook_col)
    print wrc.histogram

    import time
    sf = ps.examples.get_path("NAT.shp")
    t1 = time.time()
    queen = qf_shapefile(sf)
    wq = ps.W(queen)
    t2 = time.time()
    print "National queen: ", t2-t1
    sf = ps.examples.get_path("NAT.shp")
    t1 = time.time()
    rook = rf_shapefile(sf)
    wr = ps.W(rook)
    t2 = time.time()
    print "National rook: ", t2-t1

    t1 = time.time()
    wrps = ps.rook_from_shapefile(sf)
    t2 = time.time()
    print "PySAL rook: ", t2-t1

    t1 = time.time()
    wqps = ps.queen_from_shapefile(sf)
    t2 = time.time()
    print "PySAL queen: ", t2-t1
