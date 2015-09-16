import collections
import multiprocessing as mp
import time

import pysal


QUEEN = 1
ROOK = 2


class ContiguityWeightsLists:
    """
    Contiguity for a collection of polygons using high performance
    list, set, and dict containers
    """
    def __init__(self, collection, wttype=1):
        """
        Arguments
        =========

        collection: PySAL PolygonCollection

        wttype: int
                1: Queen
                2: Rook
        """
        self.collection = collection
        self.wttype = wttype
        self.jcontiguity()
    def jcontiguity(self):
        if self.collection.type != pysal.cg.Polygon:
            return False

        numPoly = len(self.collection)

        w = {}
        for i in range(numPoly):
            w[i] = set()

        geoms = []
        offsets = []
        c = 0  # PolyID Counter

        if self.wttype == QUEEN:
            for n in range(numPoly):
                    verts = self.collection.get(n).vertices
                    offsets += [c] * len(verts)
                    geoms += (verts)
                    c += 1

            items = collections.defaultdict(set)
            for i, vertex in enumerate(geoms):
                items[vertex].add(offsets[i])

            shared_vertices = []
            for item, location in items.iteritems():
                if len(location) > 1:
                    shared_vertices.append(location)

            for vert_set in shared_vertices:
                for v in vert_set:
                    w[v] = w[v] | vert_set
                    try:
                        w[v].remove(v)
                    except:
                        pass

        elif self.wttype == ROOK:
            for n in range(numPoly):
                verts = self.collection.get(n).vertices
                for v in range(len(verts) - 1):
                    geoms.append(tuple(sorted([verts[v], verts[v + 1]])))
                offsets += [c] * (len(verts) - 1)
                c += 1

            items = collections.defaultdict(set)
            for i, item in enumerate(geoms):
                items[item].add(offsets[i])

            shared_vertices = []
            for item, location in items.iteritems():
                if len(location) > 1:
                    shared_vertices.append(location)

            for vert_set in shared_vertices:
                for v in vert_set:
                    w[v] = w[v] | vert_set
                    try:
                        w[v].remove(v)
                    except:
                        pass
        self.w = w

def maprooklists(tpl):
    """
    shps:  a list of shps
    """
    collection = tpl[2]
    idx = tpl[1]
    c = tpl[0]  # Start counter

    geoms = []
    offsets = []
    t1 = time.time()
    for n in idx:
        verts = collection.get(n).vertices
        for v in range(len(verts) - 1):
            geoms.append(tuple(sorted([verts[v], verts[v + 1]])))
        offsets += [c] * (len(verts) - 1)
        c += 1
    t2 = time.time()
    print t2 - t1

    #Pack a dict {coordinate : set(geom{1}, geom{2}, ... geom_{n})}
    items = collections.defaultdict(set)
    for i, item in enumerate(geoms):
        items[item].add(offsets[i])

    return items

def chunks(shpidx, n, shp):
    for i in xrange(0, len(shpidx), n):
        yield (i, shpidx[i:i+n], shp)

def mergedicts(ldict):
    """
    ldict:   list of dicts to be merged
    """
    w = collections.defaultdict(set)
    items = collections.defaultdict(set)
    for d in ldict:
        for k, v in d.iteritems():
            items[k] = items[k] | v

    shared_vertices = []
    for item, location in items.iteritems():
        if len(location) > 1:
            shared_vertices.append(location)

    for vert_set in shared_vertices:
        for v in vert_set:
            w[v] = w[v] | vert_set
            try:
                w[v].remove(v)
            except:
                pass
    return w

def main():
    shpfile = '/home/jlaura/Dropbox/notebooks/pysal_weights/hexagonlattices/1024.shp'
    shp = pysal.open(shpfile)
    if shp.type != pysal.cg.Polygon:
        return False

    t1 = time.time()
    w = ContiguityWeightsLists(shp, QUEEN)
    t2 = time.time()
    print "Serial: {}".format(t2-t1)
    t1 = time.time()
    #Setup the pool
    ncores = 4
    pool = mp.Pool(processes=ncores)

    #Get the total length
    #nshps = len(pysal.open(shpfile))
    shpidx = range(len(shp))

    #Chunk the list
    partitioned_shps = list(chunks(shpidx, len(shp) / ncores, shp))
    res = pool.map(maprooklists, partitioned_shps)
    #w2 = mergedicts(res)
    t2 = time.time()
    print "Parallel: {}".format(t2 - t1)

    '''
    assert w.w == w2

    for k, v in w.w.iteritems():
        if v != w2[k]:
            print k, v, w2[k]
    '''
if __name__ == '__main__':
    main()
