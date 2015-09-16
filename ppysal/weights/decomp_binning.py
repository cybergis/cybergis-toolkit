"""
This script contains all of the binning algorithms.  They are identical in
functionality and are altered only to support the different argument passing
requirements for each type of multiprocessing call.

binning - The default, used by serial calls,iPython, apply_async
binning_map - map supports a single iterable, so the args are packed into
    a list of tuples.  The child then unpacks the iterable.
"""

import pysal as ps
import numpy as np

def binning(shps, bb, wttype="QUEEN", buckets = 8 , ids = []):
    """
    Contiguity builder for Queen (shared vertex)

    Arguments
    ---------
    shps: list of pysal.cg.Polygon objects

    bb: list of bounding box coordinates for the extent of the collection of
    shps

    wttype: string ("QUEEN" is default - only version working for now)

    buckets: int number of columns (rows) for the binning

    ids: list of ids associated with the shapes in shps

    Returns
    -------

    neighbors: dict - key is id, value is a list of neighboring ids

    """


    def bbcommon(bb, bbother):
        """
        Checks for overlaps of bounding boxes. First, east-west, then north-south.
        Element 0 is west, element 2 is east, element 3 is north, element 1 is
        south.
        All four checks must be false for chflag to be true, meaning the two
        bounding boxes do not overlap.
        """
        chflag = 0
        if not ((bbother[2] < bb[0]) or (bbother[0] > bb[2])):
            if not ((bbother[3] < bb[1]) or (bbother[1] > bb[3])):
                chflag = 1
        return chflag


    n_poly = len(shps)
    if not ids:
        ids = range(n_poly)
    DELTA = 0.0001

    bucketmin = n_poly / buckets + 2

    lengthx = ((bb[2] + DELTA - bb[0]) / bucketmin)
    lengthy = ((bb[3] + DELTA - bb[1]) / bucketmin)

    # initialize buckets
    columns = [ set() for i in range(bucketmin) ]
    rows = [ set() for i in range(bucketmin) ]


    minbox = bb[:2] * 2
    binWidth = [lengthx, lengthy] * 2

    bbcache = {}
    poly2Column = [ set() for i in range(n_poly)]
    poly2Row = [ set() for i in range(n_poly)]

    for i in range(n_poly):
        shp_i = shps[i]
        bbcache[i] = shp_i.bounding_box
        projBBox = [int((shp_i.bounding_box[:][j] -
                             minbox[j]) / binWidth[j]) for j in xrange(4)]
        for j in range(projBBox[0], projBBox[2] + 1):
            columns[j].add(i)
            poly2Column[i].add(j)
        for j in range(projBBox[1], projBBox[3] + 1):
            rows[j].add(i)
            poly2Row[i].add(j)
    w = {}
    if wttype.upper() == "QUEEN":
        # loop over polygons rather than bins
        vertCache = {}
        for polyId in xrange(n_poly):
            if polyId not in vertCache:
                vertCache[polyId] = set(shps[i].vertices)
            idRows = poly2Row[polyId]
            idCols = poly2Column[polyId]
            rowPotentialNeighbors = set()
            colPotentialNeighbors = set()
            for row in idRows:
                rowPotentialNeighbors = rowPotentialNeighbors.union(rows[row])
            for col in idCols:
                colPotentialNeighbors = colPotentialNeighbors.union(
                    columns[col])
            potentialNeighbors = rowPotentialNeighbors.intersection(
                colPotentialNeighbors)
            if polyId not in w:
                w[polyId] = set()
            for j in potentialNeighbors:
                if polyId < j:
                    if bbcommon(bbcache[polyId], bbcache[j]):
                        if j not in vertCache:
                            vertCache[j] = set(shps[j].vertices)
                        common = vertCache[polyId].intersection(vertCache[j])
                        if len(common) > 0:
                            w[polyId].add(j)
                            if j not in w:
                                w[j] = set()
                            w[j].add(polyId)
    elif wttype.upper() == "ROOK":
        # check for a shared edge
        edgeCache = {}
        # loop over polygons rather than bins
        for polyId in xrange(n_poly):
            if polyId not in edgeCache:
                iEdges = {}
                iVerts = shps[polyId].vertices
                nv = len(iVerts)
                ne = nv - 1
                for i in xrange(ne):
                    l = iVerts[i]
                    r = iVerts[i+1]
                    iEdges[(l,r)] = []
                    iEdges[(r,l)] = []
                edgeCache[polyId] = iEdges
            iEdgeSet = set(edgeCache[polyId].keys())
            idRows = poly2Row[polyId]
            idCols = poly2Column[polyId]
            rowPotentialNeighbors = set()
            colPotentialNeighbors = set()
            for row in idRows:
                 rowPotentialNeighbors = rowPotentialNeighbors.union(rows[row])
            for col in idCols:
                colPotentialNeighbors = colPotentialNeighbors.union(
                    columns[col])
            potentialNeighbors = rowPotentialNeighbors.intersection(
                colPotentialNeighbors)
            if polyId not in w:
                w[ids[polyId]] = set()
            for j in potentialNeighbors:
                if polyId < j:
                    if bbcommon(bbcache[polyId], bbcache[j]):
                        if j not in edgeCache:
                            jVerts = shps[j].vertices
                            jEdges = {}
                            nv = len(jVerts)
                            ne = nv - 1
                            for e in xrange(ne):
                                l = jVerts[e]
                                r = jVerts[e+1]
                                jEdges[(l,r)] = []
                                jEdges[(r,l)] = []
                            edgeCache[j] = jEdges
                        #for edge in edgeCache[j]:
                        if iEdgeSet.intersection(edgeCache[j].keys()):
                            w[ids[polyId]].add(ids[j])
                            if j not in w:
                                w[ids[j]] = set()
                            w[ids[j]].add(ids[polyId])
                            #break
    else:
        print "Unsupported weight type."

    neighbors = {}
    for key in w:
        if ids[key] not in neighbors:
            neighbors[ids[key]] = set( [ids[j] for j in w[key] ] )
        else:
            kid = ids[key]
            neighbors[kid] =  neighbors[kid].union(set( [ids[j] for j in w[key] ] ))
    return neighbors

def binning_map(tup_args):#shps, bb, wttype="QUEEN", buckets = 8 , ids = []):
    """
    Contiguity builder for Queen (shared vertex)

    Arguments
    ---------
    shps: list of pysal.cg.Polygon objects

    bb: list of bounding box coordinates for the extent of the collection of
    shps

    wttype: string ("QUEEN" is default - only version working for now)

    buckets: int number of columns (rows) for the binning

    ids: list of ids associated with the shapes in shps

    Returns
    -------

    neighbors: dict - key is id, value is a list of neighboring ids

    """
    bb = tup_args[1]
    wttype = "Queen"
    buckets = 8
    ids = tup_args[2]
    shps = tup_args[0]#[ids[0]:ids[-1]+1]


    def bbcommon(bb, bbother):
        """
        Checks for overlaps of bounding boxes. First, east-west, then north-south.
        Element 0 is west, element 2 is east, element 3 is north, element 1 is
        south.
        All four checks must be false for chflag to be true, meaning the two
        bounding boxes do not overlap.
        """
        chflag = 0
        if not ((bbother[2] < bb[0]) or (bbother[0] > bb[2])):
            if not ((bbother[3] < bb[1]) or (bbother[1] > bb[3])):
                chflag = 1
        return chflag


    n_poly = len(shps)
    if not ids:
        ids = range(n_poly)
    DELTA = 0.0001

    bucketmin = n_poly / buckets + 2

    lengthx = ((bb[2] + DELTA - bb[0]) / bucketmin)
    lengthy = ((bb[3] + DELTA - bb[1]) / bucketmin)

    # initialize buckets
    columns = [ set() for i in range(bucketmin) ]
    rows = [ set() for i in range(bucketmin) ]


    minbox = bb[:2] * 2
    binWidth = [lengthx, lengthy] * 2

    bbcache = {}
    poly2Column = [ set() for i in range(n_poly)]
    poly2Row = [ set() for i in range(n_poly)]

    for i in range(n_poly):
        shp_i = shps[i]
        bbcache[i] = shp_i.bounding_box
        projBBox = [int((shp_i.bounding_box[:][j] -
                             minbox[j]) / binWidth[j]) for j in xrange(4)]
        for j in range(projBBox[0], projBBox[2] + 1):
            columns[j].add(i)
            poly2Column[i].add(j)
        for j in range(projBBox[1], projBBox[3] + 1):
            rows[j].add(i)
            poly2Row[i].add(j)
    w = {}
    if wttype.upper() == "QUEEN":
        # loop over polygons rather than bins
        vertCache = {}
        for polyId in xrange(n_poly):
            if polyId not in vertCache:
                vertCache[polyId] = set(shps[i].vertices)
            idRows = poly2Row[polyId]
            idCols = poly2Column[polyId]
            rowPotentialNeighbors = set()
            colPotentialNeighbors = set()
            for row in idRows:
                rowPotentialNeighbors = rowPotentialNeighbors.union(rows[row])
            for col in idCols:
                colPotentialNeighbors = colPotentialNeighbors.union(
                    columns[col])
            potentialNeighbors = rowPotentialNeighbors.intersection(
                colPotentialNeighbors)
            if polyId not in w:
                w[polyId] = set()
            for j in potentialNeighbors:
                if polyId < j:
                    if bbcommon(bbcache[polyId], bbcache[j]):
                        if j not in vertCache:
                            vertCache[j] = set(shps[j].vertices)
                        common = vertCache[polyId].intersection(vertCache[j])
                        if len(common) > 0:
                            w[polyId].add(j)
                            if j not in w:
                                w[j] = set()
                            w[j].add(polyId)
    elif wttype.upper() == "ROOK":
        # check for a shared edge
        edgeCache = {}
        # loop over polygons rather than bins
        for polyId in xrange(n_poly):
            if polyId not in edgeCache:
                iEdges = {}
                iVerts = shps[polyId].vertices
                nv = len(iVerts)
                ne = nv - 1
                for i in xrange(ne):
                    l = iVerts[i]
                    r = iVerts[i+1]
                    iEdges[(l,r)] = []
                    iEdges[(r,l)] = []
                edgeCache[polyId] = iEdges
            iEdgeSet = set(edgeCache[polyId].keys())
            idRows = poly2Row[polyId]
            idCols = poly2Column[polyId]
            rowPotentialNeighbors = set()
            colPotentialNeighbors = set()
            for row in idRows:
                 rowPotentialNeighbors = rowPotentialNeighbors.union(rows[row])
            for col in idCols:
                colPotentialNeighbors = colPotentialNeighbors.union(
                    columns[col])
            potentialNeighbors = rowPotentialNeighbors.intersection(
                colPotentialNeighbors)
            if polyId not in w:
                w[ids[polyId]] = set()
            for j in potentialNeighbors:
                if polyId < j:
                    if bbcommon(bbcache[polyId], bbcache[j]):
                        if j not in edgeCache:
                            jVerts = shps[j].vertices
                            jEdges = {}
                            nv = len(jVerts)
                            ne = nv - 1
                            for e in xrange(ne):
                                l = jVerts[e]
                                r = jVerts[e+1]
                                jEdges[(l,r)] = []
                                jEdges[(r,l)] = []
                            edgeCache[j] = jEdges
                        #for edge in edgeCache[j]:
                        if iEdgeSet.intersection(edgeCache[j].keys()):
                            w[ids[polyId]].add(ids[j])
                            if j not in w:
                                w[ids[j]] = set()
                            w[ids[j]].add(ids[polyId])
                            #break
    else:
        print "Unsupported weight type."

    neighbors = {}
    for key in w:
        if ids[key] not in neighbors:
            neighbors[ids[key]] = set( [ids[j] for j in w[key] ] )
        else:
            kid = ids[key]
            neighbors[kid] =  neighbors[kid].union(set( [ids[j] for j in w[key] ] ))

    return neighbors
