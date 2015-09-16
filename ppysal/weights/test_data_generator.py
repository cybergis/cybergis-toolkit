import pysal as ps
import numpy as np
import scipy.spatial as ss
from random import shuffle
import itertools

#Build the lattices and save as a shapefile

lattices = [32,100,224,317]
for lattice in lattices:
    ps.weights.build_lattice_shapefile(lattice,lattice,"{}.shp".format(lattice**2))
    
polys = [1000,10000,50000,100000]

for poly in polys:
    x = range(poly)
    y = range(poly)
    
    coords = [x,y]
    
    for element in itertools.product(*coords):
        print element
    
    