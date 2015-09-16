import pysal as ps

from pysal.weights.user import build_lattice_shapefile as BLS


for k in range(10, 110, 10):
    sfName = "1000x%d.shp"%k
    BLS(1000, k, sfName)



