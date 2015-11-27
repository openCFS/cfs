#!/usr/bin/python

# import modules
import sys
sys.path.append('/home/daniel/code/cfs/debug/cfsdeps/sgpp/src/sgpp/sgopt/lib')
from pysgpp import Grid


dim = int(sys.argv[1])
level = int(sys.argv[2])

# create a two-dimensional piecewise bi-linear grid
bsplineDegree = 3
grid = Grid.createModBsplineGrid(dim,bsplineDegree)
#grid = Grid.createLinearBoundaryGrid(dim)
gridStorage = grid.getStorage()

# create regular grid
gridGen = grid.createGridGenerator()
gridGen.regular(level)
print "number of grid points:  %d, %i" % (gridStorage.size(), 16 ** 3)

fd = open('grid_points.csv', 'w')
fd.write("\n".join([gridStorage.get(i).toString() for i in xrange(gridStorage.size())]))
fd.close()
