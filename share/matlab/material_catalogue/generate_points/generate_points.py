#!/usr/bin/python

# import modules
from __future__ import print_function
import sys
#sys.path.append('/home/daniel/code/sgopt_10.0.0/lib')
import pysgpp
#from pysgpp import Grid


dim = int(sys.argv[1])
level = int(sys.argv[2])

# create a two-dimensional piecewise bi-linear grid
bsplineDegree = 3
grid = pysgpp.Grid.createModBsplineGrid(dim,bsplineDegree)
#grid = Grid.createLinearBoundaryGrid(dim)
gridStorage = grid.getStorage()

# create regular grid
gridGen = grid.createGridGenerator()
gridGen.regular(level)
print("number of grid points:  %d, %i" % (gridStorage.size(), 16 ** 3))

for i in xrange(gridStorage.size()):
  point = gridStorage.get(i)
  print(point.getCoord(0),point.getCoord(1),point.getCoord(2))

fd = open('grid_points.csv', 'w')
fd.write("\n".join([gridStorage.get(i).toString() for i in xrange(gridStorage.size())]))
fd.close()
