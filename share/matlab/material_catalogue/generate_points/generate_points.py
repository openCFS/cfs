#!/usr/bin/python

# import modules
from __future__ import print_function
import sys
sys.path.append('/home/daniel/code/sgopt_2016-03-04_166a3d9/lib')
#import pysgpp
from pysgpp import Grid

dim = int(sys.argv[1])
level = int(sys.argv[2])

bsplineDegree = 3
#grid = Grid.createModBsplineGrid(dim,bsplineDegree)
grid = Grid.createBsplineBoundaryGrid(dim,bsplineDegree)
gridStorage = grid.getStorage()

# create regular grid
gridGen = grid.getGenerator()
gridGen.regular(level)
print("number of grid points:  %d" % gridStorage.getSize())

for i in xrange(gridStorage.getSize()):
  point = gridStorage.get(i)
#  print(point.getCoord(0),point.getCoord(1),point.getCoord(2))

fd = open('grid_points.csv', 'w')
fd.write("\n".join([gridStorage.get(i).toString() for i in xrange(gridStorage.getSize())]))
fd.close()
