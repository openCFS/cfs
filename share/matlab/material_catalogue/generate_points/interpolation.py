#!/usr/bin/python

# import modules
from pysgpp import DataVector, Grid, createOperationHierarchisation, createOperationEval
from pysgpp import GridDataBase, DataMatrix

from bin.tools import readGrid, readAlpha, readDataTrivial, writeAlphaARFF, writeDataARFF

# load grid
# grid = readGrid('a.grid')
# values = readAlpha('a.arff')

data = readDataTrivial('detailed_stats_64_3D', hasclass=False)['data']
# A = DataMatrix(data.getNrows(), 3)
# p = DataVector(data.getNrows())

# for i in xrange(A.getNcols()):
#     data.getColumn(i, p)
#     p.mult(1. / 64.)
#     A.setColumn(i, p)

# p = DataVector(3)
# fd = open('material.csv', 'w')
# for i in xrange(A.getNrows()):
#     A.getRow(i, p)
#     fd.write(" ".join([str(x) for x in p.array()] + [str(data.get(i, 4))]))
#     fd.write("\n")
# fd.close()


def createDB():
    data = readDataTrivial('detailed_stats_64_3D', hasclass=False)['data']

    dim = 3
    grid = Grid.createLinearTrapezoidBoundaryGrid(dim)
    gridStorage = grid.getStorage()

    # create regular grid, level 3
    level = 6
    gridGen = grid.createGridGenerator()
    gridGen.full(level)
    print grid.getStorage().size(), data.getNrows()

    beta = DataVector(gridStorage.size())

    for k in xrange(6):
        for i in xrange(gridStorage.size()):
            gp = gridStorage.get(i)
            # calculate index in file
            ind = (1
                   + 2**(level-gp.getLevel(0))*gp.getIndex(0)*65**2
                   + 2**(level-gp.getLevel(1))*gp.getIndex(1)*65
                   + 2**(level-gp.getLevel(2))*gp.getIndex(2))
            beta[i] = data.get(ind, 3+k)

        # hierarchize
        createOperationHierarchisation(grid).doHierarchisation(beta)
        # save to DB
        db = GridDataBase(grid, beta)
        db.save("stats_64_3D_%d.db" % k)

# =============== END ========================

db = GridDataBase('stats_64_3D_0.db')

# create a two-dimensional piecewise bi-linear grid
dim = 3
grid = Grid.createLinearTrapezoidBoundaryGrid(dim)
gridStorage = grid.getStorage()
print "dimensionality:         %d" % (gridStorage.dim())

# create regular grid, level 3
level = 6
gridGen = grid.createGridGenerator()
gridGen.regular(level)
print "number of grid points:  %d" % (gridStorage.size())

# create coefficient vector
alpha = DataVector(gridStorage.size())
alpha.setAll(0.0)
print "length of alpha-vector: %d" % (len(alpha))

# set function values in alpha
for i in xrange(gridStorage.size()):
    gp = gridStorage.get(i)
    alpha[i] = db.get(gp)

# hierarchize
createOperationHierarchisation(grid).doHierarchisation(alpha)

# sanity check
opEval = createOperationEval(grid)
for i in xrange(data.getNrows()):
    p = DataVector([data.get(i, j) for j in xrange(dim)])
    p.mult(1. / 64.)
    print p, ":", opEval.eval(alpha, p), "==", data.get(i, 4)

# write grid and surpluses to file
fd = open('material.grid', 'w')
fd.write(grid.serialize())
fd.close()
writeAlphaARFF('material.alpha.arff', alpha)
