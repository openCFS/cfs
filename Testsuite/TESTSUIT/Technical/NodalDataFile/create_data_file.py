# this is for the test case NodalDataFile
# we need to know the element numbers for the named nodes inclusion_nodes
# we could use make_mesh.py from NodalForceExpression and make use of
# the stored mesh.bc. 
# Here we choose another way and assume cfs NodalDataFile -g performed.
# the -g creates the .cfs file without parsing the pdes settings and w/o FEM

from hdf5_tools import *

# get named nodes
nn = get_nodes('results_hdf5/NodalDataFile.cfs','inclusion_nodes')
nnn = len(nn)
# the corresponding coordinates
nc = get_node_coordinates('results_hdf5/NodalDataFile.cfs','inclusion_nodes')
assert nnn == len(nc)

print('# fileformat of data file is 1-based node number, nodal x-vale, nodal y-value')
print('# force is divided by number of nodes')
print('# 1:nodenumber 2:x_value 3:y_value')
for i, n in enumerate(nn):
  x, y, z = nc[i]
  print(n, round((x-.5)/nnn,15), round((y-.5)/nnn,15))
