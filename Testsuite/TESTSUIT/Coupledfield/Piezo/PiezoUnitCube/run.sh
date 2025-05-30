#!/bin/sh

# mesh
trelis -batch -nographics -nojournal UnitCube.jou

# run CFS
cfs PiezoUnitCube

# paraview results_hdf5/PiezoUnitCube.cfs 


