#!/bin/sh

# mesh
trelis -batch -nographics -nojournal Stokes.jou

# run CFS
cfs StokesBoundaryLayer2D

# paraview results_hdf5/Stokes.cfs 


