#!/bin/sh

# mesh
gmsh -3 -format msh2 cyliders.geo

# run test
cfs_master -p RotatingDisk3D.xml test
