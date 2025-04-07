#!/bin/sh

gmsh -2 -format msh2 slit.geo
gmsh -2 -format msh2 acou.geo

sed 's|meshSizeY[ ]\+=[ ]\+.\+;|meshSizeY=0.0075e-3;|g' slit.geo > slit-visco.geo
gmsh -2 -format msh2 slit-visco.geo

cfs_master -p ViscoThermalBoundaryLayer2dChannel.xml visco-thermal
/opt/cfs/build_dbg/bin/cfs -p AcousticBoundaryLayer2dChannel.xml acou-thermal
