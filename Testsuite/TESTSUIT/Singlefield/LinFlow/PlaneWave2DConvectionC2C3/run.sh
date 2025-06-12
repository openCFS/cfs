#!/bin/sh

# mesh (creates domain.msh)
gmsh -2 -format msh2 domain.geo

cfs=/opt/cfs/build_dbg/bin/cfs

for Vx in '0.0' '0.3' '0.5'; do
  job="test_M$Vx"
  echo $job
  sed 's|\(var\s\+name="Vx"\s\+value="\).\+"|\1'$Vx'"|g' PlaneWave2DConvection.xml > $job'.xml'
  $cfs $job
done
