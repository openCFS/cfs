#!/bin/sh

gmsh -3 -format msh2 Coax.geo

cfs=/opt/cfs/master/bin/cfs
for MAT in 'iso-1' 'iso-4' 'iso-9' 'aniso-19' 'aniso-49'; do
  sed 's|region name="V_1"[ ]\+material="\(.\+\)"|region name="V_1" material="'$MAT'"|' CoaxEdgeAnisotropicConductivity.xml > $MAT'.xml'
  $cfs $MAT
done
