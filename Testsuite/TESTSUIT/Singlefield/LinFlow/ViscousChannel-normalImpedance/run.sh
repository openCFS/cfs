#!/bin/sh

# mesh
trelis -batch -nographics -nojournal Thv_channel_PML.jou
# run CFS

## run cfs  
cfs_master ViscousChannel-normalImpedance
cfs_master ViscousChannel-wallEnd
