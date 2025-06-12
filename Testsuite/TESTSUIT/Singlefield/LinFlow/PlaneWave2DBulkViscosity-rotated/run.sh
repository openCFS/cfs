#!/bin/sh

# mesh
trelis -batch -nographics -nojournal AcousticChannel.jou
trelis -batch -nographics -nojournal Channel-rotated.jou
# run CFS

## run cfs  
cfs_master PlaneWave2DBulkViscosity-rotated
cfs_master AcousticChannel
