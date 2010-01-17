#!/usr/bin/env python

# This script print outs the regions and physical field, that can be found in
# result file

import vtk
import numpy
from vtk.util.misc import vtkGetDataRoot
from vtk.util import numpy_support
from vtk.util.colors import *
import sys

#################### input for the script #################### 
# the path and file name
hdf5_fn = 'PATH/FILENAME'
print hdf5_fn
################### END input for the script #################

# read file
CfsReader = vtk.vtkCFSReader()
CfsReader.SetFileName(hdf5_fn)
CfsReader.SetMultiSequenceStep(1)
CfsReader.SetTimeStep(1)
CfsReader.Update()
# start collecting the region names
numRegions = CfsReader.GetNumberOfRegionArrays()
regionNumber=-1
# region name list
regions = list()
for i in range(numRegions):
  regions.append(CfsReader.GetRegionArrayName(i))

# start collecting physical fields
physFields = list()
CfsOut = CfsReader.GetOutput()
# physical name list
arrayNames = list()
for regionNum in range(numRegions):
  CfsBlock = CfsOut.GetBlock(regionNum)
  CfsUnstructGrid = CfsBlock.SafeDownCast(CfsBlock)

  # go through the point datas
  ptData = CfsUnstructGrid.GetPointData()
  numArrays = ptData.GetNumberOfArrays()
  for i in range(numArrays):
    # only collect if it is not stored already
    if not ptData.GetArrayName(i) in arrayNames:
      arrayNames.append(ptData.GetArrayName(i))

  # go through the cell datas
  ptData = CfsUnstructGrid.GetCellData()
  numArrays = ptData.GetNumberOfArrays()
  for i in range(numArrays):
    # only collect if it is not stored already
    if not ptData.GetArrayName(i) in arrayNames:
      arrayNames.append(ptData.GetArrayName(i))


print "Regions are:"
print(regions)
print "Physical fields are:"
print(arrayNames)
