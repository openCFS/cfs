#!/usr/bin/env python

# This script extracts the result of a slice (a line in 2d and a plain in 3d)
# over time and stores it in a list of vectors. The results are stored in
# resultList.

import vtk
import numpy
from vtk.util.misc import vtkGetDataRoot
from vtk.util import numpy_support
from vtk.util.colors import *
import sys

#################### input for the script #################### 
# physical field from which to extract result
# EXAMPLE: phys_name = 'magFluxDensity'
phys_name = 'PHYSICAL FIELD'
# The region from were to extract the point
regionName = 'REGION'
# the path and file name
hdf5_fn = 'PATH/FILENAME'
# The x-, y- and z-Coordinates of the point through which the slice should go
# throguh
xCoord = 0.1
yCoord = 0.0
zCoord = 0.0
# The vector that is normal to the slice
xNormal = 1.0
yNormal = 0.0
zNormal = 0.0
print hdf5_fn
################### END input for the script #################

# This Method gives back a reference to the Data, the Array number were it can
# be found at how many degrees of freedom it has. It needs a reference to the
# reader and which kind of physical field is wanted.
def getArrayFromReader(Reader, physFieldName):
  CfsOut = Reader.GetOutput()
  CfsBlock = CfsOut.GetBlock(regionNumber)
  CfsUnstructGrid = CfsBlock.SafeDownCast(CfsBlock)
  ptData = CfsUnstructGrid.GetPointData()
  numArrays = ptData.GetNumberOfArrays()
  arrayName = ''
  dataType = ''
  for i in range(numArrays):
    arrayName = ptData.GetArrayName(i)
    # is this the correct array?
    if arrayName == physFieldName:
      # get size
      arraySize = ptData.GetArray(i).GetDataSize()
      # get dimension (1D, 2D, 3D)
      dim= arraySize/ptData.GetArray(i).GetNumberOfTuples()
      arrayNum = i
      dataType = 'pointData'
      break

  if arrayName != physFieldName:
    ptData = CfsUnstructGrid.GetCellData()
    numArrays = ptData.GetNumberOfArrays()
    arrayName = ''
    idList = vtk.vtkIdList()
    CfsUnstructGrid.GetPointCells(ptNumber,idList)
    cellNumber = idList.GetId(0)
    ptNumber = cellNumber
    # search for the correct result array
    for i in range(numArrays):
      arrayName = ptData.GetArrayName(i)
      # is this the correct array?
      if arrayName == physFieldName:
        # get size
        arraySize = ptData.GetArray(i).GetDataSize()
        # get dimension (1D, 2D, 3D)
        dim= arraySize/ptData.GetArray(i).GetNumberOfTuples()
        arrayNum = i
        dataType = 'cellData'
        break
    if arrayName != physFieldName:
      print "Error: Physical field \"",physFieldName + " \" does not exist in result file"
      exit(-1)
  return ptData, arrayNum, dim, dataType

CfsReader = vtk.vtkCFSReader()
CfsReader.SetFileName(hdf5_fn)
CfsReader.SetMultiSequenceStep(1)
CfsReader.SetTimeStep(1)
CfsReader.Update()
numRegions = CfsReader.GetNumberOfRegionArrays()
regionNumber=-1
for i in range(numRegions):
  if regionName == CfsReader.GetRegionArrayName(i):
    regionNumber = i
    break
if regionNumber == -1:
  print "could not find region:", regionName

ptData, arrayNum, dimension, dType = getArrayFromReader(CfsReader, phys_name)

timeStepRange = CfsReader.GetTimeStepRange()
numTimeSteps = timeStepRange[1]

CfsOut = CfsReader.GetOutput()
CfsBlock = CfsOut.GetBlock(regionNumber)
CfsUnstructGrid = CfsBlock.SafeDownCast(CfsBlock)

ptNumber = CfsUnstructGrid.FindPoint(xCoord, yCoord, zCoord)
foundCoord = CfsUnstructGrid.GetPoint(ptNumber) 

# The cutter uses an implicit function to perform the cutting.
# Here we define a plane, specifying its center and normal.
# Then we assign the plane to the cutter.
plane = vtk.vtkPlane()
plane.SetOrigin(xCoord, yCoord, zCoord)
plane.SetNormal(xNormal, yNormal, zNormal)
planeCut = vtk.vtkCutter()
planeCut.SetInput(CfsUnstructGrid)
planeCut.SetCutFunction(plane)
planeCut.Update()
cutMapper = vtk.vtkPolyDataMapper()
cutMapper.SetInputConnection(planeCut.GetOutputPort())
cutMapper.SetScalarRange(ptData.GetArray(arrayNum).GetRange())
cutActor = vtk.vtkActor()
cutActor.SetMapper(cutMapper)

# create point
point = vtk.vtkPointSource()
point.SetCenter(foundCoord)
#point.SetCenter(xCoord, yCoord, zCoord)
#point.SetNumberOfPoints(5)
point.SetRadius(0.0001)
# Move the point into place and create the probe filter.  For
# vtkProbeFilter, the probe point is the input, and the underlying data
# set is the source.
transL1 = vtk.vtkTransform()
transL1.Translate(0.0, 0.0, 0.0)
transL1.Scale(1, 1, 1)
transL1.RotateY(90)
tf = vtk.vtkTransformPolyDataFilter()
tf.SetInputConnection(point.GetOutputPort())
tf.SetTransform(transL1)

ball = vtk.vtkSphereSource()
ball.SetRadius(0.001)
ball.SetThetaResolution(12)
ball.SetPhiResolution(12)
balls = vtk.vtkGlyph3D()
balls.SetInputConnection(point.GetOutputPort())
balls.SetSourceConnection(ball.GetOutputPort())
mapBalls = vtk.vtkPolyDataMapper()
mapBalls.SetInputConnection(balls.GetOutputPort())
ballActor = vtk.vtkActor()
ballActor.SetMapper(mapBalls)
ballActor.GetProperty().SetColor(hot_pink)
ballActor.GetProperty().SetSpecularColor(1, 1, 1)
ballActor.GetProperty().SetSpecular(0.3)
ballActor.GetProperty().SetSpecularPower(20)
ballActor.GetProperty().SetAmbient(0.2)
ballActor.GetProperty().SetDiffuse(0.8)

contour = vtk.vtkDataSetSurfaceFilter()       
contour.SetInput(CfsUnstructGrid) 
contourMapper = vtk.vtkPolyDataMapper()
contourMapper.SetInputConnection(contour.GetOutputPort()) 
contourActor = vtk.vtkActor()                                                                           
contourActor.SetMapper(contourMapper)                                                                   
contourProp = contourActor.GetProperty()                                                                
contourProp.SetColor(0, 0, 0) 
contourProp.SetOpacity(0.3) 

## Create the RenderWindow, Renderer and both Actors
ren = vtk.vtkRenderer()
renWin = vtk.vtkRenderWindow()
renWin.AddRenderer(ren)
iren = vtk.vtkRenderWindowInteractor()
iren.SetRenderWindow(renWin)
#
## Add the actors to the renderer, set the background and size
ren.AddActor(cutActor)
ren.AddActor(contourActor)
ren.AddActor(ballActor)

ren.SetBackground(1, 1, 1)
renWin.SetSize(800, 800)

cam1 = ren.GetActiveCamera()
#cam1.SetClippingRange(11.1034, 59.5328)
#cam1.SetFocalPoint(1.71821, 0.458166, 1.3999)
#cam1.SetPosition(-2.0, -1.0, 1.0)
#cam1.SetViewUp(0.0184785, 0.479657, 0.877262)

iren.Initialize()
renWin.Render()
# visualise at the end
# iren.Start()
# search for the correct result array

resultList = list()
for tStep in range(1,numTimeSteps+1):
  print "tStep: ", tStep
  CfsReader.SetTimeStep(tStep)
  CfsReader.Update()
  CfsBlock = CfsOut.GetBlock(regionNumber)
  CfsUnstructGrid = CfsBlock.SafeDownCast(CfsBlock)
  #ptData, arrayNum, dimension, dtType = getArrayFromReader(CfsReader, phys_name)
  planeCut.SetInput(CfsUnstructGrid)
  planeCut.Update()
  cutOut = planeCut.GetOutput()
  if dType == 'pointData':
    cutData = cutOut.GetPointData()
  else:
    if dType == 'cellData':
      cutData = cutOut.GetCellData()
  resultArray = numpy_support.vtk_to_numpy(cutData.GetArray(arrayNum))
  resultList.append(resultArray)

print "Values of \"", phys_name + " \":" 
print(resultList)

iren.Start()
