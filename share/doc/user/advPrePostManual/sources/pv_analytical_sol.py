# -*- coding: utf-8 -*-
import scipy as sp
from scipy import special
 
def flatten(input, output):
    # Copy the cells etc.
    output.ShallowCopy(input)
    newPoints = vtk.vtkPoints()
    numPoints = input.GetNumberOfPoints()
    for i in range(0, numPoints):
        coord = input.GetPoint(i)
        x, y, z = coord[:3]
        x = x * 1
        y = y * 1
        z = 3 + 0.5*z
        newPoints.InsertPoint(i, x, y, z)
    output.SetPoints(newPoints)
 
def compute_analytical_sol(input, output, t):
    # Copy the cells etc.
    output.ShallowCopy(input)
    newPoints = vtk.vtkPoints()
    numPoints = input.GetNumberOfPoints()

    # Get input data array 'InputArray'
    input.GetPointData().SetActiveScalars('InputArray')
    values = input.GetPointData().GetScalars()

    # Create a new output data array 'AnalyticalSol' and add it to output
    ca = vtk.vtkFloatArray()
    ca.SetName('AnalyticalSol')
    ca.SetNumberOfComponents(1)
    ca.SetNumberOfTuples(numPoints)
    output.GetPointData().AddArray(ca)

    # Open a file for writing output
    f = open('out.txt','w')

    # Iterate over all nodes of mesh
    for i in range(0, numPoints):
        # Get input coordinates and play around with them
        coord = input.GetPoint(i)
        x, y, z = coord[:3]
        x = x * 1
        y = y * 1
        z = z * 1
        newPoints.InsertPoint(i, x, y, z)

        # Define some physical quantities
        bulkModulus = 1.156e5
        density = 1.0
        c = sqrt(bulkModulus/density)
        freq = 119.3662
        pi = 3.14159
        w = 2*pi*freq
        fac = 0.25 * complex(0,1)
        r = sqrt(x*x + y*y)

        # Get old value for current node from data array
        v_old = values.GetComponent(i,0)

        # Compute new value for current node using Hankel function from SciPy
        v_new = abs(fac * special.hankel1(0, w/c * r) )
        ca.SetValue(i, v_new )

        # Print old and new values into temp file
        print >> f, v_old, ' ', v_new
    output.SetPoints(newPoints)
    # Close temp file
    f.close()

# Print state of oneself. (Same as PrintSelf() in C++)
print self
input = self.GetInputDataObject(0, 0)
# Get current time step
t = input.GetInformation().Get(vtk.vtkDataObject.DATA_TIME_STEPS(),0)
output = self.GetOutputDataObject(0)
 
# The VTK reader for CFS++ HDF5 file will produce vtkMultiBlockDataSets
if input.IsA("vtkMultiBlockDataSet"):
    output.CopyStructure(input)
    iter = input.NewIterator()
    iter.UnRegister(None)
    iter.InitTraversal()

    # Iterate over all blocks (each block is a vtkUnstructuredGrid
    while not iter.IsDoneWithTraversal():
        curInput = iter.GetCurrentDataObject()
        curOutput = curInput.NewInstance()
        curOutput.UnRegister(None)
        output.SetDataSet(iter, curOutput)
        compute_analytical_sol(curInput, curOutput, t)
        iter.GoToNextItem();
else:
  flatten(input, output)
