#!/usr/bin/env python
import numpy as np
#from PIL import Image
# from mesh_tool import *
import sys
import math
from matplotlib import pyplot as plt

def profile(x1,x2,x, profType):
  assert(x<=1 and x >=0)
  if profType == "linear":
    return ((x2-x1) * x + x1) / 2.0
  else:
    cornerx = np.sqrt(x1**2/2.0)
    cornery = 1 - cornerx
    if x <= cornerx:
      return 0.5 - np.sqrt(x1**2 - x**2)
    if x < 1 - cornerx:
      return cornery - 0.5
    else:
      return 0.5 - np.sqrt(x1**2 - (1-x)**2)   

  
def within_profiles3d(x1,x2,y1,y2,z1,z2,coords,profType,skip_x,skip_y,skip_z):
  xm = 0.5
  ym = 0.5
  zm = 0.5
  x = coords[0]
  y = coords[1]
  z = coords[2]
  
  px = profile(x1, x2, x, profType) if not skip_x else 0.0
  py = profile(y1, y2, y, profType) if not skip_y else 0.0
  pz = profile(z1, z2, z, profType) if not skip_z else 0.0
  
  valx = (y-ym)**2 + (z-zm)**2
  valy = (x-xm)**2 + (z-zm)**2
  valz = (y-ym)**2 + (x-xm)**2
  
  if (valx <= px*px):
    return 1
  elif (valy <= py*py):
    return 2
  elif (valz <= pz*pz):
    return 3
  
  return -1

# def create_profiles_array(nx,ny,nz,x1,x2,y1,y2,z1,z2):
def create_profiles_array(args):
  res = args.res
  array = np.zeros((res,res,res))
    
  h = 1.0 / res
  
  flag = None
  
  for i in range(0,res):
    for j in range(0,res):
      for k in range(0,res):  
        coords = np.zeros(3)
        x = i * h + h / 2.0 # from array coordinates to cartesian ones (barycenter)
        y = j * h + h / 2.0
        z = k * h + h / 2.0
        coords[0] = x
        coords[1] = y
        coords[2] = z
        if args.type == 'profiles2d':
#           flag = within_profile2d(x1, x2, y1, y2, coords)
          flag = within_profiles3d(args.x1, args.x2, args.y1, args.y2, 0,0, coords, args.profile, args.skip_x, args.skip_y, args.skip_z)
        else:
          flag = within_profiles3d(args.x1, args.x2, args.y1, args.y2, args.z1, args.z2, coords, args.profile, args.skip_x, args.skip_y, args.skip_z) 
        if flag != -1:
          array[i,j,k] = flag
  return array

def visualize_structure(array, nx, ny, nz):
  import vtk
  
  # Setup a renderer, render window, and interactor
  renderer = vtk.vtkRenderer()
  renderWindow = vtk.vtkRenderWindow()
  
  renderWindow.AddRenderer(renderer);
  renderWindowInteractor = vtk.vtkRenderWindowInteractor()
  renderWindowInteractor.SetRenderWindow(renderWindow)
   
  sizex = 1.0
  sizey = 1.0
  sizez = 1.0
   
  hx = sizex / nx
  hy = sizey / ny
  hz = sizez / nz
   
  for i in range(nx):
    for j in range(ny):
      for k in range(nz):
        x = i * hx + hx / 2.0
        y = j * hy + hy / 2.0
        z = k * hz + hz / 2.0
        
        if array[i,j,k] <= 0:
          continue
          
        # create cube
        cube = vtk.vtkCubeSource()
        cube.SetCenter(x,y,z)
        cube.SetXLength(hx)
        cube.SetYLength(hy)
        cube.SetZLength(hz)
        cube.Update() 
          
        # mapper
        cubeMapper = vtk.vtkPolyDataMapper()
        cubeMapper.SetInputConnection(cube.GetOutputPort())
           
        # actor
        cubeActor = vtk.vtkActor()
        cubeActor.SetMapper(cubeMapper)
        
        if array[i,j,k] == 1:
          cubeActor.GetProperty().SetColor(255, 0, 0)
        elif array[i,j,k]  == 2:
          cubeActor.GetProperty().SetColor(0, 255, 0)
        elif array[i,j,k]  == 3:
          cubeActor.GetProperty().SetColor(0, 0, 255)
#         else:
#           print i,j,k,array[i,j,k]
        # assign actor to the renderer
        renderer.AddActor(cubeActor)

  # Render and interact
  renderWindow.Render()
  renderWindowInteractor.Start()
