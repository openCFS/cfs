#!/usr/bin/env python
import numpy as np
#from PIL import Image
# from mesh_tool import *
import sys
import math
from matplotlib import pyplot as plt

def profileLin(x1,x2,n):
  x = np.linspace(0, 1.0,n)
  return ((x2-x1) * x + x1) / 2.0
  

def profileCirc(x1,res):
  x = np.linspace(0, 1.0, res)
  r = 0.5-x1/2.0
  r2 = r**2
  cornerx = np.sqrt(r2/2.0)
  cornery = 1 - cornerx
  
  vec = np.zeros(res)
  
  for i in range(res):
    if x[i] <= cornerx:
      vec[i] = 0.5 - np.sqrt(r2 - x[i]**2)
    elif x[i] < 1 - cornerx:
      vec[i] = cornery - 0.5
    else:
      vec[i] = 0.5 - np.sqrt(r2 - (1-x[i])**2)
  
  return vec
# def profileSpline(x1,x2):

def profile(args,dir):
  if args.profile == 'linear':
    if dir == 1:
      return profileLin(args.x1, args.x2, args.res)
    if dir == 2:   
      return profileLin(args.y1, args.y2, args.res)
    if dir == 3:
      return profileLin(args.z1, args.z2, args.res)
    
  if args.profile == 'circular':
    if dir == 1:
      return profileCirc(args.x1, args.res)
    if dir == 2:   
      return profileCirc(args.y1, args.res)
    if dir == 3:
      return profileCirc(args.z1, args.res)  
# def create_profiles_array(nx,ny,nz,x1,x2,y1,y2,z1,z2):
def create_profiles_array(args):
  res = args.res
  array = np.ones((res,res,res)) * (-1)
  
  if not args.skip_x:
    vec = profile(args,1)
    draw_profile(array, vec, 1)
  if not args.skip_y:
    vec = profile(args,2)
    draw_profile(array, vec, 2)  
  if not args.skip_z:
    vec = profile(args,3)
    draw_profile(array, vec, 3)  
    
  return array

def draw_profile(array,vec,dir):
  res = array.shape[0]
  h = 1.0/res
  assert(res == len(vec))
  assert(dir >=1 and dir <=3)

  for i in range(0,res):
    for j in range(0,res):
      for k in range(0,res):  
        y = j * h + h / 2.0
        z = k * h + h / 2.0
        valx = (y-0.5)**2 + (z-0.5)**2
  
        p = vec[i]   
        if (valx <= p*p):
          if dir == 1:
            array[i,j,k] = dir
          if dir == 2:
            array[j,i,k] = dir
          if dir == 3:
            array[k,j,i] = dir  


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
