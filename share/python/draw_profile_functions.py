#!/usr/bin/env python
import numpy as np
#from PIL import Image
# from mesh_tool import *
import sys
import math
from matplotlib import pyplot as plt
from scipy import interpolate
from cmath import atan

def profileLin(x1,x2,n):
  x = np.linspace(0, 1.0,n)
  return ((x2-x1) * x + x1) / 2.0
  

def profileCirc(x1,res):
  x = np.linspace(0, 1.0, res)
  r = 0.5-x1/2.0 # radius for circular holes not stiffness
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

def f03(t):
  return (1-t)**3
def f13(t):
  return 3*t*(1-t)**2
def f23(t):
  return 3*t**2*(1-t)
def f33(t):
  return t**3

def findGradOne(vec):
  idx = 0 # idx where grad is about 1
  stop = False
  while not stop and idx < len(vec):
    if vec[idx] <= 1:
      idx += 1
    else:
      stop = True
      
  if not stop:
    idx = len(vec) - 1
    
  return idx

def profileSpline(x1,y1,res,bend):
  assert(bend <= 1 and bend >= 0)
  rx = 0.5 - x1/2.0 # radius for center (0,1)
  ry = 0.5 - y1/2.0 # radius for center (0,1)
  
  P = np.transpose(np.array([[0,1-rx],[ry*bend,1-rx],[ry,1-rx*bend],[ry,1]]))
  t = np.linspace(0, 1.0, 2*ry*res) # double over-sampling
  C = np.multiply.outer(P[:,0],f03(t)) + np.multiply.outer(P[:,1],f13(t)) + np.multiply.outer(P[:,2],f23(t)) + np.multiply.outer(P[:,3],f33(t))
  G = np.zeros(len(C[0])) # numerical finite differences
  for i in range(len(C[0])-1):
    G[i] = (C[1,i+1] - C[1,i]) / (C[0,i+1] - C[0,i])

  # interpolate for regular spacing  
  fc = interpolate.interp1d(C[0],C[1])
  fg = interpolate.interp1d(C[0],G)
  
  v = np.linspace(0,ry,ry*res)
  c = fc(v)
  g = fg(v)
  
  idx = findGradOne(g)
  
#   print "index: ",idx,v[idx],c[idx] 
  
  vec = np.zeros(res)
  
  vec[0:idx] = c[0:idx]
  vec[idx:res-idx] = c[idx] # constant value at position where grad is one
  vec[res-idx:res] = c[0:idx][::-1]
  
#   plt.plot(np.transpose(C[0,:]),np.transpose(C[1,:]))
#   plt.plot(np.transpose(P[0,:]),np.transpose(P[1,:]), 'r')
#  plt.show()
  
  return vec-0.5,c[idx]

# @param height is second return value from profileSpline
def profileSplineBisec(x1,y1,height,res,bend):
  assert(bend <= 1 and bend >= 0)
  rx = x1 / 2.0
  ry = y1 / 2.0 
    
  P = np.transpose(np.array([[0,0.5+rx],[ry*bend,0.5+rx],[0.5-ry,height],[0.5,height]]))
  t = np.linspace(0, 1.0, res) # double over-sampling
  C = np.multiply.outer(P[:,0],f03(t)) + np.multiply.outer(P[:,1],f13(t)) + np.multiply.outer(P[:,2],f23(t)) + np.multiply.outer(P[:,3],f33(t))

  # interpolate for regular spacing  
  fc = interpolate.interp1d(C[0],C[1])
  
  v = np.linspace(0,0.5,res/2)
  c = fc(v)
  
  vec = np.zeros(res)
  
  vec[0:res/2] = c
  vec[res/2:res] = c[::-1]
  
  return vec-0.5
  
# @return vector with profile or list of vectors
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
    
  if args.profile == 'spline':
    if dir == 1:
      vec1,h = profileSpline(args.x1, args.y1, args.res, args.bend)
      vec2 = profileSplineBisec(args.x1, 0.5*(args.y1+args.z1), h, args.res, args.bend)
      vec3,h = profileSpline(args.x1, args.z1, args.res, args.bend)
      return (vec1,vec2,vec3)
    if dir == 2:   
      vec1, h = profileSpline(args.y1, args.x1, args.res, args.bend)
      vec2 = profileSplineBisec(args.y1, 0.5*(args.z1+args.x1), h, args.res, args.bend)
      vec3,h = profileSpline(args.y1, args.z1, args.res, args.bend)
      return (vec1,vec2,vec3)
    if dir == 3:
      vec1, h = profileSpline(args.z1, args.y1, args.res, args.bend)
      vec2 = profileSplineBisec(args.z1, 0.5*(args.x1+args.y1), h, args.res, args.bend)
      vec3,h = profileSpline(args.z1, args.x1, args.res, args.bend)
      return (vec1,vec2,vec3)

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

## give an interpolation for the radius within 0..2pi angle for the radius vectors at index idx
def give_interpolate_radius(vec, idx, dir):
  if type(vec) == tuple:
    assert(len(vec) == 3)
    val = []
    rad = []
    
    val.append(vec[0][idx]) # 0 deg
    rad.append(0.0)
    
    val.append(vec[1][idx]) # 45 de
    rad.append(np.pi/4)
    
    val.append(vec[2][idx]) # 90
    rad.append(np.pi/2)
    
    val.append(vec[1][idx]) # 90+45
    rad.append(3./4*np.pi)
    
    val.append(vec[0][idx]) # 90+90=180
    rad.append(np.pi)
    
    val.append(vec[1][idx]) # 180+45
    rad.append(5./4*np.pi)
    
    val.append(vec[2][idx]) # 180+90=270
    rad.append(6./4*np.pi)
    
    val.append(vec[1][idx]) # 270+45
    rad.append(7./4*np.pi)
    
    val.append(vec[0][idx]) # 270+90=360=0 (full cicle) 
    rad.append(0.0)
      
    f = interpolate.interp1d(rad, val)
    return f
  else:
    assert(type(vec) <> tuple) # but a single vector
    f = interpolate.interp1d((0, 7),(vec[idx], vec[idx]))
    return f                            
                                       


# @param vec can be one vector or list of vector
def draw_profile(array,vec,dir):
  res = array.shape[0]
  h = 1.0/res
  assert(dir >=1 and dir <=3)
  
  q = []
  pa = []
#   a = []
#   ya = []
#   za = []
  
  for i in range(0,res):
    f = give_interpolate_radius(vec, i, dir)
    for j in range(0,res):
      for k in range(0,res):
        y = j * h + h / 2.0
        z = k * h + h / 2.0
        valx = (y-0.5)**2 + (z-0.5)**2
        
        phi = np.arccos(0.5*(y-0.5)/(0.5*np.sqrt(valx))) # angle between (0.5,0) and (y,z)
        if y < 0:
          phi = 2*np.pi - phi
          
        if i == 25:
          q.append(phi)
          pa.append(f(phi))
         
#           ya.append(y-0.5)
#           za.append(z-0.5)  
#           a.append(0.5*y/(0.5*np.sqrt(valx)))
#          q.append(phi)
#          v.append(valx)
        
        p = f(phi)   
        if (valx <= p*p):
          if dir == 1:
            array[i,j,k] = dir
          if dir == 2:
            array[j,i,k] = dir
          if dir == 3:
            array[k,j,i] = dir  

  return q, pa
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
