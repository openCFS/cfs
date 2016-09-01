#!/usr/bin/env python
import numpy as np
#from PIL import Image
# from mesh_tool import *
import sys
import math
from matplotlib import pyplot as plt

def profx(x1,x2,x):
  val = ((x2-x1) * x + x1) / 2.0
  return val+0.1*x
  
def profy(y1,y2,y):
  return ((y2-y1) * y + y1) / 2.0 

def profz(z1,z2,z):
  return ((z2-z1) * z + z1) / 2.0

def f03(x):
  return (1-x)**3
def f13(x):
  return 3*x*(1-x)**2;
def f23(x):
  return 3*x**2*(1-x);
def f33(x):
  return x**3

# k = 4/3.0 * (math.sqrt(2)-1);
# x = np.linspace(0, 1, 100)
# # 
# P = np.transpose(np.array([[0,0.75],[0.25*k,0.75],[0.25,0.75+0.25*k],[0.25,1]]))
# C = np.multiply.outer(P[:,0],f03(x)) + np.multiply.outer(P[:,1],f13(x)) + np.multiply.outer(P[:,2],f23(x)) + np.multiply.outer(P[:,3],f33(x))
# plt.plot(np.transpose(C[0,:]),np.transpose(C[1,:]))
# # plt.plot(np.transpose(P[0,:]),np.transpose(P[1,:]), 'r')
# 
# constx = np.linspace(0.25,0.75,100)
# consty = np.ones((100,1))
# plt.plot(constx,consty)
# 
# P = np.transpose(np.array([[0.75,1],[0.75,1-0.25*k],[1-0.25*k,0.75],[1,0.75]]))
# C = np.multiply.outer(P[:,0],f03(x)) + np.multiply.outer(P[:,1],f13(x)) + np.multiply.outer(P[:,2],f23(x)) + np.multiply.outer(P[:,3],f33(x))
# plt.plot(np.transpose(C[0,:]),np.transpose(C[1,:]))
# 
# constx = np.ones((100,1))
# consty = np.flipud(np.linspace(0.25,0.75,100))
# plt.plot(constx,consty)
# 
# C[1,:] = -C[1,:] + 1
# plt.plot(np.transpose(C[0,:]),np.transpose(C[1,:]))
# plt.xlim([0,1])
# plt.ylim([0,1])



# plt.show()
 
def part1(t):
  k = 4/3.0 * (math.sqrt(2)-1);
  if t <= 0.25:
    P = np.transpose(np.array([[0,0.75],[0.25*k,0.75],[0.25,0.75+0.25*k],[0.25,1]]))
   
#     C = np.multiply.outer(P[:,0],f03(x)) + np.multiply.outer(P[:,1],f13(x)) + np.multiply.outer(P[:,2],f23(x)) + np.multiply.outer(P[:,3],f33(x))
#     plt.plot(np.transpose(C[0,:]),np.transpose(C[1,:]))
#     plt.plot(np.transpose(P[0,:]),np.transpose(P[1,:]), 'r')
#     plt.show()
#     return np.multiply.outer(P[:,0],f03(t)) + np.multiply.outer(P[:,1],f13(t)) + np.multiply.outer(P[:,2],f23(t)) + np.multiply.outer(P[:,3],f33(t))
    return P[:,0]*f03(t) + P[:,1]*f13(t) + P[:,2]*f23(t) + P[:,3]*f33(t)
 
  if t > 0.25 and t < 0.75:
    return np.array([t,1])
  
  if t >= 0.75:
    P = np.transpose(np.array([[0.75,1],[0.75,1-0.25*k],[1-0.25*k,0.75],[1,0.75]]))
    return P[:,0]*f03(t) + P[:,1]*f13(t) + P[:,2]*f23(t) + P[:,3]*f33(t)
  
def calcRadius2d(t):
  coords = part1(t)
  print str(t) + " coords=" + str(coords) + " radius = " + str(coords[1] - 0.5) 
  return coords[1] - 0.5
# def hori1(x1,x2,y1,x):
#   c = 0.125
#   if x <= c:
#     return 0.5*x1 + (x**3*(4*c + 4*x1 - 3))/(4.0*c**3) - (x**2*(8*c + 12*x1 - 9))/(8.0*c**2) + 0.5
#     
#   if x > c and x < 0.5 + 0.5*y1 + c:
#     return 1-c
#   
#   if x >= 0.5 + 0.5*y1 + c:
#     return ((c - 1)*(c*x2 - 2*x2 - 7*c + c**2*x2 + c**2 + 2))/(2.0*c**3) - (x**2*(3*c*x2 - 6*x2 - 21*c + 8*c**2 + 6))/(2.0*c**3) - (x**3*(3*c + x2 - 1))/c**3 + (x*(3*c*x2 - 3*x2 - 12*c + 8*c**2 + 3))/c**3
#     
# def hori2(x1,x2,y2,x):
#   c = 0.125
#   if x <= c:
#     return (x**2*(8*c + 3*x1 - 3))/(2*c**2) - (x**3*(3*c + x1 - 1))/c**3 - 0.5*x1 + 0.5
#    
#   if x > c and x < 0.5 + 0.5*y2 + c:
#     return c
#    
#   if x >= 0.5 + 0.5*y2 + c:
#     return (3*c*x2 - 2*x2 - 9*c - c**3*x2 + 8*c**2 + c**3 + 2)/(2.0*c**3) + (x**2*(3*c*x2 - 6*x2 - 21*c + 8*c**2 + 6))/(2*c**3) + (x**3*(3*c + x2 - 1))/c**3 - (x*(3*c*x2 - 3*x2 - 12*c + 8*c**2 + 3))/c**3
#   
#   return 0
# 
# def vert1():
#   c = 0.125
#   return (8*c**3 - 32*c**2*x**2 - 32*c**2*x*y1 + 32*c**2*x - 8*c^2*y1**2 + 28*c**2*y1 - 20*c**2 + 24*c*x**3 + 16*c*x**2*y1 - 16*c*x**2 - 2*c*x*y1**2 + 4*c*x*y1 - 2*c*x - 2*c*y1**3 + 12*c*y1**2 - 18*c*y1 + 8*c + 4*x**3*y1 - 4*x**3 + 4*x**2*y1**2 - 8*x**2*y1 + 4*x**2 + x*y1**3 - 3*x*y1**2 + 3*x*y1 - x + y1**3 - 3*y1**2 + 3*y1 - 1)/(2*c + y1 - 1)**3

######### plotting spline #########################

  
def within_lin_profiles3d(x1,x2,y1,y2,z1,z2,coords):
  xm = 0.5
  ym = 0.5
  zm = 0.5
  x = coords[0]
  y = coords[1]
  z = coords[2]
  
  px = profx(x1, x2, x)
  py = profy(y1, y2, y)
  pz = profz(z1, z2, z)
  
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

def within_cubic_profiles2d(nx,x1,x2,y1,y2,coords):
  xm = 0.5
  ym = 0.5
  
  x = coords[0]
  y = coords[1]
  
  rad = calcRadius2d(x)
  valUp = 0.5+rad
  valLow = 0.5-rad
  
  if (y <= valUp) and (y >= valLow):
    return 255
  
#   rad = calcRadiusHor(x1, x2, y1, y2, y)
#   valUp = 0.5+0.5*rad
#   valLow = 0.5-0.5*rad
#   
#   if (x <= valUp) and (x >= valLow):
#     return 155
    
  return -1
  
  
def set_lin_profiles_array(nx,ny,nz,x1,x2,y1,y2,z1,z2):
  array = np.zeros((nx,ny,nz))
  sizex = 1.0
  sizey = 1.0
  sizez = 1.0
  
  hx = sizex / nx
  hy = sizey / ny
  hz = sizez / nz
  
  for i in range(0,nx):
    for j in range(0,ny):
      for k in range(0,nz):  
        coords = np.zeros(3)
        x = i * hx + hx / 2.0
        y = j * hy + hy / 2.0
        z = k * hz + hz / 2.0
        coords[0] = x
        coords[1] = y
        coords[2] = z
        if nz == 1:
#           flag = within_profile2d(x1, x2, y1, y2, coords)
          flag = within_lin_profiles3d(x1, x2, y1, y2, 0,0, coords)
        else:
          flag = within_lin_profiles3d(x1, x2, y1, y2, z1, z2, coords) 
        if flag != -1:
          array[i,j,k] = flag
  return array

def set_cubic_profiles_array(nx,ny,nz,x1,x2,y1,y2,z1,z2):
  array = np.zeros((nx,ny,nz))
  sizex = 1.0
  sizey = 1.0
  sizez = 1.0
  
  hx = sizex / nx
  hy = sizey / ny
  hz = sizez / nz
  
  for i in range(0,nx):
    for j in range(0,ny):
      for k in range(0,nz):  
        coords = np.zeros(3)
        x = i * hx + hx / 2.0
        y = j * hy + hy / 2.0
        z = k * hz + hz / 2.0
        coords[0] = x
        coords[1] = y
        coords[2] = z
        if nz == 1:
#           flag = within_profile2d(x1, x2, y1, y2, coords)
          flag = within_cubic_profiles2d(nx,x1, x2, y1, y2, coords)
#         else:
#           flag = within_cubic_profiles3d(x1, x2, y1, y2, z1, z2, coords) 
        if flag != -1:
          array[i,j,k] = flag
  return array

def visualize_structure(array, nx, ny, nz):
  from vtk import *
  
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
