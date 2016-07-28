#!/usr/bin/env python
import numpy as np
#from PIL import Image
from mesh_tool import *
import sys

def profx(x1,x2,x):
  return ((x2-x1) * x + x1) / 2.0
  
def profy(y1,y2,y):
  return ((y2-y1) * y + y1) / 2.0 

def profz(z1,z2,z):
  return ((z2-z1) * z + z1) / 2.0

def within_profile3d(x1,x2,y1,y2,z1,z2,coords):
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
  
def set_profile_array(nx,ny,nz,x1,x2,y1,y2,z1,z2):
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
          flag = within_profile3d(x1, x2, y1, y2, 0,0, coords)
        else:
          flag = within_profile3d(x1, x2, y1, y2, z1, z2, coords) 
        if flag != -1:
          array[i,j,k] = flag
  return array
  
  
