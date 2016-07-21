#!/usr/bin/env python
import numpy as np
#from PIL import Image
from mesh_tool import *
import sys

def bar1(x1,x2,x):
  return ((x2-x1) * x + x1) / 2.0
  
def bar2(y1,y2,y):
  return ((y2-y1) * y + y1) / 2.0 

def bar3(z1,z2,z):
  return ((z2-z1) * z + z1) / 2.0

def within_profile(x1,x2,y1,y2,coords):
  xm = 0.5
  ym = 0.5
  x = coords[0]
  y = coords[1]
  profile1 = bar1(x1, x2, x)
  profile2 = bar2(y1, y2, y)
  if ((y <= ym + profile1) and (y >= ym - profile1)) or ((x <= xm + profile2)  and (x >= xm - profile2)):
    return True
  return False

def set_structure_array(nx,ny,x1,x2,y1,y2):
  array = np.zeros((nx,ny))
  sizex = 1.0
  sizey = 1.0
  hx = sizex / nx
  hy = sizex / ny
  for j in range(0,nx):
    for i in range(0,ny):
      coords = np.zeros(2)
      x = j * hx
      y = i * hy
      coords[0] = x
      coords[1] = y
      if within_profile(x1, x2, y1, y2, coords):
        array[i,j] = 1
  return array


x1 = 0.8
x2 = 0.3
y1 = 0.5
y2 = 0.3
nx = 1000
ny = 1000


array = set_structure_array(nx, ny, x1, x2, y1, y2)
mesh = create_2d_mesh_from_array(array)
write_gid_mesh(mesh, "/home/bvu/TEST/lufo_mesh/test.mesh")
# img = Image.fromarray(array)
# img.show()
  
  
  
