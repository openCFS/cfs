#!/usr/bin/env python

# This script generates initial density distributions as peseudo density.xml files
# This is used when doing inverse homogenization.

import libxml2
import numpy
import math
from optimization_tools import *

# we assume a unit cube (2D/3D) 
# edge discretization
divider = 5 
vol_list = ["0.6"]
dim = 2 
# what is the maximal order (1 is linar, 2 quadratic, ...)
order = 6 
# the maximal number of spheres as edge basis (1*1, 2*2, 3*3, .., or 1*1*1, 2*2*2, 3*3.3, ...)
edge = 1




   

# there shall be a predefined class somewhere, I just didn't find it
class Coordinate:
  def __init__(self):
    self.x = 0.0
    self.y = 0.0 
    self.z = 0.0 

  def __init__(self, x, y, z):
    self.x = x
    self.y = y 
    self.z = z 
    
  # i, j, k ints from 0 to div-1  
  def toCoordinate(self, i, j, k, div):
    # shift the coordinate to the center of the elements
    self.x = i / float(div) + 0.5 / div
    self.y = j / float(div) + 0.5 / div 
    self.z = k / float(div) + 0.5 / div 

  # other is also Coordinate
  def dist(self, other):
    return math.sqrt((self.x - other.x)**2 + (self.y - other.y)**2 + (self.z - other.z)**2)  
 
  def printline(self):
    print str(self.x) + ", " + str(self.y) + ", " + str(self.z)
    
  def toString(self):
    return str(self.x) + ", " + str(self.y) + ", " + str(self.z) 


# helper: return numpy.ndarray element with 2/3D tolerance
# if dim = 3 the k entry is ignored
def getNDArrayEntry(data, dim, i, j, k):
  if dim == 3:
    return data[i,j,k]
  if dim == 2:
    return data[i, j]
  raise " cannot handle dimension " + str(dim)
    
# see getNDArrayEntry(data, dim, i, j, k)
def setNDArrayEntry(data, dim, i, j, k, value):
  if dim == 3:
    data[i,j,k] = value
    return
  if dim == 2:
    data[i, j] = value
    return
  raise " cannot handle dimension " + str(dim)

# make a shpere in data with from center to radius with linear gradient from inner to outer value
# data is a array from numpy (numpy.ndarray) in 2 or 3 dim
# center is a Coordinate
# inner_value is assumed to be smaller outer_value
# returns the volume
def make_sphere(data, center, radius, inner_value, outer_value, order):
  point = Coordinate(0.0, 0.0, 0.0)
  dim = data.ndim
  gap = outer_value - inner_value
  # loop over the whole data
  for i in range(divider):
    for j in range(divider):
      for k in range(cond(dim == 2, 1, divider)):
        point.toCoordinate(i, j, k, divider)       
        d = center.dist(point)
        if d < radius:
          # linear: v = (d / radius) * gap + inner_value
          ratio = 1.0
          for r in range(1, order+1):
            ratio = ratio * (d / radius)              
          v = ratio * gap + inner_value # scale down by outer_value - inner_value and shift
          old_value = getNDArrayEntry(data, dim, i, j, k)
          # we make only make smaller otherwise we would destroy multiple spheres stuff
          if v < old_value:
            setNDArrayEntry(data, dim, i, j, k, v) 

  return data.sum() / float(data.size)


# wrap the calls to make_shere by calling it for each of the center coordinate within grid
# grid multiple centers as list of Coordinate
# return data 
def multiple_make_sphere(dim, divider, grid, radius, inner_value, outer_value, order):


  data = cond(dim == 2, numpy.ones((divider, divider)), numpy.ones((divider, divider, divider)))
  
  for i in range(len(grid)):
    point = grid[i]       
    make_sphere(data, point, radius, inner_value, outer_value, order)
  
  return data
           
# find correct radius by bisection
# vol the desired resulting vol
# spheres number of shperes
# return the data  as numpy.ndarray
def find_radius(dim, vol, spheres, order):

  # set the center coordinates
  grid = [] 
  for i in range(spheres):
    for j in range(spheres):  
      for k in range(cond(dim == 2, 1, spheres)):
        grid.append(Coordinate(0.0, 0.0, 0.0))
        grid[len(grid)-1].toCoordinate(i, j, k, spheres)
        
  lower = 0.0
  upper = 30.0
  err = upper
  
  data = 0 # placeholder
  iter = 0
  
  while iter < 30 and abs(err) > 0.01:
    mid = 0.5 * (lower + upper)

    data = multiple_make_sphere(dim, divider, grid, mid, 0.002, 1.0, order)
    act_vol  = data.sum() / float(data.size)  
    err  = vol - act_vol

    if err > 0:
      upper = mid
    else:
      lower = mid

    iter = iter + 1
    print "     act_vol=" + str(act_vol) + " err=" + str(err) + " mid=" + str(mid) + " next lower=" + str(lower) + " next upper=" + str(upper) 
  
  # we are so close that left and right data is almost the same
  print "dim=" + str(dim) + " order=" + str(order) + " spheres=" + str(spheres) + " target_vol=" + str(vol) + " result_vol=" + str(data.sum() / float(data.size))
  return data

# write the data to a density.xml file
# list of data and list of setnames  
def write_density_file(filename, data_list, setname_list):
  out = open(filename, "w")
  out.write('<?xml version="1.0"?>\n')
  out.write('<cfsErsatzMaterial>\n')
  out.write('  <header>\n')
  out.write('    <design initial="0.5" lower="1e-3" name="density" region="mech" upper="1"/>\n')
  out.write('    <transferFunction application="mech" design="density" param="1" type="simp"/>\n')
  out.write('  </header>\n')
  
  for i in range(len(data_list)):
    data    = data_list[i]
    setname = setname_list[i]
    out.write('  <set id="' + setname + '">\n')
    edge = data.shape[0] # be careful!
    counter = 1

    dim = data.ndim 

    for i in range(edge):
      for j in range(edge):
        for k in range(cond(dim == 2, 1, edge)):
           val = getNDArrayEntry(data, dim, i, j, k)
           out.write('    <element nr="' + str(counter) + '" type="density" design="' + str(val) + '"/>\n')
           counter = counter + 1       
         
    out.write('  </set>\n')

  out.write(' </cfsErsatzMaterial>\n')
  out.close()
  
# main routine  
# dim either 2 or 3
# create a density file with several sets:
# divider edge discretization of mesh
# vol target volume fraction
# max_edge create sets 1*1*1, 2*2*2, 3*3*3, ...n*n*n with n = max_edge
# max_order the density gradient is linear(1), quadratic(2), cubic(3), ...
def create_density_file(dim, divider,  vol, max_edge, max_order):
  filename = "unit_" + str(dim) + "d-div_" + str(divider) + "-vol_" + str(vol) + ".density.xml"
  
  data_list = []
  setname_list = []
   
  for o in range(1, max_order+1): 
    for e in range(1, max_edge+1):
      data = find_radius(dim, vol, e, o)
      data_list.append(data)
      setname_list.append("spheres_n_" + str(e) + "-order_" + str(o))   

  write_density_file(filename, data_list, setname_list) 
  
for vol in vol_list:  
  create_density_file(dim, divider, float(vol), edge, order)  


