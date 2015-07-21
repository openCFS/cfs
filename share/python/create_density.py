#!/usr/bin/env python

# This script generates initial density distributions as pseudo density.xml files
# This is used when doing inverse homogenization and bloch mode optimization.

import libxml2
import numpy
import math
from optimization_tools import *
import argparse


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
def make_sphere(dim, center, radius, inner_value, outer_value, order, invert):
  
  point = Coordinate(0.0, 0.0, 0.0)
  
  data = numpy.ones((divider, divider)) if dim == 2 else numpy.ones((divider, divider, divider))
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

  if invert:
   data = outer_value + inner_value - data # don't make a 0 where we had 1  

  assert(numpy.amax(data) <= outer_value)
  assert(numpy.amax(data) >= inner_value)

  return data

          
# find correct radius by bisection
# vol the desired resulting vol
# return the data  as numpy.ndarray
def find_radius(dim, vol, order, invert):

  # set the center coordinates
  center = Coordinate(0.5, 0.5, 0.5)
        
  lower = 0.0
  upper = 30 # 1.4
  err = upper
  
  data = 0 # placeholder
  iter = 0
  
  while iter < 30 and abs(err) > 1e-12:
    mid = 0.5 * (lower + upper)
    
    data = make_sphere(dim, center, mid, 0.002, 1.0, order, invert)
    act_vol = data.sum() / float(data.size)
    
    err  = vol - act_vol
 
    if (not invert and err > 0) or (invert and err < 0):
      upper = mid
    else:
      lower = mid

    iter = iter + 1
    #print "     act_vol=" + str(act_vol) + " err=" + str(err) + " mid=" + str(mid) + " next lower=" + str(lower) + " next upper=" + str(upper) 
  
  # we are so close that left and right data is almost the same
  print "dim=" + str(dim) + " order=" + str(order) + " target_vol=" + str(vol) + " result_vol=" + str(data.sum() / float(data.size)) + ' min=' + str(numpy.amin(data)) + ' max=' + str(numpy.amax(data)) 
  return data

def cross(dim, vol, res):
  assert(dim == 2)
  
  # v = 2*h - h^2 
  h = 1.0-numpy.sqrt(1-vol) # the other solution is > 1
  assert(h >= 1/res and h <= 1)
  s = h * res

  print 'cross bar thickness is ' + str(h * 100) + "% which is makes " + str(int(s)) + " cells"

  data = numpy.ones((res, res)) * 1e-3 # violate exact volume
  
  start = int(res/2. - s/2. + 0.5)
  end   = int(res/2. + s/2. + 0.5)
  
  print start
  print end

  
  data[:,start : end   ] = 1
  data[  start : end, :] = 1
  
  
  
  return data

parser = argparse.ArgumentParser()
parser.add_argument("--res", help="edge discretization of length 1m", type=int, required = True )
parser.add_argument('--vol', help="volume fraction of full domain or ball only", type=float, default=0.5)
parser.add_argument('--dim', help="square (2) or cube (3)", type=int, default=2)
parser.add_argument('--order', help="order of generated shperes. Lower numbers are smoother", type=int, default=6)
parser.add_argument('--invert', help="invert to solid inside", action='store_true')
parser.add_argument('--cross', help="make a simple cross", action='store_true')
parser.add_argument('--ball', help="account vol only on the inner ball with diameter 1.0", action='store_true')

# parser.add_argument('--elem_nr', help="for debug purpose only (rotation). Ignore vol, dim, order, invert and give the design the 1-based element number", action='store_true')
parser.add_argument('--save', help="overwrite default filename")

args = parser.parse_args()

divider = args.res
  
vol = args.vol  
if args.ball:
  if args.dim == 2:
    vol *= 1.0/4.0 * numpy.pi 
  else:
    vol *= 1.0/6.0 * numpy.pi
  print "assign volume " + str(args.vol) + " to ball which restricts the total volume to " + str(vol)    
    
  
ord = ("-o_" + str(args.order)) if args.order <> 6 else ""   


data = None
filename = None

if args.cross:
  data = cross(args.dim, args.vol, args.res)
  filename = "cross_" + str(args.dim) + "d-v_" + str(args.vol) + "_" + str(divider) + ".density.xml"
else:
  data = find_radius(args.dim, vol, args.order, args.invert)
  filename = "circular_" + str(args.dim) + "d-v_" + str(args.vol) + ord  + ("-inv_" if args.invert else "_") + str(divider) + ".density.xml"

if args.save:
  filename = args.save

write_density_file(filename, data, "order_" + str(args.order) + ("_inv" if args.invert else ""))

print "generated file '" + filename + "'" 


