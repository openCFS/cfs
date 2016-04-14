#!/usr/bin/env python

# This script generates initial density distributions as pseudo density.xml files
# This is used when doing inverse homogenization and bloch mode optimization.

import libxml2
from numpy import *
import math
from optimization_tools import *
from mesh_tool import *
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

# make a shpere in data with from center to radius with gradient from inner to outer value
# data is a array from numpy (numpy.ndarray) in 2 or 3 dim
# inner_value is assumed to be smaller outer_value
# returns the volume
def make_sphere(dim, divider, radius, inner_value, outer_value, order, invert):
  
  center = Coordinate(0.5, 0.5, 0.5)
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
          if order == 'binary':
            setNDArrayEntry(data, dim, i, j, k, inner_value)  
          else:
            # linear: v = (d / radius) * gap + inner_value
            ratio = 1.0
            for r in range(1, int(order)+1):
              ratio = ratio * (d / radius)              
            v = ratio * gap + inner_value # scale down by outer_value - inner_value and shift
            old_value = getNDArrayEntry(data, dim, i, j, k)
            # we make only make smaller otherwise we would destroy multiple spheres stuff
            if v < old_value:
              setNDArrayEntry(data, dim, i, j, k, v) 

  if invert:
   data = outer_value + inner_value - data # don't make a 0 where we had 1  

  assert(numpy.amax(data) <= outer_value)
  assert(invert or numpy.amax(data) >= inner_value)

  return data

          
# find correct radius by bisection
# vol the desired resulting vol
# return the data  as numpy.ndarray
def find_radius(dim, div, vol, order, invert, lower_val):

  # set the center coordinates
  lower = 0.0
  upper = 30 # 1.4
  err = upper
  
  data = 0 # placeholder
  iter = 0
  
  while iter < 30 and abs(err) > 1e-12:
    mid = 0.5 * (lower + upper)
    
    data = make_sphere(dim, div, mid, lower_val, 1.0, order, invert)
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

def cross(dim, vol, res, lower):
  assert(dim == 2)
  
  # v = 2*h - h^2 
  h = 1.0-numpy.sqrt(1-vol) # the other solution is > 1
  assert(h >= 1/res and h <= 1)
  s = h * res

  print 'cross bar thickness is ' + str(h * 100) + "% which is makes " + str(int(s)) + " cells"

  data = numpy.ones((res, res)) * lower # violate exact volume
  
  start = int(res/2. - s/2. + 0.5)
  end   = int(res/2. + s/2. + 0.5)
  
  data[:,start : end   ] = 1
  data[  start : end, :] = 1
  
  return data


def rectangle(dim, vol, res, lower):
  assert(dim == 2)
  
  # v = 2*h - h^2 
  h = 1.0-numpy.sqrt(1-vol) # the other solution is > 1
  assert(h >= 1/res and h <= 1)
  s = h * res

  print 'cross bar thickness is ' + str(h * 100) + "% which is makes " + str(int(s)) + " cells"

  data = numpy.ones((res, res)) * lower # violate exact volume
  
  start = int(res/2. - s/2. + 0.5)
  end   = int(res/2. + s/2. + 0.5)
  
  data[:,start : end   ] = 1
  data[  start : end, :] = 1
  
  return data



def hashtag(dim, res, amplitude, thickness, speed, lower):
  
  data = numpy.ones((res, res)) if dim == 2 else numpy.ones((res,res,res)) 
  data *= lower # violate exact volume
  
  h = 1.0/res
  
  for x in range(res):
    for y in range(res):
      if dim == 2:  
        if hashtag_dist_2d(h*x, h*y, amplitude, speed) < thickness/2: 
          data[x,y] = 1
      else:    
        for z in range(res):
          if hashtag_dist_3d(h*x, h*y, h*z, amplitude, speed) < thickness/2: 
            data[x,y, z] = 1
 

  return data;

## helper for hashtag. gives for (x,y) the closests distance but only horizontally!
def hashtag_dist_2d(x, y, amplitude, speed):
  #  0.1*sin(2*x*pi+pi/2) + 0.25, 0.25, -0.1*sin(2*x*pi+pi/2) + 0.75, 0.75
  y1 = amplitude * sin(2*speed*pi*x+pi/2) + 0.25
  y2 = -amplitude * sin(2*speed*pi*x+pi/2) + 0.75
  #print "x=" + str(x) + " y=" + str(y) + " y1=" + str(y1) + " y2=" + str(y2)
  y_dist = min(abs(y-y1), abs(y-y2))
  
  inv = -1 if speed == 1 else 1 
  
  x1 = inv * amplitude * sin(2*speed*pi*y+pi/2) + 0.25
  x2 = -inv * amplitude * sin(2*speed*pi*y+pi/2) + 0.75
  x_dist = min(abs(x-x1), abs(x-x2))
  
  return min(y_dist, x_dist)

def hashtag_dist_3d(x, y, z, amplitude, speed):
 
  # we map the 2d-solution twice on the z-axis to 0.25 and 0.75 within the xy-plane
  z_dist = min(abs(z-0.25), abs(z-0.75))
  hash_dist = hashtag_dist_2d(x, y, amplitude, speed)
  xy_dist = max(z_dist, hash_dist) # max as min would show the two z-planes
  
  # we use hashtag_dist_2d and apply it on the yz-plane
  hash_dist = hashtag_dist_2d(y, z, amplitude, speed)
  x_dist = min(abs(x-0.25), abs(x-0.75))
  yz_dist = max(x_dist, hash_dist)
  
  return min(xy_dist, yz_dist)
  
parser = argparse.ArgumentParser()
parser.add_argument("--res", help="edge discretization of length 1m", type=int, required = True )
parser.add_argument('--dim', help="square (2) or cube (3)", type=int, default=2)
parser.add_argument('--lower', help="value for void material. Default 1e-3", type=float, default=1e-3)
parser.add_argument('--vol', help="volume fraction of full domain or ball only", type=float, default=0.5)
parser.add_argument('--order', help="order of generated shperes. Lower numbers are smoother. 'binary' for black and white", default="6")
parser.add_argument('--invert', help="invert to solid inside", action='store_true')
parser.add_argument('--cross', help="make a simple binary cross", action='store_true')
parser.add_argument('--rect', help="make a simple binary rectangle inclusion", action='store_true')
parser.add_argument('--hashtag', help="hashtag # based on sin-amplitude for bloch mode initial designs [0,1]", type=float)
parser.add_argument('--thickness', help="feature thickness for hashtag", type=float, default=0.1) 
parser.add_argument('--hashtag_speed', help="number of maximas, only 1,2,4, ... make sense", type=int, default=1)
parser.add_argument('--ball', help="account vol only on the inner ball with diameter 1.0", action='store_true')
parser.add_argument('--show', help="additionaly visualize the image", action='store_true')
parser.add_argument('--save', help="overwrite default filename, when it ends with an image extension the image is written")
parser.add_argument('--write_mesh', help="optionally create a sparse mesh. For more options use process_image.py", action='store_true')

# parser.add_argument('--elem_nr', help="for debug purpose only (rotation). Ignore vol, dim, order, invert and give the design the 1-based element number", action='store_true')
args = parser.parse_args()

 
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
setname = "standard"

if args.cross:
  data = cross(args.dim, vol, args.res, args.lower)
  filename = "cross_" + str(args.dim) + "d-v_" + str(args.vol) + "_" + str(args.res) + ".density.xml"
elif args.hashtag is not None: # also capture 0.0
  data = hashtag(args.dim, args.res, args.hashtag, args.thickness, args.hashtag_speed, args.lower)
  filename = "hashtag_" + str(args.dim) + "d-amp_" + str(args.hashtag) + "-th_" + str(args.thickness) + "-sp_" + str(args.hashtag_speed) + "_" + str(args.res) + ".density.xml"
else:
  data = find_radius(args.dim, args.res, vol, args.order, args.invert, args.lower)
  filename = "circular_" + str(args.dim) + "d-v_" + str(args.vol) + ("_ball" if args.ball else "") + ord  + ("-inv_" if args.invert else "_") + str(args.res) + ".density.xml"
  setname = "order_" + str(args.order) + ("_inv" if args.invert else "")

if args.save:
  filename = args.save
if filename.endswith('.png') or filename.endswith('.jpg') or filename.endswith('.jpeg') or filename.endswith('.gif') or filename.endswith('.tif'):
  get_image(data).save(filename)
  print "generated image '" + filename + "'"
elif not args.write_mesh:
  write_density_file(filename, data, setname)
  print "generated density file '" + filename + "'" 

if args.write_mesh:
  mesh = Mesh()  
  create_dense_mesh_density(data, mesh, threshold=0.5, scale=1.0, rhomin = 1e-3) # rhomin is irrelevant as we make sparse
  sparse = convert_to_sparse_mesh(mesh)
  mesh_name = filename.replace('.density.xml', '.mesh')
  write_gid_mesh(sparse, mesh_name)
  print "generated sparse mesh '" + mesh_name + "'"
  
if args.show:
  get_image(data, 800).show()

