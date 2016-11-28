#!/usr/bin/env python
import numpy as np
from numpy import outer
import sympy.solvers
#from PIL import Image
# from mesh_tool import *
import sys
import math
import vtk
from matviz_vtk import *
from matplotlib import pyplot as plt
from scipy import interpolate
from mpl_toolkits.mplot3d import Axes3D
from scipy.spatial import Delaunay
from sympy.physics.quantum.circuitplot import pyplot
from sympy import Symbol, symbols
# from basecell import calc_radius

class End_Node():
  def __init__(self):
    self.coords = (0.0,0.0,0.0)
    self.connections = []
    self.id = -1 # vtk id of point
    self.i = -1 # first index in surface map, corresponds to surface line
    self.j = -1 # second index in surface map, where on surface line 
    self.dir = -1
    
  def __init__(self,coords,id,dir,i,j):
    assert(len(coords) == 3)
    self.coords = coords
    self.id = id
    self.i = i
    self.j = j
    self.connections = [] # array of end_nodes that form a triangle together with this end node
    self.dir = dir
    
  def __str__(self):
    return "id=" + str(self.id) + " coords:" + str(self.coords) + " i=" + str(self.i) + " j=" + str(self.j) + " dir=" + str(self.dir)  
  
class Cubic_spline():
  # assume we have u_0=u_1=u_2=u_3=0 and u_4=u_5=u_6=u_7=0
  # a cubic spline is defined by its base functions and control polygon
  #f03 = (1-t)**3
  #f13 = 3*t*(1-t)**2
  #f23 = 3*t**2*(1-t)
  #f33 = t**3"
   
  bend = None
   
  #control polygon contains 4 ponts (array of lists, 1 list per point) 
  CP = None
   
  def __init__(self):
    print "Default instantiation not allowed for cubic splines!"
    sys.exit(1)
   
  # CP is numpy array with 4 coordinates; for each coordinate use a list with x- and y-component
  def __init__(self, CP = None):
    self.CP = np.transpose(CP) if CP is not None else None
    t = np.linspace(0,1,1000) # over-sampling
    C = self.eval_t(t)
    self.explicit = interpolate.interp1d(C[0,:],C[1,:])
      
  # base functions for cubic spline with 4 control points
  def f03(self,t):
    return (1-t)**3
  def f13(self,t):
    return 3*t*(1-t)**2
  def f23(self,t):
    return 3*t**2*(1-t)
  def f33(self,t):
    return t**3
    
  # evaluates spline for parameter t
  def eval_t(self,t):
    return outer(self.CP[:,0],self.f03(t)) + outer(self.CP[:,1],self.f13(t)) + outer(self.CP[:,2],self.f23(t)) + outer(self.CP[:,3],self.f33(t))
  
  # evaluates spline for given x
  # use interpolation to obtain explicit spline representation
  def eval_x(self,x):
#     assert( x.all() >= 0 and x.all() <= 1)
    # interpolate returns ndarray, need to convert it to float
    return self.explicit(x)[()]
  
  def calc_d_spline_d_t(self,t):
    return outer(3*(1-t)**2,self.CP[:,1]-self.CP[:,0]) + outer(6.0*t*(1-t), self.CP[:,2] - self.CP[:,1]) + outer(3*t**2 ,self.CP[:,3]-self.CP[:,2])
  
  def calc_param_grad_1(self):
    u = Symbol('u')
    
    dC = self.calc_d_spline_d_t(u) # dC/dt
    sol = sympy.solvers.solve(dC[0][1]-dC[0][0],u)
    t = -100
    if sol[0] > 0 and sol[0] <= 1:
      t = sol[0]
    elif  sol[1] > 0 and sol[1] <= 1:
      t = sol[1]
    else:
      print "No t found where dx/dy = 1"
     
    # conversion from t as sympy.Float to regular Python float necessary
    return float(t)
  
  def plot(self):
    t = np.linspace(0, 1, 100)
    C = self.eval_t(t)
    plt.plot(np.transpose(self.CP[0,:]),np.transpose(self.CP[1,:]))
#     t1 = self.calc_param_grad_1()
#     Ct = self.eval(t1)
    plt.plot(np.transpose(C[0,:]),np.transpose(C[1,:]))
#     plt.plot(Ct[0],Ct[1],marker='o', markersize=15)
    plt.xlim((0,0.5))
    plt.ylim((0.5,1.0))
    plt.show()
    
  def calc_coords_grad_1(self):
    t1 = self.calc_param_grad_1()
    return self.eval_t(t1)
  
  def calc_min(self):
    dC = self.calc_d_spline_d_t(u)
    
def dirToString(dir):
  assert(dir > 0 and dir < 4)
  res = 'x'
  if dir == 2:
    res = 'y'
  if dir == 3:
    res = 'z'
   
  return res

# converts angle (0 to 360) to radians and maps it to a value between 0 and pi/2.0 
def degree_to_rad_quadrant(degree):
  assert(degree >= 0.0 and degree <= 360.0)
  rad = np.pi/180. * degree
  if degree <= 90:
    return rad
  elif degree <= 180:
    return np.pi-rad
  elif degree <= 270: 
    return rad-np.pi
  else:
    return 2*np.pi-rad
  
def cartesian_to_grid_coords(x,res):
  h = 1.0 / res # assume domain is 1m x 1m x 1m
  return int((x) / h) 

def grid_to_cartesian_coords(i,res):
  h = 1.0 / res # assume domain is 1m x 1m x 1m
  return i * h + h / 2.0

# @param offset: value added to converted cartesian
def polar_to_cartesian(radius,radians,offset):
  return radius * np.cos(radians) + offset, radius * np.sin(radians) + offset

# calculates euklidian distance to origin (0,0)
# @param p: tuple with x-,y-component of point
def distance_to_center(p):
  return np.sqrt((p[0]-0.5)**2+(p[1]-0.5)**2)

# calculates angle between (0.5,1.0) and point p
def angle_to_center(p):
  if p[1] < 0.5:
    z = 1.0 - p[1]
  else:
    z = p[1]

  # in case p is exactly (0.5,0.5), angle is 0 and thus denominator is 0!
  if abs(p[0]- 0.5) <= 1e-3 and abs(p[1]- 0.5) <= 1e-3:
    return 0  
    
  return np.arccos(0.5*(z-0.5)/(0.5*np.sqrt((p[0]-0.5)**2 + (z-0.5)**2)))
 
# defines a 1D linear function
class Linear_1D():
  x1 = None
  x2 = None
  # @param x1 and x2 define the function  
  def __init__(self,x1,x2):
    self.x1 = x1
    self.x2 = x2
    
  def eval(self,x):
    return ((self.x2-self.x1) * x + self.x1) / 2.0 + 0.5
    
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

# tests if profile contains probe
# @param y: at which level of other profile do we expect the point
# e.g. for x-profile it is the y-coordinate
# @param probe: 2D coordinate of point for testing
# def contains_point(x,y,z,profile):
#   assert(x >= 0 and x <= 1)
#   assert(y >= 0 and y <= 1)
#   assert(z >= 0 and z <= 1)
#      
#   phi = angle_to_center((x,z))
#   print "(x,z): ",x,z, " angle:", degrees(phi)
#   r = calc_radius_for_quadrant(profile, y, (phi))
#   print "radius: ", r
#   val = (x-0.5)**2 + (z-0.5)**2
#     
#   print "phi:", degrees(phi)
#   print "val: ", val, " r**2: ", r**2
#      
# #   print "radius: ",r," val:",val
#      
#   if (val - r*r < 0):
#     return True
#      
#   return False

def contains_point(id_x,y,z,map):
     
  valx = (y-0.5)**2 + (z-0.5)**2
  phi = angle_to_center((y,z))
   
  r = map[int(degrees(phi)),id_x]
     
  if (valx - r*r < 1e-3 ):
    return True
     
  return False
    

# obejct defines spline in a principle plane, e.g. spline for 0 or 90 degree
class PrincipleSpline():
  
  def __init__(self):
    self.spline = None
    self.angle = -1
    self.type = "bspline"
  
  def __init__(self, x1, y1, bend, angle=0):
    rx = 0.5 - x1/2.0 # radius for center (0,1)
    ry = 0.5 - y1/2.0 # radius for center (0,1)
  
    P = np.array([[0,1-rx],[ry*bend,1-rx],[ry,1-rx*bend],[ry,1]])
    self.spline = Cubic_spline(P)
    self.angle = angle
    # coordinate where slope is 1
    self.coords_cut = self.calc_coords_grad_1()
    self.type = "bspline"
    
  # eval function in case x is one element, not a list  
  def eval_elem(self,x):
    if x <= self.coords_cut[0]:
      val = self.spline.eval_x(x)
    elif x > self.coords_cut[0] and x < 1 - self.coords_cut[0]:
      val = self.coords_cut[1]
    else: #x >= 1- coord_cut[0]
      val = self.spline.eval_x(1-x)
    
    return val
  
  # wrapper function
  #@param x: can be one argument value or list or aguments
  def eval(self,x):
    if isinstance(x, (float,int,long)):
      return self.eval_elem(x)
    
    # in case x is a list    
    res = []
    
    for i in x:
      res.append(self.eval_elem(i))

    return res
  
  def calc_coords_grad_1(self):
    return self.spline.calc_coords_grad_1()
  
# object defining spline for angles between 0,...,90 degree
class BisecSpline:
  def __init__(self):
    self.bicubic = [] # 0-th entry is spline function and 1st entry is cubic polynomial
    self.cubic = None
    self.spline = None
    self.linear = None
    self.x1 = 0
    self.y1 = 0
    self.z1 = 0
    self.type = None
    self.angle = None
  
  def __init__(self,x1,y1,z1,bend):
    self.x1 = x1
    self.y1 = y1
    self.z1 = z1
    self.bicubic = [] # 0-th entry is spline function and 1st entry is cubic polynomial
  
    assert(bend <= 1 and bend >= 0)
    self.type = None
    
    ###### case 1 : bisec spline + cubic polynomial --> bicubic ###################
    # we have a left part from a=(0,x1) which is the average of the splines x1,y1 and x1,z1
    # where the curve has grad=1 we have point b
    # the right part is from b to x=0.5 where the heigt comes from the spline y1,z1 with grad=1 at point p
    # from the point p we determine the angle phi of this bisec Profile function 
  
    # left part is same spline as orhtogonal spline up to point
    left = PrincipleSpline(x1, 0.5*(y1+z1), bend)
    b = left.coords_cut 
    self.cut = b 
    assert(b[0] <= 0.5 and b[1] >= 0.5)

    # search for point p
    right = PrincipleSpline(y1, z1, bend)
    p = right.coords_cut
  
    assert(p[0] <= 0.5 and p[1] >= 0.5)
    height = distance_to_center(p) + 0.5
  
    self.angle = angle_to_center(p)
#     self.angle = angle_to_center(p)
    
    # polynomial interpolation for right part from b to p
    lx = b[0]
    ly = b[1] 
    rx = 0.5
    A = np.array([
        [1, lx, lx**2, lx**3],
        [0, 1,  2*lx, 3*lx**2],
        [1, rx, rx**2, rx**3],
        [0, 1,  2*rx, 3*rx**2]
        ]) 
    rhs = np.array([ly, 1.0, height, 0])
    
    sol = np.linalg.solve(A, rhs)
  
    cubic = np.poly1d(sol[::-1]) # cubic polynomial
    
    #### case 2: b-spline --> bsp #############
    # if b with grad gb (approx 1) is too high for p such that the curve has a maximum within b and p we need to fallback to a b-spline from a to p
    # curve as undershoot
  
    P = np.array([[0,0.5+x1/2.0],[0.5*bend,0.5+x1/2.0],[0.5-0.5*bend,height],[0.5,height]])
    bspline = Cubic_spline(P) 
  
    #### case 3: linear --> lin ###########
    # in case undershooting for x1=0.9, y1=0.1, z1=0.1
    lin = Linear_1D(x1, x1)
    
    self.bicubic.append(left)
    self.bicubic.append(cubic)  
    self.spline = bspline
    self.linear = lin
    # to check if bicubic has under/overshooting when point p is much lower than point b
    if height >= b[1] + 1e-3:
      self.type = "bicubic"
      
    # in case function composed of b-spline and cubic function has undershoot  
    # in case b-spline has no undershoot (point p is not below bspline(x=0))
    elif height > 0.5 + x1/2.0:  
      self.type = "bSpline"
      
    # in case we have undershooting for biqua AND for spline
    else:
      self.type = "linear"
  
  def get_coords_cut(self):
    return self.cut
      
  def eval_spline(self,x):
#     assert( x.all() >= 0 and x.all() <= 1)
    # need to interpolate to assure equidistant spacing for x \in [0,1]
    res = []
    
    # make x iterable in case it's one float element and not a list
    x = np.reshape(x, np.size(x), )
    for i in x:
      assert(i >= 0 and i <=1)
      if i <= 0.5:
        val = self.spline.eval_x(i)
      else:
        val = self.spline.eval_x(1-i)
      
      assert(val is not None)
      
      res.append(val)  
    
    return res
  
  def eval_bicubic(self,x):
    # coordinate at which slope is 1
    left = self.bicubic[0]
    cubic = self.bicubic[1]
    res = []
    
    # make x iterable in case it's one float element and not a list
    x = np.reshape(x, np.size(x), )
    
    for i in x:
      assert(i >= 0 and i <=1)
      if i <= 0.5:
        if i <= left.coords_cut[0]:
          val = left.eval(i)
        else:
          val = cubic(i)
      else:
        if i <= 1.0-left.coords_cut[0]: # mirror of cubic 
         val = cubic(1-i)
        else:
         val = left.eval(1-i)
      
      res.append(val)
    return res
  
  def eval_linear(self,x):
    return self.linear.eval(x)
  
  def eval(self,x):
    if self.type == "bicubic":
      return self.eval_bicubic(x)
    elif self.type == "bSpline":
      return self.eval_spline(x)
    else: #linear case
      return self.eval_linear(x)
    
  def get_type(self):
    return self.type
  
  def plot(self):
    plt.gcf().clear()
    
    x = np.linspace(0, 1, 100)
    
    if self.type == "bicubic":
      bicubic = self.eval_bicubic(x)
      plt.plot(x,bicubic,label='bicubic',linewidth=5.0)
    if self.type == "bSpline":
      spline = self.eval_spline(x)
      plt.plot(x,spline,label='spline',linewidth=5.0)
    else:  
      linear = self.eval_linear(x)
      plt.plot(x,linear,label='linear',linewidth=5.0) 
      
    plt.legend(loc='upper left', shadow=True,prop={'size':20})
    plt.show()
    
  def plot_all(self):
    plt.gcf().clear()
    
    x = np.linspace(0, 1, 100)
    
    bicubic = self.eval_bicubic(x)
    spline = self.eval_spline(x)
    linear = self.eval_linear(x)
    
    cut = self.bicubic[0].coords_cut

    plt.plot(x,bicubic,label='bicubic',linewidth=5.0)
    plt.plot(x,spline,label='spline',linewidth=5.0)
    plt.plot(x,linear,label='linear',linewidth=5.0)
    plt.plot(cut[0],cut[1],marker='o',color='red',markersize=15) 
      
    plt.legend(loc='upper left', shadow=True,prop={'size':20})
    plt.show()  
    
  def angle(self):
    return self.angle
    
# @return vector with Profile or list of vectors
class Profile:
  def  __init__(self):
    self.bisec_angle = -1
    self.type = None
    # 0th entry: function for 0 degree; 1st entry: function for bisec; 2nd entry: function for 90 degree
    self.functions = [None] * 3
    self.direction = 0
    
  def __init__(self, args, dir):
    assert(args.profile == "linear" or args.profile == "spline")
    assert (dir == 1 or dir == 2 or dir == 3)
    self.direction = dir
    self.type = args.profile
    self.functions = [None] * 3
    if self.type == "linear":
      self.functions[0] = Linear_1D(args.x1, args.x2)
      self.functions[1] = self.splines[0]
      self.functions[2] = self.splines[0]
    else: # 'spline' case
      if dir == 1:
        self.functions[0] = PrincipleSpline(args.x1, args.y1, args.bend, 0)
        self.functions[1] = BisecSpline(args.x1, args.y1, args.z1, args.bend)
        self.functions[2] = PrincipleSpline(args.x1, args.z1, args.bend, np.pi/2.0)
      elif dir == 2:
        self.functions[0] = PrincipleSpline(args.y1, args.x1, args.bend, 0)
        self.functions[1] = BisecSpline(args.y1, args.x1, args.z1, args.bend)  
        self.functions[2] = PrincipleSpline(args.y1, args.z1, args.bend, np.pi/2.0)
      else:
        self.functions[0] = PrincipleSpline(args.z1, args.y1, args.bend, 0)
        self.functions[1] = BisecSpline(args.z1, args.y1, args.x1, args.bend)  
        self.functions[2] = PrincipleSpline(args.z1, args.x1, args.bend, np.pi/2.0)
        
      self.bisec_angle = self.functions[1].angle
      
    if args.verbose == "bisec":
      self.functions[1].plot_all()
      
# return information on profiles 
def create_profiles(args,infoXml=None):
  profiles = [None]*3 # x-,y-,z-part
  
  if not args.skip_x:
    profiles[0] = Profile(args,1)
    
  if not args.skip_y:
    profiles[1] = Profile(args,2)
    
  if not args.skip_z:
    profiles[2] = Profile(args,3)
    x = np.linspace(0, 1.0, args.res)
    
  if args.verbose == "all_profiles":
    plt.gcf().clear()
    plt.gcf().subplots_adjust(bottom=0.15)
    x = np.linspace(0, 1.0, 1000)
    
    for dir,profile in enumerate(profiles):
      if profile == None:
        continue
      plt.plot(x,profile.functions[0].eval(x),linewidth=5.0,label="dir_"+str(dir+1)+"_0")
      plt.plot(x,profile.functions[1].eval(x),linewidth=5.0,label="dir_"+str(dir+1)+"_"+str(profile.functions[1].angle[0]))
      plt.plot(x,profile.functions[2].eval(x),linewidth=5.0,label="dir_"+str(dir+1)+"_90")
    
    plt.legend(loc='upper left', shadow=True)
    plt.rcParams.update({'font.size': 18})
    #plt.ylim([0,0.5])
    plt.xlabel("x",labelpad=5)
#     plt.ylabel("radius",labelpad=20)
#     plt.savefig("profile_functions_" + str(args.stiffness) + ".png")

    plt.show()
  
  return profiles

# @param map: Profile map
# @param numLines: number of lines on surface that we want to plot
# @param dir: direction (x,y or z) as (1,2,3)
def get_surface_lines(map,numLines,dir):
  assert(dir == 1  or dir == 2 or dir == 3)
  lenx = np.size(map,1)
  
  nodes = np.zeros((numLines, lenx, 3))
  for i,grad in enumerate(np.arange(0,360,360.0/numLines)):
    radii = map[grad,:]
    rad = grad/180.0*np.pi
     
    #transformation to cartesian coordinates
    X,Y = polar_to_cartesian(radii, rad, 0.5*np.ones(lenx))
#     X = radii*np.cos(rad) + 0.5*np.ones(lenx)
#     Y = radii*np.sin(rad) + 0.5*np.ones(lenx)
    if dir == 1:
      nodes[i,:,0] = np.linspace(0.0,1.0,lenx)
      nodes[i,:,1] = X
      nodes[i,:,2] = Y
    if dir == 2:
      nodes[i,:,0] = Y
      nodes[i,:,1] = np.linspace(0.0,1.0,lenx)
      nodes[i,:,2] = X
    if dir == 3:
      nodes[i,:,0] = X
      nodes[i,:,1] = Y
      nodes[i,:,2] = np.linspace(0.0,1.0,lenx)
    
  return nodes

# for direction 'dir' with nodes 'nodes', find surface nodes on the left and right
# each surface point is also given a unique id
# @param pointId indicates the point id that already exists
def find_points_on_surface(nodes, dir, otherMap1, otherMap2, pointId):
  
  if dir == 1:
    id0 = 1
    id1 = 0
    id2 = 2
    id3 = 2
    id4 = 1
    id5 = 0
  elif dir == 2:
    id0 = 0
    id1 = 2
    id2 = 1
    id3 = 2
    id4 = 1
    id5 = 0
  else: #dir == 3
    id0 = 0
    id1 = 2
    id2 = 1
    id3 = 1
    id4 = 0
    id5 = 2
  
  assert(dir == 1 or dir == 2 or dir == 3)
  res = nodes.shape[1]
  
  nodes_ids = np.ones(nodes.shape[0:2], dtype=np.int) * (-1)
  
  for numLine,line in enumerate(nodes):
    lineLeft = [] # store points on one surface line temporarily
    lineRight = []
    for i,p in enumerate(line):
      # check if point is contained in other profiles
      if not contains_point(cartesian_to_grid_coords(p[id0], res),p[id1],p[id2],otherMap1) and not contains_point(cartesian_to_grid_coords(p[id3], res), p[id4], p[id5], otherMap2):
        nodes_ids[numLine,i] = pointId
        pointId += 1
        
  return nodes_ids, pointId
  
# list is list of lists contains surface lines and all respective points
# base used for setting right ids
# nodes is a 3d array with (x,y,z) coords for one profile
# returns array of size 2 x n with ids of end nodes in first row, second row reserved for fix_gaps()
def define_triangles(nodes_ids,nodes,cells,dir,vtkArray):
  end_nodes = []
  
  for i in range(0,nodes_ids.shape[0]):
    this_line = nodes_ids[i]
    right_line = nodes_ids[i+1 if i < nodes_ids.shape[0]-1 else 0] # wrap arround
    left_line = nodes_ids[i-1 if i > 0 else nodes_ids.shape[0]-1] # wrap arround
    for j in range(0,len(this_line)-1):
      this_id = this_line[j]
      prev_id = -1 if j == 0 else this_line[j-1]
      next_id = this_line[j+1]
      right_id = right_line[j]
      left_id = left_line[j]
      right_next_id = right_line[j+1]
      if this_id >= 0 and next_id >= 0 and right_id >= 0:
        tri = vtk.vtkTriangle()
        tri.GetPointIds().SetId(0, this_id)
        tri.GetPointIds().SetId(1, right_id)
        tri.GetPointIds().SetId(2, next_id)
        cells.InsertNextCell(tri)
        
        vtkArray.SetValue(this_id,dir)
        vtkArray.SetValue(next_id,dir)
        vtkArray.SetValue(right_id,dir)
        
        if j > 0 and (left_id < 0 or prev_id < 0):
          end_nodes.append(End_Node(nodes[i,j,:],this_id,dir,i,j))
          vtkArray.SetValue(this_id,dir+0.5)
      else:
        if this_id >= 0 and j > 0: # first point in line is not end point
          end_nodes.append(End_Node(nodes[i,j,:],this_id,dir,i,j))
          vtkArray.SetValue(this_id,dir+0.5)
          
          if right_id >= 0 and right_next_id >= 0:
            tri = vtk.vtkTriangle()
            tri.GetPointIds().SetId(0, this_id)
            tri.GetPointIds().SetId(1, right_id)
            tri.GetPointIds().SetId(2, right_next_id)
            cells.InsertNextCell(tri)
      
      if right_id >= 0 and next_id >= 0 and right_next_id >=0:
        tri = vtk.vtkTriangle()
        tri.GetPointIds().SetId(0, right_id)
        tri.GetPointIds().SetId(1, right_next_id)
        tri.GetPointIds().SetId(2, next_id)
        cells.InsertNextCell(tri)
        
      if this_id >= 0 and right_id < 0 and next_id >= 0 and right_next_id >= 0 :
        tri = vtk.vtkTriangle()
        tri.GetPointIds().SetId(0, this_id)
        tri.GetPointIds().SetId(1, right_next_id)
        tri.GetPointIds().SetId(2, next_id)
        cells.InsertNextCell(tri)
           
  return end_nodes

# calc distance between two points
def calc_distance(coords1,coords2):
  return np.sqrt((coords1[0]-coords2[0])**2 + (coords1[1]-coords2[1])**2 + (coords1[2]-coords2[2])**2)

# search for 'num' points closest points to ref; ref is an end point
def find_closest_points(ref_node,end_nodes_1,end_nodes_2,num):
  assert(num == 1 or num==2)
  if num == 1:
    n1, d1 = find_closest_point(ref_node,end_nodes_1)
    n2, d2 = find_closest_point(ref_node,end_nodes_2)
    
#     print "n1:",n1,"d1:",d1
#     print "n2:",n2," d2:",d2
    if d1 < d2:
      return n1
    else:
      return n2
  else:
    x11, x12 = find_two_closest_points(ref_node, end_nodes_1)
    x21, x22 = find_two_closest_points(ref_node, end_nodes_2)
    
    l = [x11,x12,x22,x21]
    res = sorted(l,key=lambda x: x[1])[0:2]
    
    return res[0][0],res[1][0]
  
def find_closest_point(ref_node, end_nodes):
  # iterate over all end nodes in list and compare distances to 'ref' node
  min_distance = 1e6
  min_node = None
  for node in end_nodes:
    assert(node.dir != ref_node.dir)
    
    dist = calc_distance(ref_node.coords, node.coords)
    if dist < min_distance:
      min_distance = dist
      min_node = node
  
  return min_node,min_distance     

def find_two_closest_points(ref_node, end_nodes):
  # iterate over all end nodes in list and compare distances to 'ref' node
  min_distance_1 = 1e6
  min_node_1 = None
  min_distance_2 = 1e6
  min_node_2 = None
  for node in end_nodes:
    assert(node.dir != ref_node.dir)
    
    dist = calc_distance(ref_node.coords, node.coords)
    if dist < min_distance_1:
      min_distance_2 = min_distance_1
      min_node_2 = min_node_1
      
      min_distance_1 = dist
      min_node_1 = node
      
  
  return (min_node_1,min_distance_1),(min_node_2,min_distance_2)     
  
# for each end point, check if closest point is on same line
# if true, then search for third point in other profiles' end points
# if false, then search for second and third point in other profiles' end points
def fix_end_node_gaps(this_end_nodes, other_1_end_node, other_2_end_node,cells):
  for n,node in enumerate(this_end_nodes):
      next_node = None
      if n < len(this_end_nodes)-1:
        next_node = this_end_nodes[n+1]
        if not(node.i == next_node.i and node.j+1 == next_node.j):
          next_node = None
      if next_node is None:    
        right_node = [v for v in this_end_nodes if v.i == node.i+1 and v.j == node.j]
        if right_node:
          next_node = right_node[0]
      if next_node is None and n < len(this_end_nodes)-1:
        right_next_node = [v for v in this_end_nodes if v.i == node.i+1 and v.j == node.j+1]
        if right_next_node:
          next_node = right_next_node[0]
      
      if next_node is not None and node not in next_node.connections: # consecutive neighbor on same surface line
        other_node = find_closest_points(node, other_1_end_node, other_2_end_node, 1)
        # FIXME: fix orientation of triangle
        tri = vtk.vtkTriangle()
        tri.GetPointIds().SetId(0, node.id)
        tri.GetPointIds().SetId(1, next_node.id)
        tri.GetPointIds().SetId(2, other_node.id)
        cells.InsertNextCell(tri)
        
#         print "node: ", node
#         print "next_node: ", next_node
#         print "other_node: ", other_node
        
        node.connections.append(next_node)
        node.connections.append(other_node)
        next_node.connections.append(node)
        next_node.connections.append(other_node)
        other_node.connections.append(node)
        other_node.connections.append(next_node)
      else:
        other_node_1, other_node_2 = find_closest_points(node, other_1_end_node, other_2_end_node, 2)
        # FIXME: fix orientation of triangle
        tri = vtk.vtkTriangle()
        tri.GetPointIds().SetId(0, node.id)
        tri.GetPointIds().SetId(1, other_node_1.id)
        tri.GetPointIds().SetId(2, other_node_2.id)
        cells.InsertNextCell(tri)
        
        node.connections.append(other_node_1)
        node.connections.append(other_node_2)
        other_node_1.connections.append(node)
        other_node_1.connections.append(other_node_2)
        other_node_2.connections.append(node)
        other_node_2.connections.append(other_node_1)
    
  
def create_profiles_array(args,infoXml):
  res = args.res
  array = np.ones((res,res,res)) * (-1)
  vec = None
  t = np.linspace(0, 1.0, args.res)
  
  hf = plt.figure()
  ha = hf.add_subplot(111, projection='3d')
  ha.set_xlabel('X')
  ha.set_ylabel('Y')
  ha.set_zlabel('Z')
  
  profiles = create_profiles(args,infoXml)
  
  assert(len(profiles) == 3)
  
  surfNodesXLeft = []
  surfNodesYLeft = []
  surfNodesZLeft = []
  surfNodesXRight = []
  surfNodesYRight = []
  surfNodesZRight = []
  
  id = 0 # assign an id to each surface point 
  
  if args.target == "surface_mesh":
    assert(not args.skip_x and not args.skip_y and not args.skip_z)
    
    map_x = create_profile_map(profiles[0], res)
    map_y = create_profile_map(profiles[1], res)
    map_z = create_profile_map(profiles[2], res)
    
    # first dimension of nodes : surface lines
    # second dimension of nodes: resolution of unit cube
    # third dimension: tuple with x,y,z coordinate
    nodes_1 = get_surface_lines(map_x, args.res_surf_lines, 1)
    nodes_ids_1, id = find_points_on_surface(nodes_1, 1, map_y, map_z,id)

    nodes_2 = get_surface_lines(map_y, args.res_surf_lines, 2)
    nodes_ids_2, id = find_points_on_surface(nodes_2, 2, map_x, map_z, id)
    
    nodes_3 = get_surface_lines(map_z, args.res_surf_lines, 3)
    nodes_ids_3, id = find_points_on_surface(nodes_3, 3, map_x, map_y, id)
    
    # create vtk cells and points
    cells = vtk.vtkCellArray()
    points = vtk.vtkPoints()
    points.SetNumberOfPoints(id)
    
    # set scalar info of intersection points to 1.0 to make them visible
    vtkData = vtk.vtkFloatArray()
    vtkData.SetName("intersection")
    vtkData.SetNumberOfValues(id)

    end_nodes_1 = define_triangles(nodes_ids_1,nodes_1,cells,1,vtkData)
    end_nodes_2 = define_triangles(nodes_ids_2,nodes_2,cells,2,vtkData)
    end_nodes_3 = define_triangles(nodes_ids_3,nodes_3,cells,3,vtkData)
    
    fix_end_node_gaps(end_nodes_1, end_nodes_2, end_nodes_3, cells)
    fix_end_node_gaps(end_nodes_2, end_nodes_1, end_nodes_3, cells)
    fix_end_node_gaps(end_nodes_3, end_nodes_1, end_nodes_2, cells)

    for i,line in enumerate(nodes_1):
      for j in range(len(line)):
        if nodes_ids_1[i,j] >= 0:          
          points.SetPoint(nodes_ids_1[i,j], nodes_1[i,j,0], nodes_1[i,j,1], nodes_1[i,j,2])
        if nodes_ids_2[i,j] >= 0:          
          points.SetPoint(nodes_ids_2[i,j], nodes_2[i,j,0], nodes_2[i,j,1], nodes_2[i,j,2])
        if nodes_ids_3[i,j] >= 0:          
          points.SetPoint(nodes_ids_3[i,j], nodes_3[i,j,0], nodes_3[i,j,1], nodes_3[i,j,2])    
          
    polydata = vtk.vtkPolyData()
    polydata.SetPoints(points)
    polydata.GetPointData().SetScalars(vtkData)
    
    if args.export == "surface_points":
      show_write_vtk(polydata, 1000, "surface_points.vtp")
      
    polydata.SetPolys(cells)
    
    if args.show:
      show_vtk(polydata, 1000, [], True)
      
    show_write_vtk(polydata,1000,"surface.vtp")
  else:
    for i in range(0,3):
      if profiles[i] == None:
        continue
      if args.verbose == 'profile_map' or args.export == 'radius_maps':
        create_profile_map(profiles[i], res, args.verbose,args.export == 'radius_maps', ha)
      if args.target == "volume_mesh" or args.target == "volume_vtk":
        write_profile_to_array(array, profiles[i], i+1)
      if args.target == "3dlines":
        plot_3dlines(profiles[i], res, args.res_surf_lines, i+1, ha)
  
  if args.target == '3dlines':
    plt.show()
  
  return array

# creates map with info on Profile depending on radius
# Profile contains list of tuples with vector,angle and idx where constant part begins (bisec: res/2, orthogonal: grad is 1)
def create_profile_map(profile,res,verbose=None,save=None,ha=None):
  map = np.zeros((360, res))
  h = 1.0 / res
  if save:
    out = open("profile_map_dir_"+str(profile.direction)+".txt","w")
    out.write("#i \t alpha \t radius\n")
    
  for i,x in enumerate(np.arange(0,1.0,h)):
    for alpha in range(0,360):
      map[alpha,i] = calc_radius_for_quadrant(profile, x, degree_to_rad_quadrant(alpha))
        
      if save:
        out.write(str(i) + " \t" + str(alpha) + " \t" + str(map[alpha,i]) + "\n")
        
  if verbose == 'profile_map':
    ha.set_xlabel('X')
    ha.set_ylabel('Y')
    ha.set_zlabel('Z')
    X,Y = np.meshgrid(range(res),range(360))
    ha.plot_surface(X, Y, map)
    plt.show()
    
  if save:
    out.close()
    
  return map

def plot_3dlines(profile,res,numLines,dir,ha):
  Z = np.linspace(0,2*np.pi,360)
  
  nodes = []
  
  map = create_profile_map(profile, res)
        
#   for ii in range(0,res,10):
#     radii = map[:,ii]
#     X = radii*np.cos(Z)+.5*np.ones(np.size(radii))
#     Y = radii*np.sin(Z)+.5*np.ones(np.size(radii))
#     if dir == 1:
#       ha.plot(ii*np.ones(np.size(X))/res,X,Y,'r')
#     if dir == 2:
#       ha.plot(Y,ii*np.ones(np.size(X))/res,X,'g')
#     if dir == 3:
#       ha.plot(X,Y,ii*np.ones(np.size(X))/res,'b')
      
  nodes = get_surface_lines(map, numLines, dir)
  color = None
  if dir == 1:
    color = "red"
  elif dir == 2:
    color = "blue"
  else:
    assert(dir == 3)
    color = "green"  
  for i in range(numLines):
    ha.plot(nodes[i,:,0],nodes[i,:,1],nodes[i,:,2],color=color)
    
# calculates interpolated radius for one quadrant (0 to pi/2)
# @param profile: contains three profile functions (for 0, phi (bisec) and 90 degree)
# @param x: parameter for function evaluation
# @param rad: radians for evaluation
def calc_radius_for_quadrant(profile,x,rad):
  assert(rad >= 0 and rad <= np.pi/2.0)
  
  funcs = profile.functions
  
  phi = funcs[1].angle
  assert(phi >= 0 and phi <= np.pi/2.0)
  if rad <= phi:
    return  (1 - 1.0/phi * rad) * funcs[0].eval(x) + 1.0/phi * rad * funcs[1].eval(x) - 0.5
  else : # rad <= np.pi/2.0
    fact = (rad-phi) / (np.pi/2.0-phi)
    return  (1-fact) * funcs[1].eval(x) + fact * funcs[2].eval(x) - 0.5

# rasterize profile functions
def write_profile_to_array(array,profile,dir):
  res = array.shape[0]
  assert(dir >=1 and dir <=3)
  
  map = create_profile_map(profile, res)
#   plt.savefig("spline_bend_" + str(bend) + ".png")
      
  for i in range(0,res):
    for j in range(0,res):
      for k in range(0,res):
        y = grid_to_cartesian_coords(j, res)
        z = grid_to_cartesian_coords(k, res)
        valx = (y-0.5)**2 + (z-0.5)**2
        
        phi = angle_to_center((z,y))
        
        p = map[int(phi/np.pi*180),i]
        if (valx <= p*p):
          if dir == 1:
            array[i,j,k] = dir
          if dir == 2:
            array[j,i,k] = dir
          if dir == 3:
            array[k,j,i] = dir
