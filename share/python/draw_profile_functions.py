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
  coords = None
  connections = None
  id = -1 # vtk id of point
  i = -1 # first index in surface map, corresponds to surface line
  j = -1 # second index in surface map, where on surface line 
  dir = -1 # to which structure does this end node belongs to (x,y,z)(1,2,3)?
    
  def __init__(self,coords,id,dir,i,j):
    assert(len(coords) == 3)
    self.coords = coords
    self.id = id
    self.i = i
    self.j = j
    self.connections = set() # array of end_nodes that form a triangle together with this end node
    self.dir = dir
    
  def __str__(self):
    return "id=" + str(self.id) + " coords:" + str(self.coords) + " i=" + str(self.i) + " j=" + str(self.j) + " dir=" + str(self.dir)
  
class Marching_Triangle():
  # describes active edge (by two points) where we append next triangle
  # the two points always consist of one point from one profile and and second point from another profile
  edge = None
  # stores upcoming decisions
  # 0th entry is always the next node we decide to create next triangle with
  # 1st entry is node after next node we decide to create next triangle with
  # ...
  next_node = None
  # before deciding which one is the next node, we search in neighborhood for
  # possible candidates and a eval the possibly resulting triangles with a defined metric
  other_candidate = None
  vertices = None
    
  # define active edge where edge_node_2 is assumed to be next_node in different profile than edge_node_1  
  def __init__(self,vertices,edge_node_1,edge_node_2,alt_cand):
    self.edge = []
    self.next_node = []
    self.other_candidate = alt_cand
    # array of three nodes
    self.vertices = vertices
    self.edge.append(edge_node_1)
    self.edge.append(edge_node_2)
    if len(self.next_node) > 0:
      print self.next_node[0]
    assert(len(self.next_node) == 0)
    self.next_node.append(edge_node_2)
    
  def __str__(self):
    return "vertices: " + str(self.vertices[0].id) + " "  + str(self.vertices[1].id) + " " + str(self.vertices[2].id)
        
class Cubic_spline():
  # assume we have u_0=u_1=u_2=u_3=0 and u_4=u_5=u_6=u_7=0
  # a cubic spline is defined by its base functions and control polygon
  #f03 = (1-t)**3
  #f13 = 3*t*(1-t)**2
  #f23 = 3*t**2*(1-t)
  #f33 = t**3
   
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
        add_triangle(this_id,right_id,next_id,cells)
        
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
            add_triangle(this_id,right_id,right_next_id,cells)
      
      if right_id >= 0 and next_id >= 0 and right_next_id >=0:
        add_triangle(right_id,right_next_id,next_id,cells)
        
      if this_id >= 0 and right_id < 0 and next_id >= 0 and right_next_id >= 0:
        add_triangle(this_id,right_next_id,next_id,cells)
           
  return end_nodes

def add_triangle(id1,id2,id3,cells):
  tri = vtk.vtkTriangle()
  tri.GetPointIds().SetId(0, id1)
  tri.GetPointIds().SetId(1, id2)
  tri.GetPointIds().SetId(2, id3)
  cells.InsertNextCell(tri)

# calc distance between two points
def calc_distance(p1,p2):
  return np.sqrt((p1[0]-p2[0])**2 + (p1[1]-p2[1])**2 + (p1[2]-p2[2])**2)

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
    
#     fix_end_node_gaps(end_nodes_1, end_nodes_3, end_nodes_2, cells)
#     fix_end_node_gaps(end_nodes_2, end_nodes_1, end_nodes_3, cells)
#     fix_end_node_gaps(end_nodes_3, end_nodes_1, end_nodes_2, cells)

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


def get_end_node(id,end_nodes):
  for node in end_nodes:
    if node.id == id:
      return node
  
  print "node with id ", id, " not found in list"
  return None
  
# for a given list of of end_nodes
# check if any of these nodes lie in triangle 
# by calculating barycentric coordinates for each node
# if factors for the barycentric coordinates or not between 0 and 1, we are outside the triangle
def triangle_contains_any_node(vertices,end_nodes):
  for node in end_nodes:
    
    # if we are one of the vertices skip
    if node.id == vertices[0].id or node.id == vertices[1].id or node.id == vertices[2].id:
      continue
    
    v1 = vertices[1].coords - vertices[0].coords
    v2 = vertices[2].coords - vertices[0].coords
    v3 = node.coords - vertices[0].coords
    
    dot11 = np.dot(v1,v1)
    dot12 = np.dot(v1,v2)
    dot13 = np.dot(v1,v3)
    dot22 = np.dot(v2,v2)
    dot23 = np.dot(v2,v3)
    dot33 = np.dot(v3,v3)
    
    denom = dot11 * dot22 - dot12 * dot12
    u = (dot22 * dot13 - dot12 * dot23) / denom
    v = (dot11 * dot23 - dot12 * dot13) / denom
    
#     print "triangle: ", vertices[0].id, ",",vertices[1].id, ",",vertices[2].id, " probe:", node.id, "   u: ", u, " v: ", v
    
    if (u >= 0) and (v >= 0) and (u + v < 1):
      return True 

  return False
# triangle quality is measured by the aspect ratio and whether triangle contains an end node
# we penalize the quality if triangle contains a neighbor end node             
def calc_triangle_quality(node1,end_nodes_1,node2,end_nodes_2,node3,end_nodes_3):
  neighbors_1 = give_neighbor_end_nodes(node1, end_nodes_1)
  neighbors_2 = give_neighbor_end_nodes(node2, end_nodes_2)
  neighbors_3 = give_neighbor_end_nodes(node3, end_nodes_3)

  ratio = 1
  
#   contained = None
#   contained = triangle_contains_any_node([node1,node2,node3],neighbors_1+neighbors_2+neighbors_3)  
#     print "[",node1.id,",",node2.id,",",node3.id, "] containes neighbors? ",contained
        
    
  if triangle_contains_any_node([node1,node2,node3],neighbors_1+neighbors_2+neighbors_3):
    ratio = 1000
    
  ratio *= calc_triangle_ratio(node1.coords, node2.coords, node3.coords)
  
  return ratio

# @param: vertices defining triangle
# for evaluating quality of triangle, we take the ration of the exradius to twice the inradius
def calc_triangle_ratio(v1,v2,v3):
  d1 = calc_distance(v1, v2)
  d2 = calc_distance(v1, v3)
  d3 = calc_distance(v2, v3)
  
  #aspect ratio: circumradius/2*inradius
  aspect_ratio = d1*d2*d3 / ( (d2+d3-d1) * (d1-d2+d3) * (d1+d2-d3))
  
  assert(aspect_ratio >= 1)
  
  return aspect_ratio
    
#   if aspect_ratio < threshold:
#     return True
#  
#   return False

def find_closest_point(ref_node,next_node,end_nodes):
  # iterate over all end nodes in list and compare distances to 'ref' node
  min_distance = 1e6
  min_node = None
  ref_coords = ref_node.coords + 0.5 * (next_node.coords - ref_node.coords)
   
  for node in end_nodes:
    assert(node.dir != ref_node.dir)
    dist = calc_distance(ref_coords, node.coords)
    if dist < min_distance:
      min_distance = dist
      min_node = node
  
  return min_node,min_distance 

def find_three_closest_points(ref_node,next_node,end_nodes):
  
  # iterate over all end nodes in list and compare distances to 'ref' node
  min_distance_1 = 1e6
  min_node_1 = None
  min_distance_2 = 1e6
  min_node_2 = None
  min_distance_3 = 1e6
  min_node_3 = None
  for node in end_nodes:
    
    if node.id == ref_node.id or node.id == next_node.id or node.coords[2] <= 0.5+1e-3:
      continue
    
    dist = calc_distance(ref_coords, node.coords)
    if dist < min_distance_1:
      min_distance_2 = min_distance_1
      min_node_2 = min_node_1
     
      min_distance_1 = dist
      min_node_1 = node
    elif dist < min_distance_2:
      min_distance_3 = min_distance_2
      min_node_3 = min_node_2
     
      min_distance_2 = dist
      min_node_2 = node
    elif dist < min_distance_3:
      min_distance_3 = dist
      min_node_3 = node
     
  return (min_node_1,min_distance_1),(min_node_2,min_distance_2),(min_node_3,min_distance_3)

# calculates midpoint on line between p1 and p2
def calc_midpoint(p1,p2):
  return 0.5*np.asarray([p1[0]+p2[0],p1[1]+p2[1],p1[2]+p2[2]])

def give_neighbor_end_nodes(node,end_nodes):
  # find first triangle
  candidates = []
  left_line = node.i-1 if node.i > 0 else len(end_nodes)-1
  right_line = node.i+1 if node.i < len(end_nodes)-1 else 0
  
  candidates.append((node.i,node.j-1)) # previous node on line
  candidates.append((node.i,node.j+1))
  candidates.append((left_line,node.j))
  candidates.append((left_line,node.j-1))
  candidates.append((left_line,node.j+1))
  candidates.append((right_line,node.j))
  candidates.append((right_line,node.j-1))
  candidates.append((right_line,node.j+1))
  
  res = []
  for test in candidates:
    nexts = [v for v in end_nodes if v.i == test[0] and v.j == test[1]]
    if nexts:
      res.append(nexts[0])
      
  return res

# returns next end node in end nodes list
# we assume that from our current relative position, the next end node is
# preferable in the front right direction
# if no node is found than check for left previous direction
# to check relative direction, we the previous end node in same profile
# first check      
def give_next_end_node(node,end_nodes):
  # we define front in relative coordinates by means of grid coordinates i and j
  diff_i = 0
  diff_j = 0

  # find node in connections that is on the same profile
  # get node objects
  previous = None
  for n in node.connections:
    v = get_end_node(n, end_nodes)
    if get_end_node(n, end_nodes) is not None:  
      previous = v
  
  if previous:    
    diff_i = np.sign(node.i - previous.i)
    diff_j = np.sign(node.j - previous.j)
    
    print "node ", node.id, " previous: ", previous.id, " dir_i: ", diff_i, " dir_j: ", diff_j
  
  # find first triangle
  candidates = []
  candidates.append((node.i,node.j+1)) # next node on line
  candidates.append((node.i,node.j-1))
  left_line = node.i-1 if node.i > 0 else len(end_nodes)-1
  right_line = node.i+1 if node.i < len(end_nodes)-1 else 0 
  candidates.append((right_line,node.j)) # same node on right line   
  candidates.append((right_line,node.j-1)) # previous node on right line
  candidates.append((right_line,node.j+1)) # next node on right line
  candidates.append((left_line,node.j))
  candidates.append((left_line,node.j+1))
  candidates.append((left_line,node.j-1))
  
  for test in candidates:
    nexts = [v for v in end_nodes if v.i == test[0] and v.j == test[1]]
    if nexts:
      for v in nexts:
        if v.id not in node.connections:
          return v
  
  assert(False)    
# search for 'num' points closest points to ref; ref is an end point
def find_closest_points(ref_node,next_node,end_nodes_1,end_nodes_2,num):
  assert(num == 1 or num==3)
  if num == 1:
    n1, d1 = find_closest_point(ref_node,next_node,end_nodes_1)
    n2, d2 = find_closest_point(ref_node,next_node,end_nodes_2)
   
    if d1 < d2:
      return n1
    else:
      return n2
  elif num == 3:
    x11, x12, x13 = find_three_closest_points(ref_node,next_node,end_nodes_1)
    x21, x22, x23 = find_three_closest_points(ref_node,next_node,end_nodes_2)
    
    l = [x11,x12,x13,x21,x22,x23]
    res = sorted(l,key=lambda x: x[1])[0:3]
  
    return res[0][0],res[1][0],res[2][0]

# for a given edge (list with 2 points), find a neighbor end node
# such that the resulting triangle has the best quality
def find_best_next_end_node(edge,end_nodes_1,end_nodes_2):
  others = [None] * 3
  others[0], others[1], others[2] = find_closest_points(edge[0],edge[1],end_nodes_1,end_nodes_2,3)
  # calc 3 aspect ratios for triangles with 3 possible nodes in 'others'
  ratios = np.zeros(3)
  
  for i,other in enumerate(others):
    assert(edge[0].id != edge[1].id)
    assert(edge[0].id != other.id)
    ratios[i] = calc_triangle_ratio(edge[0].coords,edge[1].coords,other.coords)
    # check for most regular triangle
    
  best_neighbor = others[np.argmin(ratios)]
  
  return best_neighbor

def update_connections(triangles,vertices):
  res = [None] * 3
  res[0] = set([v for v in vertices[0].connections if not (v==vertices[1].id or v==vertices[2].id)])
  res[1] = set([v for v in vertices[1].connections if not (v==vertices[0].id or v==vertices[2].id)])
  res[2] = set([v for v in vertices[2].connections if not (v==vertices[0].id or v==vertices[1].id)])
  
  triangles[-1].connections = res
  triangles[-1].edge[0].connections = set([v for v in triangles[-1].edge[0].connections if not (v==vertices[0].id or v==vertices[1].id or v==vertices[2].id)])
  triangles[-1].edge[1].connections = set([v for v in triangles[-1].edge[1].connections if not (v==vertices[0].id or v==vertices[1].id or v==vertices[2].id)])

def check_next_triangle(triangles,end_nodes_1,end_nodes_2,next_cand=None):
  
  # take active edge from last triangle
  active_edge = triangles[-1].edge
  a = active_edge[0]
  b = active_edge[1]
  
  if a.id == 1227 and b.id == 9167:
    return False
  
  print "\na: ", a.id, " connections: ", a.connections
  print "b: ", b.id, " connections: ", b.connections 
  
  a_nodes = end_nodes_1 if end_nodes_1[0].dir == a.dir else end_nodes_2
  b_nodes = end_nodes_1 if end_nodes_1[0].dir == b.dir else end_nodes_2
  
  cand_1_nodes = end_nodes_1 if end_nodes_1[0].dir == a.dir else end_nodes_2
  next_cand_1 = give_next_end_node(a, cand_1_nodes)
  ratio_1 = calc_triangle_quality(a,a_nodes,b,b_nodes,next_cand_1,cand_1_nodes) if next_cand_1 is not None else 1e6
  next_cand_2 = None
  cand_2_nodes = None
  ratio_2 = 1e6
  
  if not next_cand:
    cand_2_nodes = end_nodes_1 if end_nodes_1[0].dir == b.dir else end_nodes_2
    next_cand_2 = give_next_end_node(b, cand_2_nodes)
    ratio_2 = calc_triangle_quality(a,a_nodes,b,b_nodes,next_cand_2,cand_2_nodes) if next_cand_2 is not None else 1e6
    
  next = next_cand_1 if ratio_1 < ratio_2 else next_cand_2
  alt = next_cand_2 if ratio_1 < ratio_2 else next_cand_1
  
  if next_cand and next_cand.id <> next.id:
    print "next_cand.id: ", next_cand.id, " next.id: ", next.id
    asssert(next_cand.id == next.id)
    
  print next_cand_1, ratio_1, (" <- " if next_cand_1.id == next.id else " ")
  if next_cand_2:
    print next_cand_2, ratio_2, (" <- " if next_cand_2.id == next.id else " ")
    
  # if at least one candidate has good quality
  if min(ratio_1,ratio_2) < 25:
    #set connections
    a.connections.add(next.id)
    b.connections.add(next.id)
    next.connections.add(a.id)
    next.connections.add(b.id)
    triangles.append(Marching_Triangle([a,b,next],a if a.dir != next.dir else b,next,alt))
    
    print "created triangle with edge: (", triangles[-1].edge[0].id, ",", triangles[-1].edge[1].id, ") next: ", next.id
    if alt is not None:
      print "alternative: ", alt.id 
  else: # if both candidates are bad, we have go back to previous triangle an try alternative candidate
    vertices = triangles[-1].vertices
    print "triangle: ", vertices[0].id, ",",vertices[1].id, ",",vertices[2].id, " is bad -> remove"

    triangles.pop()
    update_connections(triangles, vertices)
    
    # remove all created triangles 
    while triangles[-1].other_candidate is None:
      vertices = triangles[-1].vertices
      print "remove triangle: ", vertices[0].id, ",",vertices[1].id, ",",vertices[2].id
      
      triangles.pop()
      
      # update connections 
      update_connections(triangles, vertices)
    
    vertices = triangles[-1].vertices  
    print "previous triangle: ", vertices[0].id, ",",vertices[1].id, ",",vertices[2].id, " cand: ", triangles[-1].other_candidate.id   
    next_cand = triangles[-1].other_candidate
    triangles[-1].other_candidate = None
    check_next_triangle(triangles, end_nodes_1, end_nodes_2, next_cand)
        
  return True  
# given 3 lists with end nodes of all 3 profiles
# using a marching algorithm, we connect theses nodes step by step to triangles
# start with one given triangle defined by end nodes of two profiles (assume each profile has a different color'
# we call the edge that connects two differently colored nodes an 'active' edge
# when a triangle was defined, continue with active edge and search for possible candidates (for new triangle)
# in the neighborhood of the two edge vertices              
def fix_end_node_gaps(this_end_nodes, other_1_end_node, other_2_end_node,cells):
  triangles = [] # saves all triangles found by marching
  
  # start from x profile and 0 degree and take first end node you can find
  node = [v for v in this_end_nodes if v.i == 0]
  assert(node)
  node = node[0]
  
  next = give_next_end_node(node,this_end_nodes)
  
  other = find_closest_points(node,next,other_1_end_node,other_2_end_node,1)
  other_nodes = other_1_end_node if other_1_end_node[0].dir == other.dir else other_2_end_node
  
  next.connections.add(other.id)
  next.connections.add(node.id)
  node.connections.add(other.id)
  node.connections.add(next.id)
  other.connections.add(node.id)
  other.connections.add(next.id)
  
  triangles.append(Marching_Triangle([node,next,other],next,other,node))
  
  initial = triangles[0]
  run = True
  while (len(triangles) == 1 or node not in triangles[-1].vertices) and len(triangles) > 0 and run:
    
    run = check_next_triangle(triangles, this_end_nodes, other_nodes)
#     vertices = [active_edge[0],active_edge[1],next_node]    
#     # find active edge
#     if next_node.dir != active_edge[0].id:
#       triangles.append(Marching_Triangle(vertices,active_edge[0],next_node))
#     elif next_node.dir != active_edge[1].id:
#       triangles.append(Marching_Triangle(vertices,active_edge[1],next_node))
#     else:
#       print "Ohoh...."
#       
#     print triangles[-1].vertices[0].id,triangles[-1].vertices[1].id,triangles[-1].vertices[2].id
#     print triangles[-1].vertices[0].coords,triangles[-1].vertices[1].coords,triangles[-1].vertices[2].coords
#       
#     if len(triangles) > 3:
#       if triangles[-1].is_identical(triangles[-3]):
#         break  
        
  for triangle in triangles:
    verts = triangle.vertices
    add_triangle(verts[0].id, verts[1].id, verts[2].id, cells)
      
