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
import timeit



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

# calculates euklidian distance to origin (0,0)
# @param p: tuple with x-,y-component of point
def distance_to_center(p):
  return np.sqrt((p[0]-0.5)**2+(p[1]-0.5)**2)

# calculates angle between (0.5,0.5) and point p
def angle_to_center(p):
  phi = np.arccos(0.5*(p[1]-0.5)/(0.5*np.sqrt((p[0]-0.5)**2 + (p[1]-0.5)**2)))
  if p[0] < 0.5:
    phi = 2.0*np.pi - phi
    
  return phi
 
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
  phi = angle_to_center((y,z)) # angle between (0.5,0.5) and (y,z)
  
  if y < 0.5:
    phi = 2.0 * np.pi - phi
          
  r = map[int(phi/np.pi*180),id_x]
  
  if (valx <= r*r):
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
    if x < self.coords_cut[0]:
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
    assert(b[0] <= 0.5 and b[1] >= 0.5)

    # search for point p
    right = PrincipleSpline(y1, z1, bend)
    p = right.coords_cut
  
    assert(p[0] <= 0.5 and p[1] >= 0.5)
    height = distance_to_center(p) + 0.5
  
    self.angle = 2*np.pi-angle_to_center(p)
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
    
    bicubic = self.eval_bicubic(x)
    spline = self.eval_spline(x)
    linear = self.eval_linear(x)
    
    plt.plot(x,bicubic,label='bicubic',linewidth=5.0)
    plt.plot(x,spline,label='spline',linewidth=5.0)
    plt.plot(x,linear,label='linear',linewidth=5.0) 
      
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
    self.dir = dir
    self.type = args.type
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
        self.functions[1] = BisecSpline(args.y1, args.z1, args.x1, args.bend)  
        self.functions[2] = PrincipleSpline(args.y1, args.z1, args.bend, np.pi/2.0)
      else:
        self.functions[0] = PrincipleSpline(args.z1, args.y1, args.bend, 0)
        self.functions[1] = BisecSpline(args.z1, args.x1, args.y1, args.bend)  
        self.functions[2] = PrincipleSpline(args.z1, args.x1, args.bend, np.pi/2.0)
        
      self.bisec_angle = self.functions[1].angle
      
    if args.verbose == "bisec":
      self.functions[1].plot()
      
#       t = np.linspace(0, 1.0, args.res)
#       plt.plot(t,vec1,label='x1y1', linewidth=5.0)
#       plt.plot(t,vec2,label='bisec', linewidth=5.0)
#       plt.plot(t,vec3,label='x1z1', linewidth=5.0)
#       plt.legend(loc='upper left', shadow=True)
#       plt.rcParams.update({'font.size': 18})
#       plt.savefig("3splines.png")
#       plt.show()
    
      
# return information on profiles 
def create_profiles(args,infoXml=None):
  profiles = [None]*3 # x-,y-,z-part
  
  if not args.skip_x:
    profiles[0] = Profile(args,1)
    
  if not args.skip_y:
    profiles[1] = Profile(args,2)
    
#     plt.gcf().clear()
#     plt.gcf().subplots_adjust(bottom=0.15)
#     x = np.linspace(0, 1.0, args.res)
#     plt.plot(x,profiles[1].functions[0].eval(x),label="y_0")
#     plt.plot(x,profiles[1].functions[2].eval(x),label="y_90")

  if not args.skip_z:
    profiles[2] = Profile(args,3)
    x = np.linspace(0, 1.0, args.res)
#     plt.plot(x,profiles[2].functions[0].eval(x),label="z_0")
#     plt.plot(x,profiles[2].functions[2].eval(x),label="z_90")
#     plt.legend(loc='upper left', shadow=True)
#     plt.show()
    
  if args.verbose == "all_profiles":
    plt.gcf().clear()
    plt.gcf().subplots_adjust(bottom=0.15)
    x = np.linspace(0, 1.0, 1000)
    
    for dir,profile in enumerate(profiles):
      if profile == None:
        continue
#       plt.plot(x,profile.functions[0].eval(x),label=str(profile.functions[0].angle),linewidth=5.0,label="dir_"+str(i+1)+"_0")
      plt.plot(x,profile.functions[0].eval(x),linewidth=5.0,label="dir_"+str(dir+1)+"_0")
      plt.plot(x,profile.functions[1].eval(x),linewidth=5.0,label="dir_"+str(dir+1)+"_"+str(profile.functions[1].angle[0]))
#       plt.plot(x,profile.functions[1].eval(x),label=str(profile.bisec_angle))
#       plt.plot(x,profile.functions[1].eval(x),label=str(profile.bisec_angle))
      plt.plot(x,profile.functions[2].eval(x),linewidth=5.0,label="dir_"+str(dir+1)+"_90")
#       plt.plot(x,profile.functions[2].eval(x),label=str(profile.functions[2].angle))
    
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
    X = radii*np.cos(rad) + 0.5*np.ones(lenx)
    Y = radii*np.sin(rad) + 0.5*np.ones(lenx)
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
def find_points_on_surface(nodes,dir,otherMap1,otherMap2, pointId):
  nodesLeft = []
  nodesRight = []
  pointId += 1
  
  if dir == 1:
    id0 = 0
    id1 = 2
    id2 = 0
    id3 = 1
  elif dir == 2:
    id0 = 1
    id1 = 2
    id2 = 0
    id3 = 1
  else: #dir == 3
    id0 = 1
    id1 = 2
    id2 = 0
    id3 = 2
        
  assert(dir == 1 or dir == 2 or dir == 3)
  res = nodes.shape[1]
  for numLine,prof in enumerate(nodes):
    lineLeft = [] # store points on one surface line temporally
    lineRight = []
    for i,p in enumerate(prof):
      if not contains_point(i,p[id0],p[id1],otherMap1) and not contains_point(i, p[id2], p[id3], otherMap2):
        if p[dir-1] < 0.5: # using plane through origin (0.5,0.5) to decide if surface point is on left or right side of structure
          lineLeft.append((p,pointId))
        else:
          lineRight.append((p,pointId))
        pointId += 1
        
    nodesLeft.append(lineLeft)
    nodesRight.append(lineRight[::-1])   
    
  return nodesLeft, nodesRight, pointId
# list is list of lists contains surface lines and all respective points
# base used for setting right ids
def define_triangles(list,cells):
  prevLine = None
  for i,line in enumerate(list):
    if i < len(list) - 1:
      prevLine = list[i+1]
    else:
      prevLine = list[0] # last line in list
      
    for j in range(0,len(line),1): # get coordinates of each point
      if (j+1 < np.size(line, 0)) and (j < np.size(prevLine, 0)):
        tri = vtk.vtkTriangle()
        tri.GetPointIds().SetId(0, line[j][1]) 
        tri.GetPointIds().SetId(1, prevLine[j][1])
        tri.GetPointIds().SetId(2, line[j+1][1])
        cells.InsertNextCell(tri)
         
#         print "triangle:",line[j][1],prevLine[j][1],line[j+1][1]
#         print "point",line[j][1], ": ",line[j][0][0],line[j][0][1],line[j][0][2]
#         print "point",prevLine[j][1], ": ",prevLine[j][0][0],prevLine[j][0][1],prevLine[j][0][2]
#         print "point",line[j+1][1], ": ",line[j+1][0][0],line[j+1][0][1],line[j+1][0][2]
      if (j+1 < np.size(line, 0)) and (j+1 < np.size(prevLine, 0)):
        tri = vtk.vtkTriangle()
        tri.GetPointIds().SetId(0, line[j+1][1])
        tri.GetPointIds().SetId(1, prevLine[j][1])
        tri.GetPointIds().SetId(2, prevLine[j+1][1])
        cells.InsertNextCell(tri)
    
#         print "triangle:",line[j+1][1],prevLine[j][1],prevLine[j+1][1]
       
#         print "point",line[j+1][1], ": ",line[j+1][0][0],line[j+1][0][1],line[j+1][0][2]
#         print "point",prevLine[j][1], ": ",prevLine[j][0][0],prevLine[j][0][1],prevLine[j][0][2]
#         print "point",prevLine[j+1][1], ": ",prevLine[j+1][0][0],prevLine[j+1][0][1],prevLine[j+1][0][2]

#     if i == 9:
#       break  
  return cells
  
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
      map_x = create_profile_map(profiles[0], args.res) if profiles[0] <> None else None
      map_y = create_profile_map(profiles[1], args.res) if profiles[1] <> None else None
      map_z = create_profile_map(profiles[2], args.res) if profiles[2] <> None else None
      
      idx = -1 # enumerating all surface nodes in order to create triangles
      
      # first dimension of nodes : surface lines
      # second dimension of nodes: resolution of unit cube
      # third dimension: tuple with x,y,z coordinate
      nodes = get_surface_lines(map_x, args.res_surf_lines, 1)
      
      #ha.scatter(nodes[:,:,0],nodes[:,:,1],nodes[:,:,2])
      
      #ha.scatter(nodes[0:4,:,0],nodes[0:4,:,1],nodes[0:4,:,2],color='red')
      
      surfNodesXLeft, surfNodesXRight, idx = find_points_on_surface(nodes, 1, map_y, map_z, idx)
      
#       for i,line in enumerate(surfNodesXRight):
#         for tuple in line:
#           ha.scatter(tuple[0][0],tuple[0][1],tuple[0][2])
#           ha.text(tuple[0][0],tuple[0][1],tuple[0][2],tuple[1])
#         if i == 10:
#           break
#       plt.show()
      
      out = open("surfNodesXLeft.coords","w")
      for i,line in enumerate(surfNodesXLeft):
        for tuple in line:
          out.write("id=" + str(tuple[1]) + ": x=" + str(tuple[0][0]) + " y=" + str(tuple[0][1]) + " z=" + str(tuple[0][2]) + "\n")
      out.close()
      nodes = get_surface_lines(map_y, args.res_surf_lines, 2)
      surfNodesYLeft, surfNodesYRight, idx = find_points_on_surface(nodes, 2, map_x, map_z, idx)  
            
      nodes = get_surface_lines(map_z, args.res_surf_lines, 3)
      surfNodesZLeft, surfNodesZRight, idx = find_points_on_surface(nodes, 3, map_x, map_y, idx)
      
      surfNodes = [surfNodesXLeft, surfNodesXRight, surfNodesYLeft, surfNodesYRight, surfNodesZLeft, surfNodesZRight]
      
      if args.show:
        for list in surfNodes:
          for line in list:
            for tuple in line: # tuple consists of on array with 3 coordinate components and 1 point id
              ha.scatter(tuple[0][0],tuple[0][1],tuple[0][2])
            break

        plt.show()
      
      # create vtk cells and points
      cells = vtk.vtkCellArray()
      points = vtk.vtkPoints()
      points.SetNumberOfPoints(idx)
      
      cells = define_triangles(surfNodesXLeft,cells)
      cells = define_triangles(surfNodesXRight,cells)
      cells = define_triangles(surfNodesYLeft,cells)
      cells = define_triangles(surfNodesYRight,cells)
      cells = define_triangles(surfNodesZLeft,cells)
      cells = define_triangles(surfNodesZRight,cells)
      
      for part in surfNodes:
        for line in part:
          for tuple in line:
            p = tuple[0]
            # set point with respectived ids as set in find_points_on_surface()
            points.SetPoint(tuple[1], (tuple[0][0], tuple[0][1], tuple[0][2]))
          
#       out = open("vtk_point_ids.txt","w")
#       for i in range(points.GetNumberOfPoints()):
#         p = points.GetPoint(i)
#         out.write("id=" + str(i) + ": x=" + str(p[0]) + " y=" + str(p[1]) + " z=" + str(p[2]) + "\n")
#       
#       out.close()  
            
      polydata = vtk.vtkPolyData()
      polydata.SetPoints(points)
      polydata.SetPolys(cells)
      
      if args.show:
        show_vtk(polydata, 1000, [], True)
        
      show_write_vtk(polydata,1000,"surface.vtp")
  else:
    for i in range(0,3):
      if profiles[i] == None:
        continue
      if args.verbose == 'profile_map':
        create_profile_map(profiles[i], res, args.verbose, ha)
      if args.target == "volume_mesh" or args.target == "volume_vtk":
        write_profile_to_array(array, profiles[i], i+1)
      if args.target == "3dlines":
        plot_3dlines(profiles[i], res, args.res_surf_lines, i+1, ha)
  
  if args.target == '3dlines':
    plt.show()
  
  return array

# creates map with info on Profile depending on radius
# Profile contains list of tuples with vector,angle and idx where constant part begins (bisec: res/2, orthogonal: grad is 1)
def create_profile_map(profile,res,verbose=None,ha=None):
  map = np.zeros((360, res))
  h = 1.0 / res
  for i,x in enumerate(np.arange(0,1.0,h)):
    
    for alpha in range(0,360):
      rad = np.pi/180. * alpha
      if alpha <= 90:
        map[alpha,i] = calc_radius_for_quadrant(profile, x, rad)
      elif alpha <= 180:
        map[alpha,i] = calc_radius_for_quadrant(profile, x, np.pi-rad)
      elif alpha <= 270: 
        map[alpha,i] = calc_radius_for_quadrant(profile, x, rad-np.pi)
      else: # 270 <= alpha <= 360
        map[alpha,i] = calc_radius_for_quadrant(profile, x, 2*np.pi-rad)

#       if i == 1:
#         print i,alpha,map[alpha,i]
#         print x+h/2.0,rad,calc_radius_for_quadrant(profile, x+h/2.0, rad)
#         sys.exit()
#         print i," \t",alpha," \t",map[alpha,i]
#         
#     if i == 1:
#       sys.exit()
  if verbose == 'profile_map':
    ha.set_xlabel('X')
    ha.set_ylabel('Y')
    ha.set_zlabel('Z')
    X,Y = np.meshgrid(range(res),range(360))
    ha.plot_surface(X, Y, map)
    plt.show()
    
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
  for i in range(numLines):
    ha.plot(nodes[i,:,0],nodes[i,:,1],nodes[i,:,2],'r')
    
#   elif rad <= np.pi:
#     return calc_radius_for_quadrant(profile, phi, x, np.pi-rad)
#   elif rad <= 3*np.pi/2:
#     return calc_radius_for_quadrant(profile, phi, x, rad-np.pi)
#   else: #rad between 3pi/2 and 2*pi
#     return calc_radius_for_quadrant(profile, phi, x, 2*np.pi-rad)

# calculates interpolated radius for one quadrant (0 to pi/2)
# @param profile: contains three profile functions (for 0, phi (bisec) and 90 degree)
# @param x: parameter for function evaluation
# @param rad: radiant for evaluation
def calc_radius_for_quadrant(profile,x,rad):
  assert(rad >= 0 and rad <= np.pi/2.0)
  
  funcs = profile.functions
  phi = funcs[1].angle
  if x == 0.3:
    1.0/phi * rad * funcs[1].eval(x)
  if rad <= phi:
    return  (1 - 1.0/phi * rad) * funcs[0].eval(x) + 1.0/phi * rad * funcs[1].eval(x) - 0.5
  else : # rad <= np.pi/2.0
    fact = (rad-phi) / (np.pi/2.0-phi)
    return  (1-fact) * funcs[1].eval(x) + fact * funcs[2].eval(x) - 0.5
## give an interpolation for the radius within 0..2pi angle for the radius vectors at index idx
def give_interpolate_radius(vec, idx):
  if type(vec) == tuple:
    assert(len(vec) == 3)
    val = []
    rad = []
    
    val.append(vec[0][0][idx]) # 0 deg
    assert(vec[0][1] == 0.0)
    rad.append(0.0) 
    
    val.append(vec[1][0][idx]) # 45 de
    rad.append(vec[1][1])
    
    val.append(vec[2][0][idx]) # 90
    assert(vec[2][1] == np.pi/2)
    rad.append(np.pi/2)
    
    val.append(vec[1][0][idx]) # 90+45
    rad.append(np.pi-vec[1][1]) #np.pi/2+phi)
    
    val.append(vec[0][0][idx]) # 90+90=180
    rad.append(np.pi)
    
    val.append(vec[1][0][idx]) # 180+45
    rad.append(np.pi+vec[1][1])
    
    val.append(vec[2][0][idx]) # 180+90=270
    rad.append(6.0/4.0*np.pi)
    
    val.append(vec[1][0][idx]) # 270+45
    rad.append(2*np.pi-vec[1][1])
    
    val.append(vec[0][0][idx]) # 270+90=360=0 (full circle) 
    rad.append(2*np.pi)
      
    f = interpolate.interp1d(rad, val)
    return f
  else:
    assert(type(vec) <> tuple) # but a single vector
    f = interpolate.interp1d((0, 7),(vec[idx], vec[idx])) # 7 is close to 2*pi
    return f                            
                                       

# @param vec can be one vector or list of vector
def write_profile_to_array(array,profile,dir):
  res = array.shape[0]
  h = 1.0 / res
  assert(dir >=1 and dir <=3)
  
  map = create_profile_map(profile, res)
#   plt.savefig("spline_bend_" + str(bend) + ".png")
      
  for i in range(0,res):
    for j in range(0,res):
      for k in range(0,res):
        y = j * h + h / 2.0
        z = k * h + h / 2.0
        valx = (y-0.5)**2 + (z-0.5)**2
        
        phi = angle_to_center((y,z))
        
        p = map[int(phi/np.pi*180),i]
        if (valx <= p*p):
          if dir == 1:
            array[i,j,k] = dir
          if dir == 2:
            array[j,i,k] = dir
          if dir == 3:
            array[k,j,i] = dir
