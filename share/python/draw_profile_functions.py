#!/usr/bin/env python
import numpy as np
from numpy import outer
import sympy.solvers
from sympy import Symbol, symbols
import sys
import math
try:
  import vtk
  from vtk.util.numpy_support import vtk_to_numpy
  from matviz_vtk import *
except:
  print("WARNING: failed to load vtk!")  

import matplotlib
matplotlib.use('tkagg')
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy import interpolate, spatial

try:
  from meshpy.tet import MeshInfo, build
except:
  print("Failed to load meshpy - need it for tetrahedralized mesh")

res = 1
res_surf_lines = 1
# file object to log activities while fixing surface gaps
logger = None
# file object for info xml
infoXml = None
interpolation = None

class End_Node():
  coords = None
  id = -1 # vtk id of point
  i = -1 # first index in surface map, corresponds to surface line
  j = -1 # second index in surface map, where on surface line 
  dir = -1 # to which structure does this end node belongs to (x,y,z)(1,2,3)?
  next = None # store next end node in clock-wise direction when we have to step over a valley
  neighbors = None # list with coordinates of all (inner + end) nodes
    
  def __init__(self,coords,id,dir,i,j,neighbors):
    assert(len(coords) == 3)
    self.coords = coords
    self.id = id
    self.i = i
    self.j = j
    self.dir = dir
    self.neighbors = neighbors
    
  def __str__(self):
#     return str(self.id) + "   " + str(self.coords[0]) + " " + str(self.coords[1]) + " " + str(self.coords[2]) + " " + str(self.dir)
    return "id=" + str(self.id) + " coords:" + str(self.coords) + " i=" + str(self.i) + " j=" + str(self.j) + " dir=" + str(self.dir)
  
  def __eq__(self,other):
    return self.id == other.id
  
class Marching_Triangle():
  # describes active edge (by two points) where we append next triangle
  # the two points always consist of one point from one profile and and second point from another profile
  edge = None
  # stores upcoming decisions
  # 0th entry is always the next node we decide to create next triangle with
  # 1st entry is node after next node we decide to create next triangle with
  # ...
  next = None
  # before deciding which one is the next node, we search in neighborhood for
  # possible candidates and a eval the possibly resulting triangles with a defined metric
  other_candidates = []
  vertices = None
    
  # define active edge where edge_node_2 is assumed to be next_node in different profile than edge_node_1  
  def __init__(self,vertices,edge_node_1,edge_node_2,alt_cand):
    self.edge = []
    self.other_candidates = alt_cand
    # array of three nodes
    self.vertices = vertices
    self.edge.append(edge_node_1)
    self.edge.append(edge_node_2)
    self.next = edge_node_2
    if edge_node_1.dir == edge_node_2.dir:
      print("edge_node_1:",edge_node_1.dir)
      print("edge_node_2:",edge_node_2.dir)
    assert(edge_node_1.dir != edge_node_2.dir)
    
  def __str__(self):
    return "vertices: " + str(self.vertices[0].id) + " "  + str(self.vertices[1].id) + " " + str(self.vertices[2].id)
  
  def is_triangle(self,vert1,vert2,vert3):  
    if vert1 not in self.vertices or vert2 not in self.vertices or vert3 not in self.vertices:
      return False
    else:
      return True
    
class Cubic_spline():
  # assume we have u_0=u_1=u_2=u_3=0 and u_4=u_5=u_6=u_7=0
  # a cubic spline is defined by its base functions and control polygon
  #f03 = (1-t)**3
  #f13 = 3*t*(1-t)**2
  #f23 = 3*t**2*(1-t)
  #f33 = t**3
   
  t_1 = None   
  #control polygon contains 4 ponts (array of lists, 1 list per point) 
  CP = None
   
  # CP is numpy array with 4 coordinates; for each coordinate use a list with x- and y-component
  def __init__(self, CP = None):
    self.CP = np.transpose(CP) if CP is not None else None
    t = np.linspace(0,1,1000) # over-sampling
    C = self.eval_t(t)
    self.explicit = interpolate.interp1d(C[0,:],C[1,:])
    
    if infoXml is not None:
      infoXml.write('      <controlPolygon>\n')
      cp = self.CP
      infoXml.write('        <P1 x="' + str(cp[0][0]) + '" y="' + str(cp[1][0]) + '"/>\n')
      infoXml.write('        <P2 x="' + str(cp[0][1]) + '" y="' + str(cp[1][1]) + '"/>\n')
      infoXml.write('        <P3 x="' + str(cp[0][2]) + '" y="' + str(cp[1][2]) + '"/>\n')
      infoXml.write('        <P4 x="' + str(cp[0][3]) + '" y="' + str(cp[1][3]) + '"/>\n')
      infoXml.write('      </controlPolygon>\n')
#       coords_cut = self.calc_coords_grad_1()
#       infoXml.write('      <gradient dx/dy="1" t="' + str(self.calc_param_grad_1()) + '" x="' + str(coords_cut[0]) + '" y="' + str(coords_cut[1]) + '"/>\n')
      
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
    ret = self.explicit(x)[()]
    
    return ret
  
  def calc_d_spline_d_t(self,t):
    return outer(3*(1-t)**2,self.CP[:,1]-self.CP[:,0]) + outer(6.0*t*(1-t), self.CP[:,2] - self.CP[:,1]) + outer(3*t**2 ,self.CP[:,3]-self.CP[:,2])
  
  def calc_param_grad_1(self):
    if self.t_1 is not None:
      return self.t_1
    
    u = Symbol('u')
    
    dC = self.calc_d_spline_d_t(u) # dC/dt
    sol = sympy.solvers.solve(dC[0][1]-dC[0][0],u)
    t = -100
    if sol[0] > 0 and sol[0] <= 1:
      t = sol[0]
    elif  sol[1] > 0 and sol[1] <= 1:
      t = sol[1]
    else:
      print("No t found where dx/dy = 1 ",sol)
    
    # conversion from t as sympy.Float to regular Python float necessary  
    self.t_1 = float(t)
    return self.t_1
  
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
    ret = self.eval_t(t1) # array of array with 1 elem, due to outer product
    return (ret[0][0],ret[1][0])
  
  def calc_min(self):
    dC = self.calc_d_spline_d_t(u)
    
def dirToString(dir):
  assert(dir >= 0 and dir <= 2)
  res = 'x'
  if dir == 1:
    res = 'y'
  if dir == 2:
    res = 'z'
   
  return res

# checks if logging should be used and adds a line break
def log(text,linebreak=True):
  if logger:
    lb = "\n" if linebreak else ""
    logger.write(text + lb)
    logger.flush()

# returns the two directions of the plane whose normal shows in profile direction
# e.g. we have x profile, then we live in the z-y plane -> return 2,1
def give_normal_plane_axes(profile_dir):
  if profile_dir == 0: # x profile
    return 2,1 # z,y plane
  elif profile_dir == 1: # y profile
    return 0,2 # x,z plane
  else: # z profile
    return 1,0 # y,x plane

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
  return i * h + h/2.0

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

class Heaviside():
  beta = 0
  eta = 0
  x1 = 0
  height = 1.0
  
  def __init__(self,beta,eta,x1,height):
    self.beta = beta
    self.eta = eta
    self.x1 = x1
    self.height = height
  
  def eval(self,x):
    # y0
    h0 = 0.5 + self.x1 / 2.0
    # ymax
    h1 = self.height
    g0 = calc_tanh(self.beta, self.eta, 0)
    g1 = calc_tanh(self.beta, self.eta, 1.0)
    a = (h1-h0) / (g1 - g0)
    b = h0 - a * g0
    return a * calc_tanh(self.beta, self.eta, 2*x) + b

# shifted tanh, returns values are between 0 and 1 for x \in [0,0.5]
def calc_tanh(beta,eta,x):
  return 1.0 - 1.0 / (np.exp(2.0*beta*(x-eta)) + 1)
   
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

# for a given profile and point p, check if p lies on or inside profile
def contains_point(p,profile):
  assert(p[0] >= 0.0 and p[0] <= 1.0)
  assert(p[1] >= 0.0 and p[1] <= 1.0)
  assert(p[2] >= 0.0 and p[2] <= 1.0)
  
  # depending on direction, estimate plane in which we perform check
  if profile.direction == 0:# check for x profile -> z-y plane
    x = p[2]
    y = p[0]
    z = p[1]
  elif profile.direction == 1: #check for y profile -> x-z plane
    x = p[0]
    y = p[1]
    z = p[2]
  else: # dir == 2; check for z profile --> y-x plane
    x = p[1]
    y = p[2]
    z = p[0]
    
  phi = angle_to_center((x,z))
  r = calc_radius_for_quadrant(profile, y, phi)
  val = (x-0.5)**2 + (z-0.5)**2
  
  if (val - r*r <= 1e-3):
    return True
  
  return False
    

# obejct defines spline in a principle plane, e.g. spline for 0 or 90 degree
class PrincipleSpline():
  
  spline = None
  angle = -1
  coords_cut = None 
  
  def __init__(self, x1, y1, bend, angle=0):
    rx = 0.5 - x1/2.0 # radius for center (0,1)
    ry = 0.5 - y1/2.0 # radius for center (0,1)
    
    if infoXml is not None:
      infoXml.write('    <bspline degree="' + str(degrees(angle)) + '" rad1="' + str(x1) + '" rad2="' + str(y1) + '" bend="' + str(bend) + '">\n')
    
    self.angle = angle
    P = np.array([[0,1-rx],[ry*bend,1-rx],[ry,1-rx*bend],[ry,1]])
    self.spline = Cubic_spline(P)
    # coordinate where slope is 1
    self.coords_cut = self.calc_coords_grad_1()
    
    if infoXml is not None:
      infoXml.write('    </bspline>\n')
    
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
    if isinstance(x, (float,int)):
      return self.eval_elem(x)
    
    # in case x is a list    
    ret = []
    
    for i in x:
      ret.append(self.eval_elem(i))
    
    if type(ret) == np.float64:
      return float(ret)
    elif len(ret) == 1:
      return float(ret[0])
    else:
      return ret
  
  def calc_coords_grad_1(self):
    return self.spline.calc_coords_grad_1()
    
# object defining spline for angles between 0,...,90 degree
class BisecSpline:
  def __init__(self):
    self.bicubic = [] # 0-th entry is spline function and 1st entry is cubic polynomial
    self.cubic = None
    self.spline = None
    self.linear = None
    self.heaviside = None
    self.x1 = 0
    self.y1 = 0
    self.z1 = 0
    self.type = None
    self.angle = None
  
  def __init__(self,x1,y1,z1,bend,beta,eta,force=None):
    self.x1 = x1
    self.y1 = y1
    self.z1 = z1
    self.bicubic = [] # 0-th entry is spline function and 1st entry is cubic polynomial
  
    assert(bend <= 1 and bend >= 0)
    self.type = None
    
    if infoXml:
      infoXml.write('    <bisectionFunction>\n')
      infoXml.write('      <bicubic>\n')
#       infoXml.write('      <bicubic coeff0="' + str(sol[0]) + '" coeff1="' + str(sol[1]) + '" coeff2="' + str(sol[2]) + '" coeff3="' + str(sol[3]) +'"/>\n')
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
    height = p[1]#distance_to_center(p) + 0.5
  
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
    
    if infoXml:
      infoXml.write('        <polynomial order="cubic" a0="' + str(sol[0]) + '" a1="' + str(sol[1]) + '" a2="' + str(sol[3]) +'"/>\n')
      infoXml.write('      </bicubic>\n')
      infoXml.write('      <bspline rad1 = "' + str(x1) + '" rad2="' + str(height) + '" bend="' + str(bend) + '">\n')
    
    #### case 2: b-spline --> bsp #############
    # if b with grad gb (approx 1) is too high for p such that the curve has a maximum within b and p we need to fallback to a b-spline from a to p
    # curve as undershoot
    
    P = np.array([[0,0.5+x1/2.0],[0.5*bend,0.5+x1/2.0],[0.5-0.5*bend,height],[0.5,height]])
    bspline = Cubic_spline(P)
    
    if infoXml:
      infoXml.write('      </bspline>\n')
      infoXml.write('      <linear x1="' + str(x1) + '"/>\n')
    
    #### case 3: linear --> lin ###########
    # in case undershooting for x1=0.9, y1=0.1, z1=0.1
    lin = Linear_1D(x1, x1)
    
    self.bicubic.append(left)
    self.bicubic.append(cubic)
    self.spline = bspline
    self.linear = lin
    
    #### case 4: heaviside --> heavi###########
    self.heaviside = Heaviside(beta, eta, x1, height)
    
    if force:
      self.type = force
    else:
      # to check if bicubic has under/overshooting when point p is much lower than point b
      # sample points and check for maximum value
      t = np.linspace(0, 1, 100)
      samples = self.eval_bicubic(t)
      if max(samples) < p[1]:
        self.type = "bicubic"
      # in case function composed of b-spline and cubic function has undershoot  
      # in case b-spline has no undershoot (point p is not below bspline(x=0))
      elif height > 0.5 + x1/2.0:  
        self.type = "bspline"
      # in case we have undershooting for biqua AND for spline
      else:
        self.type = "linear"
      
    if infoXml:
      infoXml.write('      <selection type="' + str(self.type) + '" angle="' + str(degrees(self.angle)) + '"/>\n')
      infoXml.write('    </bisectionFunction>\n\n')  
      
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
  
  def eval_heaviside(self,x):
    # make x iterable in case it's one float element and not a list
    x = np.reshape(x, np.size(x), )
    res = []
    for i in x:
      assert(i >= 0 and i <=1)
      if i <= 0.5:
        res.append(self.heaviside.eval(i))
      else:
        res.append(self.heaviside.eval(1.0-i))
    
    return res
  
  def eval(self,x):
    ret = None
    if self.type == "bicubic":
      ret =  self.eval_bicubic(x)
    elif self.type == "bspline":
      ret = self.eval_spline(x)
    elif self.type == "heaviside":
      ret = self.eval_heaviside(x)
    else: #linear case
      ret = self.eval_linear(x)
    
    if type(ret) == np.float64 or type(ret) == float:
      return float(ret)
    elif len(ret) == 1:
      return float(ret[0])
    else:
      return ret
    
  def get_type(self):
    return self.type
  
  def plot(self):
    plt.gcf().clear()
    
    x = np.linspace(0, 1, 100)
    
    if self.type == "bicubic":
      bicubic = self.eval_bicubic(x)
      plt.plot(x,bicubic,label='bicubic',linewidth=5.0)
    if self.type == "bspline":
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
    heavi = self.eval_heaviside(x)
    
    cut = self.bicubic[0].coords_cut
    
    assert(self.type == 'bicubic' or self.type == 'bspline' or self.type == 'linear' or self.type == 'heaviside')

    plt.ylim((0.5,1.0))
    bc_label = 'bicubic*' if self.type == 'bicubic' else 'bicubic'
    sp_label = 'bspline*' if self.type == 'bspline' else 'bspline'
    lin_label = 'linear*' if self.type == 'linear' else 'linear' 
    hv_label = 'heaviside*' if self.type == 'heaviside' else 'heaviside'
    plt.plot(x,bicubic,label=bc_label,linewidth=5.0)
    plt.plot(x,spline,label=sp_label,linewidth=5.0)
    plt.plot(x,linear,label=lin_label,linewidth=5.0)
    plt.plot(x,heavi,label=hv_label,linewidth=5.0)
    plt.plot(cut[0],cut[1],marker='o',color='red',markersize=15) 
    plt.legend(loc='upper left', shadow=True,prop={'size':20})
    plt.show()  
    
  def angle(self):
    return self.angle
    
# @return vector with Profile or list of vectors
class Profile:
  def __init__(self, args, dir):
    assert (dir == 0 or dir == 1 or dir == 2)
    self.bisec_angle = -1
    self.direction = dir
    # 0th entry: function for 0 degree; 1st entry: function for bisec; 2nd entry: function for 90 degree
    self.functions = [None] * 3
    # depending on profile, store the radii of the two boundary circles
    self.radius_left = 0
    self.radius_right = 0

    if infoXml:  
      infoXml.write('  <profile dir="' + str(dir) + '">\n')
      
    if dir == 0:
      self.functions[0] = PrincipleSpline(args.x1, args.y1, args.bend, 0)
      self.functions[1] = BisecSpline(args.x1, args.y1, args.z1, args.bend,args.beta,args.eta,args.force_bisec)
      self.functions[2] = PrincipleSpline(args.x1, args.z1, args.bend, np.pi/2.0)
      self.radius_left = args.x1 / 2.0
      self.radius_right = args.x1 / 2.0
    elif dir == 1:
      self.functions[0] = PrincipleSpline(args.y1, args.x1, args.bend, 0)
      self.functions[1] = BisecSpline(args.y1, args.x1, args.z1, args.bend,args.beta,args.eta,args.force_bisec)  
      self.functions[2] = PrincipleSpline(args.y1, args.z1, args.bend, np.pi/2.0)
      self.radius_left = args.y1 / 2.0
      self.radius_right = args.y1 / 2.0
    else: # dir == 2
      self.functions[0] = PrincipleSpline(args.z1, args.y1, args.bend, 0)
      self.functions[1] = BisecSpline(args.z1, args.y1, args.x1, args.bend,args.beta,args.eta,args.force_bisec)  
      self.functions[2] = PrincipleSpline(args.z1, args.x1, args.bend, np.pi/2.0)
      self.radius_left = args.z1 / 2.0
      self.radius_right = args.z1 / 2.0
      
    self.bisec_angle = self.functions[1].angle
    
    if infoXml:  
      infoXml.write('  </profile>\n\n')
  
    if args.verbose == "bisec":
      self.functions[1].plot_all()
      
# return information on profiles 
def create_profiles(args,infoXml=None):
  profiles = [None]*3 # x-,y-,z-part
  
  if not args.skip_x:
    profiles[0] = Profile(args,0)
    
  if not args.skip_y:
    profiles[1] = Profile(args,1)
    
  if not args.skip_z:
    profiles[2] = Profile(args,2)
    x = np.linspace(0, 1.0, args.res)
    
  if args.verbose == "all_splines" or args.verbose == "all_profiles":
    plt.gcf().clear()
    plt.gcf().subplots_adjust(bottom=0.15)
    x = np.linspace(0, 1.0, 1000)
    
    for dir,profile in enumerate(profiles):
      if profile == None:
        continue
      plt.plot(x,profile.functions[0].eval(x),linewidth=5.0,label="dir_"+str(dir+1)+"_0")
      if args.verbose == "all_profiles":
        plt.plot(x,profile.functions[1].eval(x),linewidth=5.0,label="dir_"+str(dir+1)+"_"+str(profile.functions[1].angle))
      plt.plot(x,profile.functions[2].eval(x),linewidth=5.0,label="dir_"+str(dir+1)+"_90")
    
    plt.ylim((0.5,1.0))
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
  assert(dir == 0  or dir == 1 or dir == 2)
  lenx = np.size(map,1)
  
  nodes = np.zeros((numLines, lenx, 3))
  for i,grad in enumerate(np.arange(0,360,360.0/numLines)):
    radii = map[int(grad),:]
    rad = grad/180.0*np.pi
     
    #transformation to cartesian coordinates
    X,Y = polar_to_cartesian(radii, rad, 0.5*np.ones(lenx))
    idx_y, idx_x = give_normal_plane_axes(dir)
    idx_z = [v for v in (0,1,2) if v != idx_x and v != idx_y]
    assert(idx_z != idx_x and idx_z != idx_y)
    nodes[i,:,idx_x] = X
    nodes[i,:,idx_y] = Y
    nodes[i,:,idx_z] = np.linspace(0.0,1.0,lenx)
      
  return nodes

# for direction 'dir' with nodes 'nodes', find surface nodes on the left and right
# each surface point is also given a unique id
# min_distance defines smallest distance between a surface node and a neighbor inner node that we allow
# if the distance between these two is smaller, we omit the surface node
def find_points_on_surface(nodes, profile, otherProfile1, otherProfile2, pointId, min_distance=1.0/res):
  dir = profile.direction
  assert(dir == 0 or dir == 1 or dir == 2)
  res = nodes.shape[1]
  
  nodes_ids = np.ones(nodes.shape[0:2], dtype=np.int) * (-1)
  
  # for each line, check if a point is contained in two other profiles
  # if not -> point is a surface point
  # if yes, omit point 
  # if point p is surface point, check if next point p+1 is also a surface point
  # if not, check distance from this inner point (p+1) to previous surface point (p)
  # if distance is too small (e.g. < 0.5* grid spacing (h)), then omit p
  
  for numLine,line in enumerate(nodes):
    for i,p in enumerate(line):
      # detect crossing border from inner to outer:
      if i > 0 and nodes_ids[numLine,i-1] == -1:
        # it is possible that we already performed this check on
        if contains_point(line[i-1], otherProfile1) or contains_point(line[i-1], otherProfile2):
          intersect = find_intersection_point(line[i-1], p, profile, otherProfile1, otherProfile2)
          if calc_distance(p, intersect) < min_distance:
            nodes_ids[numLine,i] = -1
            continue

      # check if point is contained in other profiles
      # assume we start with a point that is a surface point
      if not contains_point(p, otherProfile1) and not contains_point(p, otherProfile2):  
        nodes_ids[numLine,i] = pointId
        pointId += 1
      else: # point is not on surface
        assert(True)
        if nodes_ids[numLine,i-1] > -1:
          intersect = find_intersection_point(line[i-1], p, profile, otherProfile1, otherProfile2)
          # calc distance to intersection point
          if calc_distance(line[i-1], intersect) < min_distance:
            # if too small, don't draw previous surface point
            nodes_ids[numLine,i-1] = -1
            pointId -=1
            
  return nodes_ids, pointId

# use bisection algorithm to find intersection point between two profiles
# left and right are tuples/lists with x,y,z coordinates
# assume:
# left is a surface point on profile
# right is a neighbor (same profile) of left but not a surface point
def find_intersection_point(left, right, profile, otherProfile1, otherProfile2):
  dir = profile.direction
  phi = -1.0
  lower = -1.0
  upper = -1.0
  
  # as both left and right lie on same surface line, the angle to (0.5,0.5) is the same
  # find out on which coordinate component we have to perform bisection
  dir1,dir2 = give_normal_plane_axes(dir)
  # convert to degrees and modulo 45 degrees to stay in first quadrant
  phi = degrees(angle_to_center((left[dir1],left[dir2]))) % 45
#   print "phii: ", phi
  phi = radians(phi)
#   print "phi in rad: ", phi
  lower = left[dir]
  upper = right[dir]
    
  return bisection(lower,upper,phi,profile, otherProfile1, otherProfile2)

# for a given interval [lower,upper] and profile angle:
# find a point that is on the intersection of two/three profiles
# a point fulfills this requirement if the interval is small
# and lower is a surface point whereas upper is an inner point     
def bisection(lower,upper,phi,profile, otherProfile1, otherProfile2):
  dir = profile.direction
#   print "dir: ", dir
  l = lower
  u = upper
  
#   print "lower: ", lower, " upper: ", upper
  
  midpoint = -1.0
  midpoint_node = np.zeros(3)
  while abs(u-l) > 1e-4:
    midpoint = 0.5 * (l + u)
    # get coordinates of node with midpoint as one coordinate component (depends on profile direction)
    plane_coordinates = polar_to_cartesian(calc_radius_for_quadrant(profile, midpoint, phi), phi, 0.5)
    # direction of axes of plane normal to profile.direction
    # e.g. we have x profile --> dir1=2, dir2=1 (z,y plane)
    dir1, dir2 =  give_normal_plane_axes(profile.direction)
    assert(dir != dir1 and dir != dir2)
    # assign plane coordinates to correct 3d coordinate components
    midpoint_node[dir] = midpoint
    midpoint_node[dir1] = plane_coordinates[1]
    midpoint_node[dir2] = plane_coordinates[0]
     
    otherDir1 = otherProfile1.direction
    otherDir2 = otherProfile2.direction
    
    # midpoint is inner, take interval from [lower,midpoint]
    if contains_point(midpoint_node, otherProfile1) or contains_point(midpoint_node, otherProfile2):
      u = midpoint
    else:
    # midpoint is on surface, take interval from [midpoint,upper]
      l = midpoint
      
  return midpoint_node
  
# creates triangles between end nodes of same profile where
# we have e.g. a valley with 1 or 2 nodes
# def postprocess_end_nodes(end_nodes,all_nodes_ids,cells):
#   global tris
#   delete_ids = []
#   # a valley with one end node in between:
#   # this node (i,j) the one after next node and (i,j+2)
#   # and the left(right) next one (i-1,j+1) exists, but not next node
#   for node in end_nodes:
#     # don't check on the boundary
#     if node.j > np.size(all_nodes_ids,1)-3 or node.i > np.size(all_nodes_ids,0) - 2:
#        continue
#     # check if the one after the next one exist
#     # nn is node after next node 
#     n = get_end_node_by_grid_coords(node.i, node.j+1, end_nodes)[0]
#     n_inner = True if all_nodes_ids[node.i,node.j+1] >= 0 else False
#     nn = get_end_node_by_grid_coords(node.i,node.j+2,end_nodes)[0]
#     nn_inner = True if all_nodes_ids[node.i,node.j+2] >= 0 else False
#     ln,ln_idx = get_end_node_by_grid_coords(node.i-1,node.j+1,end_nodes)
#     rn,rn_idx = get_end_node_by_grid_coords(node.i+1,node.j+1,end_nodes)
#     
#     if nn and not (n or n_inner) and (ln or rn):
#       # both can be none but if both or not none, something went wrong
#       assert(bool(ln) is not bool(rn) or (ln == None and rn == None))
#       
#       # check which one is in the valley
#       print("found a valley")
#       next = ln if ln is not None else rn
#       next_idx = ln_idx if ln is not None else rn_idx
#       
#       add_triangle(node.id, next.id, nn.id, cells)
#         
#       delete_ids.append(next.id)
#       node.next = nn
#       nn.next = node
#       
#       continue
#     
#     # check valleys that contains 2 end nodes
#     # node after next next node (node after nn)
#     nnn                   = get_end_node_by_grid_coords(node.i, node.j+3, end_nodes)[0]
#     rn,rnn_idx   = get_end_node_by_grid_coords(node.i+1,node.j+1,end_nodes)
#     rnn,rnn_idx = get_end_node_by_grid_coords(node.i+1,node.j+2,end_nodes)
#     ln,ln_idx       = get_end_node_by_grid_coords(node.i-1,node.j+1,end_nodes)
#     lnn, lnn_idx  = get_end_node_by_grid_coords(node.i-1,node.j+2,end_nodes)
#     
#     # next node the node after next node should not be of type end node
#     if nnn and not (n or n_inner)  and not ( nn or nn_inner) and (ln and lnn or rn and rnn):
#       print("found valley containing 2 nodes")
#       # check which one is first in the valley
#       next = ln if ln is not None else rn
#       next_idx = ln_idx if ln is not None else rn_idx
#       # check which one is first in the valley
#       next_2 = lnn if lnn else rnn
#       next_2_idx = lnn_idx if lnn else rnn_idx
#       
#       add_triangle(node.id, next.id,next_2.id, cells)
#       add_triangle(node.id, next_2.id,nnn.id, cells)
#         
#       node.next = nnn
#       nnn.next = node
#       
#       delete_ids.append(next.id)
#       delete_ids.append(next_2.id)
#       
#   # remove all end nodes that were in valleys
#   end_nodes[:] = [v for v in end_nodes if v.id not in delete_ids]
  
# list is list of lists contains surface lines and all respective points
# base used for setting right ids
# nodes is a 3d array with (x,y,z) coords for one profile
# returns array of size 2 x n with end nodes objects for fix_gaps()
def define_triangles(nodes_ids,nodes,cells,dir,vtkArray):
  end_nodes = []
  global tris
  
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
      # divide rectangle into 2 triangles, this one is the first triangle
      if this_id >= 0 and next_id >= 0 and right_id >= 0:
        add_triangle(this_id,right_id,next_id,cells)
        
        vtkArray.SetValue(this_id,dir)
        vtkArray.SetValue(next_id,dir)
        vtkArray.SetValue(right_id,dir)
        
        if j > 0 and (left_id < 0 or prev_id < 0):
          # next, left next, left, left previous, previous, right previous, right, right next
          left_next_id = left_line[j+1]
          left_prev_id = left_line[j-1]
          right_prev_id = right_line[j-1]
          
          neighbors = give_all_neighbor_coords(i, j, nodes, nodes_ids)
          
          end_nodes.append(End_Node(nodes[i,j,:],this_id,dir,i,j,neighbors))
          vtkArray.SetValue(this_id,dir+0.5)
      else:
        if this_id >= 0 and j > 0: # first point in line is not end point
          
          neighbors = give_all_neighbor_coords(i, j, nodes, nodes_ids)
          end_nodes.append(End_Node(nodes[i,j,:],this_id,dir,i,j,neighbors))
          vtkArray.SetValue(this_id,dir+0.5)
          
          # fill gaps that form a triangle
          if right_id >= 0 and right_next_id >= 0:
            add_triangle(this_id,right_id,right_next_id,cells)
      
      # divide rectangle into 2 triangles, this one is the second triangle 
      if right_id >= 0 and next_id >= 0 and right_next_id >=0:
        add_triangle(right_id,right_next_id,next_id,cells)
      
      # fill gaps that form a triangle    
      if this_id >= 0 and right_id < 0 and next_id >= 0 and right_next_id >= 0:
        add_triangle(this_id,right_next_id,next_id,cells)
           
  return end_nodes

# returs a list with indices associated with locations of neighbors in the grid
# e.g. next neighbor has indices (i,j+1), left next (i-1,j+1)
def give_all_neighbor_coords(i,j,nodes,nodes_ids):
  right_line_idx = i+1 if i < nodes_ids.shape[0]-1 else 0
  left_line_idx = i-1 if i > 0 else nodes_ids.shape[0]-1
  # indices of neighborhood
  idx = [(i,j+1),(left_line_idx,j+1),(left_line_idx,j),(left_line_idx,j-1),(i,j-1),(right_line_idx,j-1),(right_line_idx,j),(right_line_idx,j+1)]
  neighbors = []
  
  for v in idx:
    # skip the ones that don't exist
    if v[0] < 0 or v[0] >= nodes_ids.shape[0] or v[1] < 0 or v[1] >= nodes_ids.shape[1]:
      continue
    
    # if v is a surface node id
    if nodes_ids[v[0],v[1]] != -1:
      # for each neighbor, store its coordinates and id
      neighbors.append((nodes[v[0],v[1],:],nodes_ids[v[0],v[1]]))
      
  return neighbors
  
def add_triangle(id1,id2,id3,cells):
  tri = vtk.vtkTriangle()
  tri.GetPointIds().SetId(0, id1)
  tri.GetPointIds().SetId(1, id2)
  tri.GetPointIds().SetId(2, id3)
  cells.InsertNextCell(tri)

# calc distance between two points
def calc_distance(p1,p2):
  return np.sqrt((p1[0]-p2[0])**2 + (p1[1]-p2[1])**2 + (p1[2]-p2[2])**2)

def generate_basecell(args,info,log):
  global res, res_surf_lines, interpolation, logger
  res = args.res
  res_surf_lines = args.res_surf_lines
  interpolation = args.interpolation
  
  logger = log
  global infoXml
  infoXml = info
  
  array = np.ones((res,res,res)) * (-1)
  vec = None
  t = np.linspace(0, 1.0, args.res)
  
  profiles = create_profiles(args,infoXml)
  
  hf = None
  
  if args.verbose != "off":
    hf = plt.figure()
    ha = hf.add_subplot(111, projection='3d')
    ha.set_xlabel('X')
    ha.set_ylabel('Y')
    ha.set_zlabel('Z')
  
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
    nodes_1 = get_surface_lines(map_x, args.res_surf_lines, 0)
    nodes_ids_1, id = find_points_on_surface(nodes_1, profiles[0], profiles[1], profiles[2],id)

    nodes_2 = get_surface_lines(map_y, args.res_surf_lines, 1)
    nodes_ids_2, id = find_points_on_surface(nodes_2, profiles[1], profiles[0], profiles[2], id)
    
    nodes_3 = get_surface_lines(map_z, args.res_surf_lines, 2)
    nodes_ids_3, id = find_points_on_surface(nodes_3, profiles[2], profiles[0], profiles[1], id)
    
    # create vtk cells and points
    cells = vtk.vtkCellArray()
    
    # set scalar info of intersection points to 1.0 to make them visible
    vtkData = vtk.vtkFloatArray()
    vtkData.SetName("intersection")
    vtkData.SetNumberOfValues(id)

    end_nodes_1 = define_triangles(nodes_ids_1,nodes_1,cells,0,vtkData)
    end_nodes_2 = define_triangles(nodes_ids_2,nodes_2,cells,1,vtkData)
    end_nodes_3 = define_triangles(nodes_ids_3,nodes_3,cells,2,vtkData)
    
    dump_end_nodes(end_nodes_1+end_nodes_2+end_nodes_3)
    
    # creates triangles between end nodes of same profile where
    # we have e.g. a valley with 1 or 2 nodes
#     postprocess_end_nodes(end_nodes_1,nodes_ids_1,cells)
#     postprocess_end_nodes(end_nodes_2,nodes_ids_2,cells)
#     postprocess_end_nodes(end_nodes_3,nodes_ids_3,cells)
    
    points = vtk.vtkPoints()
    points.SetNumberOfPoints(id)
    for i,line in enumerate(nodes_1):
      for j in range(len(line)):
        if nodes_ids_1[i,j] >= 0:          
          points.SetPoint(nodes_ids_1[i,j], nodes_1[i,j,0], nodes_1[i,j,1], nodes_1[i,j,2])
        if nodes_ids_2[i,j] >= 0:          
          points.SetPoint(nodes_ids_2[i,j], nodes_2[i,j,0], nodes_2[i,j,1], nodes_2[i,j,2])
        if nodes_ids_3[i,j] >= 0:          
          points.SetPoint(nodes_ids_3[i,j], nodes_3[i,j,0], nodes_3[i,j,1], nodes_3[i,j,2])

        
    id = triangulate_boundary_circles(profiles[0],nodes_ids_1,id,points,cells,vtkData)
    id = triangulate_boundary_circles(profiles[1],nodes_ids_2,id,points,cells,vtkData)
    id = triangulate_boundary_circles(profiles[2],nodes_ids_3,id,points,cells,vtkData)
    
    if not args.skip_surface_gaps:  
      fix_profile_intersection_gaps(profiles,end_nodes_1+end_nodes_2+end_nodes_3, cells)
    
#     for i in range(cells.GetNumberOfCells()):
#       idList = vtk.vtkIdList()
#       cells.GetNextCell(idList)
#     ps, cs = read_vtk("surface.vtp")
#     surface_to_volume_mesh(ps,cs)

    polydata = vtk.vtkPolyData()
    polydata.SetPoints(points)
    polydata.GetPointData().SetScalars(vtkData)
    
    if args.export == "surface_points":
      show_write_vtk(polydata, 1000, "surface_points.vtp")
      
    polydata.SetPolys(cells)
    
    if args.show:
      show_vtk(polydata, 1000, [], True)
    
    if not args.skip_surface_gaps:  
      stlName = args.save if args.save.endswith(".stl") else args.save + ".stl"  
      write_stl(polydata,stlName)  
     
    if args.save_vtp:  
      show_write_vtk(polydata,1000,args.save+".vtp")
  else:
    for i in range(0,3):
      if profiles[i] == None:
        continue
      if args.verbose == 'profile_map' or args.export == 'radius_maps' or args.verbose == "polar_plot":
        create_profile_map(profiles[i], res, args.verbose,args.export == 'radius_maps', ha)
      if args.verbose == 'interpolation':
        y = []
        rad = np.linspace(0,pi/2.0,500)
        for r in rad:
          y.append(calc_radius_for_quadrant(profiles[i], 0.5, r))
        plt.gcf().clear()
        plt.plot(rad,y,linewidth=5.0)
        plt.show()
      if args.target == "volume_mesh" or args.target == "volume_vtk":
        # if basecell is symmetric, calculate only 1/8 and mirror the rest
        symmetric = True if args.x1 == args.x2 and args.y1 == args.y2 and args.z1 == args.z2 else False 
        write_profile_to_array(array, profiles[i], i, symmetric)
      if args.target == "3dlines":
        plot_3dlines(profiles[i], res, args.res_surf_lines, i, ha)
        
  if args.target == '3dlines':
    plt.show()
  
  return array

# creates map with info on profile depending on radius
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
        
  if verbose == "polar_plot":
    plt.gcf().clear()
    ax = plt.axes(polar=True)
    theta = np.linspace(0, 2.0*np.pi,360)
    plt.plot(theta,map[:,0],linewidth=5.0)
    plt.plot(theta,map[:,int(res/4)],linewidth=5.0)
    plt.plot(theta,map[:,int(res/2)],linewidth=5.0)
    plt.rcParams.update({'font.size': 18})
    plt.show()
        
  if verbose == 'profile_map':
    ha.set_xlabel('X')
    ha.set_ylabel('Y')
    ha.set_zlabel('Z')
    X,Y = np.meshgrid(list(range(res)),list(range(360)))
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
  if dir == 0:
    color = "red"
  elif dir == 1:
    color = "blue"
  else:
    assert(dir == 2)
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
  phi = float(funcs[1].angle)  # bisec angle
  assert(phi >= 0 and phi <= np.pi/2.0)
  val = None
  # interpolation is global variable
  assert(interpolation == "linear" or interpolation == "heaviside")
  if interpolation == "linear":
    val = calc_radius_linear(funcs,phi,x,rad)
  else: #heaviside
    val = calc_radius_heaviside(funcs,phi,x,rad)
  
  assert(val is not None)
  return val - 0.5

# rasterize profile functions
# if basecell is symmetric, rasterize only 1/8 of structure
# and mirror the rest
def write_profile_to_array(array,profile,dir,symmetric):
  res = array.shape[0]
  assert(dir >=0 and dir <=2)
  
  map = create_profile_map(profile, res)
  
  bound = res if not symmetric else int(res/2)
  
  for i in range(0,bound):
    for j in range(0,bound):
      for k in range(0,bound):
        y = grid_to_cartesian_coords(j, res)
        z = grid_to_cartesian_coords(k, res)
        valx = (y-0.5)**2 + (z-0.5)**2
        
        phi = angle_to_center((z,y))
        
        p = map[int(phi/np.pi*180),i]
        if (valx <= p*p):
          if dir == 0:
            array[i,j,k] = dir
          if dir == 1:
            array[j,i,k] = dir
          if dir == 2:
            array[k,j,i] = dir
  
  if symmetric:
    # mirror octant
    array[0:bound,0:bound,bound:res] = array[0:bound,0:bound,0:bound][:,:,::-1]
    array[0:bound,bound:res,0:bound] = array[0:bound,0:bound,0:bound][:,::-1,:]
    array[0:bound,bound:res,bound:res] = array[0:bound,0:bound,0:bound][:,::-1,::-1]
    array[bound:res,0:bound,0:bound] = array[0:bound,0:bound,0:bound][::-1,:,:]
    array[bound:res,0:bound,bound:res] = array[0:bound,0:bound,0:bound][::-1,:,::-1]
    array[bound:res,bound:res,0:bound] = array[0:bound,0:bound,0:bound][::-1,::-1,:]
    array[bound:res,bound:res,bound:res] = array[0:bound,0:bound,0:bound][::-1,::-1,::-1] 
  
# find end node object in list with matching grid coordinates
def get_end_node_by_grid_coords(i,j,end_nodes):
  for idx,node in enumerate(end_nodes):
    if node.i == i and node.j == j:
      return node, idx
  
  return None, -1

# vertices contains 3 end nodes that form a triangle
# for each end node, check if any of its neighbors is contained in triangle
# spanned by these 3 end nodes    
def triangle_contains_any_node(vertices):
  for vertice in vertices:
    for n in vertice.neighbors:
      # if we are one of the vertices skip
      if n[1] == vertices[0].id or n[1] == vertices[1].id or n[1] == vertices[2].id:
        continue
      
#       log("check if triangle (" + str(vertices[0].id) + "," + str(vertices[1].id) + "," + str(vertices[2].id) + " contains node " + str(n[1]))
      if triangle_contains_point(vertices,n[0]):
        return True    
  
  return False

# check if projection of point onto triangle plane lies within triangle
# @param vertices(coordinates) defining the triangle
# @param point(coordinates) for testing
# basic idea: express location of point via barycentric coordinates (using vertices)
# if barycentric coordinates(s,t) >= 0 and sum <= 1, triangle contains point
# calculation of barycentric coordinates:
# - define line going from origin to point in parametric form
# - define plane spanned by triangle's vertices
# intersect line and plane to find projection (and parameters)  
def triangle_contains_point(vertices,point):
  origin=np.array((0.5,0.5,0.5))
  # line equation from origin to point in parametric form                             
  # P = origin + k * d
  # plane equation in parametric form
  # P = A + s*u + t*v

  # setting line equation equals plane equation in order 
  # to find parameters for projected point onto plane
  # origin - A = s*u + t*v - k*d
  
  coords = []
  for v in vertices:
    coords.append(v.coords)
  
  A = np.array(coords[0])
  B = np.array(coords[1])
  C = np.array(coords[2])
  
  u = B - A
  v = C - A
  d = point - origin
  
  # intersect line and plane to find projected point (and respective parameters)
  mat = np.transpose(np.asarray((u,v,-d), dtype=float))
  rhs = origin - A 
  
  sol = linsolve_3x3(mat,rhs)
  # sol[0] -> s, sol[1] -> t, sol[2] -> d
#   print("point: " + str(point) + " s: " + str(sol[0]) + " t:" + str(sol[1]) + " d: " + str(sol[2]))
  if sol[0] > 0.0 and sol[1] > 0.0 and sol[0]+sol[1] < 1.0:
    log("projected point inside triangle (" + str(vertices[0].id) + "," + str(vertices[1].id) + ","  + str(vertices[2].id) + ")")
    log("point: " + str(point) + " s: " + str(sol[0]) + " t:" + str(sol[1]) + " sum: " + str(sol[0]+sol[1]))
    return True
  else:
    return False
  
# applying Cramer's rule
# assume A is a numpy array, b a list
def linsolve_3x3(A,b):
  assert(isinstance(A,np.ndarray))
  res = [None] * 3
  
  det = np.linalg.det(A)
  
  assert(det > 0.0 or det < 0.0)
  
  A1 = np.zeros((3,3))
  A1[:,0] = b 
  A1[:,1] = A[:,1]
  A1[:,2] = A[:,2]
  
  A2 = np.zeros((3,3))
  A2[:,0] = A[:,0]
  A2[:,1] = b
  A2[:,2] = A[:,2]
  
  A3 = np.zeros((3,3))
  A3[:,0] = A[:,0] 
  A3[:,1] = A[:,1]
  A3[:,2] = b
  
  # x1 = det(A1)/det(A)
  res[0] = np.linalg.det(A1) / det
  # x2 = det(A2)/det(A)
  res[1] = np.linalg.det(A2) / det
  # x3 = det(A3)/det(A) 
  res[2] = np.linalg.det(A3) / det
  
  return res
    
# @param: vertices defining triangle (coordinates)
# for evaluating quality of triangle, we take the ration of the exradius to twice the inradius
def calc_triangle_ratio(v1,v2,v3):
  d1 = calc_distance(v1, v2) # a
  d2 = calc_distance(v1, v3) # b
  d3 = calc_distance(v2, v3) # c
  
  aspect_ratio = 0
  denom = (d2+d3-d1) * (d1-d2+d3) * (d1+d2-d3)
  if not close(denom,0,1e-10):
    aspect_ratio = d1*d2*d3 / denom
  else: # all three points form a line
    aspect_ratio = 1e6
 
  assert(aspect_ratio >= 1)
  
  return aspect_ratio

# for a given end node, find another end node in 'end_nodes' list
# which is the nearest (euclidean distance) one    
def find_closest_point(ref_node,end_nodes):
  # iterate over all end nodes in list and compare distances to 'ref' node
  min_distance = 1e6
  min_node = None
   
  for node in end_nodes:
    dist = calc_distance(ref_node.coords, node.coords)
    if dist < min_distance:
      min_distance = dist
      min_node = node
  
  log("closest node to id=" + str(ref_node.id) + ": id=" + str(min_node.id))
  return min_node 

# returns list with end nodes in neighborhood with certain radius
# search for neighborhood is performed via kd-tree query  
# @param 'node' of interest
# @param list with all 'end_nodes' 
# @param ball 'radius' for neighborhood
def give_next_end_nodes_in_ball(tree,node,end_nodes):
  radius=6.0/res
  # get indices of nearest neighbors within ball radius
  indices = tree.query_ball_point(node.coords,radius)
  candidates = []
  # extract nearest end nodes
  for i in indices:
    candidates.append(end_nodes[i])
  
  assert(len(candidates) > 0)
  return candidates

# returns nearest end node to 'node' which is on a different profile than 'node'
# @param kd-tree which we are using to find neighbors within a ball
# @param node for which we are searching a neighbor for
# @param end_nodes list for search
def give_next_neighbor_in_other_profile(tree,node,end_nodes):
  candidates = give_next_end_nodes_in_ball(tree, node, end_nodes)
  next = None
  distance = 1e6
  # search for end node on other profile with smallest distance to 'node'
  for c in candidates:
    if c.dir != node.dir:
      d = calc_distance(c.coords, node.coords)
      if d < distance:
        distance = d
        next = c
  
  assert(next is not None)  
  assert(next.dir != node.dir)
  
  return next

# similar to give_next_neighbor_in_other_profile, but here
# we are looking for a neighbor on a given list of end_nodes defining a circle
# (in the same profile!)
# @param kd-tree which we are using to find neighbors within a ball
# @param node for which we are searching a neighbor for
# @param list with search end nodes defining a circle
def give_next_neighbor_on_circle(tree,node,end_nodes):
  candidates = give_next_end_nodes_in_ball(tree, node, end_nodes)
  next = None
  distance = 1e6
  # search for end node on other profile with smallest distance to 'node'
  for c in candidates:
    d = calc_distance(c.coords, node.coords)
    if d < distance and c != node:
      distance = d
      next = c
  
  assert(next is not None)  
  
  return next

# from a list of candidates, choose the third triangle vertice such that the resulting triangle has the smallest aspect_ratio
def give_best_next_neighbor(history,triangles,candidates, vert1, vert2, quality_bound):
  if len(candidates) == 1 and (candidates[0].id == vert1.id or candidates[0].id == vert2.id):
    return candidates[0]
  best_ratio = 1e9
  best_neighbor = None
  
  log("possible candidates: ",linebreak=False)
  assert(candidates is not None)
  for alt in candidates:
    log(str(alt.id) + " ",linebreak=False)  
  log("")
  
  feasible_cands = give_feasible_candidates(history,triangles,candidates,vert1,vert2)
  
  for node in feasible_cands:
    ratio = calc_triangle_ratio(vert1.coords, vert2.coords, node.coords)
    log("triangle (" + str(vert1.id) + "," + str(vert2.id) + "," + str(node.id) + ") ratio=" + str(ratio))
    if ratio < best_ratio:
      best_ratio = ratio
      best_neighbor = node
  
  assert(quality_bound is not None)
  if quality_bound and best_ratio > quality_bound:
    return None
      
  assert(best_neighbor is not None)
  return best_neighbor

def give_next_end_node_on_circle(node,end_nodes,max_i):
#   assert(list_contains_end_nodes(node,end_nodes))
  if len(end_nodes) == 1:
    return node
  
  # find direct neighbors
  next = [v for v in end_nodes if v.i == (node.i+1)%(max_i+1) and v.j == node.j]
  
  assert(len(next) == 1)
 
  return next[0]

def give_best_next_end_node_on_circle(triangles,outer_node,inner_node,outer_end_nodes,inner_end_nodes,first,second):
    
  max_i_outer = get_max_grid_coords(outer_end_nodes)
  max_i_inner = get_max_grid_coords(inner_end_nodes) 
#   log("max_i_outer:" + str(max_i_outer))
#   log("max_i_inner:" + str(max_i_inner))
  
  next_outer = give_next_end_node_on_circle(outer_node,outer_end_nodes,max_i_outer)
  next_inner = give_next_end_node_on_circle(inner_node,inner_end_nodes,max_i_inner)
  
  if len(triangles) > 2:
    if outer_node == first:
      next_outer = first
    if inner_node == second:
      next_inner = second
  
  log("next_outer:" + str(next_outer.id))
  log("next_inner:" + str(next_inner.id))
  
  rat_outer = calc_triangle_ratio(outer_node.coords,inner_node.coords,next_outer.coords)
  rat_inner = calc_triangle_ratio(outer_node.coords,inner_node.coords,next_inner.coords)
  
    
  next = next_outer if rat_outer < rat_inner else next_inner
  log("next: " + str(next.id) + " ratio: " + str(min(rat_outer,rat_inner)))
  
  return next
  
# for a given node on outer circle, return closest node on inner circle
def give_next_end_node_on_inner_circle(outer_node,inner_end_nodes):
  assert(not list_contains_end_nodes(outer_node, inner_end_nodes))
  return find_closest_point(outer_node, inner_end_nodes)

# if we don't have alternatives (other candidates for next end node):
# take care of triangle if detected a bad triangle (aspect ratio too big, contains projection of neighbors, ...)
# corresponds to going one step backwards in algorithm
# remove latest triangle and check second best candidate    
def handle_bad_triangle(history,triangles,end_nodes,tree,quality_bound):
  next = None
  # if we don't have alternatives for neighbor end nodes
#   if alternatives is None or len(alternatives) == 0:
  while next is None:
    # removes triangles as long as we don't have candidates for these triangles
    if not pop_triangles(triangles):
      return False

    tri = triangles[-1] 
    alternatives = tri.other_candidates
    a = tri.edge[0]
    b = tri.edge[1]
    log("active edge: " + str(a.id) + "," + str(b.id))
    tri.next = None
    next = give_best_next_neighbor(history,triangles,alternatives, a, b, quality_bound)
    if next is not None:
      # next is good, so it is not a candidate anymore
      tri.other_candidates = [v for v in tri.other_candidates if v.id != next.id]
      # reset all information with new next node  
      tri.edge[0] = a if next.dir!=a.dir else b
      assert(tri.edge[0].dir != next.dir) 
      tri.edge[1] = next
      tri.next = next
      tri.vertices[2] = next
      best_ratio = calc_triangle_ratio(a.coords, b.coords, next.coords)
      log("set new triangle " + str(tri.vertices[0].id) + "," + str(tri.vertices[1].id) + "," + str(tri.vertices[2].id) + " with ratio " + str(best_ratio))
      log("set active edge to (" + str(tri.edge[0].id) + "," + str(tri.edge[1].id) + ")")
#       triangles.append(tri)
      return True
    else:
      triangles.append(tri)
      # best neighbor not found, so candidates are invalid
      triangles[-1].other_candidates = []
  return False  

def add_good_triangle(triangles, a, b, next_cands, next):
  
# decide whether a or b is member based on triangle ratios
# and not profile directions as a.dir might be equal b.dir
  edge_0 = a if a.dir != next.dir else b
  alternatives = [v for v in next_cands if v.id != next.id]
  triangles.append(Marching_Triangle([a, b, next], edge_0, next, alternatives))
  best_ratio = calc_triangle_ratio(a.coords, b.coords, next.coords)
  assert(best_ratio < 30)
  log("next:" + str(next.id))
  log("created triangle " + str(a.id) + "," + str(b.id) + "," + str(next.id) + "  edge: " + str(triangles[-1].edge[0].id) + "," + str(triangles[-1].edge[1].id) + " with ratio " + str(best_ratio))
  log("")

def check_next_triangle(history,triangles,end_nodes,tree,quality_bound,initial_edge=None):
  # take active edge from last triangle
  active_edge = triangles[-1].edge
  a = active_edge[0]
  b = active_edge[1]
  
#   if a.id == 2042 and b.id == 520:
#     return False
  
  # don't compare last element in triangles to avoid self-comparison
  if edge_already_connected(history+triangles[:-1],active_edge):
    return False

  log("\na: " + str(a.id))
  log("b: " + str(b.id))
  
  next_cands = give_next_end_nodes_in_ball(tree, a, end_nodes)
  next = give_best_next_neighbor(history,triangles,next_cands,a,b,quality_bound)
  if next is not None:
    add_good_triangle(triangles, a, b, next_cands, next)
    return True
  else: # if both candidates are bad, we have go back to previous triangle an try alternative candidate
    # dummy, will be removed in handle_bad_triangle    
    triangles.append(Marching_Triangle([a,b,b],a,b,[]))   
    # returns False if no alternatives left at all  
    return handle_bad_triangle(history,triangles,end_nodes,tree,quality_bound)
      
# pop triangles that have no alternative candidates    
def pop_triangles(triangles):
  log("pt: number of triangles: " + str(len(triangles)))
  log("pt: number of candidates: " + str(len(triangles[-1].other_candidates))) 

  while len(triangles[-1].other_candidates) == 0:
    vertices = triangles[-1].vertices
    
    if len(triangles) == 1:
      return False
    
    log("number of triangles: " + str(len(triangles)))
    log("remove triangle: " + str(vertices[0].id) + "," + str(vertices[1].id) + "," + str(vertices[2].id))
    triangles.pop()
    
  if len(triangles) == 1:
    return False
  vertices = triangles[-1].vertices
  if triangles[-1].edge[0] == triangles[-1].next:
    triangles[-1].edge[0] = [v for v in triangles[-1].vertices if v.id != triangles[-1].next.id and v.id != triangles[-1].edge[1].id][0]
  else: #triangles[-1].edge[1] == triangles[-1].next
    triangles[-1].edge[1] = [v for v in triangles[-1].vertices if v.id != triangles[-1].next.id and v.id != triangles[-1].edge[0].id][0]
      
  log("previous triangle: " + str(vertices[0].id) + "," + str(vertices[1].id) + "," + str(vertices[2].id))
  log("set edge to " + str(triangles[-1].edge[0].id) + "," + str(triangles[-1].edge[1].id))
  log("with other_candidates: ",linebreak=False)
  for cand in triangles[-1].other_candidates:
    log(str(cand.id) + " ",linebreak=False)
  log("")
    
  return True
 
# given 3 lists with end nodes of all 3 profiles
# using a marching algorithm, we connect theses nodes step by step to triangles
# start with one given triangle defined by end nodes of two profiles (assume each profile has a different color'
# we call the edge that connects two differently colored nodes an 'active' edge
# when a triangle was defined, continue with active edge and search for possible candidates (for new triangle)
# in the neighborhood of the two edge vertices              
def fix_profile_intersection_gaps(profiles,end_nodes,cells):
  # find profiles with biggest radius
  start_dir = 0
  start_radius = 0
  for p in profiles:
    if p.radius_left > start_radius:
      start_dir = p.direction
      start_radius = p.radius_left
  
  # save all valid created triangles
  # we need to store it to make sure no triangle is defined more than once
  triangles_history = []    
  # start from profile with biggest radius
  nodes = end_nodes
#   nodes = [v for v in end_nodes if v.dir == start_dir]
  for n in nodes:
    log(str(n))
  assert(nodes)
  tree = build_tree(end_nodes)
  count = 0
  for n in nodes:
    start_node = n
    
    next = give_next_neighbor_in_other_profile(tree,start_node,end_nodes)
    
    log("\nstarting with " + str(start_node.id) + "," + str(next.id))
    log("start: " + str(start_node.id) + " " + str(start_node.coords))
    log("next: " + str(next.id) + " " + str(next.coords))
    
    out, out2 = start_triangulation(triangles_history,start_node,next,end_nodes,tree,cells)
    if out2:
      break
    if not out:
      log("skip " + str(n.id))
    
#     count += 1
#     if count == 20:
#       break 
    
  for tri in triangles_history:
    verts = tri.vertices
    add_triangle(verts[0].id, verts[1].id, verts[2].id, cells)
        
  check_duplicated_triangles(triangles_history)    
        
# generate points on the bounding circle (left and right) of a profile
# radius is the current radius of a circle where we already have the boundary points
# step is the step size we use to reduce the radius 
# the reduced radius is then used to sample end nodes
# when introducing new nodes in the mesh, we also have to assign them an unique id
# id+1 is the current number of points in the surface mesh
def generate_end_nodes_in_circle(profile,arc_length,radius,points,vtkData,id,set_id=True):
  global res
  dir = profile.direction
  
  # number of points change with circle radius
  tmp_res = 0
  # in case step size is so big such that radius -= step becomes negative
  if radius < 0:
    radius = 1e-4
  
  step_angle = degrees(arc_length / radius)
  
  log("arc length:" + str(arc_length))
  log("radius: " + str(radius))
  log("step_angle: " + str(step_angle))
  
  nodes_left = []
  nodes_right = []
  
  for line,angle in enumerate(np.arange(0,360,step_angle)):
    if angle > 0 and (360.0 - angle) < step_angle/2.0:
      break
    tmp_res += 1
    x,y = polar_to_cartesian(radius, radians(angle), 0.5)
    coords_left = np.zeros(3)
    coords_right = np.zeros(3)
    idx_y, idx_x = give_normal_plane_axes(dir)
    idx_z = [v for v in (0,1,2) if v != idx_x and v != idx_y][0]
    assert(idx_z != idx_x and idx_z != idx_y)
    coords_left[idx_x] = x
    coords_left[idx_y] = y
    coords_left[idx_z] = 0.0
    
    coords_right[idx_x] = x
    coords_right[idx_y] = y
    coords_right[idx_z] = 1.0 
    
    # assign any other direction to end node as triangulation routine assumes two different directions  
    nodes_left.append(End_Node(coords_left,id,dir,line,0,[]))
    if set_id:
      id += 1  
      vtk_id = points.InsertNextPoint(coords_left)
      
      vtkData.InsertValue(vtk_id,-1)
      if vtk_id != nodes_left[-1].id:
        print(("vtk_id" + str(vtk_id) + " id: " + str(nodes_left[-1].id)))
      assert(vtk_id == nodes_left[-1].id)  
    
    # assign any other direction to end node as triangulation routine assumes two different directions
    nodes_right.append(End_Node(coords_right,id,dir,line,res-1,[]))
    if set_id:
      id += 1  
      vtk_id = points.InsertNextPoint(coords_right)
      vtkData.InsertValue(vtk_id,-1)
      if vtk_id != nodes_right[-1].id:
        print(("vtk_id" + str(vtk_id) + " id: " + str(nodes_right[-1].id)))
      assert(vtk_id == nodes_right[-1].id)   
    # if radius is too small, we only have one point left
    if radius < 1e-3:
      break    
  
  return nodes_left, nodes_right, id, tmp_res

# meshes the left and right boundary circles of a profile with triangles
# for this, we need to add additional points on smaller circles within starting circle
# we need nodes_ids to determine ids of points on starting circle 
def triangulate_boundary_circles(profile,nodes_ids,id,points,cells,vtkData):
  radius = profile.radius_left
  arc_length = radius * radians(360.0/np.size(nodes_ids,0))
  outer_points_left, outer_points_right,id, blah = generate_end_nodes_in_circle(profile,arc_length,radius,points,vtkData,id,False)
  set_correct_point_ids(outer_points_left, outer_points_right,nodes_ids)
  
  inner_points_left = []
  inner_points_right = []
  triangles = []
  step = 2.0/res
  tmp_res = 0
  
  # create points on circles lying in the same plane as original circle
  while len(outer_points_left) > 1 and len(outer_points_right) > 1:
    radius -= step
    inner_points_left, inner_points_right, id, tmp_res = generate_end_nodes_in_circle(profile,arc_length,radius,points,vtkData,id)
    
    # previous_points_* are the points on the outer circle
    # left: x=0, right: x = 1
    for lists in [(outer_points_left,inner_points_left),(outer_points_right,inner_points_right)]:
      # set fictional profile affiliation to distinguish which
      # nodes are on outer circle and which ones are on inner circle
      # in fact, they are all in the same profile
      for p in lists[0]:
        p.dir = 5
      for p in lists[1]:
        p.dir = 6
        
      start = lists[0][0]
      next = give_next_end_node_on_inner_circle(start,lists[1])
      log("\nstart:" + str(start.id))
      log("next:" + str(next.id))
      start_boundary_circle_triangulation(start,next,lists[0],lists[1],cells)
     
    outer_points_left = inner_points_left
    outer_points_right = inner_points_right
      
  return id

# this is called after first call of generate_end_nodes_in_circle
# as that method returns end nodes on the outmost circles, but with wrong ids
# here we set the right (original) ids     
def set_correct_point_ids(nodes_left, nodes_right, nodes_ids):
  for node in nodes_left:
    assert(nodes_ids[node.i,node.j] > -1)
    node.id = nodes_ids[node.i,node.j]

  for node in nodes_right:
    assert(nodes_ids[node.i,node.j] > -1)
    node.id = nodes_ids[node.i,node.j]  
    
# checks if node is in list    
def list_contains_end_nodes(node,list):
  for n in list:
    if n.id == node.id:
      return True
  return False

# starts triangulation algorithm for given initial triangle
# all three nodes must be of type End_Node
# next_node is on same 'circle/direction' as start_node
# other_node is on other 'circle/direction' as start_node
# @param history: list of triangles that we don't want to modify anymore,
# need this list to make sure no triangle is defined more than once
# @param start: starting node
# @param next: node that lies on other profile than start
# @param end_nodes: list with all end_nodes
# @param kd-tree: needed for neigborhood search
# @param cells: vtk list storing all valid triangles
def start_triangulation(history,start,next,end_nodes,tree,cells):
  triangles = [] # saves all triangles found by marching
  initial_edge = [start,next]
  
  other_cands = give_next_end_nodes_in_ball(tree, start, end_nodes)
  other = give_best_next_neighbor(history,triangles,other_cands,start,next,100)
  
  # if first triangle already exists, go back to loop in fix_profile_intersection_gaps
  for tri in history:
    #if other is None or edge_already_connected(history, (start,next)):
    if other is None:
      #log("triangle " + str(start.id) + "," + str(next.id) + "," + str(other.id) + " already exists.")
      log("edge (" + str(start.id) + "," + str(next.id) + ") already exists")
      return False, False

  assert(other is not None)
  
  log("other: " + str(other))

  add_good_triangle(triangles, start, next, [], other)
  
  quality_bound = 5  
  end = False
  stop = False
  while not end and quality_bound < 30: 
    # 3 possibilities: - False, nothing to go back we need to increase quality bound
    # - True: Process next triangle
    # - True: Reached the end
    failed = not check_next_triangle(history,triangles, end_nodes, tree, quality_bound, initial_edge)
    # check if we have just created a triangle where new active edge is 
    # identical to active edge of very first triangle --> we are finished!
    if failed:
      # in this case, we want to output the so far valid triangles for debugging purposes
      if len(triangles) > 1:
        stop = False
        break
      # Start again with more tolerant quality bound
      quality_bound *= 1.5
      log("increased quality_bound to " + str(quality_bound))
    else:
      edge = triangles[-1].edge
      # if we have reached the very first edge, we're done
      if len(triangles) > 2 and (edge[0].id == initial_edge[0].id and edge[1].id == initial_edge[1].id) or (edge[0].id == initial_edge[1].id and edge[1].id == initial_edge[0].id):
        log("filled")
        end = True
      # find out if we have reached an edge which is already part of another triangle
      # if this is the case, we are done for this given starting node
      # this is the case when we have gaps going through all 3 profiles 
      # (e.g x1=y1=z1)
      # triangles[:-1]:don't compare edge with its own triangle
      if len(triangles) > 2 and edge_already_connected(history+triangles[:-1],edge):
        log("filled")
#         count += 1
#         print("count",count)
#         if count == 2:
#           history += triangles
#           return False,True
        break 
      
  if quality_bound >= 30:
    Exception("Quality bound exceeded limit of 30!")
    
  history += triangles
  
  return True, stop  

# connects points on an outer (larger radius) circle with 
# points on an inner circle(smaller radius) by using triangles; 
# similar to start_triangulation, but simpler as trianglulation
# here is straightforward and we don't need to check quality bounds
# @param first is first end node of triangulation and lies on outer circle
# @param second is second end node of triangulation and lies on inner circle
def start_boundary_circle_triangulation(first,second,outer_end_nodes,inner_end_nodes,cells):
  triangles = [] # saves all triangles found by marching
  
  end = False
  # we know that first lies on outer circle and next on inner circle
  outer_dir = first.dir 
  inner_dir = second.dir
  outer_node = first
  inner_node = second
  count = 0
  next = None # third vertex of triangle
  while not end: 
    log("\nouter:" + str(outer_node.id))
    log("inner:" + str(inner_node.id))
    
    next = give_best_next_end_node_on_circle(triangles,outer_node, inner_node, outer_end_nodes, inner_end_nodes,first,second)
    assert(next is not None)
    log("next:" + str(next.id))
    # find node that is first vertex of next triangle
    # next is second vertex of next triangle
    add_good_triangle(triangles, outer_node, inner_node, [], next)
    
    edge = triangles[-1].edge
    if len(triangles) > 2 and (edge[0].id == first.id and edge[1].id == second.id) or (edge[0].id == second.id and edge[1].id == first.id):
      log("filled")
      end = True
      
    outer_node = next if next.dir == outer_dir else edge[0]
    inner_node = next if next.dir == inner_dir else edge[0]
    
    assert(next.dir != edge[0].dir)
    
  for triangle in triangles:
    verts = triangle.vertices
    add_triangle(verts[0].id, verts[1].id, verts[2].id, cells)
  
  return True

def get_max_grid_coords(end_nodes):
  max_i = 0
  for n in end_nodes:
    if n.i > max_i:
      max_i = n.i
  
  return max_i

def read_vtk(filename):
  reader = vtk.vtkXMLPolyDataReader()
  reader.SetFileName(filename)
  reader.Update()
  polydata = reader.GetOutput()
  points = vtk_to_numpy(polydata.GetPoints().GetData())
  cells = vtk_to_numpy(polydata.GetPolys().GetData())
  
  points = points.tolist()
    
  print(points)
  print(cells) 

  return points,cells

def surface_to_volume_mesh(ps,cs):
  mesh_info = MeshInfo()
  
  mesh_info.set_points(ps)
  mesh_info.set_facets(tris)
  
  mesh = build(mesh_info)
  print("Mesh Points:")
  for i, p in enumerate(mesh.points):
      print(str(i) + "," + str(p) + "\n")
  print("Point numbers in tetrahedra:")
  for i, t in enumerate(mesh.elements):
      print(str(i) + "," + str(t) + "\n")
  mesh.write_vtk("test.vtk")

# helper function for calc_radius_for_quadrant
# return linear interpolation between principal spline and bisec
# @param funcs: array with 3 entries: spline1, bisec, spline2
# @param rad: angle in radians
# @param x: cartesian x-coordinate
def calc_radius_linear(funcs,phi,x,rad):
  assert(phi >= 0 and phi <= np.pi/2.0)
  if rad <= phi:
    alpha = 1.0/phi * rad # scale section between 0 and 1
    return  (1 - alpha) * funcs[0].eval(x) + alpha * funcs[1].eval(x)
  else : # rad <= np.pi/2.0
    alpha = (rad-phi) / (np.pi/2.0-phi)  # scale section between 0 and 1
    return  (1-alpha) * funcs[1].eval(x) + alpha * funcs[2].eval(x)

def close(x,ref,eps=1e-6):
  return abs(x-ref) < eps

# helper function for calc_radius_for_quadrant
# return heaviside interpolation between principal spline and bisec
# see calc_radius_linear for params
def calc_radius_heaviside(funcs,phi,x,rad):
  # get beta and eta for heaviside function from bisec
  beta = funcs[1].heaviside.beta
  eta = funcs[1].heaviside.eta
    
  eps = 1e-5 #for numerical comparison
    
  # a + c * tanh(...) in [0,1]
  assert(calc_tanh(beta, eta, 1) <= 1)
  assert(calc_tanh(beta, eta, 0) >= 0)
  c = 1.0 / (calc_tanh(beta, eta, 1) - calc_tanh(beta, eta, 0))
  a = - c * calc_tanh(beta, eta, 0)
  if abs(a+c*calc_tanh(beta, eta, 0)) > eps:
    print(a+c*calc_tanh(beta, eta, 0))
  assert(abs(a+c*calc_tanh(beta, eta, 0)) < eps) 
  assert(abs(a+c*calc_tanh(beta, eta, 1)) >= 1-eps)
    
  if rad <= phi:
    alpha = 1.0/phi * rad # scale section between 0 and 1
    scale = a+c*calc_tanh(beta, eta, 1-alpha)
    assert(close(a+c*calc_tanh(beta, eta, 0),0))
    assert(close(a+c*calc_tanh(beta, eta, 1),1))
    assert(scale >= -eps and scale <= 1.0 + eps)
    assert(funcs[0].eval(x) >= 0.5-eps and funcs[1].eval(x) >= 0.5-eps)
    
    return scale * funcs[0].eval(x) + (1-scale) * funcs[1].eval(x)
  else:
    alpha = (rad-phi) / (np.pi/2.0-phi)  # scale section between 0 and 1
    scale = a+c*calc_tanh(beta, eta, alpha)
    
    assert(scale >= -eps and scale <= 1.0 + eps)
    assert(funcs[1].eval(x) >= 0.5-eps and funcs[2].eval(x) >= 0.5-eps)
    return (1-scale) * funcs[1].eval(x) + scale * funcs[2].eval(x)

# helper function for building a kd-tree
# @param list containing end nodes objects
# @returns scipy.spatial.kdtree object
def build_tree(end_nodes):
  # copy all end node coordinates to one 2d array
  points = np.ones((len(end_nodes),3)) * (-1.0)
  for idx,n in enumerate(end_nodes):
    points[idx] = n.coords
  # build tree
  return spatial.cKDTree(points)#

# calculates barycenter/centroid of a triangle
# @param coordinates of the three vertices defining that triangle
def calc_triangle_barycenter(vert1,vert2,vert3):
  return 1.0/3.0 * (np.array(vert1)+np.array(vert2)+np.array(vert3))

# returns 3 points
# 1st point lies on median between triangle's barycenter and vert1 and distance from vert1 to that point is eps
# 2nd point lies on median between triangle's barycenter and vert2 and distance from vert1 to that point is eps
# 3rd point lies on median between triangle's barycenter and vert3 and distance from vert1 to that point is eps
def calc_points_on_triangle_medians(vert1,vert2,vert3,eps=1e-3):
  centroid = calc_triangle_barycenter(vert1, vert2, vert3)
  
  # directional vector from vert1 to centroid 
  d1 = centroid - vert1 
  p1 = vert1 + eps * d1
  # directional vector from vert2 to centroid
  d2 = centroid - vert2 
  p2 = vert2 + eps * d2
  # directional vector from vert3 to centroid
  d3 = centroid - vert3 
  p3 = vert3 + eps * d3
  
  return p1,p2,p3

# check if triangle defined by its 3 vertices overlaps with
# other triangles in list
# @param list of triangles that are o.k. so far
# @param vertices(end nodes) defining triangle that we want to check
# a triangle overlaps a second one if three of its inner points,
# projected on plane of second triangle, lie inside second triangle
# the three inner points are chosen such that they lie on the
# three medians of the triangle and are close to the three vertices
# (distance controlled by eps)
def triangle_overlap_others(triangles,vertices):
  assert(len(vertices) == 3)
  # points for testing
  sample1, sample2, sample3 = calc_points_on_triangle_medians(vertices[0].coords, vertices[1].coords, vertices[2].coords,eps=0.01)
  bary1 = calc_triangle_barycenter(vertices[0].coords, vertices[1].coords, vertices[2].coords)
  # don't check last triangle as it might be the one we want to modify
  for tri in triangles:
    if tri.next is None:
      log("skip overlap check for (" + str(tri.vertices[0].id) + "," + str(tri.vertices[1].id) + "," + str(tri.vertices[2].id) + ") and (" + str(vertices[0].id) + "," + str(vertices[1].id) + "," + str(vertices[2].id) + ")")
      continue
    
    # if triangles are too far away from each other, we don't care about overlapping projections
    bary2 = calc_triangle_barycenter(tri.vertices[0].coords, tri.vertices[1].coords, tri.vertices[2].coords)
    if (calc_distance(bary1, bary2) > 0.3):
#       log("distance = " + str(calc_distance(bary1, bary2)))
      continue
    
#     log("check if triangle (" + str(vertices[0].id) + "," + str(vertices[1].id) + "," + str(vertices[2].id) + ") overlap other with (",linebreak=False)
#     for vert in tri.vertices:
#       log(str(vert.id) + ",", linebreak=False)
#     log("")
    
    if triangle_contains_point(tri.vertices,sample1) or triangle_contains_point(tri.vertices,sample2) or triangle_contains_point(tri.vertices,sample3):
      return True
    
  # no overlapping found  
  return False

def dump_end_nodes(list):
  out = open("end_nodes.txt","w")
  for n in list:
    out.write(str(n) + "\n")
  
  out.close()

def check_duplicated_triangles(triangles):
  for i,this in enumerate(triangles):
    for j,other in enumerate(triangles):
      if i == j:
        continue
      if this.is_triangle(other.vertices[0],other.vertices[1],other.vertices[2]):
        raise Exception("found duplicated triangle (" + str(this.vertices[0].id) + "," + str(this.vertices[1].id) + "," + str(this.vertices[2].id) + ")")

# returns a list of end nodes that can form feasible triangles with vert1 and vert2
# @param history is a list with Marching Triangles that should not be modified anymore
# @param triangles is list with current Marchin Triangles that we are working on, 
#        e.g. new triangles are added or bad ones are deleted
# a candidate 'cand' is feasible when all following criteria are fulfilled:
# - if cand is on same profile as vert1 or vert2, it must be a direct neighbor
# - triangle defined by cand, vert1, vert2 must not already exist and must not overlap other already created triangles
# - vert1, vert2 and cand must differ from each other (not the same id) 
def give_feasible_candidates(history,triangles,candidates,vert1,vert2):
  cands = []
  for node in candidates:
    if node.id == vert1.id or node.id == vert2.id:
      continue
    # if not direct neighbor; take care of modulo res as i = 0 and i = res-1 are also neighbors!
    if node.dir == vert1.dir and (abs(node.i - vert1.i) > 1 and not abs(node.i - vert1.i) == res-1 or abs(node.j - vert1.j) > 1 and not abs(node.j - vert1.j) == res-1):
      continue
    if node.dir == vert2.dir and (abs(node.i - vert2.i) > 1 and not abs(node.i - vert2.i) == res-1 or abs(node.j - vert2.j) > 1 and not abs(node.j - vert2.j) == res-1):
      continue
    if len(triangles) > 0 and triangles[-1].is_triangle(vert1,vert2,node):
      continue
    if node.id == vert1.id or node.id == vert2.id:
      continue
    if triangle_overlap_others(history+triangles,(vert1,vert2,node)):
      continue
    # found a feasible candidate
    cands.append(node)
    log("feasible candidate for " + str(vert1.id) + "," + str(vert2.id) + ": " + str(node.id))
    
  return cands
  
# find out if we have reached an edge which is already part of another triangle
# if this is the case, we are done for this given starting node
# this is the case when we have gaps going through all 3 profiles 
# (e.g x1=y1=z1)
def edge_already_connected(triangles,edge):
  for tri in triangles:
    # copy ids in lists
    ref = [tri.vertices[0].id,tri.vertices[1].id,tri.vertices[2].id]
    comp = [edge[0].id,edge[1].id]
    # check if all items in comp are also included in ref
    if all([item in ref for item in comp]):
      log("edge " + str(comp) + " already connected to triangle " + str(ref))
      
      return True
    
  return False