#!/usr/bin/env python
import numpy as np
from numpy import outer
import sympy.solvers
import sys
import math
import vtk
from vtk.util.numpy_support import vtk_to_numpy
from matviz_vtk import *
import matplotlib
matplotlib.use('tkagg')
from matplotlib import pyplot as plt
from scipy import interpolate
from mpl_toolkits.mplot3d import Axes3D
from sympy import Symbol, symbols
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

class End_Node():
  coords = None
  connections = None
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
    self.connections = set() # array of end_nodes that form a triangle together with this end node
    self.dir = dir
    self.neighbors = neighbors
    
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
#     if len(self.next_node) > 0:
#       print self.next_node[0]
#     assert(len(self.next_node) == 0)
#     self.next_node.append(edge_node_2)
    self.next_node = edge_node_2
    
  def __str__(self):
    return "vertices: " + str(self.vertices[0].id) + " "  + str(self.vertices[1].id) + " " + str(self.vertices[2].id)
        
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
    
    global infoXml
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
    return self.explicit(x)[()]
  
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
      print("No t found where dx/dy = 1")
    
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
    return self.eval_t(t1)
  
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
    g1 = calc_tanh(self.beta, self.eta, 0.5)
    a = (h1-h0) / (g1 - g0)
    b = h0 - a * g0
    return a * calc_tanh(self.beta, self.eta, x) + b

# shifted tanh, returns values are between 0 and 1 for x \in [0,0.5]
def calc_tanh(beta,eta,x):
  return 1.0 - 1.0 / (np.exp(2.0*beta*(2*x-eta)) + 1)
   
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

# def contains_point(id_x,y,z,map):
#      
#   valx = (y-0.5)**2 + (z-0.5)**2
#   phi = angle_to_center((y,z))
#    
#   r = map[int(degrees(phi)),id_x]
#      
#   if (valx - r*r < 1e-3 ):
#     return True
#      
#   return False

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
#   print "(x,z): ",x,z, " angle:", degrees(phi)
  r = calc_radius_for_quadrant(profile, y, phi)
  val = (x-0.5)**2 + (z-0.5)**2
  
#   print "p: ", p, " angle: ", phi, " radius: ", r
  
#   print "radius: ",r," val:",val
  
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
    
    global infoXml
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
    self.heavi = None
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
    
    global infoXml
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
    rhs = np.array([ly[0], 1.0, height[0], 0])
    
    sol = np.linalg.solve(A, rhs)
    
    cubic = np.poly1d(sol[::-1]) # cubic polynomial
    
    if infoXml:
      infoXml.write('        <polynomial order="cubic" a0="' + str(sol[0]) + '" a1="' + str(sol[1]) + '" a2="' + str(sol[3]) +'"/>\n')
      infoXml.write('      </bicubic>\n')
      infoXml.write('      <bspline rad1 = "' + str(x1) + '" rad2="' + str(height[0]) + '" bend="' + str(bend) + '">\n')
    
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
    self.heaviside = Heaviside(beta, eta, x1, height[0])
#     print("f(0.5):",self.heaviside.eval(0.5))
#     print("bicubic(0.5):",self.eval_bicubic(0.5))
#     print("bspline(0.5):",self.eval_spline(0.5))
    
    if force:
      self.type = force
    else:
      # to check if bicubic has under/overshooting when point p is much lower than point b
      if p[1] >= b[1] + 1e-3:
        self.type = "bicubic"
      # in case function composed of b-spline and cubic function has undershoot  
      # in case b-spline has no undershoot (point p is not below bspline(x=0))
      elif height > 0.5 + x1/2.0:  
        self.type = "bspline"
      # in case we have undershooting for biqua AND for spline
      else:
        self.type = "linear"
      
    if infoXml:
      infoXml.write('      <selection type="' + str(self.type) + '" angle="' + str(degrees(self.angle[0])) + '"/>\n')
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
    if self.type == "bicubic":
      return self.eval_bicubic(x)
    elif self.type == "bspline":
      return self.eval_spline(x)
    elif self.type == "heaviside":
      return self.eval_heaviside(x)
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
    assert(args.profile == "linear" or args.profile == "spline")
    assert (dir == 0 or dir == 1 or dir == 2)
    self.bisec_angle = -1
    self.direction = dir
    self.type = args.profile
    # 0th entry: function for 0 degree; 1st entry: function for bisec; 2nd entry: function for 90 degree
    self.functions = [None] * 3
    # depending on profile, store the radii of the two boundary circles
    self.radius_left = 0
    self.radius_right = 0

    global infoXml    
    if infoXml:  
      infoXml.write('  <profile dir="' + str(dir) + '">\n')
      
    if self.type == "linear":
      self.functions[0] = Linear_1D(args.x1, args.x2)
      self.functions[1] = self.splines[0]
      self.functions[2] = self.splines[0]
    else: # 'spline' case
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
  global res, res_surf_lines 
  res = args.res
  res_surf_lines = args.res_surf_lines
  profiles = [None]*3 # x-,y-,z-part
  
  if not args.skip_x:
    profiles[0] = Profile(args,0)
    
  if not args.skip_y:
    profiles[1] = Profile(args,1)
    
  if not args.skip_z:
    profiles[2] = Profile(args,2)
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
# def find_points_on_surface(nodes, dir, otherMap1, otherMap2, pointId):
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
#       print "p: ", p
      # detect crossing border from inner to outer:
      if i > 0 and nodes_ids[numLine,i-1] == -1:
        # it is possible that we already performed this check on
        if contains_point(line[i-1], otherProfile1) or contains_point(line[i-1], otherProfile2):
          intersect = find_intersection_point(line[i-1], p, profile, otherProfile1, otherProfile2)
          if calc_distance(p, intersect) < min_distance:
#             print "omit ", p
            nodes_ids[numLine,i] = -1
            continue

      # check if point is contained in other profiles
      # assume we start with a point that is a surface point
      if not contains_point(p, otherProfile1) and not contains_point(p, otherProfile2):  
        nodes_ids[numLine,i] = pointId
        pointId += 1
      else: # point is not on surface
        assert(True)
#         print "id: ", pointId-1
        # detect border from outer to inner
        if nodes_ids[numLine,i-1] > -1:
          intersect = find_intersection_point(line[i-1], p, profile, otherProfile1, otherProfile2)
#           print "found intersection point ", intersect, " for profile ", profile.direction
#           print "distance of last node to intersection: ", calc_distance(line[i-1], intersect)  
          # calc distance to intersection point
          if calc_distance(line[i-1], intersect) < min_distance:
#             print "removed point ", line[i-1]
            # if too small, don't draw previous surface point
            nodes_ids[numLine,i-1] = -1
            pointId -=1
#             print "omit ", line[i-1]
            
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
#     print "midpoint: ", midpoint
#     print "coords_grad_1: ", profile.functions[0].calc_coords_grad_1()
#     print "phi: ", degrees(phi), " radius: ", calc_radius_for_quadrant(profile, midpoint, phi)[0]
    # get coordinates of node with midpoint as one coordinate component (depends on profile direction)
    plane_coordinates = polar_to_cartesian(calc_radius_for_quadrant(profile, midpoint, phi)[0], phi, 0.5)
#     print "plane_coordinates: ", plane_coordinates
    # direction of axes of plane normal to profile.direction
    # e.g. we have x profile --> dir1=2, dir2=1 (z,y plane)
    dir1, dir2 =  give_normal_plane_axes(profile.direction)
    assert(dir != dir1 and dir != dir2)
    # assign plane coordinates to correct 3d coordinate components
    midpoint_node[dir] = midpoint
    midpoint_node[dir1] = plane_coordinates[1]
    midpoint_node[dir2] = plane_coordinates[0]
     
#     print "midpoint_node: ", midpoint_node
    
    otherDir1 = otherProfile1.direction
    otherDir2 = otherProfile2.direction
    
    # midpoint is inner, take interval from [lower,midpoint]
    if contains_point(midpoint_node, otherProfile1) or contains_point(midpoint_node, otherProfile2):
#       print "inner"
      u = midpoint
    else:
#       print "outer"
    # midpoint is on surface, take interval from [midpoint,upper]
      l = midpoint
      
#     print "continue with interval [", l, u, "], length:", abs(u-l)
    
  return midpoint_node
  
# creates triangles between end nodes of same profile where
# we have e.g. a valley with 1 or 2 nodes
def postprocess_end_nodes(end_nodes,all_nodes_ids,cells):
  global tris
  delete_ids = []
  # a valley with one end node in between:
  # this node (i,j) the one after next node and (i,j+2)
  # and the left(right) next one (i-1,j+1) exists, but not next node
  for node in end_nodes:
    # don't check on the boundary
    if node.j > np.size(all_nodes_ids,1)-3 or node.i > np.size(all_nodes_ids,0) - 2:
       continue
    # check if the one after the next one exist
    # nn is node after next node 
    n = get_end_node_by_grid_coords(node.i, node.j+1, end_nodes)[0]
    n_inner = True if all_nodes_ids[node.i,node.j+1] >= 0 else False
    nn = get_end_node_by_grid_coords(node.i,node.j+2,end_nodes)[0]
    nn_inner = True if all_nodes_ids[node.i,node.j+2] >= 0 else False
    ln,ln_idx = get_end_node_by_grid_coords(node.i-1,node.j+1,end_nodes)
    rn,rn_idx = get_end_node_by_grid_coords(node.i+1,node.j+1,end_nodes)
    
    if nn and not (n or n_inner) and (ln or rn):
      # both can be none but if both or not none, something went wrong
      assert(bool(ln) is not bool(rn) or (ln == None and rn == None))
      
      # check which one is in the valley
      next = ln if ln is not None else rn
      next_idx = ln_idx if ln is not None else rn_idx
      
      add_triangle(node.id, next.id, nn.id, cells)
        
      delete_ids.append(next.id)
      node.next = nn
      nn.next = node
      
      continue
    
    # check valleys that contains 2 end nodes
    # node after next next node (node after nn)
    nnn                   = get_end_node_by_grid_coords(node.i, node.j+3, end_nodes)[0]
    rn,rnn_idx   = get_end_node_by_grid_coords(node.i+1,node.j+1,end_nodes)
    rnn,rnn_idx = get_end_node_by_grid_coords(node.i+1,node.j+2,end_nodes)
    ln,ln_idx       = get_end_node_by_grid_coords(node.i-1,node.j+1,end_nodes)
    lnn, lnn_idx  = get_end_node_by_grid_coords(node.i-1,node.j+2,end_nodes)
    
    # next node the node after next node should not be of type end node
    if nnn and not (n or n_inner)  and not ( nn or nn_inner) and (ln and lnn or rn and rnn):
      # check which one is first in the valley
      next = ln if ln is not None else rn
      next_idx = ln_idx if ln is not None else rn_idx
      # check which one is first in the valley
      next_2 = lnn if lnn else rnn
      next_2_idx = lnn_idx if lnn else rnn_idx
      
      add_triangle(node.id, next.id,next_2.id, cells)
      add_triangle(node.id, next_2.id,nnn.id, cells)
        
      node.next = nnn
      nnn.next = node
      
      delete_ids.append(next.id)
      delete_ids.append(next_2.id)
      
  # remove all end nodes that were in valleys
  end_nodes[:] = [v for v in end_nodes if v.id not in delete_ids]
  
# list is list of lists contains surface lines and all respective points
# base used for setting right ids
# nodes is a 3d array with (x,y,z) coords for one profile
# returns array of size 2 x n with ids of end nodes in first row, second row reserved for fix_gaps()
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
          
          if right_id >= 0 and right_next_id >= 0:
            add_triangle(this_id,right_id,right_next_id,cells)
      
      if right_id >= 0 and next_id >= 0 and right_next_id >=0:
        add_triangle(right_id,right_next_id,next_id,cells)
        
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

def create_profiles_array(args,info,log):
  global res
  logger = log
  infoXml = info
  
  res = args.res
  array = np.ones((res,res,res)) * (-1)
  vec = None
  t = np.linspace(0, 1.0, args.res)
  
  profiles = create_profiles(args,infoXml)
  
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
     
    # creates triangles between end nodes of same profile where
    # we have e.g. a valley with 1 or 2 nodes
    postprocess_end_nodes(end_nodes_1,nodes_ids_1,cells)
    postprocess_end_nodes(end_nodes_2,nodes_ids_2,cells)
    postprocess_end_nodes(end_nodes_3,nodes_ids_3,cells)
    
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
          
    # ha is 3dplot object
    id = triangulate_boundary_circles(profiles[0],nodes_ids_1,id,points,cells,vtkData)
    id = triangulate_boundary_circles(profiles[1],nodes_ids_2,id,points,cells,vtkData)
    id = triangulate_boundary_circles(profiles[2],nodes_ids_3,id,points,cells,vtkData)
    print("number of points:",id)
#     plt.show()

    fix_profile_intersection_gaps(end_nodes_1, end_nodes_3, cells)
    fix_profile_intersection_gaps(end_nodes_2, end_nodes_3, cells)
    
    print("number of cells:",cells.GetNumberOfCells())
    global tris
    
    tris = []
    for i in range(cells.GetNumberOfCells()):
      idList = vtk.vtkIdList()
      cells.GetNextCell(idList)
    
    ps, cs = read_vtk("surface.vtp")
    surface_to_volume_mesh(ps,cs)

    polydata = vtk.vtkPolyData()
    polydata.SetPoints(points)
    polydata.GetPointData().SetScalars(vtkData)
    
    if args.export == "surface_points":
      show_write_vtk(polydata, 1000, "surface_points.vtp")
      
    polydata.SetPolys(cells)
    
    if args.show:
      show_vtk(polydata, 1000, [], True)
      
    write_stl(polydata)  
      
#     show_write_vtk(polydata,1000,"surface.vtp")
  else:
    for i in range(0,3):
      if profiles[i] == None:
        continue
      if args.verbose == 'profile_map' or args.export == 'radius_maps':
        create_profile_map(profiles[i], res, args.verbose,args.export == 'radius_maps', ha)
      if args.target == "volume_mesh" or args.target == "volume_vtk":
        write_profile_to_array(array, profiles[i], i)
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
  assert(dir >=0 and dir <=2)
  
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
          if dir == 0:
            array[i,j,k] = dir
          if dir == 1:
            array[j,i,k] = dir
          if dir == 2:
            array[k,j,i] = dir


# find end node object in list with matching id
def get_end_node_by_id(id,end_nodes):
  for node in end_nodes:
    if node.id == id:
      return node
  
  return None

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
  center = (0.5,0.5,0.5)
  for vertice in vertices:
    for n in vertice.neighbors:
      # if we are one of the vertices skip
      if n[1] == vertices[0].id or n[1] == vertices[1].id or n[1] == vertices[2].id:
        continue
    
      b = center - vertices[0].coords
      right1 = vertices[1].coords-vertices[0].coords
      right2 = vertices[2].coords-vertices[0].coords
      right3 = center-n[0]
      A = np.transpose(np.asarray((right1,right2,right3), dtype=float))
      
      sol = linsolve_3x3(A,b)
      
      # sol[0] -> k, sol[0] -> l, sol[0] -> s
      if sol[0] >= 0 and sol[1] >= 0 and sol[0]+sol[1] <= 1:
        if logger:
          logger.write("line through node " + str(n[1]) + " intersects  with triangle (" + str(vertices[0].id) + "," + str(vertices[1].id) + "," + str(vertices[2].id) + ")\n")
        return True    
  
  return False 
# applying Cramer's rule
# assume A is a numpy array, b a list
def linsolve_3x3(A,b):
  assert(isinstance(A,np.ndarray))
  res = [None] * 3
  
  det = np.linalg.det(A)
  
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
    
# triangle quality is measured by the aspect ratio and whether triangle contains an end node
# we penalize the quality if triangle contains a neighbor end node             
def calc_triangle_quality(node1,end_nodes_1,node2,end_nodes_2,node3,end_nodes_3):
  neighbors_1 = give_neighbor_end_nodes(node1, end_nodes_1)
  neighbors_2 = give_neighbor_end_nodes(node2, end_nodes_2)
  neighbors_3 = give_neighbor_end_nodes(node3, end_nodes_3)

  ratio = 1
  
  # in case we test a triangle where one of the 3 points appears two times
  # -> triangle not valid 
  if node1.id == node2.id or node1.id == node3.id or node2.id == node3.id:
    return 1e6
        
  if triangle_contains_any_node([node1,node2,node3]):
    ratio = 1000
    
  ratio *= calc_triangle_ratio(node1.coords, node2.coords, node3.coords)
  
  return ratio

# @param: vertices defining triangle (coordinates)
# for evaluating quality of triangle, we take the ration of the exradius to twice the inradius
def calc_triangle_ratio(v1,v2,v3):
  d1 = calc_distance(v1, v2) # a
  d2 = calc_distance(v1, v3) # b
  d3 = calc_distance(v2, v3) # c
  
  inradius = 0.5* np.sqrt((d2+d3-d1)*(d3+d1-d2)*(d1+d2-d3) / (d1+d2+d3))
  exradius = d1*d2*d3/np.sqrt((d1+d2+d3)*(d2+d3-d1)*(d3+d1-d2)*(d1+d2-d3))
  
  if logger:
    logger.write("inradius " + str(inradius) + "\n")
    logger.write("exradius " + str(exradius) + "\n")
  
  aspect_ratio = d1*d2*d3 / ( (d2+d3-d1) * (d1-d2+d3) * (d1+d2-d3))
  
#   assert(aspect_ratio >= 1)
  if aspect_ratio <= 1:
    return 1e6
  
  return aspect_ratio
    
def find_closest_point(ref_node,next_node,end_nodes):
  # iterate over all end nodes in list and compare distances to 'ref' node
  min_distance = 1e6
  min_node = None
  ref_coords = ref_node.coords + 0.5 * (next_node.coords - ref_node.coords)
   
  for node in end_nodes:
#     assert(node.dir != ref_node.dir)
    dist = calc_distance(ref_coords, node.coords)
    if dist < min_distance:
      min_distance = dist
      min_node = node
  
  if logger:
    logger.write("closest node to id=" + str(ref_node.id) + ": id=" + str(min_node.id) + "\n")
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
# dir and other_dir helps determining walking direction in case node is a starting node
# res: in case resolution of grid on which end_nodes live on has different resolution then basecell grid
# --> this is the case whenn triangulating boundary circles, there we have to create points inside the circle with different resolution
def give_next_end_node(node,end_nodes,dir=None,other_dir=None,initial_edge=None):
  # in case we have reached one vertice of the very first triangle found by algorithm
#   if initial_edge and (node.id == initial_edge[0].id or node.id == initial_edge[1].id):
#     return [node]
  diff_i = None
  diff_j = None
  
  if dir is not None and other_dir is not None:
    # we define front in relative coordinates by means of grid coordinates i and j
    if dir == 0 and other_dir == 2 or dir == 2 and other_dir == 0:
      if node.dir == 0:
        if node.j < res/2:
          diff_i = 1
          diff_j = 0
        else:
          diff_i = 0
          diff_j = -1
          
      elif node.dir == 2:
        if node.coords[0] < 0.5:
          diff_i = 0
          diff_j = 1
        else:
          diff_i = 0
          diff_j = -1
    elif dir == 1 and other_dir == 2 or dir == 2 and other_dir == 1:
      if node.dir == 1:
        if node.j < res/2:
          diff_i = 1
          diff_j = 0
        else:  
          diff_i = -1
          diff_j = 0
      elif node.dir == 2:
        if node.coords[0] < 0.5:
          diff_i = 1
          diff_j = 0
        else:
          diff_i = 1
          diff_j = 0
    else: # not handled yet
      #assert(False)
      diff_i = 1
      diff_j = 0
      
  # find node in connections that is on the same profile
  # get node objects
  previous = None
  for n in node.connections:
    v = get_end_node_by_id(n, end_nodes)
    if get_end_node_by_id(n, end_nodes) is not None:  
      previous = v
      # if node is starting node, it cannot be node previous to this node
      if initial_edge and (n != initial_edge[0].id or n != initial_edge[1].id):
        break
    
  if previous:
    diff_i = np.sign(node.i - previous.i)
    diff_j = np.sign(node.j - previous.j)
    # if node.i is res-1 and previous.i is, diff_i is 1
    # but this is wrong, as previous.i = 0 only results from modulo the resolution
    # in this case, diff_i should be -1!
    if previous.i == 0 and node.i == res_surf_lines-1:
      diff_i = -1
    if node.i == 0 and previous.i == res_surf_lines-1:
      diff_i = 1
    if logger:
      logger.write("node " + str(node.id) + " i,j: " + str(node.i) + "," + str(node.j) + " previous: " + str(previous.id) + " dir_i: " + str(diff_i) + " dir_j: " + str(diff_j) + "\n")
  else:
    if logger:
      logger.write("first node " + str(node.id) + " i,j: " + str(node.i) + "," +  str(node.j) + " dir_i: " + str(diff_i) + " dir_j: " + str(diff_j) + "\n")
  
  candidates = []
  i = node.i
  j = node.j
  
  left_line = node.i-1 if node.i > 0 else res_surf_lines-1
  right_line = node.i+1 if node.i < res_surf_lines-1 else 0
  
  #candidates.append(())
  if diff_j == 0:
    if diff_i == 1:
      candidates.append((i,j-1))
      candidates.append((right_line,j-1))#candidates.append((i+1,j-1))
      candidates.append((right_line,j))#candidates.append((i+1,j))
      candidates.append((right_line,j+1))#candidates.append((i+1,j+1))
      candidates.append((i,j+1))
    else:
      candidates.append((i,j+1))
      candidates.append((left_line,j+1))#candidates.append((i-1,j+1))
      candidates.append((left_line,j))#candidates.append((i-1,j))
      candidates.append((left_line,j-1))#candidates.append((i-1,j-1))
      candidates.append((i,j-1))
  else:
    if diff_i == 0: 
      if diff_j == 1:
        candidates.append((right_line,j))
        candidates.append((right_line,j+1))
        candidates.append((i,j+1))
        candidates.append((left_line,j+1))
        candidates.append((left_line,j))
      else: # diff_j = -1
        candidates.append((left_line,j))
        candidates.append((left_line,j-1))
        candidates.append((i,j-1))
        candidates.append((right_line,j-1))
        candidates.append((right_line,j))
    else: # all diagonal cases
      if diff_i == 1 and diff_j == 1:
        candidates.append((right_line,j-1))
        candidates.append((right_line,j))
        candidates.append((right_line,j+1))
        candidates.append((i,j+1))
        candidates.append((left_line,j+1))
        candidates.append((left_line,j))
      elif diff_i == 1 and diff_j == -1:
        candidates.append((left_line,j))
        candidates.append((left_line,j-1))
        candidates.append((i,j-1))
        candidates.append((right_line,j-1))
        candidates.append((right_line,j))
        candidates.append((right_line,j+1))
      elif diff_i == -1 and diff_j == 1:
        candidates.append((right_line,j))
        candidates.append((right_line,j+1))
        candidates.append((i,j+1))
        candidates.append((left_line,j+1))
        candidates.append((left_line,j))
        candidates.append((left_line,j-1))
      else: # diff_i == -1 and diff_j == -1
        candidates.append((left_line,j+1))
        candidates.append((left_line,j))
        candidates.append((left_line,j-1))
        candidates.append((i,j-1))
        candidates.append((right_line,j-1))
        candidates.append((right_line,j))   

  result = []  
    
  for test in candidates:
    nexts = [v for v in end_nodes if v.i == test[0] and v.j == test[1]]
    if nexts:
      for v in nexts:
        if logger:
          logger.write("v: " + str(v) + "\n")
        if v.id not in node.connections:
          result.append(v)
        if initial_edge and (v.id == initial_edge[0].id or v.id == initial_edge[1].id):
          result.append(node)
          
  if initial_edge and (node.id == initial_edge[0].id or node.id == initial_edge[1].id):
    result.append(node)
  
  if logger:
    logger.write("possible next nodes of " + str(node.id) +" :")  
    for r in result:
      logger.write(str(r.id) + " ")
    logger.write("\n")
  
  return result

# similar to give_next_end_node but much simpler
def give_next_end_node_on_circle(node,end_nodes,initial_edge,max_i,max_j):
  if len(end_nodes) == 1:
    return [node]
#   assert(initial_edge is not None)
  assert(list_contains_end_nodes(end_nodes, node))
  result = []
  next = [v for v in end_nodes if v.i == (node.i+1)%max_i and v.j == node.j]
  
  for v in next:
    if v.id not in node.connections:
      result.append(v)
    if initial_edge and (v.id == initial_edge[0].id or v.id == initial_edge[1].id):
      result.append(node)      
  if initial_edge and (node.id == initial_edge[0].id or node.id == initial_edge[1].id):
    result.append(node)
      
  assert(len(result)>0)
  
  if logger:
    if next:
      logger.write("this: " + str(node.id) + " i:" + str(node.i) + " j:" + str(node.j) + " nexts on circle: ")
      for n in result:
        logger.write(str(n.id) + " ")
      logger.write("\n")
        
  return result
  

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

# remove connections to vertices (array with 2 nodes ids)
def update_connections(triangles,vertices):
#   res = [None] * 3
  vertices[0].connections = set([v for v in vertices[0].connections if not (v==vertices[1].id or v==vertices[2].id)])
  vertices[1].connections = set([v for v in vertices[1].connections if not (v==vertices[0].id or v==vertices[2].id)])
  vertices[2].connections = set([v for v in vertices[2].connections if not (v==vertices[0].id or v==vertices[1].id)])
  
#   print "removed connection from ", vertices[0].id, " to ", vertices[1].id, " and ", vertices[2].id
#   print "removed connection from ", vertices[1].id, " to ", vertices[0].id, " and ", vertices[2].id
#   print "removed connection from ", vertices[2].id, " to ", vertices[0].id, " and ", vertices[1].id
#   
#   triangles[-1].connections = res
  # remove connections to previous triangle vertices
  triangles[-1].edge[0].connections = set([v for v in triangles[-1].edge[0].connections if not (v==vertices[0].id or v==vertices[1].id or v==vertices[2].id)])
  triangles[-1].edge[1].connections = set([v for v in triangles[-1].edge[1].connections if not (v==vertices[0].id or v==vertices[1].id or v==vertices[2].id)])
  triangles[-1].vertices[0].connections = set([v for v in triangles[-1].vertices[0].connections if not (v==vertices[0].id or v==vertices[1].id or v==vertices[2].id)])
  triangles[-1].vertices[1].connections = set([v for v in triangles[-1].vertices[1].connections if not (v==vertices[0].id or v==vertices[1].id or v==vertices[2].id)])
  triangles[-1].vertices[2].connections = set([v for v in triangles[-1].vertices[2].connections if not (v==vertices[0].id or v==vertices[1].id or v==vertices[2].id)])
  
#   print "removed connections to ", vertices[0].id, vertices[1].id, vertices[2].id, " from triangle with vertices ", triangles[-1]

# from a list of candidates, choose the third triangle vertice such that the resulting triangle has the smallest aspect_ratio
def give_best_next_neighbor(candidates, vert1, vert2):
  if len(candidates) == 1 and (candidates[0].id == vert1.id or candidates[0].id == vert2.id):
    return candidates[0]
  best_ratio = 1e9
  best_neighbor = None
  
  for node in candidates:
    if node.id == vert1.id or node.id == vert2.id:
      continue
    ratio = calc_triangle_ratio(vert1.coords, vert2.coords, node.coords)
    if logger:
      logger.write("node: " + str(node.id) + " ratio: " + str(ratio) + "\n")
    if ratio < best_ratio:
      best_ratio = ratio
      best_neighbor = node
      
  assert(best_neighbor is not None)
  return best_neighbor

def check_next_triangle_with_cand(cand,triangles,end_nodes_1,end_nodes_2,initial_edge,alternatives):
   # take active edge from last triangle
  active_edge = triangles[-1].edge
  a = active_edge[0]
  b = active_edge[1]
  
  # end node list that contains node a
  a_nodes = end_nodes_1 if list_contains_end_nodes(end_nodes_1, a) else end_nodes_2
  b_nodes = end_nodes_1 if list_contains_end_nodes(end_nodes_1, b) else end_nodes_2
  
  cand_nodes = end_nodes_1 if end_nodes_1[0].dir == cand.dir else end_nodes_2
  ratio = calc_triangle_quality(a, a_nodes, b, b_nodes, cand, cand_nodes)
  
  if logger:
    logger.write("\na: " + str(a.id) + " connections: " + str(a.connections) + "\n")
    logger.write("b: " + str(b.id) + " connections: " + str(b.connections) + "\n")
    logger.write("given next node: " + str(cand.id) + " ratio: " + str(ratio) + "\n")
  
  if ratio < 20:
    vertices = [a,b,cand]
    #set connections
    a.connections.add(cand.id)
    b.connections.add(cand.id)
    cand.connections.add(a.id)
    cand.connections.add(b.id)
    
    if logger:
      logger.write("added connection from " + str(a.id) + " to " + str(cand.id) + "\n")
      logger.write("added connection from " + str(b.id) + " to " + str(cand.id) + "\n")
      logger.write("added connection from " + str(cand.id) + " to " + str(a.id) + "\n")
      logger.write("added connection from " + str(cand.id) + " to " + str(b.id) + "\n")
    
    triangles[-1].next_node = cand
    triangles[-1].vertices = vertices
    triangles[-1].edge[0] = a if a.dir != cand.dir else b
    triangles[-1].edge[1] = cand 
    
    if logger:
      logger.write("created triangle with edge: (" + str(triangles[-1].edge[0].id) + "," + str(triangles[-1].edge[1].id) + ") next: " + str(cand.id) + "\n")
      logger.write("alternatives: ")
      if alternatives: 
        for alt in alternatives:
          logger.write(str(alt.id) + " ")
      logger.write("\n")     
    check_next_triangle(triangles,end_nodes_1,end_nodes_2,initial_edge=initial_edge)
  else:
    if logger:
      logger.write("bad triangle " + str(a.id) + "," + str(b.id) + "," + str(cand.id) + "\n") 
    handle_bad_triangle(triangles,end_nodes_1,end_nodes_2,initial_edge,alternatives) 

# if we don't have alternatives (other candidates for next end node):
# take care of triangle if detected a bad triangle (aspect ratio too big, contains projection of neighbors, ...)
# corresponds to going one step backwards in algorithm
# if there are alternatives, remove connections to del_node and check new triangle with one alternative    
#def handle_bad_triangle(triangles,end_nodes_1,end_nodes_2,initial_edge,alternatives=None,del_node=None):
def handle_bad_triangle(triangles,end_nodes_1,end_nodes_2,initial_edge,alternatives=None):
  next_cand = None
  # if we don't have alternatives for neighbor end nodes
  if alternatives is None or len(alternatives) == 0:
    pop_triangles(triangles)
    next_cand = triangles[-1].other_candidates[-1]
    triangles[-1].other_candidates.pop(-1)
    remove_connection_to_node(triangles[-1].next_node,triangles[-1])
    triangles[-1].edge = triangles[-2].edge
    if logger:
      logger.write("set active edge to (" + str(triangles[-1].edge[0].id) + "," + str(triangles[-1].edge[1].id) + ")\n")
      logger.write("check triangle " + str(triangles[-1].edge[0].id) + "," + str(triangles[-1].edge[1].id) + " with " + str(next_cand.id) + "\n")
  else:
    # get best neighbor from list with alternatives
    next_cand = give_best_next_neighbor(alternatives, triangles[-1].edge[0], triangles[-1].edge[1])
    if next_cand:
      # remove chosen next node from alternatives
      alternatives = [v for v in alternatives if v.id != next_cand.id]
    
#     remove_connection_to_node(del_node,)
  # we need to revert the active edge as we're going on step back but want to keep info 
  # on already checked and not checked neighbor candidates
  # remove connection to previous next node
  
  check_next_triangle_with_cand(next_cand,triangles, end_nodes_1, end_nodes_2, initial_edge,alternatives)
    
def check_next_triangle(triangles,end_nodes_1,end_nodes_2,initial_edge=None):
  # take active edge from last triangle
  active_edge = triangles[-1].edge
  a = active_edge[0]
  b = active_edge[1]
  
  max_i_1 = 0
  max_j_1 = 0
  max_i_2 = 0
  max_j_2 = 0
  # determines if we are meshing a circular area
  on_circle = False
  if end_nodes_1[0].dir == end_nodes_2[0].dir:
    on_circle = True 
    # calc max value for i and j (surface grid coordinates)
    for n in end_nodes_1:
      if n.i > max_i_1:
        max_i_1 = n.i
      if n.j > max_j_1:  
        max_j_1 = n.j
    for n in end_nodes_2:
      if n.i > max_i_2:
        max_i_2 = n.i
      if n.j > max_j_2:  
        max_j_2 = n.j 
  
  max_i_1 += 1
  max_j_1 += 1
  max_i_2 += 1
  max_j_2 += 1      
           
#   if a.id == 7418 or b.id == 7418:
#     return False
  if logger:
    logger.write("\na: " + str(a.id) + " connections: " + str(a.connections) + "\n")
    logger.write("b: " + str(b.id) + " connections: " + str(b.connections) + "\n")
  
  a_nodes = end_nodes_1 if list_contains_end_nodes(end_nodes_1, a) else end_nodes_2
  b_nodes = end_nodes_1 if list_contains_end_nodes(end_nodes_1, b) else end_nodes_2 
  
  ratio_1 = 1e6
  ratio_2 = 1e6
  next_cand_1 = None
  next_cand_2 = None
  cand_2_nodes = None
  candidates_2 = None
  
  if list_contains_end_nodes(end_nodes_1, a):
    cand_1_nodes = end_nodes_1
    max_i = max_i_1
    max_j = max_j_1  
  else:
    cand_1_nodes = end_nodes_2
    max_i = max_i_2
    max_j = max_j_2
  if on_circle:
    candidates_1 = give_next_end_node_on_circle(a, cand_1_nodes,initial_edge,max_i,max_j)
  else:
    candidates_1 = give_next_end_node(a, cand_1_nodes, end_nodes_1[0].dir, end_nodes_2[0].dir, initial_edge=initial_edge)

  # in this case we stepped over a valley
  if len(candidates_1) == 0:
    if a.next is None:
      if logger:
        logger.write("could not find next neighbor of " + str(a.id) + "\n")
      handle_bad_triangle(triangles,end_nodes_1,end_nodes_2,initial_edge)
    # check if next neighbor is already saved in a node
    next_cand_1 = a.next
  else:
    # test all possible triangles with neighbors and take the one with the smallest aspect ratio
    next_cand_1 = give_best_next_neighbor(candidates_1, a, b)
    ratio_1 = calc_triangle_quality(a,a_nodes,b,b_nodes,next_cand_1,cand_1_nodes)

  if list_contains_end_nodes(end_nodes_1, b):
    cand_2_nodes = end_nodes_1
    max_i = max_i_1
    max_j = max_j_1  
  else:
    cand_2_nodes = end_nodes_2
    max_i = max_i_2
    max_j = max_j_2

  if on_circle:
    candidates_2 = give_next_end_node_on_circle(b, cand_2_nodes,initial_edge,max_i,max_j)
  else:  
    candidates_2 = give_next_end_node(b, cand_2_nodes,end_nodes_1[0].dir, end_nodes_2[0].dir,initial_edge=initial_edge)
  
  # in this case we stepped over a valley
  if len(candidates_2) == 0:
    if b.next is None:
      if logger:
        logger.write("could not find next neighbor of " + str(b.id) + "\n")
      handle_bad_triangle(triangles,end_nodes_1,end_nodes_2,initial_edge)
    # check if next neighbor is already saved in a node
    next_cand_2 = b.next
  else:
    next_cand_2 = give_best_next_neighbor(candidates_2, a, b)
    ratio_2 = calc_triangle_quality(a,a_nodes,b,b_nodes,next_cand_2,cand_2_nodes) if next_cand_2 is not None else 1e6
  
  next = next_cand_1 if ratio_1 < ratio_2 else next_cand_2
  alt = next_cand_2 if ratio_1 < ratio_2 else next_cand_1
  
  if next_cand_1 is None:
    print("next_cand_1 is None: a=" + str(a.id) + " b=" + str(b.id)  + "\n")
  if next is None:
    print("next is None: a=" + str(a.id) + " b=" + str(b.id) + "\n")
    print("next_cand_1: " + str(next_cand_1) + "\n")
    print("next_cand_2: " + str(next_cand_2) + "\n")
  
  assert(next_cand_1 is not None)
  assert(next is not None)  
  
  if logger:
    logger.write(str(next_cand_1) + " ratio=" + str(ratio_1) + (" <- " if next_cand_1.id == next.id else " ") + "\n")
    if next_cand_2:
      logger.write(str(next_cand_2) + " ratio=" + str(ratio_2) + (" <- " if next_cand_2.id == next.id else " ") + "\n")
    
  alternatives = candidates_1[:]+candidates_2[:]
  # remove chosen next node from alternatives
  alternatives = [v for v in alternatives if v.id != next.id]
    
  # if at least one candidate has good quality
  if min(ratio_1,ratio_2) < 20:
    #set connections
    a.connections.add(next.id)
    b.connections.add(next.id)
    next.connections.add(a.id)
    next.connections.add(b.id)
    
    #triangles.append(Marching_Triangle([a,b,next],a if a.dir != next.dir else b,next,alternatives))
    # decide whether a or b is member based on triangle ratios
    # and not profile directions as a.dir might be equal b.dir 
    edge_0 = a if not list_contains_end_nodes(a_nodes, next) else b
    
    triangles.append(Marching_Triangle([a,b,next],edge_0,next,alternatives))
    
    if logger: 
      logger.write("added connection from " + str(a.id) + " to " + str(next.id) + "\n")
      logger.write("added connection from " + str(b.id) + " to " + str(next.id) + "\n")
      logger.write("added connection from " + str(next.id) + " to " + str(a.id) + "\n")
      logger.write("added connection from " + str(next.id) + " to " + str(b.id) + "\n")
      logger.write("created triangle with edge: (" + str(triangles[-1].edge[0].id) + "," + str(triangles[-1].edge[1].id) + ") next: " + str(next.id) + "\n")
      logger.write("alternatives: ")
      if alternatives is not None:
        for alt in alternatives:
          logger.write(str(alt.id) + " ")
      logger.write("\n")
          
  else: # if both candidates are bad, we have go back to previous triangle an try alternative candidate
    if logger:
      logger.write("bad triangle " + str(a.id) + "," + str(b.id) + "," + str(next.id) + "\n") 
      logger.write("alternatives: ")
      if alternatives is not None:
        for alt in alternatives:
          logger.write(str(alt.id) + " ")
      logger.write("\n")
           
    triangles.append(Marching_Triangle([a,b,b],a,b,alternatives))     
    handle_bad_triangle(triangles,end_nodes_1,end_nodes_2,initial_edge,alternatives)
      
  return True

# for a given Marching Triangle, remove all connections from vertices to node
def remove_connection_to_node(node,triangle):
  vertices = triangle.vertices
  for vert in vertices:
    if logger:
      logger.write("remove connection to " + str(node.id) + " from " + str(vert.id) + "\n")
    vert.connections = set([v for v in vert.connections if v != node.id])

# pop triangles that have no alternative candidates    
def pop_triangles(triangles):
  while len(triangles[-1].other_candidates) == 0:
    vertices = triangles[-1].vertices
    
    if logger:
      logger.write("remove triangle: " + str(vertices[0].id) + "," + str(vertices[1].id) + "," + str(vertices[2].id) + "\n")
    
    triangles.pop()
    
    # update connections 
    update_connections(triangles, vertices)
  
  vertices = triangles[-1].vertices
  if logger:
    logger.write("previous triangle: " + str(vertices[0].id) + "," + str(vertices[1].id) + "," + str(vertices[2].id) + "\n")
    logger.write("with other_candidates: ")
    for cand in triangles[-1].other_candidates:
      logger.write(str(cand.id) + " ")
    
# given 3 lists with end nodes of all 3 profiles
# using a marching algorithm, we connect theses nodes step by step to triangles
# start with one given triangle defined by end nodes of two profiles (assume each profile has a different color'
# we call the edge that connects two differently colored nodes an 'active' edge
# when a triangle was defined, continue with active edge and search for possible candidates (for new triangle)
# in the neighborhood of the two edge vertices              
def fix_profile_intersection_gaps(this_end_nodes,other_end_nodes,cells):
  # start from x profile and 0 degree and take first end node you can find
  nodes = [v for v in this_end_nodes if v.i == 0]
  assert(nodes)
  for n in nodes:
    start_node = n
    
    next = give_next_end_node(start_node,this_end_nodes,this_end_nodes[0].dir,other_end_nodes[0].dir)[0]
  
    other = find_closest_point(start_node, next, other_end_nodes)[0]
    
    start_triangulation(start_node,next,other,this_end_nodes,other_end_nodes,cells)  
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
  
  if logger:
    logger.write("arc length:" + str(arc_length) + "\n")
    logger.write("radius: " + str(radius) + "\n")
    logger.write("step_angle: " + str(step_angle) + "\n")
  
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
  global res, res_surf_lines, logger
  radius = profile.radius_left
  arc_length = radius * radians(360.0/np.size(nodes_ids,0))
  previous_points_left, previous_points_right,id, blah = generate_end_nodes_in_circle(profile,arc_length,radius,points,vtkData,id,False)
  set_correct_point_ids(previous_points_left, previous_points_right,nodes_ids)
  
  points_left = []
  points_right = []
  triangles = []
  step = 4.0/res
  tmp_res = 0
  
  # create points on circles lying in the same plane as original circle
  while len(previous_points_left) > 1 and len(previous_points_right) > 1:
    radius -= step
    points_left, points_right, id, tmp_res = generate_end_nodes_in_circle(profile,arc_length,radius,points,vtkData,id)
    
    # create first triangle
    for lists in [(previous_points_left,points_left),(previous_points_right,points_right)]:
      start_node = lists[0][0]
      next = give_next_end_node_on_circle(start_node, lists[0], None, res_surf_lines, res)[0]
      other = find_closest_point(start_node, next, lists[1])[0]
#       if logger:
#         logger.write("starting with " + str(start_node.id) + "," + str(next.id) + "," + str(other.id) + "\n")
#         logger.write("initial edge: " + str(start_node.id) + "," + str(other.id) + "\n")
#         logger.write("start_node" + str(start_node.id) + "\n")
#         logger.write("next" + str(next.id) + "\n")
#         logger.write("other" + str(other.id) + "\n")
      start_triangulation(start_node, next, other,lists[0],lists[1],cells)
             
    previous_points_left = points_left
    previous_points_right = points_right
      
    clear_all_connections(previous_points_left)
    clear_all_connections(previous_points_right)
    
  return id
# for given end node list, remove all connections of end nodes
def clear_all_connections(list):
  for n in list:
    n.connections = set()    
    
def set_correct_point_ids(nodes_left, nodes_right, nodes_ids):
  for node in nodes_left:
    assert(nodes_ids[node.i,node.j] > -1)
    node.id = nodes_ids[node.i,node.j]

  for node in nodes_right:
    assert(nodes_ids[node.i,node.j] > -1)
    node.id = nodes_ids[node.i,node.j]  
    
# checks if node is in list    
def list_contains_end_nodes(list,node):
  for n in list:
    if n.id == node.id:
      return True
  return False

# starts triangulation algorithm for given initial triangle
# all three nodes must be of type End_Node
# next_node is on same 'circle/direction' as start_node
# other_node is on other 'circle/direction' as start_node
def start_triangulation(start,next,other,this_end_nodes,other_end_nodes,cells):
  triangles = [] # saves all triangles found by marching
  global tris
  
  if logger: 
    logger.write("\nstarting with " + str(start.id) + " dir = " + str(this_end_nodes[0].dir) + " other_dir= " + str(other_end_nodes[0].dir) + "\n")
    logger.write("other: " + str(other) + "\n")

  next.connections.add(other.id)
  next.connections.add(start.id)
  start.connections.add(other.id)
  start.connections.add(next.id)
  other.connections.add(start.id)
  other.connections.add(next.id)
  
  triangles.append(Marching_Triangle([start,next,other],next,other,[start]))
  initial_edge = [start,other]
  if logger:
    logger.write("created triangle with edge: (" + str(next.id) + "," + str(other.id) + ") next: " + str(next.id) + "\n")
    logger.write("initial edge: " + str(start.id) + "," + str(other.id) + "\n")
    
  run = True
  end = False
  while len(triangles) > 0 and run and not end:
    run = check_next_triangle(triangles, this_end_nodes, other_end_nodes, initial_edge)
    # check if we have just created a triangle where new active edge is 
    # identical to active edge of very first triangle --> we are finished!
    edge = triangles[-1].edge
    if len(triangles) > 2 and (edge[0].id == initial_edge[0].id and edge[1].id == initial_edge[1].id) or (edge[0].id == initial_edge[1].id and edge[1].id == initial_edge[0].id):
      if logger:
        logger.write("filled \n")
      end = True 
      
  for triangle in triangles:
    verts = triangle.vertices
    add_triangle(verts[0].id, verts[1].id, verts[2].id, cells)
    
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
#   print("num points:", len(pp))
#   print("num cells:", len(cc))
#   print(cc)

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
