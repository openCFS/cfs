#coding=utf-8
from enum import Enum
import math
import sys

from matplotlib import pyplot as plt
import matplotlib
from numpy import outer
from scipy import interpolate
from sympy import Symbol, symbols
import sympy.solvers

from mpl_toolkits.mplot3d import Axes3D
import numpy as np

#matplotlib.use('tkagg')
  
try:
  from skimage import measure
except:
  print("Warning: Failed to load skimag - need it for Marching cubes.")
  
try:
  import vtk
  from vtk.util.numpy_support import vtk_to_numpy
  import matviz_vtk
except:
  print("WARNING: failed to load vtk!")  

res = 1
res_surf_lines = 1
# file object for info xml
infoXml = None
interpolation = None
symmetric = False # basecell symmetric(x1=x2=y1=y2=z1=z2)?

class Cubic_spline():
  # assume we have u_0=u_1=u_2=u_3=0 and u_4=u_5=u_6=u_7=13
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
    sol = sympy.solvers.solve(dC[0][1]-dC[0][0],u) # dx=dy
    t = -100
    if sol[0] > 0 and sol[0] <= 1:
      t = sol[0]
    elif  sol[1] > 0 and sol[1] <= 1:
      t = sol[1]
    else:
      print("No t found where dx/dy = 1 ",sol)
      return None
    # conversion from t as sympy.Float to regular Python float necessary  
    self.t_1 = float(t)
    return self.t_1
  
  # left indicates if we are plotting spline for left branch of basecell bar
  # if not left, then we're at right branch -> compute spline as left branch
  # and mirror the result
  def plot(self,left=True):
    matplotlib.rcParams.update({'font.size': 20})
    matplotlib.rcParams['lines.linewidth'] = 5
    if left:
      t = np.linspace(0, 1, 100)
      C = self.eval_t(t)
      plt.plot(np.transpose(self.CP[0,:]),np.transpose(self.CP[1,:]),marker='o', markersize=15,color='green')
      t1 = self.calc_param_grad_1()
      Ct = self.eval_t(t1)
      plt.plot(np.transpose(C[0,:]),np.transpose(C[1,:]),color='black')
      plt.plot(Ct[0],Ct[1],marker='o', markersize=15)
    else:
      t = np.linspace(0, 1, 100)
      C = self.eval_t(t)
      plt.plot(np.transpose(1-self.CP[0,:]),np.transpose(self.CP[1,:]),marker='o', markersize=15,color='green')
      t1 = self.calc_param_grad_1()
      Ct = self.eval_t(t1)
      plt.plot(1-np.transpose(C[0,:]),np.transpose(C[1,:]),color='black')
      plt.plot(1-Ct[0],Ct[1],marker='o', markersize=15)
    plt.xlim((0,0.4))
      
    plt.ylim((0.6,1.0))
    plt.show()
    
  def calc_coords_grad_1(self):
    t1 = self.calc_param_grad_1()
    if t1 is None:
      return None
    ret = self.eval_t(t1) # array of array with 1 elem, due to outer product
    return np.array((ret[0][0],ret[1][0]))
  
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
  
def cartesian_to_grid_coord(x,res,eps=1e-6):
  h = 1.0 / res # assume domain is 1m x 1m x 1m
  return int(round((x-h/2.0) / h+eps)) 

def grid_to_cartesian_coords(voxel,res,h=None):
  assert(len(voxel) == 3)
  if h is None:
    h = []
    if type(res) is int:
      for i in range(len(voxel)):
        h.append(1.0 / res)
    else:
      for r in res:
        h.append(1.0 / res) # assume domain is 1m x 1m x 1m
  
  p = np.zeros(len(voxel))  
  for i,v in enumerate(voxel):
    p[i] = v * h[i] + h[i]/2.0
     
  return p 

# @param offset: value added to converted cartesian
def polar_to_cartesian(radius,radians,offset=0.5):
  return radius * np.cos(radians) + offset, radius * np.sin(radians) + offset

# convert cartesian coordinates of p=(px,py) to polar ones
# center is the origin of coordinate system e.g. (0.5,0.5)
def cartesian_to_polar(px,py,center):
  return np.sqrt((px - center[0])**2 + (py - center[1])**2)

# calculates euklidian distance to origin (0,0)
# @param p: tuple with x-,y-component of point
def distance_to_center(p):
  return np.sqrt((p[0]-0.5)**2+(p[1]-0.5)**2)

# calculates angle between (0.5,0.5) and point p in radians
def angle_to_center(p):
  """
  >>> print(round(angle_to_center(np.array([1.0,0.0])),5))
  5.49779
  >>> print(round(angle_to_center(np.array([0.5,1.0])),5))
  1.5708
  >>> print(round(angle_to_center(np.array([0.0,0.0])),5))
  3.92699
  >>> print(round(angle_to_center(np.array([0.5,-1.0])),5))
  4.71239
  """
  x = p[0] - 0.5
  y = p[1] - 0.5
  
  val = math.atan2(y,x) # returns value between pi and -pi
  if y < 0:
    val += 2.0*np.pi
    
  return val
 
# defines a 1D linear function
class Linear_1D():
  x1 = None
  x2 = None
  # @param x1 and x2 define the function  
  def __init__(self,x1,x2):
    self.x1 = x1
    self.x2 = x2
    
  # f(0) = x1/2.0 + 0.5
  # f(0.5) = x2/2.0 + 0.5  
  def eval(self,x):
    return (self.x2-self.x1) * 2.0 * x + self.x1

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
  assert(x >= -1e-6 and x <= 1.0 + 1e-6)
  return 1.0 - 1.0 / (np.exp(2.0*beta*(x-eta)) + 1)

# obejct defines spline in a principle plane, e.g. spline for 0 or 90 degree
class PrincipleSpline():
  
  spline = None
  angle = -1
  coords_cut = None
  left = True
  x1 = None # assume spline lives in 2d plane
  y1 = None 
  
  def __init__(self, x1, y1, bend, angle=0, left_flag=True):
    rx = 0.5 - x1/2.0 # radius for center (0,1)
    ry = 0.5 - y1/2.0 # radius for center (0,1)
    self.x1 = x1
    self.y1 = y1
    self.left = left_flag
    
    if infoXml is not None:
      infoXml.write('    <bspline degree="' + str(np.degrees(angle)) + '" rad1="' + str(x1) + '" rad2="' + str(y1) + '" bend="' + str(bend) + '">\n')
    
    self.angle = angle
    P = np.array([[0,1-rx],[ry*bend,1-rx],[ry,1-rx*bend],[ry,1]])
    
    self.spline = Cubic_spline(P)
    # coordinate where slope is 1
    self.coords_cut = self.calc_coords_grad_1()
    
    if self.coords_cut is None:
      print("ERROR: Cannot create spline x1=",x1," y1=",y1)
      sys.exit()
      
    if not self.left:
      self.coords_cut = [1.0 - self.coords_cut[0],self.coords_cut[1]] 
    
    if infoXml is not None:
      infoXml.write('    </bspline>\n')
    
  # eval function in case x is one element, not a list  
  def eval_elem(self,x):
    if self.left:
      if x <= self.coords_cut[0]:
        val = self.spline.eval_x(x)
      elif x > self.coords_cut[0] and x <= 0.5:
        val = self.coords_cut[1]
      else: #x >= 1- coord_cut[0]
        val = 0.5
        #raise Exception("Spline not defined for x=" + str(x))
    else:
      if x >= self.coords_cut[0]:
        val = self.spline.eval_x(1-x) # mirror left part to get right part
      elif x < self.coords_cut[0] and x >= 0.5:
        val = self.coords_cut[1]
      else:
        val = 0.5    
    return val
  
  # wrapper function
  # @param x: can be one argument value or list or aguments
  def eval(self,x):
    if isinstance(x, (np.float32,np.float64,float,int)):
      ret = self.eval_elem(x)
      return ret
    
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
    self.left = None
  
  # @param force: bisec curve to be enforced
  # @param interpolation: interpolation type between splines (linear or heaviside)
  def __init__(self,x1,y1,z1,bend,beta,eta,interpolation,force=None,left=True,x2=-1):
    self.x1 = x1
    self.y1 = y1
    self.z1 = z1
    self.bicubic = [] # 0-th entry is spline function and 1st entry is cubic polynomial
    self.left = left
  
    assert(bend <= 1 and bend >= 0)
    assert(interpolation == "linear" or interpolation == "heaviside")
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
  
    # search for point p
    right = PrincipleSpline(y1, z1, bend)
    p = right.coords_cut
    p[0] = 1.0 - p[0] # move p to 1st quadrant to calculate right angle
    height = distance_to_center(p) + 0.5
    self.angle = angle_to_center(p)
    
    # left part is same spline as orhtogonal spline up to point b
    # don't need eta, as it is set to 0.5 in this case
    left_spline = PrincipleSpline(x1, weighted_by_angle(self.angle,interpolation,y1,z1,beta), bend)
    b = left_spline.coords_cut 
    self.cut = b 
    assert(b[0] <= 0.5 and b[1] >= 0.5)
     
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
    lin_x2 = min(x1,x2)/2.0+0.5
    if height > lin_x2:
      lin_x2 = height
    lin = Linear_1D(x1/2.0+0.5,lin_x2)
    
    self.bicubic.append(left_spline)
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
      t = None 
      if self.left:
        t = np.linspace(0, 0.5, 100)
      else:
        t = np.linspace(0.5, 1, 100)
        
      assert(t is not None)
      samples = self.eval_bicubic(t)
      if max(samples) <= height+1e-3:
        self.type = "bicubic"
      # in case function composed of b-spline and cubic function has undershoot  
      # in case b-spline has no undershoot (point p is not below bspline(x=0))
      elif height > 0.5 + x1/2.0:  
        self.type = "bspline"
      # in case we have undershooting for biqua AND for spline
      else:
        self.type = "linear"
    
    if infoXml:
      infoXml.write('      <selection type="' + str(self.type) + '" angle="' + str(np.degrees(self.angle)) + '"/>\n')
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
      if self.left:
        assert(max(x) <= 0.5 + 1e-6)
        val = self.spline.eval_x(i)
      else:
        assert(max(x) >= 0.5 and max(x) <= 1.0)
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
    val = None
    if self.left:
      for i in x:
        if i <= left.coords_cut[0]:
          val = left.eval(i)
        else:
          val = cubic(i)
        res.append(val)  
    else: #right
#       assert(max(x) >= 0.5)
      for i in x:
        if i <= 1.0-left.coords_cut[0]: # mirror of cubic 
         val = cubic(1-i)
        else:
         val = left.eval(1-i)
        res.append(val)
    
    return res
  
  def eval_linear(self,x):
    if self.left:
      return self.linear.eval(x)
    else:
      return self.linear.eval(1-x)
  
  def eval_heaviside(self,x):
    # make x iterable in case it's one float element and not a list
    x = np.reshape(x, np.size(x), )
    res = []
    for i in x:
      assert(i >= 0 and i <=1)
      if self.left:
        assert(max(x) <= 0.5)
        res.append(self.heaviside.eval(i))
      else:
        assert(max(x) >= 0.5 and max(x) <= 1.0)
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
    plt.rcParams.update({'font.size': 28})
    plt.figure(figsize=(10,10))
    plt.rc('axes', linewidth=3.5)

    
    x = None 
    if self.left:
      x = np.linspace(0, 0.5, 100)
    else:
      x = np.linspace(0.5, 1, 100)
    
    bicubic = self.eval_bicubic(x)
    spline = self.eval_spline(x)
    linear = self.eval_linear(x)
    heavi = self.eval_heaviside(x)
    
    cut = self.bicubic[0].coords_cut
    
    assert(self.type == 'bicubic' or self.type == 'bspline' or self.type == 'linear' or self.type == 'heaviside')

    if self.left:
      plt.xlim((0,0.5))
    else:
      plt.xlim((0.5,1.0))
    
    plt.ylim((0.5,1.0))
    plt.xlabel("x")
    plt.ylabel("r")
    bc_label = 'bicubic*' if self.type == 'bicubic' else 'bicubic'
    sp_label = 'bspline*' if self.type == 'bspline' else 'bspline'
    lin_label = 'linear*' if self.type == 'linear' else 'linear' 
    hv_label = 'heaviside*' if self.type == 'heaviside' else 'heaviside'
    plt.plot(x,bicubic,label=bc_label,linewidth=5.0)
    plt.plot(x,spline,label=sp_label,linewidth=5.0)
    plt.plot(x,linear,label=lin_label,linewidth=5.0)
    plt.plot(x,heavi,label=hv_label,linewidth=5.0)
    if self.left:
      plt.plot(cut[0],cut[1],marker='o',color='black',markersize=15)
    else:
      plt.plot(1-cut[0],cut[1],marker='o',color='black',markersize=15) 
    plt.legend(loc='upper left', shadow=True,prop={'size':20})
#     plt.show()  
    plt.savefig('test.png', dpi=800)
    
  def angle(self):
    return self.angle
    
# @return vector with Profile or list of vectors
class Profile:
  def __init__(self, args, dir):
    assert (dir == 0 or dir == 1 or dir == 2)
    self.bisec_angle = -1
    self.direction = dir
    self.splines_left = [None] * 4
    self.bisecs_left = [None] * 4
    self.splines_right = [None] * 4
    self.bisecs_right = [None] * 4
    # depending on profile, store the radii of the two boundary circles
    self.radius_left = 0
    self.radius_right = 0

    if infoXml:
      infoXml.write('  <profile dir="' + str(dir) + '">\n')
      
    plot_dir = None if not args.plot_bisec else args.plot_bisec[0]
    assert(plot_dir is None or plot_dir == "x" or plot_dir == "y" or plot_dir == "z")
    plot_bisec = -1 if not args.plot_bisec else int(args.plot_bisec[1])
    assert(plot_bisec <= 7) 
     
    # right hand side coordinate system
    # major axis points away, first minor axis points rightwards, second minor points upward (s. paraview)
    # major x, minors: z,y
    # major y, minors: x,z
    # major z, minors: y,x  
    if dir == 0:
      self.splines_left[0] = PrincipleSpline(args.x1, args.z2, args.bend, 0)
      self.splines_left[1] = PrincipleSpline(args.x1, args.y2, args.bend, np.pi/2.0)
      self.splines_left[2] = PrincipleSpline(args.x1, args.z1, args.bend, np.pi)
      self.splines_left[3] = PrincipleSpline(args.x1, args.y1, args.bend, 1.5*np.pi)
      
#       self.splines_left[0].spline.plot()
      
      self.bisecs_left[0] = BisecSpline(args.x1, args.z2, args.y2, args.bend,args.beta,args.eta,args.interpolation,args.force_bisec,x2=args.x2)
      self.bisecs_left[1] = BisecSpline(args.x1, args.y2, args.z1, args.bend,args.beta,args.eta,args.interpolation,args.force_bisec,x2=args.x2)
      self.bisecs_left[2] = BisecSpline(args.x1, args.z1, args.y1, args.bend,args.beta,args.eta,args.interpolation,args.force_bisec,x2=args.x2)
      self.bisecs_left[3] = BisecSpline(args.x1, args.y1, args.z2, args.bend,args.beta,args.eta,args.interpolation,args.force_bisec,x2=args.x2)
      
      self.splines_right[0] = PrincipleSpline(args.x2, args.z1, args.bend, 0, False)
      self.splines_right[1] = PrincipleSpline(args.x2, args.y2, args.bend, np.pi/2.0, False)
      self.splines_right[2] = PrincipleSpline(args.x2, args.z1, args.bend, np.pi, False)
      self.splines_right[3] = PrincipleSpline(args.x2, args.y1, args.bend, 1.5*np.pi, False)
       
      self.bisecs_right[0] = BisecSpline(args.x2, args.z2, args.y2, args.bend,args.beta,args.eta,args.interpolation,args.force_bisec,left=False,x2=args.x1)
      self.bisecs_right[1] = BisecSpline(args.x2, args.y2, args.z1, args.bend,args.beta,args.eta,args.interpolation,args.force_bisec,left=False,x2=args.x1)
      self.bisecs_right[2] = BisecSpline(args.x2, args.z1, args.y1, args.bend,args.beta,args.eta,args.interpolation,args.force_bisec,left=False,x2=args.x1)
      self.bisecs_right[3] = BisecSpline(args.x2, args.y1, args.z2, args.bend,args.beta,args.eta,args.interpolation,args.force_bisec,left=False,x2=args.x1)
      
      self.radius_left = args.x1 / 2.0
      self.radius_right = args.x2 / 2.0
    elif dir == 1:
      self.splines_left[0] = PrincipleSpline(args.y1, args.x2, args.bend, 0)
      self.splines_left[1] = PrincipleSpline(args.y1, args.z2, args.bend, np.pi/2.0)
      self.splines_left[2] = PrincipleSpline(args.y1, args.x1, args.bend, np.pi)
      self.splines_left[3] = PrincipleSpline(args.y1, args.z1, args.bend, 1.5*np.pi)
      
      self.bisecs_left[0] = BisecSpline(args.y1, args.x2, args.z2, args.bend,args.beta,args.eta,args.interpolation,args.force_bisec,x2=args.y2)
      self.bisecs_left[1] = BisecSpline(args.y1, args.z2, args.x1, args.bend,args.beta,args.eta,args.interpolation,args.force_bisec,x2=args.y2)
      self.bisecs_left[2] = BisecSpline(args.y1, args.x1, args.z1, args.bend,args.beta,args.eta,args.interpolation,args.force_bisec,x2=args.y2)
      self.bisecs_left[3] = BisecSpline(args.y1, args.z1, args.x2, args.bend,args.beta,args.eta,args.interpolation,args.force_bisec,x2=args.y2)
      
      self.splines_right[0] = PrincipleSpline(args.y2, args.x2, args.bend, 0, False)
      self.splines_right[1] = PrincipleSpline(args.y2, args.z2, args.bend, np.pi/2.0, False)
      self.splines_right[2] = PrincipleSpline(args.y2, args.x1, args.bend, np.pi, False)
      self.splines_right[3] = PrincipleSpline(args.y2, args.z1, args.bend, 1.5*np.pi, False)
       
      self.bisecs_right[0] = BisecSpline(args.y2, args.x2, args.z2, args.bend,args.beta,args.eta,args.interpolation,args.force_bisec,left=False,x2=args.y1)
      self.bisecs_right[1] = BisecSpline(args.y2, args.z2, args.x1, args.bend,args.beta,args.eta,args.interpolation,args.force_bisec,left=False,x2=args.y1)
      self.bisecs_right[2] = BisecSpline(args.y2, args.x1, args.z1, args.bend,args.beta,args.eta,args.interpolation,args.force_bisec,left=False,x2=args.y1)
      self.bisecs_right[3] = BisecSpline(args.y2, args.z1, args.x2, args.bend,args.beta,args.eta,args.interpolation,args.force_bisec,left=False,x2=args.y1)
      
      self.radius_left = args.y1 / 2.0
      self.radius_right = args.y2 / 2.0
    else: # dir == 2
      self.splines_left[0] = PrincipleSpline(args.z1, args.y2, args.bend, 0)
      self.splines_left[1] = PrincipleSpline(args.z1, args.x2, args.bend, np.pi/2.0)
      self.splines_left[2] = PrincipleSpline(args.z1, args.y1, args.bend, np.pi)
      self.splines_left[3] = PrincipleSpline(args.z1, args.x1, args.bend, 1.5*np.pi)
      
      self.bisecs_left[0] = BisecSpline(args.z1, args.y2, args.x2, args.bend,args.beta,args.eta,args.interpolation,args.force_bisec,x2=args.z2)
      self.bisecs_left[1] = BisecSpline(args.z1, args.x2, args.y1, args.bend,args.beta,args.eta,args.interpolation,args.force_bisec,x2=args.z2)
      self.bisecs_left[2] = BisecSpline(args.z1, args.y1, args.x1, args.bend,args.beta,args.eta,args.interpolation,args.force_bisec,x2=args.z2)
      self.bisecs_left[3] = BisecSpline(args.z1, args.x1, args.y2, args.bend,args.beta,args.eta,args.interpolation,args.force_bisec,x2=args.z2)
      
      self.splines_right[0] = PrincipleSpline(args.z2, args.y2, args.bend, 0, False)
      self.splines_right[1] = PrincipleSpline(args.z2, args.x2, args.bend, np.pi/2.0, False)
      self.splines_right[2] = PrincipleSpline(args.z2, args.y1, args.bend, np.pi, False)
      self.splines_right[3] = PrincipleSpline(args.z2, args.x1, args.bend, 1.5*np.pi, False)
      
      self.bisecs_right[0] = BisecSpline(args.z2, args.y2, args.x2, args.bend,args.beta,args.eta,args.interpolation,args.force_bisec,left=False,x2=args.z1)
      self.bisecs_right[1] = BisecSpline(args.z2, args.x2, args.y1, args.bend,args.beta,args.eta,args.interpolation,args.force_bisec,left=False,x2=args.z1)
      self.bisecs_right[2] = BisecSpline(args.z2, args.y1, args.x1, args.bend,args.beta,args.eta,args.interpolation,args.force_bisec,left=False,x2=args.z1)
      self.bisecs_right[3] = BisecSpline(args.z2, args.x1, args.y2, args.bend,args.beta,args.eta,args.interpolation,args.force_bisec,left=False,x2=args.z1)
      
      self.radius_left = args.z1 / 2.0
      self.radius_right = args.z2 / 2.0
      
    if plot_dir and plot_dir == dirToString(dir):
      if plot_bisec <= 3: # left side
        self.bisecs_left[plot_bisec].plot_all()
      else: # right side
        self.bisecs_right[plot_bisec-4].plot_all()
      plt.show()
    
    if infoXml:  
      infoXml.write('  </profile>\n\n')
  
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
    
  if args.verbose == "all_splines" or args.verbose == "all_bisecs":
    figs = []
    for i in range(0,3):
      f = plt.figure(i,figsize=(9, 15))
      f.suptitle("dir "+str(i), fontsize=20)
      figs.append(f)

    x_left = np.linspace(0, 0.5, 1000)
    x_right = np.linspace(0.5, 1.0, 1000)
    
    for dir,profile in enumerate(profiles):
      if profile == None:
        continue
      
      count = 411 # need this for add_suplot
      for i in range(0,4):
        sub1 = figs[dir].add_subplot(count)
        sub1.set_ylim((0.5,1.0))
        if args.verbose == "all_splines":
          sub1.set_title(str(np.degrees(profile.splines_left[i].angle)) + "°")
          sub1.plot(x_left,profile.splines_left[i].eval(x_left),linewidth=5.0)
          sub1.plot(x_right,profile.splines_right[i].eval(x_right),linewidth=5.0)
        else: # all_bisecs
          sub1.set_title(str(np.degrees(profile.bisecs_left[i].angle)) + "°")
          sub1.plot(x_left,profile.bisecs_left[i].eval(x_left),linewidth=5.0,color='navy')
          sub1.plot(x_right,profile.bisecs_right[i].eval(x_right),linewidth=5.0,color='navy')
        count += 1
      
    plt.show()
  
  
  return profiles

# @param map: Profile map
# @param numLines: number of lines on surface that we want to plot
# @param dir: direction (x,y or z) as (1,2,3)
def get_surface_lines(profile):
  nodes = np.zeros((res_surf_lines,res, 3))
  interval = np.linspace(0, 1.0, res)
  
  for i,alpha in enumerate(np.arange(0,360,360.0/res_surf_lines)):
    for j,x in enumerate(interval):
      nodes[i,j] = get_surface_point_candidate(profile,np.radians(alpha),x)
  
  return nodes

# @param alpha in radians
def get_surface_point_candidate(profile,alpha,x):
  assert(alpha >= 0 and alpha <= 2.0*np.pi)
 
  radius = calc_radius(profile, x, alpha)
  px,py = polar_to_cartesian(radius, alpha)
 
  point = np.zeros(3)
  major_dir = profile.direction
  minor_dir_1, minor_dir_2 = give_normal_plane_axes(major_dir)
  point[major_dir] = x
  point[minor_dir_1] = px
  point[minor_dir_2] = py
 
  return point

# calc distance between two points
def calc_distance(p1,p2):
  return np.linalg.norm(np.asarray(p1)-np.asarray(p2))

# returns points and cells describing basecell 
def generate_basecell(args,info):
  global res, res_surf_lines, interpolation
  res = args.res
  res_surf_lines = args.res_surf_lines
  interpolation = args.interpolation
  
  global infoXml
  infoXml = info
  
  array = np.full((res,res,res),args.lower,dtype=float)
  if args.multiple_regions:
    array = np.ones((res,res,res),dtype=np.int) * (-1)
  
  profiles = create_profiles(args,infoXml)
  
  hf = None
  polydata = None
  new_surf_points = None
  faces = None
  bp_lists = None
  
  if args.verbose != "off" or args.target == "3dlines":
    hf = plt.figure()
    ha = hf.add_subplot(111, projection='3d')
    ha.set_xlabel('X')
    ha.set_ylabel('Y')
    ha.set_zlabel('Z')
  
  assert(len(profiles) == 3)
  
  id = 0 # assign an id to each surface point
  
  points = None
  cells = None
  for i in range(0,3):
    if profiles[i] == None:
      continue
    if args.verbose == "profile_map" or args.export == "radius_maps" or args.verbose == "polar_plot":
      create_profile_map(profiles[i], res, args.verbose,args.export == "radius_map", ha)
    if args.verbose == "interpolation":
      y1 = []
      y2 = []
      rad = np.linspace(0,pi/2.0,500)
      for r in rad:
        y1.append(calc_radius(profiles[i], 0.495, r))
        y2.append(calc_radius(profiles[i], 0.51, r))
        
      plt.gcf().clear()
      fig = plt.figure(1)
      sub1 = fig.add_subplot(211)
      sub1.set_title("x=0.1")
      sub1.plot(rad,y1,linewidth=5.0)
      sub2 = fig.add_subplot(212)
      sub2.set_title("x=0.9")
      sub2.plot(rad,y2,linewidth=5.0)
      
      plt.show()
      
    if args.target == "volume_mesh" or args.target == "surface_mesh" or args.target == "marching_cubes" or args.target == "image":
      global symmetric
      # if basecell is symmetric, calculate only 1/8 and mirror the rest
      symmetric = True if args.x1 == args.x2 and args.y1 == args.y2 and args.z1 == args.z2 else False
      symmetric = False
      write_profile_to_array(array, profiles[i],args.multiple_regions)
    if args.target == "3dlines":
      if args.save_vtp: #write 3 .vtp files
        points = vtk.vtkPoints()
        lines = vtk.vtkCellArray()
        polygon = vtk.vtkPolyData()
        
        points,lines = write_polylines_to_vtk(profiles[i],res,args.res_surf_lines,points,lines)
        polygon.SetPoints(points)
        polygon.SetLines(lines)
        matviz_vtk.show_write_vtk(polygon,1000,"3dlines_"+str(i)+".vtp")
      else:  
        plot_3dlines(profiles[i], res, args.res_surf_lines, i, ha)
  
  if args.target == "volume_mesh":
#     h = 1.0/res
#     x = np.arange(0,1+h,h)
#     from pyevtk.hl import gridToVTK
#     gridToVTK("array",x,x,x,cellData={"array":array})

    points, cells = voxels_to_points_and_cells(array,multRegions=args.multiple_regions)
      
  if args.target == "surface_mesh" or args.target == "marching_cubes":
    ############################ new surface mesh approach ####################
    # Use marching cubes to obtain the surface mesh of voxelized structure
    # marching_cubes expect float values (not double)
    h = np.float32(1.0/args.res)
    # coords of vertices lie in [0,1-h]
    import time
#     start = time.time()
    verts = []
    faces = []
    normals = []
     
    points, cells, normals, values = measure.marching_cubes_lewiner(array,spacing=(h,h,h),allow_degenerate=False,step_size=1)
    
    # marching_cubes returns float values
    points = np.asarray(points)
    # scale structure to [0,1]^3
    # moves structure to [0,1]^3
    points += (h/2.0,h/2.0,h/2.0)
  
    # extract points on the boundary circles
    # each entry contains a list representing one boundary face of the base cell
    # points: list with 3d coords
    if args.target == "surface_mesh" or args.target == "marching_cubes":
      points, bp_lists = adjust_and_extract_boundary_points(profiles, points)
      if args.target == "surface_mesh" :
        points, cells = mesh_boundary_circles(points, cells, bp_lists, args.bc_flags)
      
#       points, cells, _ = pymesh.remove_duplicated_vertices_raw(points,cells)
#       marching_cubes.write_vtp(points, cells, (h,h,h), "test.vtp")
      
      if args.simplify:  
        marching_cubes.write_vtp(points,cells,(h,h,h),name="mc_lewiner.vtp")
        points, cells = collapse_short_edges(points, cells)
      
  if args.target == '3dlines' and not args.save_vtp:
    plt.show()
    
  return array, points, cells

# creates map with info on profile depending on radius
# Profile contains list of tuples with vector,angle and idx where constant part begins (bisec: res/2, orthogonal: grad is 1)
def create_profile_map(profile,res,verbose=None,save=None,ha=None):
  map = np.zeros((360, res))
  h = 1.0 / res
  if save:
    out = open("profile_map_dir_"+str(profile.direction)+".txt","w")
    out.write("#i \t alpha \t radius\n")
  
  if symmetric:  
    for i,x in enumerate(np.arange(0,1.0,h)):
      for alpha in range(0,360):
        map[alpha,i] = calc_radius_for_quadrant(profile, x, degree_to_rad_quadrant(alpha))
        if save:
          out.write(str(i) + " \t" + str(alpha) + " \t" + str(map[alpha,i]) + "\n")
  else:
    for i,x in enumerate(np.arange(0,1.0,h)):
      for alpha in range(0,360):
        map[alpha,i] = calc_radius(profile, x, np.radians(alpha))
            
  if verbose == "polar_plot":
    plt.gcf().clear()
    plt.rcParams.update({'font.size': 30})
    plt.figure(figsize=(13,13))
    ax = plt.axes(polar=True)
    #ax.set_rlabel_position(-22.5) 
#     ax.set_rticks([0.05,0.1,0.2,0.3])
    #ax.set_rticks([0,0.05, 0.1, 0.15, 0.2])
    theta = np.linspace(0, 2.0*np.pi,360)
    #plt.plot(theta,map[:,0],linewidth=5.0)
#     for i,bisec in enumerate(profile.bisecs_left):
#       phi = bisec.angle + np.pi/2.0 * i
#       plt.plot(phi,map[int(np.degrees(phi)),int(res/2)+1],'k.',color="red",markersize=20)
#     plt.plot(theta,map[:,int(res/4)],linewidth=5.0)
    dir = ""
    plt.plot(theta,map[:,5],linewidth=6.0,label=dirToString(profile.direction)+"=0.1")
    plt.plot(theta,map[:,int(res/4)],linewidth=6.0,label=dirToString(profile.direction)+"=0.25")
    plt.plot(theta,map[:,int(res/2)-2],linewidth=6.0,label=dirToString(profile.direction)+"=0.5")
    plt.legend(loc='best', bbox_to_anchor=(1.07, 1.15))
    
#     for i,bisec in enumerate(profile.bisecs_left):
#       phi = bisec.angle + np.pi/2.0 * i
#       plt.plot(phi,map[int(np.degrees(phi)),int(res/2)-1],'k.',color="red",markersize=20)
#     plt.plot(theta,map[:,int(res/2)-1],linewidth=5.0)
#     plt.show()
    plt.savefig('polar-plot_dir-'+str(profile.direction) + '.png', dpi=1000)
        
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
  nodes = []
   
  nodes = get_surface_lines(profile)
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
    
# calculates interpolated radius for one quadrant (0 to pi/2) for symmetric basecell
# @param profile: contains three profile functions (for 0, phi (bisec) and 90 degree)
# @param x: parameter for function evaluation
# @param rad: radians for evaluation
def calc_radius_for_quadrant(profile,x,rad):
  assert(rad >= 0 and rad <= np.pi/2.0)
  funcs = [profile.splines_left[0],profile.bisecs_left[0],profile.splines_left[1]]
  val = None
  # interpolation is global variable
  assert(interpolation == "linear" or interpolation == "heaviside")
  if interpolation == "linear":
    val = calc_radius_linear_symm(funcs,x,rad)
  else: #heaviside
    val = calc_radius_heaviside_symm(funcs,x,rad)
  
  assert(val is not None)
  return val - 0.5

# similar to calc_radius_for_quadrant, but for asymmetric case
def calc_radius(profile,x,rad):
  assert(rad >= 0 and rad <= 2.0*np.pi)
  splines = None
  bisecs = None
  if x <= 0.5:
    splines = profile.splines_left
    bisecs = profile.bisecs_left
  else:
    splines = profile.splines_right
    bisecs = profile.bisecs_right
  
  assert(splines is not None and bisecs is not None)
  # interpolation is global variable
  assert(interpolation == "linear" or interpolation == "heaviside")
  if interpolation == "linear":
    val = calc_radius_linear(splines,bisecs,x,rad)
  else: #heaviside
    val = calc_radius_heaviside_asymm(splines,bisecs,x,rad)
  
  assert(val is not None)
  return val - 0.5
    
# rasterize given profile function
# if basecell is symmetric, rasterize only 1/8 of structure
# and mirror the rest
# @param array: stores voxelized info on profile
# @param profile: profile of interest
# @param multiple_regions: fill array with 0,1 or 2 (profile direction), else make array binary
def write_profile_to_array(array,profile,multiple_regions):
  res = array.shape[0]
  
  bound = res if not symmetric else int(res/2)
  
  for i in range(0,bound):
    for j in range(0,bound):
      for k in range(0,bound):
        x, y, z = grid_to_cartesian_coords((i,j,k), res)
        valx = cartesian_to_polar(y, z, (0.5,0.5))
        
        phi = angle_to_center((y,z))
        
        r = calc_radius(profile, x, phi)
        # get right indices depending on profile's direction
        major = profile.direction
        minor_1, minor_2 = give_normal_plane_axes(major)
        idx = [i,j,k]
        
        # comparing radii only gives voxels inside the profile, but not
        # voxel that are crossed by a profile surface line;
        # thus, check also if the projection of the voxel center  lies 
        # within the voxel bounds --> better: check for many point samples
        # inside the voxel, but this is too costly
#         projection = calc_projection(profile, center)
        if (valx-r <= 1e-6):
          array[idx[major],idx[minor_1],idx[minor_2]] = major if multiple_regions else 1

  if symmetric:
    # mirror octant
    array[0:bound,0:bound,bound:res] = array[0:bound,0:bound,0:bound][:,:,::-1]
    array[0:bound,bound:res,0:bound] = array[0:bound,0:bound,0:bound][:,::-1,:]
    array[0:bound,bound:res,bound:res] = array[0:bound,0:bound,0:bound][:,::-1,::-1]
    array[bound:res,0:bound,0:bound] = array[0:bound,0:bound,0:bound][::-1,:,:]
    array[bound:res,0:bound,bound:res] = array[0:bound,0:bound,0:bound][::-1,:,::-1]
    array[bound:res,bound:res,0:bound] = array[0:bound,0:bound,0:bound][::-1,::-1,:]
    array[bound:res,bound:res,bound:res] = array[0:bound,0:bound,0:bound][::-1,::-1,::-1] 
  
# helper function for calc_radius_for_quadrant
# return linear interpolation between principal spline and bisec
# @param funcs: array with 3 entries: spline1, bisec, spline2
# @param rad: angle in radians
# @param x: cartesian x-coordinate
def calc_radius_linear_symm(funcs,phi,x,rad):
  phi = float(funcs[1].angle)  # bisec angle
  assert(phi >= 0 and phi <= np.pi/2.0)
  if rad <= phi:
    alpha = 1.0/phi * rad # scale section between 0 and 1
    return  (1 - alpha) * funcs[0].eval(x) + alpha * funcs[1].eval(x)
  else : # rad <= np.pi/2.0
    alpha = (rad-phi) / (np.pi/2.0-phi)  # scale section between 0 and 1
    return  (1-alpha) * funcs[1].eval(x) + alpha * funcs[2].eval(x)

# helper function for calc_radius_for_quadrant
# return heaviside interpolation between principal spline and bisec
# see calc_radius_linear for params
def calc_radius_heaviside_symm(funcs,x,rad):
  phi = float(funcs[1].angle)  # bisec angle
  
  return calc_radius_heaviside((funcs[0],funcs[2]), funcs[1], x, rad,phi)
    
# returns two splines, one bisec function and angle offset depending on angle rad
# called by calc_radius_linear and calc_radius_heaviside
# offset is required as bisec angle is always between 0 and 90 degrees
def get_splines_and_bisec(splines,bisecs,rad):
  assert(rad >= 0)
  assert(len(splines) == 4 and len(bisecs) == 4)
  spline_1 = None
  spline_2 = None
  s1_idx = -1 # rigth splines/bisecs array indices 
  s2_idx = -1 # depends on angle rad
  offset = -1 # bisec angle is always between 0 and 90 degree; need to add offset
  # find out which bisection angle to use
  if rad <= np.pi/2.0: # between 0 and 90 degree
    s1_idx = 0
    s2_idx = 1
    offset = 0
    assert(np.isclose(splines[s1_idx].angle,0,1e-6))
    assert(np.isclose(splines[s2_idx].angle,np.pi/2.0,1e-6))
  elif rad <= np.pi: # between 90 and 180 degree
    s1_idx = 1
    s2_idx = 2
    offset = np.pi/2.0
    assert(np.isclose(splines[s1_idx].angle,np.pi/2.0,1e-6))
    assert(np.isclose(splines[s2_idx].angle,np.pi,1e-6))
  elif rad <= 1.5*np.pi: # between 180 and 270 degree
    s1_idx = 2
    s2_idx = 3
    offset = np.pi
    assert(np.isclose(splines[s1_idx].angle,np.pi,1e-6))
    assert(np.isclose(splines[s2_idx].angle,1.5*np.pi,1e-6))
  else: # between 270 and 360 degree
    s1_idx = 3
    s2_idx = 0
    offset = 1.5*np.pi
    assert(np.isclose(splines[s1_idx].angle,1.5*np.pi,1e-6))
    assert(np.isclose(splines[s2_idx].angle,0,1e-6))
    
  assert(s1_idx >= 0 and s1_idx < 4)
  assert(s2_idx >= 0 and s2_idx < 4)
  spline_1 = splines[s1_idx]
  spline_2 = splines[s2_idx]
  bisec = bisecs[s1_idx]
  
  return spline_1, spline_2, bisec, offset
# calculates profile value (asymmetric case) using linear interpolation
# between splines and bisecs
# @param splines: all 4 splines of certain bar we want to draw
# @param bisecs: all 4 bisecs of certain bar we want to draw
# decision whether we're at left or right part of cell was done by 
# function calling this function
def calc_radius_linear(splines,bisecs,x,rad):
  spline_1, spline_2, bisec, offset = get_splines_and_bisec(splines,bisecs,rad)
  
  # bisection angle
  phi = bisec.angle + offset
  
  if rad <= phi:
    # project interval [offset,phi] onto [0,1]
    alpha = (rad - offset) / (phi - offset)
    return  (1 - alpha) * spline_1.eval(x) + alpha * bisec.eval(x)
  else:
    # project interval [phi,offset+pi/2] onto [0,1] 
    alpha = (rad - phi) / (offset+np.pi/2.0 - phi)
    return  (1 - alpha) * bisec.eval(x) + alpha * spline_2.eval(x)

# @param alpha: interpolation argument lying in range [0,1], project it if outside range
# @param beta and eta: arguments of tanh
def calc_heaviside_interpolation(alpha,beta,eta):
  assert(alpha >= 0 and alpha <= 1)
  eps = 1e-5 #for numerical comparison
  
  # a + c * tanh(...) in [0,1]
  assert(calc_tanh(beta, eta, 1) <= 1)
  assert(calc_tanh(beta, eta, 0) >= 0)
  
  c = 1.0 / (calc_tanh(beta, eta, 1) - calc_tanh(beta, eta, 0))
  a = -c * calc_tanh(beta, eta, 0)
  
  if abs(a+c*calc_tanh(beta, eta, 0)) > eps:
    print("CHI:",a+c*calc_tanh(beta, eta, 0))
  assert(abs(a+c*calc_tanh(beta, eta, 0)) < eps) 
  assert(abs(a+c*calc_tanh(beta, eta, 1)) >= 1-eps)
  
  scale = a+c*calc_tanh(beta, eta, alpha)
  assert(scale >= -eps and scale <= 1.0 + eps)
  
  return scale
# same as calc_radius_linear but with heaviside interpolation
# between splines and bisecs
def calc_radius_heaviside(splines,bisec,x,rad,phi,offset=0):
  # get beta and eta for heaviside function from bisec
  beta = bisec.heaviside.beta
  eta = bisec.heaviside.eta
    
  if rad <= phi:
    # project interval [offset,phi] onto [0,1]
    alpha = (rad - offset) / (phi - offset)
    scale = calc_heaviside_interpolation(1-alpha, beta, eta)
    return scale * splines[0].eval(x) + (1-scale) * bisec.eval(x)
  else:
    # project interval [phi,offset+pi/2] onto [0,1] 
    alpha = (rad - phi) / (offset+np.pi/2.0 - phi)
    scale = calc_heaviside_interpolation(alpha, beta, eta)
    return (1-scale) * bisec.eval(x) + scale * splines[1].eval(x)

def calc_radius_heaviside_asymm(splines,bisecs,x,rad):  
  spline_1, spline_2, bisec, offset = get_splines_and_bisec(splines,bisecs,rad)
  # bisection angle in radians
  phi = bisec.angle + offset
  
  return calc_radius_heaviside((spline_1,spline_2),bisec,x,rad,phi,offset)
    
def write_polylines_to_vtk(profile,res,numLines,points,lines):
  nodes = get_surface_lines(profile)
  for i in range(numLines):
    vtk_ids = []
    for j in range(res):
      vtk_ids.append(points.InsertNextPoint(nodes[i,j,0], nodes[i,j,1], nodes[i,j,2]))
    lines.InsertNextCell(res) # specify number of following nodes
    for j in range(res):
      lines.InsertCellPoint(vtk_ids[j])
  
  return points,lines

# @param angle in radians between 0 and pi/2
# @param interpolation: either "linear" or "heaviside", passed by basecell.py
# @param y1,z1: we want to interpolate between these two values
# @para beta param of tanh(), in case of heaviside interpolation
# returns weighted radius for spline part of bicubic bisec
def weighted_by_angle(angle,interpolation,y1,z1,beta=None):
  """
  # python doctest
  >>> weighted_by_angle(np.pi/4,"heaviside",0.2,0.3,10) != 0.2
  True
  >>> print(round(weighted_by_angle(0,"heaviside",0.2,0.3,10),2))
  0.2
  >>> print(round(weighted_by_angle(np.pi/2,"heaviside",0.2,0.3,10),2))
  0.3
  """ 
  assert(angle >= 0 and angle <= np.pi/2.0)
  w = (angle-np.pi/2.0) / (-np.pi/2.0) # linear interpolation 
  if interpolation == "heaviside":
    assert(beta is not None)
    alpha = angle / (np.pi/2.0) # project angle in [0,pi/2] onto [0,1]
    # set eta = 0.5 as we want weight = 0.5 for an angle of 0.5
    w = 1 - calc_heaviside_interpolation(alpha, beta, 0.5)
  
  assert(w >= 0 and w <= 1)  
  return w*y1+(1-w)*z1


# adjust all points on boundary circles by using information from radius
# extract these adjusted points into 6 lists (xleft,xright,yleft,yright,zleft,zright)
# -> need this for meshing the circles
def adjust_and_extract_boundary_points(profiles,points):
  xmin = min(points, key=lambda t: t[0])[0]
  xmax = max(points, key=lambda t: t[0])[0]
  ymin = min(points, key=lambda t: t[1])[1]
  ymax = max(points, key=lambda t: t[1])[1]
  zmin = min(points, key=lambda t: t[2])[2]
  zmax = max(points, key=lambda t: t[2])[2]
  
#   print("dimensions:",xmin,ymin,zmin,xmax,ymax,zmax)
  
  vertices = [None] * len(points)
  lists = [[] for i in range(6)]
  major_dir = -1
  for id in range(len(points)):
    p = points[id]
    # x
    if np.isclose(p[0],xmin,1e-6):
      # move x to 0
      points[id][0] = 0
      lists[0].append((points[id],id))
      major_dir = 0
    elif np.isclose(p[0],xmax,1e-6):
      # move x to 1
      points[id][0] = 1.0
      lists[3].append((points[id],id))
      major_dir = 0
    # y
    elif np.isclose(p[1],ymin,1e-6):
      # move y to 0
      points[id][1] = 0
      lists[1].append((points[id],id))
      major_dir = 1
    elif np.isclose(p[1],ymax,1e-6):
      # move y to 1
      points[id][1] = 1.0
      lists[4].append((points[id],id))
      major_dir = 1  
    # z
    elif np.isclose(p[2],zmin,1e-6):
      # move z to 0
      points[id][2] = 0
      lists[2].append((points[id],id))
      major_dir = 2
    elif np.isclose(p[2],zmax,1e-6):
      # move z to 1
      points[id][2] = 1.0
      lists[5].append((points[id],id))
      major_dir = 2
    else:
      vertices[id] = points[id]
      continue
    
    assert(major_dir > -1 and major_dir < 3)
    minor1, minor2 = give_normal_plane_axes(major_dir)
    phi = angle_to_center((p[minor1],p[minor2]))
    vertices[id] = radius_to_3d_coords(profiles[major_dir], p[major_dir], phi)
  
  for l in lists:
    for i in range(len(l)):
      l[i] = (vertices[l[i][1]],l[i][1])
  
  for i,l in enumerate(lists):
    major = i%3
    minor1, minor2 = give_normal_plane_axes(major)
    assert(len(l) > 0)
    # sort points in cyclic order
    # here we live in [0,1]^3
    cell_center = np.asarray([0.5,0.5,0.5])
    l.sort(key=lambda c:math.atan2(np.asarray(c[0][minor1])-cell_center[0], np.asarray(c[0][minor2])-cell_center[1]))
    
  return vertices, lists
 
# returns list of coordinates for given boundary 'bound'
# using order {0:"x_left", 1:"y_left", 2:"z_left", 3:"x_right", 4:"y_right", 5:"z_right"}
# @param bpoints: six lists with boundary points for each face of the bounding box
# order: min_x,min_y,min_z,max_x,max_y,max_z
# @param cc_3d: basecell center, has also to be projected onto right plane
def extract_2d_bc_boundary_coords(bpoints,bound):
  result = []
  # depending on boundary flag, determine which coordinate component to compare
  # order: xmin, ymin, zmin, ymax, ymax, zmax
  # bound -> major_dir: 0->0, 1->1, 2->2, 3->0, 4->1, 5->2 
  major_dir = bound%3
  minor_1, minor_2 = give_normal_plane_axes(major_dir)
  
  for p in bpoints[bound]:
#     print("bc_bounds[bound]:",bc_bounds,bound,bc_bounds[bound])
    result.append( ((p[0][minor_1],p[0][minor_2]),p[1]) )
  
  return result  

# @param flags: 6 flags indicating which boundary circle to mesh and which not
# order of flags: xmin,xmax,ymin,ymax,zmin,zmax
# @param bp_lists: 6 lists with boundary points for each bounding box face
# @param points: list with 3d coords
def mesh_boundary_circles(points,cells,bpoints,flags=None):
  print("meshed boundary circles")
  assert(len(bpoints) == 6)
   
  lists_2d = []
  points = list(points)
  cells = list(cells)
  
  # if not set, mesh everything
  if flags is None:
    flags = [True] * 6
  # give basecell boundary circles to be meshed
  bounds_to_mesh = np.where(flags)[0]
  for b in bounds_to_mesh:
    bcoords_2d = extract_2d_bc_boundary_coords(bpoints,b)
    points, cells = mesh_basecell_boundary(points,cells,bcoords_2d,b)
  
  return points, cells    

# @returns 3d cartesian coordinates
# @param profile of interest
# @param point of evaluation of profile, e.g. x=0.5 means y=0.5 if we have y-profile
def radius_to_3d_coords(profile,x,phi):
  r = calc_radius(profile, x, phi)
  px,py = polar_to_cartesian(r,phi)
  point = np.ones(3)*(-1)
  major_dir = profile.direction
  minor_1, minor_2 = give_normal_plane_axes(major_dir)
  point[major_dir] = x
  point[minor_1] = px
  point[minor_2] = py
 
  return point

def voxels_to_points_and_cells(array,multRegions=False,wx=1.0,wy=1.0,wz=1.0):
  resx, resy, resz = array.shape
  
  hx = wx / resx
  hy = wy / resx
  hz = wz / resx
  
  centers = []
  
  for i in range(0,resx):
    for j in range(0,resy):
      for k in range(0,resz):
        x = i * hx + hx/2.0
        y = j * hy + hy/2.0 
        z = k * hz + hz/2.0
        
        if (array[i,j,k] >= 0 and multRegions) or(array[i,j,k] > 0 and not multRegions):
          if i > 0 and j > 0 and k > 0 and i < resx-1 and j < resy-1 and k < resz - 1:
            if multRegions:
              if array[i-1,j,k] < 0 or array[i+1,j,k] < 0 or array[i,j-1,k] < 0 or array[i,j+1,k] < 0 or array[i,j,k-1] < 0 or array[i,j,k+1] < 0:
                centers.append([x,y,z])
            else:      
              if array[i-1,j,k] == 0 or array[i+1,j,k] == 0 or array[i,j-1,k] == 0 or array[i,j+1,k] == 0 or array[i,j,k-1] == 0 or array[i,j,k+1] == 0:
                centers.append([x,y,z])
          else:
            centers.append([x,y,z])       
        
  # dummy objects, create_centered_bars needs them
  vtk_points = vtk.vtkPoints()
  vtk_cells = vtk.vtkCellArray()
  # simple python lists
  points, cells = matviz_vtk.create_centered_bars(vtk_cells,vtk_points,centers,[hx,hy,hz])
  
  return points, cells

# here we only deal with 3 dimensions
# @param points: list of all created 3d point coords
# @param bound: which of the 6 boundary faces to mesh
# order: xmin,xmax,ymin,ymax,zmin,zmax
def mesh_basecell_boundary(points,cells,coords_2d,bound):
  # here we live in [0,1]^3
  cell_center = np.asarray([0.5,0.5,0.5])
  # need this for mapping between planar and space coordinate
  major_dir = bound%3
  minor_dir_1, minor_dir_2 = give_normal_plane_axes(major_dir)
  # sort points in circle order
  coords_2d.sort(key=lambda c:math.atan2(np.asarray(c[0][0])-cell_center[0], np.asarray(c[0][1])-cell_center[1]))

  test = [ [np.float64(elem[0][0]),np.float64(elem[0][1])] for elem in coords_2d]

  import triangle

  tri = dict(vertices=np.asarray(test))
  mesh = triangle.triangulate(tri, 'q32.5a0.1C')
#   import triangle.plot as plot
#   import matplotlib.pyplot as plt  
#   triangle.plot.compare(plt, tri, mesh)
#   plt.show()
  mesh_points = mesh['vertices']
  mesh_tris = mesh['triangles']
  
  # up to len(l), l and mesh_points have the same ordering of points
  # map from local mesh_points point ids to global ones
  map = [-1] * len(mesh_points)
  assert(len(mesh_points) > 0)
  for i in range(0,len(coords_2d)):
#     print(coords_2d[i][1])
#     print(map[i])
    map[i] = int(coords_2d[i][1])
    assert(type(coords_2d[i][1]) is int)
    
  ######### mapping back to 3d ########################
  # 0,1,2 -> 0.0  3,4,5 -> 1.0
  comp = 0.0 if 0 <= bound <= 2 else 1.0
  new_points = []
  next_id = len(points)
  # map from 2d point to 3d point
  for i in range (len(coords_2d),len(mesh_points)):
    new_p = np.zeros(3)
    new_p[major_dir] = comp
    new_p[minor_dir_1] = mesh_points[i][0]
    new_p[minor_dir_2] = mesh_points[i][1]
    # id of new pointcoords_2d
    new_points.append(new_p)
    assert(type(next_id) is int)
    map[i] = next_id
    next_id +=1
  
  points.extend(new_points)
  # use lookup table to set new triangles from meshed boundary circle  
  for tri in mesh_tris:
    cells.append((map[tri[0]], map[tri[1]], map[tri[2]]))
    assert(type(map[tri[0]]) is int)
    assert(type(map[tri[1]]) is int)
    assert(type(map[tri[2]]) is int)
#   import matplotlib
#   matplotlib.use('tkagg')
#   from matplotlib import pyplot as plt
#   plt.gcf()
#   labels = []
#   for i in range(len(points)):
#     labels.append(i)
#   points = np.asanyarray(points)
#   plt.plot(points[:,minor_dir_1],points[:,minor_dir_2],'o')
#   for i, label in enumerate(labels):
#     plt.text(points[i,minor_dir_1],points[i,minor_dir_2],labels[i])
#   plt.triplot(points[:,minor_dir_1],points[:,minor_dir_2],cells)
#   plt.show()
    
  return points,cells


def cartesian_to_voxel_coords(point,minx,miny,minz,hx,hy,hz,log=False):
  """
    Returns voxel coordinates (i,j,k) for given cartesian coords (x,y,z)
    @param point described in cartesian coords
    @param minx: smallest x value of domain
    @param miny: smallest y value of domain
    @param minz: smallest z value of domain
    @param hx,hy,hz: lattice spacing
    
    >>> cartesian_to_voxel_coords((0.0,0.0,0.0),0.0,0.0,0.0,1,1,1)
    (0, 0, 0)
    
    >>> cartesian_to_voxel_coords((0.4,0.4,0.4),0.0,0.0,0.0,0.1,0.1,0.1)
    (4, 4, 4)
    
    >>> cartesian_to_voxel_coords((0.4,0.4,0.4),-2.0,5.0,0.4,0.1,0.1,0.1)
    (23, -45, 0)
  """
  if log:
    print("point:",point)
  i = int((point[0]-minx) / hx-1e-6)
  j = int((point[1]-miny) / hy-1e-6)
  k = int((point[2]-minz) / hz-1e-6)
  
  return i,j,k

def voxel_to_cartesian_coords(voxel,lbounds,h):
  """ 
    Returns cartesian coordinates (x,y,z) for given voxel coords (i,j,k)
    @param voxel: (i,j,k)
    @param lbounds: lower bounds (3 values) of domain
    @param h: lattice spacings for the 3 coordinate directions
    
    >>>   
  """
  
  assert(len(lbounds) == len(h))
  x = voxel[0] * h[0] + lbounds[0] + 1e-6
  y = voxel[1] * h[1] + lbounds[1] + 1e-6
  z = voxel[2] * h[2] + lbounds[2] + 1e-6
  
  return np.array([x,y,z])

def calc_edge_lengths(mesh):
  import pymesh
  minEdge = 9999
  av = 0
  maxEdge = -9999
  for i,f in enumerate(mesh.faces):
    # a triangle faces has 3 edges
    v0 = mesh.vertices[f[0]]
    v1 = mesh.vertices[f[1]]
    v2 = mesh.vertices[f[2]]
   
    tmp = min(np.linalg.norm(v1-v0),np.linalg.norm(v2-v0),np.linalg.norm(v2-v1))
    if minEdge > tmp:
      minEdge = tmp
    tmp = max(np.linalg.norm(v1-v0),np.linalg.norm(v2-v0),np.linalg.norm(v2-v1))
    if maxEdge < tmp:
      maxEdge = tmp

    av += np.linalg.norm(v1-v0)
    av += np.linalg.norm(v2-v0)
    av += np.linalg.norm(v2-v1)

  av = av / 3.0 / len(mesh.faces)
   
  print("min,max,av:",minEdge,maxEdge,av)
  return minEdge, maxEdge, av

# use pymesh to collapse short edges and afterwards to repair obtuse triangles
def collapse_short_edges(verts,faces,abs_thresh=None):
  import pymesh
  mesh = pymesh.form_mesh(np.asarray(verts),np.asarray(faces))
  minl,maxnl,avl = calc_edge_lengths(mesh)
  t = abs_thresh
  if abs_thresh is None:
    t = 0.8*avl
  mesh, info = pymesh.collapse_short_edges(mesh, abs_threshold=t,preserve_feature=True)
  print("info:",info)
  mesh, info = pymesh.remove_obtuse_triangles(mesh,130)
  print("info:",info)
  return mesh.vertices, mesh.faces
