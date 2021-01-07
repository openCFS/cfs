#!/usr/bin/env python
# this ia a tool for spaghetti visualization and calculation.
# code segments are based on LuFoVI/code/spaghetti.py by Daniel and Eloise 
import numpy as np
from numpy.linalg import norm 
import os

# for automatic differentiation
import autograd
import autograd.numpy as agnp

# interactive
if __name__ == '__main__':
  import optimization_tools as ot
  import matplotlib
  # necessary for remote execution, even when only saved: http://stackoverflow.com/questions/2801882/generating-a-png-with-matplotlib-when-display-is-undefined
  #matplotlib.use('Agg')
  #matplotlib.use('tkagg')
  import matplotlib.pyplot as plt
  import matplotlib.colors as colors
  import matplotlib.cm as cmx
  from matplotlib.path import Path
  import matplotlib.patches as patches
  import argparse


# to conveniently access global values
class Global:
  # the default value are vor easy debugging via import spaghetti
  def __init__(self):
    self.shapes = []         # array of Spaghetti 
    self.rhomin = 1e-6     
    self.radius = .05        # up to now constant for all spaghetti for arc
    self.boundary = 'poly'   # up to now only 'poly'
    self.transition = 0.015  # paramter for boundary: 2*h
    self.combine= 'max'      # up to now only 'max'
    self.n = [30,30,1]        # [nx, ny, nz]
    self.opts = {} # there is an optional dictionary in xml in the python element
    
    # this a arrays of size (nx+1,ny+1) with 
    self.idx_field = None # entries are int vectors of size shapes. -2 far outside, -1 far inside, >= 0 nearest part by shape 
    self.dist_field = None # entries are nodal closest distance by shape
    
glob = Global()    
    
## This functions are called from openCFS via SpaghettiDesign.cc

# called from SpaghettiDesign.cc constructor
# @radius a constant in meter, e.g. .15
# @boundary only 'poly'
# @transition transition zone which is 2*h
# @combine how to combine? currently only 'max'
# @nx, ny, nz are for rho. Currently nx == ny and nz == 1
def cfs_init(rhomin, radius, boundary, transition, combine, nx, ny, nz, dict):
  glob.rhomin = rhomin
  glob.radius = radius
  glob.boundary = boundary
  glob.transition = transition
  glob.combine = combine
  assert nx == ny and nz == 1   
  glob.n = [nx, ny, 1]
  glob.opts = dict
  print('cfs_init called: ', glob.n, dict) 
  
  
## set spaghetti. Initially create it, otherwise update
# @param id the 0-based index of the spaghetti, corresponds index within glob.shapes
def cfs_set_spaghetti(id, px, py, qx, qy, a_list, p):    
  assert id < len(glob.shapes) + 1 # we may add another spaghetti
  
  P = [px,py]  
  Q = [qx,qy]
  
  if id >= len(glob.shapes):
    base = sum([len(s.var()) for s in glob.shapes])
    # def __init__(self,id,base,radius,P,Q,a,p):
    glob.shapes.append(Spaghetti(id, base, glob.radius, P, Q, a_list, p))
    print('cfs_set_spaghetti: create ', glob.shapes[-1])
  else:
    glob.shapes[id].set(P, Q, a_list, p)
    print('cfs_set_spaghetti: update ', glob.shapes[id])
  
## give back the density as 1D numpy arry with the current spaghetti setting
def cfs_map_to_density():
  print('cfs_map_to_density: called')
  
  assert glob.n[0] == glob.n[1]
  nx = glob.n[0]
  
  if glob.idx_field is None or glob.dist_field is None:
    glob.idx_field, glob.dist_field = create_idx_field() 

  var = glob.shapes[0].var()

  rho = np.ones((nx, nx))  
  for j in range(nx):
    for i in range(nx):
      rho[i,j] = integrate_rho(var,i, j, ad=False) 

  return rho.reshape(np.prod(glob.n), order='F')
  
# as we cannot create a numpy array in C (it should work but fails in reality) we get it here.
# it shall have the size of rho as a 1D array  
def cfs_get_drho_vector():  
  print('cfs_get_drho_vector: returns zero vector of size ', np.prod(glob.n))
  assert np.prod(glob.n) >= 1

  return np.zeros((np.prod(glob.n)))
  
## get the gradient of a density based function w.r.t. all shape variables.
# we calculate it for all and don't skip fixed variables
# @param drho is a 1D np.array with d func/d rho - this is to be handled via chainrule. 
# @return d func/ d shape  
def cfs_get_gradient(drho):
  n  = glob.n
  nx = n[0]
  # see density() for comment on layout
  assert len(drho) == np.prod(n)

  if glob.idx_field is None or glob.dist_field is None:
    glob.idx_field, glob.dist_field = create_idx_field() 



  print('cfs_get_gradient: return trivial vector only - Implement it!')

  total = sum([len(s.var()) for s in glob.shapes])
  return np.ones(total)    
      
# own autograd compatible norm as agnp.linalg.norm seems not to work
def agnorm(X):
  sum = 0
  for e in X:
    sum += e * e
    
  return agnp.sqrt(sum)

# find the mininimum of (value, index) tuples, where only value is compared and value might contain None
def idx_min(a,b):  
  if a[0] == None:
    return b
  else:
    if b[0] == None or a[0] < b[0]:
      return a
    else:
      return b  
  

# convert string to integer dof: x= 0 or int to string
def dof(val):
  if val == 0:
    return 'x'
  if val == 1:
    return 'y'
  if val == 2:
    return 'z'
  if val == 'x':
    return 0
  if val == 'y':
    return 1
  if val == 'z':
    return 2
  assert(False)

# minimal and maximal are vectors.
def create_figure(res, minimal, maximal):

  dpi_x = res / 100.0 

  fig = matplotlib.pyplot.figure(dpi=100, figsize=(dpi_x, dpi_x))
  ax = fig.add_subplot(111)
  ax.set_xlim(min(0,minimal[0]), max(1,maximal[0]))
  ax.set_ylim(min(0,minimal[1]), max(1,maximal[1]))
  return fig, ax

def dump_shapes(shapes):
  for s in shapes:
     print(s)   

colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
# transforms ids from 0 to 6 to color codes 'b' to 'k'. Only for matplotlib, not for vtk!
def matplotlib_color_coder(id):
  return colors[id % len(colors)]
    
class Spaghetti: 
  # @param id starting from 0
  # @param base sum of all len(var()) of all previous shapes. Starts with 0
  def __init__(self,id,base,radius,P,Q,a,p):
    self.id = id
    self.base = base
    self.color = matplotlib_color_coder(id)
    self.radius = radius  
    self.set(P,Q,a,p)

    # derivate of fast_dist_ad by automatic differentiation
    self.grad_fast_dist = autograd.grad(self.fast_dist_ad,argnum=0)

  
  # does not only se the variables but also computes the internal helpers (C, ...) and deletes glob index fields
  def set(self, P,Q,a,p):
    
    glob.idx_field = None
    glob.dist_field = None
    
    assert len(P) == 2 and len(Q) == 2  
    self.P = np.array(P) # start coordinate
    self.Q = np.array(Q) # end coordinate
    self.a = a # vector of normals, might be empty
    self.p = p # profile is full structure thickness
    self.w = p/2. # half structure
    
    # helper
    self.n = len(a) + 1 # number of segments
  
    self.U = self.Q - self.P # base line p -> q
    self.U0 = self.U / norm(self.U)
  
    self.V = np.array((-self.U[1],self.U[0])) # normal to u
    self.V0 = self.V / norm(self.V)
    
    self.E = [] # endpoints within p-q line but without p and q
    self.H_int = [] # real summits without p and q
    for i in range(1,self.n):
      self.E.append(self.P + i * self.U/self.n)
      self.H_int.append(self.E[-1] + self.a[i-1] * self.V0)
    assert len(self.E) == len(self.a)

    # h = p + h_int + q where the summits are surrounded by p and q as "artificial" summits
    self.H = []
    self.H.append(self.P) # use the numpy arrays
    self.H.extend(self.H_int)
    self.H.append(self.Q)

    # list of segments as tuple of start and endpoint (P,H_*,Q), for 1 segement (no internal summit) this is P->Q
    # note there is also list straight segment parts L from P->K_1s, K_1e->K_2s, ...->Q
    self.segs = []
    for i in range(len(self.H)-1):
      self.segs.append((self.H[i], self.H[i+1]))

    # list of segment lines as vectors T
    self.T = []
    for seg in self.segs:
      self.T.append(seg[1] - seg[0])

    # list of normals M for each segment
    self.M = []
    for T in self.T:
      self.M.append([-T[1], T[0]] / norm(T))

    # For the arcs the center points C_* 
    # the idea is to find the angle alpha between two segments, have vector B at half the angle
    # now scale radius/sin(alpha/2) along B0 to find C 
    self.C = []
    # note there is also list straight segment parts L from P->K_1s, K_1e->K_2s, ...->Q. For one seg this is P->Q
    self.L = []
    self.L.append([self.P]) # we will append K1
    # list of arc is only meant as a helper for plotting the arcs
    self.arc = [] 
    r = self.radius
    for i in range(len(self.segs)-1):
      H = self.segs[i][1] # summit point is end of first segment
      v1 = -self.T[i]     # from H to start of first segment 
      v2 = self.T[i+1]    # from H to end of second segment
      M1 = self.M[i]   # normal of first segment
      M2 = self.M[i+1] # normal of second segment  
      # get angle between two segments, standard form for angle between two vectors
      # cos(alpha) = v1*v2/ (||v1|| * ||v2||)
      cosa = np.dot(v1,v2)/(norm(v1) * norm(v2))
      # If v1 and v2 are aligned, it might happen that numerically abs(cosa)>1.
      # Eloise did rounding for 12 digits, however this kills AD with wrong results w.o. warning!
      cosa = 1 if cosa > 1 else cosa
      cosa = -1 if cosa < -1 else cosa

      alpha = np.arccos(cosa)
      B = (M1 + M2)/2
      B0 = B/norm(B)
      
      # curvature is negative. center is right of spline in direction p2-p1 -> C = H - r/np.sin(alpha/2) * B0, K1=C + r*M1
      # curvature is positive. center is left of spline in direction p2-p1 -> C = H + r/np.sin(alpha/2) * B0, K1=C - r*M1
      nc = -1.0 if np.dot(M1, v2) < 0 else 1.0
      
      # from function of hypotenuse of a right-angled triangle
      C = H + nc*(r/np.sin(alpha/2) * B0)
      K1 = C - nc*(r*M1)
      K2 = C - nc*(r*M2)
    
      self.C.append(C)
      self.L[-1].append(K1) # finish last (initial) L
      self.L.append([K2]) # to have added K1 for next arc in next iteration or Q after the loop
      
      # only as helper for plotting the stuff.
      self.arc.append((C,cosa, K1, K2, nc < 0))
      
    self.L[-1].append(self.Q)
    assert(len(self.L[-1]) == 2)
    
    assert len(self.C) == len(self.H_int)
    assert len(self.C) == len(self.L)-1 # one segment more than arcs
    assert len(self.L) == len(self.T)   # L is part of the T vector

  # give optimzation variables as defined array such that we can easily differentiate
  # @return p_x,p_y,q_x,q_y,p,a_1,...a_n-1
  def var(self):
    r = [self.P[0],self.P[1],self.Q[0],self.Q[1],self.p]
    for a in self.a:
      r.append(a)

    return r  

  # helper which gives the name of the variables returned by var()
  def varnames(self):
    r = ['px','py','qx','qy','p']
    for i in range(len(self.a)):
      r.append('a' + str(i+1))

    return r


  # search closest distance to all segments and circs from point X
  # @where if 's' only segments, if 'c' only end points and arcs, None for all
  # @return distance or None when where = 's' and outside, second value part index: 0..n-1 = segment. n...2*n-2=arc, 2n-1=P,2n=Q, -1 = None value 
  def dist(self, X, where = None):

    w = self.w # profile, we are negative inside and 0 at the boundary
    n = self.n # number of segments, one arc less
    
    minimal = (None,-1)

    # first n segments 0...n-1
    if not where or where == 's':
      for i in range(len(self.L)):    
        g = self.L[i][0] # start point of the real line segment which is P or K of left arc
        h = self.L[i][1] # end point of the real line segment which is Q or K of right arc
        M = self.M[i]
        # we are in range if (X-h) @ (h-g) <= 0 and (X-g) @ (h-g) >= 0.
        if (X-h) @ (h-g) <= 0 and (X-g) @ (h-g) >= 0:  # one can precompute h-g
          # (X-g) @ M is distance to segment. Positive on M side, negative on the other side
          minimal = idx_min(minimal, (abs((X-g) @ M) - w,i)) # don't forget abs!
    
    # end points and arcs
    if not where or where == 'c':
      # n-1 arcs n ... 2n-2
      r = self.radius
      for i, C in enumerate(self.C):
        # check if we are in the code (base is the center)
        v1 = self.T[i] 
        v2 = -self.T[i+1]
        XC = X - C
    
        if v1 @ XC >= 0 and v2 @ XC >=0: # we are in the code  
          d = abs(norm(XC) - r) - w
          minimal = idx_min(minimal, (d, n+i))
      
      assert minimal[1] <= 2*n-2
      
      # distance to start point
      minimal = idx_min(minimal,(norm(X-self.P)-w,2*n-1)) 
      # distance to end point
      minimal = idx_min(minimal, (norm(X-self.Q)-w,2*n)) # faster would be root after min() 

    return minimal 

  # gives the closest distance by part idx - see fast_dist_ad for documentation
  # this version tries to be fast by using cached data
  # @param idx is the param, -1 for far inside and -2 for for outside. 
  # @return +/-1 3*glob.transition if far inside or far outside
  def fast_dist(self,X,idx):

    # this comes from create_idx_field which is -1/-2 if poly() does not need to be applied    
    if idx == -1:
      return -3 * glob.transition
    if idx == -2:
      return 3 * glob.transition
    
    # constants
    n = self.n
    w = self.w

    P = self.P
    Q = self.Q
    
    # endpoints are the most easy ones
    if idx == 2*n-1:
      return norm(X-P) - w 
    if idx == 2*n:
      return norm(X-Q) - w 

    # the index of our segment or arc
    i = idx if idx < n else idx - n
    assert i >= 0 and i < 2*n 
    H = self.H[i]
    M = self.M[i]

    # segments: abs((X-g) @ M) - w but instad of g we can use H
    if idx < n:
      return abs((X-H) @ M) - w
    
    assert i >= 0 and i <= n-2
    r = self.radius
    C = self.C[i]
    return abs(norm(X-C) - r) - w

  # gives the closest distance by part idx - see second return value from dist()!
  # does not check if the part is applicable to the coordinate (e.g. if we are for arc within the cone).
  # is meant as basis for automatic differentiation, not fast as it used no precompiled
  # param var the variables which are p_x,p_y,q_x,q_y,a_1,...a_n-1,w -> var()
  # param X the coordinate as array, 
  # param idx  0..n-1 = segment. n...2*n-2=arc, 2n-1=P,2n=Q,
  # only the distance itself 
  def fast_dist_ad(self,var,X,idx):
    # main changes for autograd
    # - use autograd.numpy (agnp)
    # - any array for computation needs to be agnp.array, e.g. V = [-U[1],U[0]] cannot be differentiated

    if idx == -1:
      return -3 * glob.transition
    if idx == -2:
      return 3 * glob.transition

    # design variables
    assert len(var) == 4 + 1 + len(self.a) 
   
    P = agnp.array(var[0:2])
    Q = agnp.array(var[2:4])
    p = var[4]
    w = p/2.0
    a = var[5:]

    # constants
    n = len(a)+1
    assert n == self.n

    # endpoints are the most easy ones
    if idx == 2*n-1:
      return agnorm(X-P) - w 
    if idx == 2*n:
      return agnorm(X-Q) - w 

    # the index of our segment or arc
    i = idx if idx < n else idx - n
    assert i >= 0 and i < n 
    U = Q - P
    V0 = agnp.array([-U[1],U[0]]) / agnorm(U) # normal to U and normalized
    assert len(self.a) == n-1
 
    H_s = P if i == 0   else P + i/n * U + a[i-1] * V0    # summit of begin of segment
    H_e = Q if i >= n-1 else P + (i+1)/n * U + a[i] * V0  # summit of end of segment, which is Q for last segment
      
    T = H_e - H_s # segment as vector   
    M = agnp.array([-T[1], T[0]]) / agnorm(T) # normal

    # segments: abs((X-g) @ M) - w but instad of g we can use H_e
    if idx < n:
      return agnp.abs(agnp.dot((X-H_e), M)) - w
    
    # in the arc case, we need two segmens and such three summits. we now use s(start=old s), c(center=old e), f(final, new)
    assert i >= 0 and i <= n-2
    r = self.radius
    H_c = H_e
    H_f = Q if i == n-2 else P + (i+2)/n * U + a[i+1] * V0
    
    v1 = H_s - H_c
    v2 = H_f - H_c
    M1 = M
    M2 = agnp.array([-v2[1],v2[0]]) / agnorm(v2)
    
    cosa = agnp.dot(v1,v2)/(agnorm(v1) * agnorm(v2))
    cosa =  1 if cosa >  1 else cosa
    cosa = -1 if cosa < -1 else cosa

    alpha = agnp.arccos(cosa)

    B = (M1 + M2)/2
    B0 = B/agnorm(B)
      
    nc = -1.0 if agnp.dot(M1, v2) < 0 else 1.0
      
    C = H_c + nc*(r/agnp.sin(alpha/2) * B0)
   
    # arcs: abs(norm(XC) - r) - w
    return abs(agnorm(X-C) - r) - w
   
  # give gradient via autograd for given X and idx for all parameters
  # @return array of size var
  def ddist(self,var, X,idx): 
    t = self.grad_fast_dist(var, X,idx)
    assert len(var) == len(t)

    # t is a tuple of unshaped np.arrays - these are strange to extract    
    J = [t[i][()] for i in range(len(t))] 
    
    return J  
   
  # give string with summary information
  def summary(self):
    return "id=" + str(self.id) + " n=" + str(self.n)

  # print the shape info
  def __str__(self):
    return "id=" + str(self.id) + " P=" + str(self.P) + " Q=" + str(self.Q) + " a=" + str(self.a) \
         + " p=" + str(self.p) + " radius=" + str(self.radius) + " color=" + self.color  

  # shape info for a given index
  def to_string(self, idx):
    return "shape=" + str(self.id) + " color=" + str(self.color) + " idx=" + str(idx) + " a=" + str(self.a[idx])
   

   
   
# returns the polynomial boundary modelling. 
# @param dist positive inside, negative outside, boundary at 0
# @param transition is 2*h in feature mapping paper
# @return if outside |transition| -> rhomin or 1 
def poly(dist):   
   phi = -dist # positive inside, negative outside
   h = glob.transition/2.0
   rm = glob.rhomin
   if phi <= -h:
     return rm
   if phi >= h:
     return 1.0
   return 3.0/4.0*(1.0 - rm) * (phi/h - phi**3/(3*h**3)) + .5 * (1+rm)   
   
# returns the nodal density value, is ad_differentiable
# the fast is not for performant calculation (it is NOT) but for having the idx given
# @param the variables - make it robust to handle all noodles
# @param X the coordinate
# @param idx the precomouted part index (see Spaghetti.dist()) - make a vector of indices
def fast_rho_ad(var, X, idx):   
  s = glob.shapes[0]
  
  val = s.fast_dist_ad(var,X,idx)
  return poly(val)
   
def fast_rho(X, idx):   
  s = glob.shapes[0]
  
  val = s.fast_dist(X,idx)
  return poly(val)


# differentiate fast_rho_ad w.r.t. to the var vector -> gives a vector
grad_fast_rho = autograd.grad(fast_rho_ad,argnum=0) 
   
# ad differentiate fast_rho_ad, does some convience for the return type
def drho(var, X, idx):   
  t = grad_fast_rho(var, X,idx)
  assert len(var) == len(t)

  # t is a tuple of unshaped np.arrays - these are strange to extract    
  J = [t[i][()] for i in range(len(t))] 
    
  return J 
   
   
# create a idx field for fast access where >= 0 for closest part, -1 for inside and -2 for too far outside
# uses glob   
# @return fields of size (n+1)^dim, to capture the nodes of the density elements. 
#         The first field contains index vectors (len(shapes)), the second a vector of rhos (len(shapes)) 
def create_idx_field():
  
  shapes = glob.shapes
  n = glob.n[0]
  N = n+1
  
  h = glob.transition *.55 # is save to add a little buffer
  
  x_ls = np.linspace(0,1,N)
  
  idx = np.ones((N,N,len(shapes)),dtype=int) 
  dist  = np.ones((N,N,len(shapes)))

  for j, y in enumerate(x_ls):
    for i, x in enumerate(x_ls):
      X = [x,y]
      for si, s in enumerate(shapes):
        d, k = s.dist(X)
        dist[i,j,si] = d
        idx[i,j,si]  = k if d > -h and d < h else (-1 if d < h else -2)   

  return idx, dist
      

# integrate a rho for element with indices i and j (however the ordering is)
# uses glob.idx_field
# @param if ad fast_dist_ad evaluated, otherwise glob.dist_field is used 
def integrate_rho(var, i, j, ad = False):   
  rho = -1.0
      
  idx_field = glob.idx_field
  s = glob.shapes[0]
  
  # we take the indices as we need them for higher order integration 
  idx1 = idx_field[i,j][0]
  idx2 = idx_field[i+1,j][0]
  idx3 = idx_field[i,j+1][0]
  idx4 = idx_field[i+1,j+1][0]
 
  dx = 1/(glob.n[0]+1)
 
  if idx1 == idx2 == idx3 == idx4 == -1:
    rho = 1.0
  elif idx1 == idx2 == idx3 == idx4 == -2:
    rho = glob.rhomin              
  else:
    if ad:
      x = i * dx
      y = j * dx

      rho = .25 * (poly(s.fast_dist_ad(var,[x,y], idx1)) 
                 + poly(s.fast_dist_ad(var,[x+dx,y], idx2))       
                 + poly(s.fast_dist_ad(var,[x,y+dx], idx3))
                 + poly(s.fast_dist_ad(var,[x+dx,y+dx], idx4)))
    else:
      dist_field = glob.dist_field    
     
      rho = .25 * (poly(dist_field[i,j]) 
                 + poly(dist_field[i+1,j])
                 + poly(dist_field[i,j+1]) 
                 + poly(dist_field[i+1,j+1]))    

  return rho

 
# generates a density map for a unit square. 
# this is a trivial implementation, serving for reference which whall be deleted in near future     
def density(nx):
  if len(glob.shapes) != 1: 
    print("Warnung: density(nx) only implemented for first shape")    

  s = glob.shapes[0]
     
  # the serial element list in cfs is row wise orderd with lower left first and upper right last
  
  # when we operate with optimization_tools with arrays of density, they are column wise ordered
  # with first row in array is left column in image space from bottom  up
  # wich is first colum in array is lower row in image space from left to right
  # stupid enough a simple reshape() does not do the arrangement we want in optimization_tools
  
  # however, a numpy.reorder(order='F') does the job and we need again Fortran reordering to 
  # go back from linear to matrix. 
  #
  # That's simply because optimization_tools was written without deeper thinking in the beginning :(   
     
  nodal_val = np.ones((nx+1,nx+1))
  nodal_idx = np.ones((nx+1,nx+1),dtype=int) * -1
  
  for i,x in enumerate(np.linspace(0,1,nx+1)):
    for j,y in enumerate(np.linspace(0,1,nx+1)): 
      X = [x,y]
      val,idx = s.dist(X)
      nodal_val[i,j] = poly(val)
      nodal_idx[i,j] = idx
      
  rho = np.ones((nx,nx))     
  for i in range(0,nx):
    for j in range(0,nx): 
      rho[i,j] = .25 * (nodal_val[i,j] + nodal_val[i+1,j] + nodal_val[i,j+1] + nodal_val[i+1,j+1]) 
  
  return rho    
      
# reads 2D and returns list of Shaghettis
# @param radius if given overwrites the value from the xml header
def read_xml(filename, set = None, radius = None):
 
  xml = ot.open_xml(filename)
 
  shapes = []
  sq = 'last()' if set == None else '@id="' + str(set) + '"'

  if not radius:
    radius = float(ot.xpath(xml, '//header/spaghetti/@radius')) 

  while True: # exit with break
    idx = len(shapes)
    base = '//set[' + sq + ']/shapeParamElement[@shape="' + str(idx) + '"]'

    test = xml.xpath(base) 
    if len(test) == 0:
      break

    # our xpath helper in optimization tools expects exactly one hit
    Px = float(ot.xpath(xml, base + '[@dof="x"][@tip="start"]/@design'))
    Py = float(ot.xpath(xml, base + '[@dof="y"][@tip="start"]/@design'))
    Qx = float(ot.xpath(xml, base + '[@dof="x"][@tip="end"]/@design'))                   
    Qy = float(ot.xpath(xml, base + '[@dof="y"][@tip="end"]/@design'))
    # with of noodle is 2*w -> don't confuse with P
    p  = float(ot.xpath(xml, base + '[@type="profile"]/@design'))
    a = []
    last = -1
    list = xml.xpath(base + '[@type="normal"]') 
    for normal in list: 
      nr = int(normal.get('nr'))
      des = float(normal.get('design'))
      if not nr > last:
        raise('numbering for normal in shape ' + str(idx) + ' seems out of order: nr=' + str(nr))
      a.append(des)
      last = nr 

    base = sum([len(s.var()) for s in shapes])
    noodle = Spaghetti(id=idx, base=base, radius=radius, P=(Px,Py), Q=(Qx,Qy), a=a, p=p)
    shapes.append(noodle)  
    print('# read noodle', noodle)
      
  return shapes   

    
# creates a matplotlib figure     
# return fig
def plot_data(res, shapes, detailed):

  # could respect non-unit regions and out of bounds movement
  minimal = [0,0]
  maximal = [1,1] 
  
  fig, sub = create_figure(res, minimal, maximal)
  
  for s in shapes:
    
    if detailed:
      # plot tangent lines with extended summits
      for seg in s.segs:
        p1 = seg[0]
        p2 = seg[1]
        l = plt.Line2D((p1[0],p2[0]),(p1[1],p2[1]), color=s.color)
        sub.add_line(l)
  
      # start and endpoint
      fig.gca().add_artist(plt.Circle(s.P, 0.01, color = s.color))
      fig.gca().add_artist(plt.Circle(s.Q, 0.01, color = s.color))
      sub.add_line(plt.Line2D((s.P[0],s.Q[0]),(s.P[1],s.Q[1]), color= 'gray'))
      
      for E in s.E:
        fig.gca().add_artist(plt.Circle(E, 0.005, color = 'black'))
  
      for H in s.H_int: # the outer H which is P and Q is already in L
        fig.gca().add_artist(plt.Circle(H, 0.005, color = 'red'))
        
      for C in s.C:  
        fig.gca().add_artist(plt.Circle(C, 0.005, color = 'blue'))
  
      for L in s.L:  
        fig.gca().add_artist(plt.Circle(L[0], 0.005, color = 'gray'))
        fig.gca().add_artist(plt.Circle(L[1], 0.005, color = 'gray'))
      
      # plot normals
      assert len(s.T) == len(s.M)
      for i, T in enumerate(s.T):
        M = s.M[i]
        p1 = s.H[i] + .5*T
        p2 = p1 + .1*M 
        
        sub.add_line(plt.Line2D((p1[0],p2[0]),(p1[1],p2[1]), color= 'red'))

    # plot rectangles
    for i, L in enumerate(s.L):
      L1 = L[0]
      L2 = L[1]
      M = s.M[i]
      w = s.w
      
      verts = [L1-w*M, L2-w*M, L2+w*M, L1+w*M, L1-w*M,]
      codes = [Path.MOVETO, Path.LINETO, Path.LINETO, Path.LINETO,Path.CLOSEPOLY,]
      path = Path(verts, codes)
      patch = patches.PathPatch(path,facecolor=s.color, lw=1, alpha=.3)
      sub.add_patch(patch)

    # plot half circles for start and end of noodle 
    M = s.M[0] # normal of the first segment tells us where to draw the radius
    angle = np.arctan2(M[1], M[0])*180/np.pi
    sub.add_patch(patches.Arc(s.P, 2*s.w, 2*s.w, theta1=angle, theta2=angle-180, edgecolor=s.color, lw=1))
    M = s.M[-1]
    angle = np.arctan2(M[1], M[0])*180/np.pi
    sub.add_patch(patches.Arc(s.Q, 2*s.w, 2*s.w, theta1=angle+180, theta2=angle, edgecolor=s.color, lw=1))
  
    # plot arcs
    r = 2*s.radius
    for C,cosa, K1, K2, pc in s.arc:
      v2 = K2-C
      cosa = min(max(cosa,-1),1) # can be numerically out of bounds, e.g. with a=0 -> cosa=-1.0000000000000002
      alpha = 180-np.arccos(cosa)*180/np.pi # opening angle between two segments
      beta = np.arctan2(v2[1], v2[0])*180/np.pi # Orientation of arc defined by C->K2

      gamma1 = beta if pc else beta-alpha
      gamma2 = beta+alpha if pc else beta
      
      sub.add_patch(patches.Arc(C, r-2*s.w, r-2*s.w, theta1=gamma1, theta2=gamma2, edgecolor=s.color, lw=1))
      sub.add_patch(patches.Arc(C, r+2*s.w, r+2*s.w, theta1=gamma1, theta2=gamma2, edgecolor=s.color, lw=1))
      if detailed:
        sub.add_patch(patches.Arc(C, r, r, theta1=gamma1, theta2=gamma2, edgecolor=s.color))

  return fig


# write distance values and that stuff
def write_vtk(name,N, detailed, derivative):
  from pyevtk.hl import gridToVTK
  
  shapes = glob.shapes
  
  none_distance = 1.1
  
  x_ls = np.linspace(0,1,N)
  
  dist_segs    = np.ones((N,N,1)) * none_distance # can be None
  dist_circ    = np.zeros((N,N,1)) # defined everywhere
  dist_idx     = np.zeros((N,N,1),dtype=int) # part indices 
  dist_fast_ad = np.zeros((N,N,1)) # using fast_dist_ad based on part idx
  dist         = np.zeros((N,N,1))
  rho          = np.zeros((N,N,1)) # nodal application of boundary function
  rho_ad       = np.zeros((N,N,1)) # use rho_ad() for differentiable nodal density function
  
  total = sum([len(s.var()) for s in shapes])
  assert total == shapes[-1].base + len(shapes[-1].var()) 
  ddist         = np.zeros((total,N,N,1)) # ad for fast_dist_ad
  drho_ad       = np.zeros((total,N,N,1)) # ad for fast_rho_ad 
  
  for s in shapes:
    var = s.var()
    for j, y in enumerate(x_ls):
      for i, x in enumerate(x_ls):
        X = [x,y]
        val, idx = s.dist([x,y])
        dist[i,j,0]     = val
        dist_idx[i,j,0] = idx
        rho[i,j,0]      = poly(val)
        if detailed:
          dist_fast_ad[i,j,0] = s.fast_dist_ad(var,X,idx)
          segs = s.dist(X, 's')[0]
          if(segs): # can be None
            dist_segs[i,j,0] = segs 
          dist_circ[i,j,0] = s.dist(X, 'c')[0]
          rho_ad[i,j,0] = fast_rho_ad(var,X,idx)
        if derivative:
          g = drho(var,X,idx)
          for e in range(len(s.var())):
            drho_ad[s.base+e,i,j,0] = g[e]

        if derivative and detailed:
          g = s.ddist(var,X,idx)
          for e in range(len(s.var())):
            ddist[s.base+e,i,j,0] = g[e]

  pd={"dist": dist}
  pd["dist_idx"]  = dist_idx
  pd["rho"]       = rho
  if detailed:
    pd["fast_dist_ad"]  = dist_fast_ad
    pd["dist_segs"] = dist_segs
    pd["dist_circ"] = dist_circ
    pd["rho_ad"] = rho_ad
  if derivative:
    for s in shapes:
      for e, n in enumerate(s.varnames()):
        pd["d_rho / d_s_" + str(s.id) + '_' + n] = drho_ad[s.base+e]

  if derivative and detailed:
    for s in shapes:
      for e, n in enumerate(s.varnames()):
        pd["d_dist(s_" + str(s.id) + ") / d_" + n] = ddist[s.base+e]


  gridToVTK(name, x_ls,x_ls,np.zeros(1) , pointData=pd)
  print("# wrote '" + name + ".vtr'") 

# plot the distance function on a horizontal line crossing H_1 for the first shape
def lineplot(res):
  s = glob.shapes[0]
  assert len(s.H_int) >= 1 
  y = s.H_int[0][1]
  h = 1./(res-1)
  
  dis = np.ones(res) # dist
  fdisad = np.ones(res) # fast_dist_ad
  idx = np.ones(res,dtype=int)
  ad = np.ones(res)
  var = np.copy(s.var())
  for i, x in enumerate(np.linspace(0,1,res)):
    X = [x,y]
    v,k = s.dist(X) 
    dis[i] = v
    fdisad[i] = s.fast_dist_ad(var,X,k)
    idx[i] = k
    ad[i] = s.ddist(var,X,k)[5]
    
  # finite difference for a  
  org = s.a
  ap = np.copy(s.a)
  epsa = 1e-5
  ap[0] += epsa
  s.set(s.P, s.Q, ap, s.p)  
  dv = s.var()
  dap   = np.ones(res) # using dist
  dapad = np.ones(res) # using fast_dist_ad
  for i, x in enumerate(np.linspace(0,1,res)):
    X = [x,y]
    dap[i] = s.dist(X)[0]
    dapad[i] = s.fast_dist_ad(dv,X,idx[i])
  s.set(s.P, s.Q, org, s.p) # reset   
  
  print('# lineplot of distance for height ',y, ' var=',var)
  print('# set xlabel "x"; set ylabel "distance at y=' + str(y) + '"; set y2label "d distance / dx"; set ytics nomirror; set y2tics nomirror')
  print('# plot "line.dat" u 1:2 w l t "dist", "line.dat" u 1:3 w l axis x1y2 t "d dist / dx", 0')
  print("#(1) x \t(2) dist \t(3) fast_dist_ad \t(4) d dist / dx \t(5) d dist/ da \t(6) d fast_dist_ad / da \t(7)  d dist/da (AD) \t(8) idx ")
  for i, x in enumerate(np.linspace(0,1,res)):
    v = dis[i]
    vf = fdisad[i]
    prev = dis[i-1] if i > 0 else dis[i]
    prfd = dapad[i-1] if i > 0 else dapad[i]
    da = (dap[i]-v) / epsa
    dda = (dapad[i]-v) / epsa
    print(x, v, vf, (v-prev)/h, da, dda, ad[i], idx[i])
    #print(str(x) + ' \t' + str(v) + ' \t' + str((v-prev)/h) + ' \t' + str(idx[i]))






# __name__ is 'spaghetti' if imported or '__main__' if run as commandline tool
if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument("input", help="a .density.xml")
  parser.add_argument("--radius", help="overwrite value from .density.file", type=float)
  parser.add_argument("--set", help="set within a .density.file", type=int)
  parser.add_argument('--save', help="save the image to the given name with the given format. Might be png, pdf, eps, vtp")
  parser.add_argument('--hide', help="suppress technical details for spaghetti", action='store_true')
  parser.add_argument('--rhomin', help="void outside the feature", type=float, default=1e-6)
  parser.add_argument('--transition', help="size of the transition zone (2*h)", type=float, default=.1)
  parser.add_argument('--density', help="write a density.xml to the given filename with density_res")
  parser.add_argument('--density_res', help="resolution for density",type=int, default=60)
  parser.add_argument('--vtk', help="write vtk file for given name (w/o extenstion)")
  parser.add_argument('--vtk_res', help="resolution for vtk export", type=int, default=300)
  parser.add_argument('--vtk_detailed', help="additional vtk output", action='store_true')
  parser.add_argument('--vtk_sens', help="additional sensitvity output via vtk", action='store_true')
  parser.add_argument('--lineplot', help="plots the distance value for the horizontal line crossing H1 in the given res", type=int)
  parser.add_argument('--noshow', help="don't show the image", action='store_true')  

  args = parser.parse_args()
  
  if not os.path.exists(args.input):
    print("error: cannot find '" + args.input + "'")
    os.sys.exit()
  
  shapes = read_xml(args.input, args.set, args.radius)
  
  glob.shapes = shapes
  glob.rhomin = args.rhomin
  glob.transition = args.transition

  if args.lineplot:
    if not args.noshow:
      print("error: use lineplot with --noshow")
      os.sys.exit()
    lineplot(args.lineplot)

  if args.density:
    rho = density(args.density_res)
    ot.write_density_file(args.density,rho)

  if args.vtk:
    write_vtk(args.vtk, args.vtk_res, args.vtk_detailed, args.vtk_sens)

  fig = plot_data(800,shapes,not args.hide)
  if args.save:
    print("write '" + args.save + "'")
    fig.savefig(args.save)
  if not args.noshow:
    fig.show()
    input("Press Enter to terminate.")
  
else:
  # either a manual import shapghetti in python or a call from openCFS via SpaghettiDesign.cc
  # in the later case the cfs_init(), ... functions need to be provided
  
  f = 'line.density.xml'
  X = [0.5,0.6]
  #shapes = read_xml(f)
  #s = shapes[0]
  #print(s.dist([0.1,0.4]))
  #print(s.dist([0.3,0.4]))
  #var = s.var()
  
  #print(s.dist(X),s.ddist(var,X,2))
