#!/usr/bin/env python
# this ia a tool for spaghetti visualization and calculation.
# code segments are based on LuFoVI/code/spaghetti.py by Daniel and Eloise 
import numpy as np
from numpy.linalg import norm 
import os
from pyevtk.hl import gridToVTK
from itertools import product

# for automatic differentiation
import autograd
import autograd.numpy as agnp
from pygments.unistring import combine
from py._builtin import enumerate

# interactive
if __name__ == '__main__':
  import optimization_tools as ot
  import matplotlib
  # necessary for remote execution, even when only saved: http://stackoverflow.com/questions/2801882/generating-a-png-with-matplotlib-when-display-is-undefined
  #matplotlib.use('Agg')
#  matplotlib.use('tkagg')
  import matplotlib.pyplot as plt
  import matplotlib.colors as colors
  import matplotlib.cm as cmx
  from matplotlib.path import Path
  import matplotlib.patches as patches
  import argparse


# to conveniently access global values
class Global:
  # the default value are for easy debugging via import spaghetti
  def __init__(self):
    self.shapes = []         # array of Spaghetti 
    self.rhomin = 1e-30
    self.radius = 0.25       # up to now constant for all spaghetti for arc
    self.boundary = 'poly'   # up to now only 'poly' and 'linear'
    self.transition = 0.05   # paramter for boundary: 2*h
    self.combine= 'max'      # up to now only 'max'
    self.n = [10,10,1]       # [nx, ny, nz]
    self.opts = {} # there is an optional dictionary in xml in the python element
    
    # from opts
    self.silent = False # suppress command line output
    self.order = 6
    
    # these are arrays of size (nx+1,ny+1) with 
    self.idx_field = None # entries are int vectors of size shapes. -2 far outside, -1 far inside, >= 0 nearest part by shape 
    self.dist_field = None # entries are nodal closest distance by shape
    
    # this is an array of size (nx,ny,sum([len(s.var())) with the shape gradients at each element
    self.grad_field = None
    
    # this are fields of field elements of n-size with data to given to cfs via <result>.
    # this are gradient information but also arbitrary data can be added to this dictionary
    # the string key is reported to cfs via cfs_info_field_keys
    self.info_field = {}
    
  # step size as numpy array
  def dx(self):
    return 1/np.array(self.n) # not 1/n+1 !
  
  # total number of variables (not cached)
  def total(self):
    return sum([len(s.var()) for s in self.shapes])  
  
  # give optimzation variables of all spaghetti as defined array such that we can easily differentiate
  # @return p_x,p_y,q_x,q_y,p,a_1,...a_n-1,p_x,p_y,q_x,q_y,p,a_1,...a_n-1,...
  def var_all(self):
    vars = []
    for i in range(len(self.shapes)):
      vars.extend(self.shapes[i].var())
    return vars
    
glob = Global()    
    
## This functions are called from openCFS via SpaghettiDesign.cc

# called from SpaghettiDesign.cc constructor
# @radius a constant in meter, e.g. .15
# @boundary 'poly' or 'linear'
# @transition transition zone which is 2*h
# @combine how to combine? 'max', 'KS' (Kreisselmeier-Steinhauser), 'p-norm'
# @nx, ny, nz are for rho. Currently nx == ny and nz == 1
# @design tupel with design names as strings, usually only 'density'
# @dict dictionary transparently given from the xml file to python
def cfs_init(rhomin, radius, boundary, transition, combine, nx, ny, nz, design, dict):
  # non-zero value avoids divide by 0 in autograd. Seems to also work with 0 though
  glob.rhomin = 1e-30
  glob.radius = radius
  glob.boundary = boundary
  glob.transition = transition
  glob.combine = combine
  assert nx == ny and nz == 1   
  glob.n = [nx, ny, 1]
  glob.design = design
  glob.opts = dict
  
  if not 'silent' in dict:
    raise("expect 'silent' in dict")
  glob.silent = dict['silent'] == '1' 
  if not 'order' in dict:
    raise("expect 'order' in dict")
  glob.order = int(dict['order']) 
  glob.gradvtk = False
  if 'gradvtk' in dict:
    if str(dict['gradvtk']) != 'false' and str(dict['gradvtk']) != '0':
      glob.gradvtk = True
  if 'p' in dict:
    glob.p = float(dict['p'])
  else:
    glob.p = 8
  
  if not glob.silent:
    print('cfs_init designs:',glob.design)
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
    if not glob.silent:
      print('cfs_set_spaghetti: create ', glob.shapes[-1])
  else:
    glob.shapes[id].set(P, Q, a_list, p)
    if not glob.silent:
      print('cfs_set_spaghetti: update ', glob.shapes[id])
  
## give back the density as 1D numpy arry with the current spaghetti setting
def cfs_map_to_design():
  if not glob.silent:
    print('cfs_map_to_design: called for designs',glob.design)
  
  assert glob.n[0] == glob.n[1]
  nx = glob.n[0]
  
  if glob.idx_field is None or glob.dist_field is None:
    glob.idx_field, glob.dist_field, glob.idx_field_shapes_only = create_idx_field() 

  var = glob.var_all()

  rho = np.ones((nx, nx))
  angles = np.ones((nx, nx))
  for j in range(nx):
    for i in range(nx):
      tt = combine_designs(var, i, j, ad=False)
      rho[i, j] = tt[0]
      if len(glob.design) > 1:
        angles[i, j] = tt[1]

  if len(glob.design) > 1:
    if glob.design[0] == 'density':
      return np.concatenate((rho.reshape((np.prod(glob.n)), order='F'),angles.reshape((np.prod(glob.n)), order='F')), axis=0)
    elif glob.design[0] == 'rotAngle':
      return np.concatenate((angles.reshape((np.prod(glob.n)), order='F'),rho.reshape((np.prod(glob.n)), order='F')), axis=0)
    else:
      print('Warning: design designations not implemented!')
      assert(False)
    #return np.array([rho.reshape(np.prod(glob.n), order='F'), angles.reshape(np.prod(glob.n), order='F')])
    #return np.insert(angles.reshape(np.prod(glob.n), order='F'), np.arange(np.prod(glob.n)), rho.reshape(np.prod(glob.n), order='F'))
  else:
    return rho.reshape(np.prod(glob.n), order='F')
  
  
## get the gradient of a density based function w.r.t. all shape variables.
# we calculate it for all and don't skip fixed variables.
# The element wise drho_vec * d_shape can be queried via cfs_get_field_info
# @param drho is a 1D np.array with d func/d rho - this is to be handled via chainrule.
# @param label the name of the current drho_vec's function to allow for cfs_get_field_info  
# @return d func/ d shape  
def cfs_get_gradient(ddes_vec, label):
  if not glob.silent:
    print('cfs_get_gradient',label,ddes_vec) 
  n  = glob.n
  nx = n[0]
  # see density() for comment on layout
  assert len(ddes_vec) == np.prod(n)*len(glob.design)

  # prepare data, reset in Spaghetti.set()
  if glob.idx_field is None or glob.dist_field is None:
    glob.idx_field, glob.dist_field, glob.idx_field_shapes_only = create_idx_field()
  


  var = glob.var_all()
  var = np.array(var)
  
  # prepare gradient, reset in Spaghetti.set()
  if glob.grad_field is None:
    glob.grad_field = np.empty((nx,nx,len(glob.design),glob.total()))
    # d_rho / d_shape
    df =  autograd.jacobian(combine_designs,argnum=0)
    for i in range(nx):
      for j in range(nx):
        for k in range(len(glob.design)):
          tt = df(var,i, j,ad=True) # autograd gives no nice vector
          #print(str(tt))
          #print(tt.shape)
          if glob.design[k] == 'density':
            glob.grad_field[i,j,k] = np.array(tt[0]) # make a handleable vector
          elif glob.design[k] == 'rotAngle':
            glob.grad_field[i,j,k] = np.array(tt[1]) # make a handleable vector
          else:
            print('design ' + str(glob.design[k]) + ' not permitted!')
            assert(False)
          # the components of t are nonzero only where the boundary is gray

  sens = np.zeros(len(var))   # return vector 
  # drho_vec represents a 2D matrix, glob.grad_field is three dimensional with the last dimension the variable
  for snum, s in enumerate(glob.shapes):
    var = glob.shapes[snum].var()
    varnames = glob.shapes[snum].varnames()   
    for i in range(len(var)):
      for k in range(len(glob.design)):
         ds = glob.grad_field[:,:,k,s.base+i].reshape(np.prod(glob.n), order='F') # the shape sensitivity as vector
         # sens_field = element wise the sensitivity of the function times the shape variables.
         # for optmization we need only the sum, we store it here for optional debug output via cfs
         sens_field = ddes_vec[k*np.prod(glob.n):(k+1)*np.prod(glob.n)] * ds   
         # '/' is no valid name for the cfs hdf5 writer but it works with gridToVTK?!
         if k == 0:
           glob.info_field['d_' + label + '_' + glob.design[k] + '_by_d_s' + str(snum) + '_' + varnames[i]] = np.concatenate((sens_field, np.zeros((np.prod(glob.n)))))
         elif k == 1:
           glob.info_field['d_' + label + '_' + glob.design[k] + '_by_d_s' + str(snum) + '_' + varnames[i]] = np.concatenate((np.zeros((np.prod(glob.n))),sens_field))
         sens[s.base+i] += sens_field.sum() # eqals for i in range(nx): for j in range(nx): sens +=  drho[i,j] * glob.grad_field[i,j]

  # make a field
  if glob.gradvtk:
    from pyevtk.hl import gridToVTK
    drho = ddes_vec.reshape((n[0],n[1]),order='F') # see density()
    total = glob.total()
    shapes = glob.shapes
    drho_ad = np.zeros((total,nx,nx,1))
    d_c_d_shapes_loc = np.zeros((total,nx,nx,1))
    drho_ar = np.expand_dims(drho,axis=2)
    pd={"d_c / d_rho": drho_ar}
    for i in range(nx):
      for j in range(nx):
        drho_ad[:,i,j,0] =  glob.grad_field[i,j]
        d_c_d_shapes_loc[:,i,j,0] =  drho[i,j] * glob.grad_field[i,j]
    for s in shapes:
      for e, n in enumerate(s.varnames()):
        pd["d_rho / d_s" + str(s.id) + '_' + n] = drho_ad[s.base+e]
        pd["d_c / d_s" + str(s.id) + '_' + n] = d_c_d_shapes_loc[s.base+e]
    x_ls = np.linspace(0,1,nx+1)
    name = str(glob.n[0]) + 'x' + str(glob.n[1]) + '-' + glob.boundary + '-tr_' + str(glob.transition) + '-rm_' + str(glob.rhomin) + '-order_' + str(glob.order)
    gridToVTK(name, x_ls,x_ls,np.zeros(1) , cellData=pd)
    if not glob.silent:
      print("# wrote '" + name + ".vtr'")

  if not glob.silent:
    print('cfs_get_gradient',label,'->',sens)

  return sens

# as we cannot create a numpy array in C (it should work but fails in reality) we get it here.
# it shall have the size of rho as a 1D array  
def cfs_get_drho_vector():
  
  size = len(glob.design) * np.prod(glob.n)
  assert size >= 1
    
  if not glob.silent:  
    print('cfs_get_drho_vector: returns zero vector of size ', size)

  return np.zeros(size)

# give all keys in glob.info_field as list of strings
def cfs_info_field_keys():
  return list(glob.info_field.keys())
  
# get info filed data by key. Return a numpy.array
def cfs_get_info_field(key):
  if not key in glob.info_field:
    print("requested key '" + key + "' by cfs_get_info_field not valid list ", cfs_info_field_keys)
  return glob.info_field[key]
  
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

    # derivative of fast_dist_ad by automatic differentiation, used for output, not for cfs
    self.grad_fast_dist = autograd.grad(self.fast_dist_ad,argnum=0)

  
  # does not only set the variables but also computes the internal helpers (C, ...) and deletes glob index fields
  def set(self, P,Q,a,p):
    
    glob.idx_field = None
    glob.idx_field_shapes_only = None
    glob.dist_field = None
    glob.grad_field = None
    
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
    # Careful: segments are 0-based here and 1-based in paper!
    self.T = []
    for seg in self.segs:
      self.T.append(seg[1] - seg[0])

    # list of normals M for each segment
    # Careful: normals are 0-based here and 1-based in paper!
    self.M = []
    for T in self.T:
      self.M.append([-T[1], T[0]] / norm(T))

    # For the arcs the center points C^i
    # the idea is to find the angle alpha between two segments, have vector B at half the angle
    # now scale radius/sin(alpha/2) along B0 to find C 
    self.C = []
    self.C.append(self.P) # add P as C^0
    # note there is also list straight segment parts L from P->K_1s, K_1e->K_2s, ...->Q. For one seg this is P->Q
    # L is exclusively used in visualization
    self.L = []
    self.L.append([self.P]) # we will append K1
    # list of arc is only meant as a helper for plotting the arcs
    self.arc = [] 
    r = self.radius
    for i in range(len(self.segs)-1):
      H = self.segs[i][1] # summit point is end of first segment
      v1 = -self.T[i]     # from H to start of first segment 
      v2 = self.T[i+1]    # from H to end of second segment
      v1 = v1 / agnorm(v1)
      v2 = v2 / agnorm(v2)

      # get angle between two segments, standard form for angle between two vectors
      # cos(alpha) = v1*v2/ (||v1|| * ||v2||)
      # If a=0 -> alpha = 0 and numerics goes crazy for the AD case -> see solution there
      cosa = agnp.dot(v1,v2)
      cosa = np.clip(cosa,-1,1)
      alpha = agnp.arccos(cosa)
      if agnorm(v1+v2) > 1e-12:
        scaling = r/agnp.sin(alpha/2) # would be r/sin(0)
        B = v1 + v2
        B0 = B/agnorm(B)
      else:
        B0 = agnp.array([-v1[1], v1[0]])
        scaling = r

      # from function of hypotenuse of a right-angled triangle
      C = H + scaling * B0

      # projection onto segment
      K1 = H+v1*np.dot(v1,C-H)
      K2 = H+v2*np.dot(v2,C-H)
    
      self.C.append(C)
      self.L[-1].append(K1) # finish last (initial) L
      self.L.append([K2]) # to have added K1 for next arc in next iteration or Q after the loop
      
      # only as helper for plotting the stuff.
      self.arc.append((C,cosa, K1, K2))
      
    self.C.append(self.Q) # append Q as C^n
    self.L[-1].append(self.Q)
    assert(len(self.L[-1]) == 2)
    
    assert len(self.C)-2 == len(self.H_int)
    assert len(self.C)-2 == len(self.L)-1 # one segment more than arcs
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
  # @what 'all' for normal use, 'distance' only returns distance and 'index' only returns index 
  def dist(self, X, where = None, what = 'all'):

    w = self.w # profile, we are negative inside and 0 at the boundary
    n = self.n # number of segments, one arc less
    
    minimal = (None,-1)

    # first n segments 0...n-1
    if not where or where == 's':
      for i in range(len(self.T)):
        M = self.M[i]
        C0 = self.C[i] # point on the line perpendicular to the start of the straight segment i (P for first segment)
        C1 = self.C[i+1] # point on the line perpendicular to the end of the straight segment i (Q for last segment)
        t = self.T[i]
        H = self.H[i]
        # we are in range if (X-C^i+1) @ T^i <= 0 and (X-C^i) @ T^i >= 0.
        if (X-C1) @ t <= 0 and (X-C0) @ t >= 0:
          # (X-g) @ M is distance to segment. Positive on M side, negative on the other side. Instead of g, we can use H
          minimal = idx_min(minimal, (abs((X-H) @ M) - w,i)) # don't forget abs!
    
    # arcs
    if not where or where == 'c':
      # n-1 arcs n ... 2n-2
      r = self.radius
      for i, C in enumerate(self.C[1:-1]): # arcs only are around interior C
        # check if we are in the cone (base is the center)
        v1 = self.T[i] 
        v2 = -self.T[i+1]
        XC = X - C
    
        if v1 @ XC >= 0 and v2 @ XC >=0: # we are in the cone
          
          #assert len(self.E) == len(self.C)-2
          #E = self.E[i]
          d = abs(norm(XC) - r) - w
          
        # why is this differentiation here? Seems discontinuous?!
          # if X is on the same side of E as C is and if X is furhter away from E than C is
          #if norm((X-C) + (X-E)) > norm(X-E) and norm(X-E) > norm(X-C):
          #  d = abs(norm(XC)+r) - w
          #else:
          #  d = abs(norm(XC)-r) - w
          
          
          #print('dist cone', X,C,XC,norm(XC),abs(norm(XC) - r),d)
          minimal = idx_min(minimal, (d, n+i))
      
      assert minimal[1] <= 2*n-2
      
      # distance to start point
      minimal = idx_min(minimal,(norm(X-self.P)-w,2*n-1)) 
      # distance to end point
      minimal = idx_min(minimal, (norm(X-self.Q)-w,2*n)) # faster would be root after min() 
    
    if what == 'index':
      return minimal[1]
    elif what == 'distance':
      return minimal[0]
    elif what == 'all':
      return minimal
    else:
      assert(False) 


  # for debug purpose, the part idx human reable printed
  # the parts are 1-based!
  def nice_idx(self, idx):
    n = self.n
    i = idx
    assert i >= -2 and i <= n + n-1 + 2 # n segments, n-1 arcs, P, Q
    if i == -2:
      return "outside"
    if i == -1:
      return "inside"
    if i < n:
      return "seg_" + str(i+1)
    if i < n + n-1:
      return "arc_" + str(i - n + 1)
    if i == n + n-1:
      return "P"
    if i == n + n-1 + 1:
      return "Q"
    assert False  


  # gives the closest distance by part idx - see fast_dist_ad for documentation
  # this version tries to be fast by using cached data
  # @param idx is the param, -1 for far inside and -2 for for outside. 
  # @return +/-1 3*glob.transition if far inside or far outside
  def fast_dist(self,X,idx):

    # this comes from create_idx_field which is -1/-2 if boundary() does not need to be applied    
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
    C = self.C[i+1]
    return abs(norm(X-C) - r) - w

  # gives the closest distance by part idx - see second return value from dist()!
  # does not check if the part is applicable to the coordinate (e.g. if we are for arc within the cone).
  # is meant as basis for automatic differentiation, not fast as it used no precompiled
  # param var the variables which are p_x,p_y,q_x,q_y,a_1,...a_n-1,w -> var()
  # param X the coordinate as array, 
  # param idx  0..n-1 = segment. n...2*n-2=arc, 2n-1=P,2n=Q -> see nice_idx()
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
    
    # in the arc case, we need two segments and such three summits. we now use s(start=old s), c(center=old e), f(final, new)
    assert i >= 0 and i <= n-2
    r = self.radius
    H_c = H_e
    H_f = Q if i == n-2 else P + (i+2)/n * U + a[i+1] * V0
    
    v1 = H_s - H_c
    v2 = H_f - H_c
    v1 = v1 / agnorm(v1)
    v2 = v2 / agnorm(v2)

    # the scaling is based on the angle between the v, If a=0 -> alpha = 0 and numerics goes crazy
    if agnorm(v1+v2) > 1e-12:
      B = v1 + v2
      B0 = B/agnorm(B)
      cosa = agnp.dot(v1,v2)
      assert cosa >= -.9999999999 and cosa <= .9999999999
      #cosa = agnp.clip(cosa,-1,1) 
      alpha = agnp.arccos(cosa)
      scaling = r/agnp.sin(alpha/2)
    else:
      B0 = agnp.array([-v1[1], v1[0]])
      scaling = r # this is the case when the v and B are aligned

    C = H_c + scaling * B0
   
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
   
   
   
# returns the boundary modelling. 
# @param dist positive inside, negative outside, boundary at 0
# @param transition is 2*h in feature mapping paper
# @return if outside |transition| -> rhomin or 1 
def boundary(dist):   
  phi = -dist # positive inside, negative outside
  h = glob.transition/2.0
  rm = glob.rhomin
  if phi <= -h:
    return rm
  if phi >= h:
    return 1.0
  if glob.boundary == 'linear':
    return .5*((1-rm)*phi/h+1+rm)
  elif glob.boundary == 'poly':
    return 3.0/4.0*(1.0 - rm) * (phi/h - phi**3/(3*h**3)) + .5 * (1+rm)
  else:
    print("Error: boundary type '" + glob.boundary + "' not implemented!")
    os.sys.exit()
   
# returns the nodal density value, is ad_differentiable
# the fast is not for performant calculation (it is NOT) but for having the idx given
# @param the variables - make it robust to handle all noodles
# @param X the coordinate
# @param idx the precomouted part index (see Spaghetti.dist()) - make a vector of indices
def fast_rho_ad(var, X, idx):   
  s = glob.shapes[0]
  
  val = s.fast_dist_ad(var,X,idx)
  return boundary(val)
   
def fast_rho(X, idx):   
  s = glob.shapes[0]
  
  val = s.fast_dist(X,idx)
  return boundary(val)


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
# @param discretization is by default glob.n[0]+1     
# @return fields of size (n+1)^dim, to capture the nodes of the density elements. 
#         The first field contains index vectors (len(shapes)), the second a vector of rhos (len(shapes)) 
def create_idx_field(discretization = None):
  
  shapes = glob.shapes
  N = discretization if discretization else glob.n[0]+1
  
  
  h = glob.transition *.55 # is save to add a little buffer
  
  x_ls = np.linspace(0,1,N)
  
  idx = np.ones((N,N,len(shapes)),dtype=int) 
  idx_shapes_only = np.ones((N,N,len(shapes)),dtype=int) 
  dist  = np.ones((N,N,len(shapes)))

  for j, y in enumerate(x_ls):
    for i, x in enumerate(x_ls):
      X = [x,y]
      for si, s in enumerate(shapes):
        d, k = s.dist(X)
        dist[i,j,si] = d
        idx_shapes_only[i,j,si] = k
        idx[i,j,si]  = k if d > -h and d < h else (-1 if d < h else -2)
        # print('cif: i,j,X,d,k',i,j,X,d,k)   

  return idx, dist, idx_shapes_only
      

# integrate a rho for element with indices i and j (however the ordering is)
# uses glob.idx_field and glob.idx_field_shapes_only
# @param if ad fast_dist_ad evaluated, otherwise glob.dist_field is used 
def integrate_rho(var_all, shape, i, j, ad = False):
  shape_num = shape.id
  var = var_all[shape.base:shape.base+len(shape.var())]
  idx_field = glob.idx_field
  order = glob.order
  
  # we take the indices 
  idx1 = idx_field[i,j][shape_num]
  idx2 = idx_field[i+1,j][shape_num]
  idx3 = idx_field[i,j+1][shape_num]
  idx4 = idx_field[i+1,j+1][shape_num]
 
  # we quickly deal with elements inside of or far away from single shapes
  if idx1 == idx2 == idx3 == idx4 == -1:
    return agnp.ones((order*order))
  elif idx1 == idx2 == idx3 == idx4 == -2:
    return glob.rhomin*agnp.ones((order*order))
  
  # we take non-cropped indices as we need them for higher order integration 
  idx_field = glob.idx_field_shapes_only
  idx1 = idx_field[i,j][shape_num]
  idx2 = idx_field[i+1,j][shape_num]
  idx3 = idx_field[i,j+1][shape_num]
  idx4 = idx_field[i+1,j+1][shape_num]
 
  dx = glob.dx()[0]
  deta = dx/(order-1)
  x = i * dx
  y = j * dx
  if order < 1: 
    raise('smallest order is 1 (midpoint')
  XX = [[x+k[0]*deta,y+k[1]*deta] for k in product(range(order),repeat=2)] if order >= 2 else [[x+.5*dx,y+.5*dx]]
  
  if idx1 == idx2 == idx3 == idx4:
    if ad:
      # this case is rho = agnp.array([boundary(shape.fast_dist_ad(var,X,idx1)) for X in XX])
      # but fast_dist_ad is manually inlined in order to save doing costly agnp computations for each integration point
      # design variables
      assert len(var) == 4 + 1 + len(shape.a)

      P = agnp.array(var[0:2])
      Q = agnp.array(var[2:4])
      p = var[4]
      w = p/2.0
      a = var[5:]

      # constants
      n = len(a)+1
      assert n == shape.n

      # endpoints are the most easy ones
      if idx1 == 2*n-1:
        return agnp.array([boundary(agnorm(X-P) - w) for X in XX])
      if idx1 == 2*n:
        return agnp.array([boundary(agnorm(X-Q) - w) for X in XX])

      # the index of our segment or arc
      i = idx1 if idx1 < n else idx1 - n
      assert i >= 0 and i < n
      U = Q - P
      V0 = agnp.array([-U[1],U[0]]) / agnorm(U) # normal to U and normalized
      assert len(shape.a) == n-1

      H_s = P if i == 0   else P + i/n * U + a[i-1] * V0    # summit of begin of segment
      H_e = Q if i >= n-1 else P + (i+1)/n * U + a[i] * V0  # summit of end of segment, which is Q for last segment

      T = H_e - H_s # segment as vector
      M = agnp.array([-T[1], T[0]]) / agnorm(T) # normal

      # segments: abs((X-g) @ M) - w but instad of g we can use H_e
      if idx1 < n:
        return agnp.array([boundary(agnp.abs(agnp.dot((X-H_e), M)) - w) for X in XX])

      # in the arc case, we need two segmens and such three summits. we now use s(start=old s), c(center=old e), f(final, new)
      assert i >= 0 and i <= n-2
      r = shape.radius
      H_c = H_e
      H_f = Q if i == n-2 else P + (i+2)/n * U + a[i+1] * V0

      v1 = H_s - H_c
      v2 = H_f - H_c
      v1 = v1 / agnorm(v1)
      v2 = v2 / agnorm(v2)
  
      # the scaling is based on the angle between the v, If a=0 -> alpha = 0 and numerics goes crazy
      scaling = r # this is the case when the v and B are aligned
      if agnorm(v1+v2) > 1e-12:
        B = v1 + v2
        B0 = B/agnorm(B)
        cosa = agnp.dot(v1,v2)
        assert cosa >= -.9999999999 and cosa <= .9999999999
        #cosa = agnp.clip(cosa,-1,1)
        alpha = agnp.arccos(cosa)
        scaling = r/agnp.sin(alpha/2)
      else:
        B0 = agnp.array([-v1[1], v1[0]])
        scaling = r
  
      C = H_c + scaling * B0

      # arcs: abs(norm(XC) - r) - w
      return agnp.array([boundary(abs(agnorm(X-C) - r) - w) for X in XX])

    else: # no ad-case
      rho = agnp.array([boundary(shape.fast_dist(X,idx1)) for X in XX])
  else: # not all idx same
    #vi = [shape.dist(X) for X in XX]
    #print(vi)
    if ad:
      rho = agnp.array([boundary(shape.fast_dist_ad(var,X,shape.dist(X, None, 'index'))) for X in XX])
    else:
      rho = agnp.array([boundary(shape.dist(X, None, 'distance')) for X in XX])
    
  return rho

# get material rotation angle for element with indices i and j (however the ordering is)
# uses glob.idx_field and glob.idx_field_shapes_only
# @param if ad fast_dist_ad evaluated, otherwise glob.dist_field is used 
def get_material_rotation(var_all, shape_num, i, j, ad):
  shape = glob.shapes[shape_num]
  var = var_all[shape.base:shape.base+len(shape.var())]
  dx = glob.dx()[0]
  X = [(i+0.5) * dx, (j+0.5)*dx]
  d, idx = shape.dist(X)
  
  
  P = agnp.array(var[0:2])
  Q = agnp.array(var[2:4])
  p = var[4]
  w = p/2.0
  a = var[5:]

  # constants
  n = len(a)+1
  assert n == shape.n

  # endpoints are the most easy ones
  # perpendicular to X-P or X-Q
  # vec_vert = X-P => vec = (p1-x1, x2-p2)
  # rotAngle is angle between x-axis and vec
  if idx == 2*n-1:
    XP = P-X
    vec =  agnp.array([XP[1],-XP[0]])
    return -agnp.arctan2(vec[1],vec[0])
  if idx == 2*n:
    XQ = Q-X
    vec = agnp.array([XQ[1],-XQ[0]])
    return -agnp.arctan2(vec[1],vec[0])
  
  # the index of our segment or arc
  i = idx if idx < n else idx - n
  assert i >= 0 and i < n
  U = Q - P
  V0 = agnp.array([-U[1],U[0]]) / agnorm(U) # normal to U and normalized
  
  H_s = P if i == 0   else P + i/n * U + a[i-1] * V0    # summit of begin of segment
  H_e = Q if i >= n-1 else P + (i+1)/n * U + a[i] * V0  # summit of end of segment, which is Q for last segment

  T = H_e - H_s # segment as vector
  
  # segments: parallel to segment
  if idx < n:
    vec = T
    return -agnp.arctan2(vec[1],vec[0])

  # arcs: perpendicular to X-E
  assert i >= 0 and i <= n-2
  M = agnp.array([-T[1], T[0]]) / agnorm(T) # normal
  r = shape.radius
  H_c = H_e
  H_f = Q if i == n-2 else P + (i+2)/n * U + a[i+1] * V0

  v1 = H_s - H_c
  v2 = H_f - H_c
  M1 = M
  M2 = agnp.array([-v2[1],v2[0]]) / agnorm(v2)

  B = (M1 + M2)/2
  B0 = B/agnorm(B)

  nc = -1.0 if agnp.dot(M1, v2) < 0 else 1.0

  # the scaling is based on the angle between the M, If a=0 -> alpha = 0 and numerics goes crazy
  scaling = r # this is the case when the M and B are aligned
  if agnorm(M1-M2) > 1e-10:
    cosa = agnp.dot(v1,v2)/(agnorm(v1) * agnorm(v2))
    assert cosa >= -.9999999999 and cosa <= .9999999999
    #cosa = agnp.clip(cosa,-1,1)
    alpha = agnp.arccos(cosa)
    scaling = r/agnp.sin(alpha/2)

  C = H_c + nc*(scaling * B0)
  XC = C-X
  vec = agnp.array([XC[1],-XC[0]])
  return -agnp.arctan2(vec[1],vec[0])

# compute average material rotation angle weighted by density field, if passed
def combine_angles(angles, density=None):
  # use polar coordinates
  X = agnp.cos(angles)
  Y = agnp.sin(angles)
  # scale averaging importance using density
  if density is not None:
    assert len(density) == len(angles)
    X = agnp.multiply(density, X)
    Y = agnp.multiply(density, Y)
  # get largest vector sum
  # 180° flips are allowed, as they don't change material properties
  # largest vector sum has least cancellations avoidable by 180° flips
  m = 0
  iter = product([1,-1],repeat=len(angles)-1)
  for xx in iter:
    flip = list(xx)
    flip.append((1))
    flip = agnp.array(flip)
    sumvecX = agnp.sum(agnp.multiply(flip,X))
    sumvecY = agnp.sum(agnp.multiply(flip,Y))
    # maximum of squared norm gives same vector as maximum of norm 
    n = sumvecX**2+sumvecY**2
    if n > m:
      # save maximum
      m = n
      vecX = sumvecX
      vecY = sumvecY
  return agnp.arctan2(vecY,vecX)

# integrate rho for element with indices i and j (however the ordering is) for all shapes and take (smooth) maximum
# get rotation angles and take average angle between shapes
# uses glob.idx_field
# @param if ad fast_dist_ad evaluated
def combine_designs(var, i, j, ad):
  num_shapes = len(glob.shapes)
  p = glob.p
  order = glob.order
  # var needs to be passed for autograd to work 
  rho_shapes_ip = agnp.array([integrate_rho(var, s, i, j, ad) for s in glob.shapes])
  #angles_shapes = agnp.array([get_material_rotation(var, s.id, i, j, ad) for s in glob.shapes])
  #print(directions_shapes.shape)
  if glob.combine == 'p-norm':
    rho = agnp.sum(agnp.power(agnp.sum(agnp.power(rho_shapes_ip,p),0),(1/p)))/(order*order)
  elif glob.combine == 'KS':
    rho = agnp.sum(agnp.log(agnp.sum(agnp.exp(p*rho_shapes_ip),0)))/(p*order*order)
  else:
    rho = agnp.sum(agnp.max(rho_shapes_ip,0))/(order*order)
    #rho_ip = agnp.sum(rho_shapes)/(order*order)
  rho_shapes = agnp.sum(rho_shapes_ip,1)/(order*order)
  # for angle average use only non-zero values
  sidx = agnp.nonzero(rho_shapes > 1e-25)
  if len(sidx[0]) != 0:
    angles_shapes = agnp.array([get_material_rotation(var, idx, i, j, ad) for idx in sidx[0]])
    angle = combine_angles(angles_shapes, rho_shapes[sidx[0]])
  else:
    # fallback in void, won't have influence
    angle = 0
  #X2_scaled = agnp.multiply(rho_shapes,agnp.cos(2*angles_shapes))
  #Y2_scaled = agnp.multiply(rho_shapes,agnp.sin(2*angles_shapes))
  #angle = 0.5*agnp.arctan2(agnp.sum(Y2_scaled,axis=0),agnp.sum(X2_scaled,axis=0))
  return agnp.array([rho, angle])



 
# generates a density map for a unit square. 
# this is a trivial implementation, serving for reference which whall be deleted in near future     
def density(nx):
  if len(glob.shapes) != 1: 
    print("Warnung: density(nx) only implemented for first shape")    
  assert nx == glob.n[0]

  s = glob.shapes[0]
  
  rho = cfs_map_to_design()
  
  return rho.reshape((nx,nx),order='F')
#      
#   # the serial element list in cfs is row wise orderd with lower left first and upper right last
#   
#   # when we operate with optimization_tools with arrays of density, they are column wise ordered
#   # with first row in array is left column in image space from bottom  up
#   # wich is first colum in array is lower row in image space from left to right
#   # stupid enough a simple reshape() does not do the arrangement we want in optimization_tools
#   
#   # however, a numpy.reorder(order='F') does the job and we need again Fortran reordering to 
#   # go back from linear to matrix. 
#   #
#   # That's simply because optimization_tools was written without deeper thinking in the beginning 
#   # and orientet itself ot the 99 lines Matlab code (as Matlab has a columnwise ordering)   
#      
#   nodal_val = np.ones((nx+1,nx+1))
#   nodal_idx = np.ones((nx+1,nx+1),dtype=int) * -1
#   
#   for i,x in enumerate(np.linspace(0,1,nx+1)):
#     for j,y in enumerate(np.linspace(0,1,nx+1)): 
#       X = [x,y]
#       val,idx = s.dist(X)
#       nodal_val[i,j] = boundary(val)
#       nodal_idx[i,j] = idx
#       
#   rho = np.ones((nx,nx))     
#   for i in range(0,nx):
#     for j in range(0,nx): 
#       rho[i,j] = .25 * (nodal_val[i,j] + nodal_val[i+1,j] + nodal_val[i,j+1] + nodal_val[i+1,j+1]) 
#   
#   return rho    
      
# reads 2D and returns list of Shaghettis
# @param radius if given overwrites the value from the xml header
def read_xml(filename, set = None, radius = None, cfseval = False):
 
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
    if cfseval:
      cfs_set_spaghetti(idx, Px, Py, Qx, Qy, a, p)
    print('# read noodle', noodle)
      
  return shapes   

    
# creates a matplotlib figure     
# return fig
def plot_data(res, shapes, detail):

  # could respect non-unit regions and out of bounds movement
  minimal = [0,0]
  maximal = [1,1] 
  
  fig, sub = create_figure(res, minimal, maximal)
  
  # for subscripts and superscripts using detail = 2
  SUB = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")
  SUP = str.maketrans("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹")
  
  for s in shapes:
    
    if detail > 0:
      # plot tangent lines with extended summits
      for num, seg in enumerate(s.segs):
        p1 = seg[0]
        p2 = seg[1]
        l = plt.Line2D((p1[0],p2[0]),(p1[1],p2[1]), color=s.color)
        sub.add_line(l)
        if detail > 2:
          p3 = 0.5*(p1+p2)
          v = p2-p1
          n = np.array((-v[0],v[1]))
          p4 = np.array((p3[0]+0.03*n[0], p3[1]-0.03*n[1]+0.01))
          angle = np.arctan2(v[1],v[0])*180/np.pi
          trans_angle = plt.gca().transData.transform_angles(np.array((angle,)),p4.reshape((1, 2)))[0]
          t = plt.text(p4[0], p4[1], '$t^'+str(num+1)+'$', fontsize=16, rotation=angle, rotation_mode='anchor', color=s.color)

      # plot rectangles
      for num, L in enumerate(s.L):
        L1 = L[0]
        L2 = L[1]
        M = s.M[num]
        w = s.w
        
        verts = [L1-w*M, L2-w*M, L2+w*M, L1+w*M, L1-w*M,]
        codes = [Path.MOVETO, Path.LINETO, Path.LINETO, Path.LINETO,Path.CLOSEPOLY,]
        path = Path(verts, codes)
        col = s.color if detail < 2 else [1,1,1]
        alph = .3 if detail < 2 else 1
        filll = True if detail < 2 else False
        patch = patches.PathPatch(path,facecolor=col, edgecolor=s.color, lw=1, alpha=alph, fill=filll)
        sub.add_patch(patch)
        if detail > 1 and num == 0:
          p1 = 0.8*L1+0.2*L2
          sign = np.sign(s.a[0]) if (len(s.a) > 0 and s.a[0] != 0) else 1
          p2 = p1+sign*w*M
          sub.add_line(plt.Line2D((p1[0],p2[0]),(p1[1],p2[1]), color= 'red'))
          plt.annotate('p/2', 0.7*p2+0.3*p1, fontsize=14, xytext=(-5,6), textcoords='offset points', color = 'red')

      # start and endpoint
      fig.gca().add_artist(plt.Circle(s.P, 0.01, color = s.color))
      fig.gca().add_artist(plt.Circle(s.Q, 0.01, color = s.color))
      sub.add_line(plt.Line2D((s.P[0],s.Q[0]),(s.P[1],s.Q[1]), color= 'gray'))
      if detail == 2:
        plt.annotate('$P$', s.P, fontsize=16, xytext=(3,3), textcoords='offset points', color = s.color)
        plt.annotate('$Q$', s.Q, fontsize=16, xytext=(3,3), textcoords='offset points', color = s.color)
      if detail > 2:
        plt.annotate('$P=E^0=H^0=C^0$', s.P, fontsize=16, xytext=(3,3), textcoords='offset points', color = s.color)
        plt.annotate('$Q=E^{m}=H^{m}=C^m$', s.Q, fontsize=16, xytext=(3,3), textcoords='offset points', color = s.color)
        plt.annotate('$u$', 0.5*(s.P+s.Q), fontsize=16, xytext=(1,1), textcoords='offset points', color = 'gray')

      for num, E in enumerate(s.E):
        fig.gca().add_artist(plt.Circle(E, 0.005, color = 'black'))
        if detail > 1:
          H = s.H[num+1]
          sub.add_line(plt.Line2D((E[0],H[0]),(E[1],H[1]), color= 'green'))
          plt.annotate('$a_'+str(num+1)+'$', .5*(H+E), fontsize=16, xytext=(3,-6), textcoords='offset points', color = 'green')
          vec = s.P-s.Q
          gamma2 = 180/np.pi*np.arctan2(vec[1],vec[0])
          if s.a[num] > 0:
            gamma1 = gamma2-90
          else:
            gamma1 = gamma2
            gamma2 = gamma1+90
          sub.add_patch(patches.Arc(E, 0.05, 0.05, theta1=gamma1, theta2=gamma2, edgecolor='green', lw=1))
          angle = np.pi/180*.5*(gamma1+gamma2)
          center = E + 0.0125*np.array((np.cos(angle),np.sin(angle)))
          fig.gca().add_artist(plt.Circle(center, 0.002, color = 'green'))
        if detail > 2:
          plt.annotate('$E^'+str(num+1)+'$', E, fontsize=16, xytext=(1,1), textcoords='offset points', color = 'black')
          if num == 0:
            nvec = vec/norm(vec)
            p1 = E-0.1*nvec
            vec = s.Q-s.P
            nvec = vec/norm(vec)
            normal = np.array((-nvec[1],nvec[0]))
            plt.arrow(p1[0],p1[1], 0.1*nvec[0],0.1*nvec[1], head_width=0.01, color='black')
            leg = 'u0'
            plt.annotate(leg.translate(SUP), p1+.06*nvec, fontsize=16, xytext=(1,0), textcoords='offset points', color = 'black')
            plt.arrow(p1[0],p1[1], 0.1*normal[0],0.1*normal[1], head_width=0.01, color='black')
            leg = 'v0'
            plt.annotate(leg.translate(SUP), p1+.04*normal, fontsize=16, xytext=(-16,1), textcoords='offset points', color = 'black')

      for num, H in enumerate(s.H_int): # the outer H which is P and Q is already in L
        fig.gca().add_artist(plt.Circle(H, 0.005, color = 'red'))
        if detail > 2:
          plt.annotate('$H^'+str(num+1)+'$', H, fontsize=16, xytext=(1,1), textcoords='offset points', color = 'red')
        
      for num, C in enumerate(s.C[1:-1]): # arcs are only around interior C
        fig.gca().add_artist(plt.Circle(C, 0.005, color = 'blue'))
        if detail > 2:
          plt.annotate('$C^'+str(num+1)+'$', C, fontsize=16, xytext=(1,1), textcoords='offset points', color = 'blue')
        if detail > 1:
          L=s.L[num+1][0]
          vec = C-L
          w = s.w
          L = L+w*vec/norm(vec)
          sub.add_line(plt.Line2D((L[0],C[0]),(L[1],C[1]), color='gray', lw=1))
          L=s.L[num][1]
          vec = C-L
          w = s.w
          L = L+w*vec/norm(vec)
          sub.add_line(plt.Line2D((L[0],C[0]),(L[1],C[1]), color='gray', lw=1))
        if detail > 3:
          L=s.L[num][1]
          vec = C-L
          sub.add_line(plt.Line2D((L[0],C[0]),(L[1],C[1]), color='dodgerblue', lw=1, linestyle='-'))
          plt.annotate('r', L+.5*vec, fontsize=16, xytext=(0,0), textcoords='offset points', color = 'dodgerblue')
          H = s.H[num+1]
          sub.add_line(plt.Line2D((L[0],H[0]),(L[1],H[1]), color='dodgerblue', lw=1.5, linestyle='-'))
          sub.add_line(plt.Line2D((H[0],C[0]),(H[1],C[1]), color='dodgerblue', lw=1, linestyle='--'))
          v1 = C-H
          a1 = np.arctan2(v1[1], v1[0])*180/np.pi # Orientation of arc alpha/2
          v2 = s.H[num] - H
          a2 = np.arctan2(v2[1], v2[0])*180/np.pi # Orientation of arc alpha/2  
          alpha2 = np.min((np.abs(a1-a2),360.-np.abs(a1-a2)))

          sub.add_patch(patches.Arc(H, 0.1, 0.1, theta1=a1, theta2=a1+alpha2, edgecolor='dodgerblue', lw=1))
          angle = a1+.5*(alpha2)
          center = H + 0.05*np.array((np.cos(angle),np.sin(angle)))
          leg = r'$\alpha$'+str(num+1)
          plt.annotate(leg.translate(SUP), center, fontsize=16, xytext=(0,0), textcoords='offset points', color = 'dodgerblue')
  
      for L in s.L:  
        fig.gca().add_artist(plt.Circle(L[0], 0.005, color = 'gray'))
        fig.gca().add_artist(plt.Circle(L[1], 0.005, color = 'gray'))
      
      # plot normals
      if detail < 2:
        assert len(s.T) == len(s.M)
        for i, T in enumerate(s.T):
          M = s.M[i]
          p1 = s.H[i] + .5*T
          p2 = p1 + .1*M 
          
          sub.add_line(plt.Line2D((p1[0],p2[0]),(p1[1],p2[1]), color= 'red'))

    # plot half circles for start and end of noodle 
    M = s.M[0] # normal of the first segment tells us where to draw the radius
    angle = np.arctan2(M[1], M[0])*180/np.pi
    sub.add_patch(patches.Arc(s.P, 2*s.w, 2*s.w, theta1=angle, theta2=angle-180, edgecolor=s.color, lw=1))
    M = s.M[-1]
    angle = np.arctan2(M[1], M[0])*180/np.pi
    sub.add_patch(patches.Arc(s.Q, 2*s.w, 2*s.w, theta1=angle+180, theta2=angle, edgecolor=s.color, lw=1))
  
    # plot arcs
    r = 2*s.radius
    for C,cosa, K1, K2 in s.arc:
      v1 = K1-C
      v2 = K2-C
      cosa = min(max(cosa,-1),1) # can be numerically out of bounds, e.g. with a=0 -> cosa=-1.0000000000000002
      alpha = np.arctan2(v1[1], v1[0])*180/np.pi # Orientation of arc defined by C->K2
      beta = np.arctan2(v2[1], v2[0])*180/np.pi # Orientation of arc defined by C->K2

      gamma1 = min(alpha,beta)
      gamma2 = max(alpha,beta)
      
      if detail == 2:
        sub.add_patch(patches.Arc(C, r, r, theta1=0, theta2=360, edgecolor='gray', lw=1))
      sub.add_patch(patches.Arc(C, r-2*s.w, r-2*s.w, theta1=gamma1, theta2=gamma2, edgecolor=s.color, lw=1))
      sub.add_patch(patches.Arc(C, r+2*s.w, r+2*s.w, theta1=gamma1, theta2=gamma2, edgecolor=s.color, lw=1))
      if detail > 0:
        sub.add_patch(patches.Arc(C, r, r, theta1=gamma1, theta2=gamma2, edgecolor=s.color, lw=1.5))

  return fig


# write distance values and that stuff
def write_vtk(name,N, detailed, derivative):
  from pyevtk.hl import gridToVTK
  
  shapes = glob.shapes
  
  none_distance = 1.1
  
  x_ls = np.linspace(0,1,N)
  
  # this is generally not optimized but for understanding and validation
  dist_segs    = np.ones((N,N,1)) * none_distance # can be None
  dist_circ    = np.zeros((N,N,1)) # defined everywhere
  dist_idx     = np.zeros((N,N,1),dtype=int) # part indices 
  dist_fast_ad = np.zeros((N,N,1)) # using fast_dist_ad based on part idx
  dist         = np.zeros((N,N,1))
  rho          = np.zeros((N,N,1)) # nodal application of boundary function
  rho_ad       = np.zeros((N,N,1)) # use rho_ad() for differentiable nodal density function
  
  # this is using optimized code meant for cfs usage
  assert len(glob.shapes) == 1
  field_idx    = create_idx_field(N) # has many shapes!

  field_dist_ad = np.zeros((N,N,1))
  
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
        rho[i,j,0]      = boundary(val)
        if detailed:
          dist_fast_ad[i,j,0]  = s.fast_dist_ad(var,X,idx)
          field_dist_ad[i,j,0] = s.fast_dist_ad(var,X,field_idx[i,j,0])
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
    pd["field_idx"] = field_idx
    pd["field_dist"] = field_dist_ad
  if derivative:
    for s in shapes:
      for e, n in enumerate(s.varnames()):
        pd["d_rho / d_s" + str(s.id) + '_' + n] = drho_ad[s.base+e]

  if derivative and detailed:
    for s in shapes:
      for e, n in enumerate(s.varnames()):
        pd["d_dist(s" + str(s.id) + ") / d_" + n] = ddist[s.base+e]

  gridToVTK(name, x_ls,x_ls,np.zeros(1) , pointData=pd)
  print("# wrote '" + name + ".vtr'") 

# plot the distance function on a horizontal line crossing H_1 for the first shape
def lineplot(res):
  s = glob.shapes[0]
  assert len(s.H_int) >= 1 
  #y = s.H_int[0][1]
  x = .5
  h = 1./(res-1)
  
  dis = np.ones(res) # dist
  fdisad = np.ones(res) # fast_dist_ad
  idx = np.ones(res,dtype=int)
  ad = np.ones(res)
  var = np.copy(s.var())
  for i, y in enumerate(np.linspace(0,1,res)):
    X = [x,y]
    v,k = s.dist(X) 
    dis[i] = v
    fdisad[i] = s.fast_dist_ad(var,X,k)
    idx[i] = k
    t = s.ddist(var,X,k)
    #print(X,s.nice_idx(k),t)
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
  parser.add_argument('--detail', help="level of technical details for spaghetti plot", choices=[0, 1, 2, 3, 4], default=1, type=int)
  parser.add_argument('--rhomin', help="void outside the feature", type=float, default=1e-6)
  parser.add_argument('--transition', help="size of the transition zone (2*h)", type=float, default=.1)
  parser.add_argument('--boundary', help="type of boundary modelling ('poly' or 'linear')", choices=['poly', 'linear'], default='poly')
  parser.add_argument('--combine', help="type of (smooth) maximum function for combination of shapes", choices=['max', 'KS', 'p-norm'], default='max')
  parser.add_argument('--order', help="number of integration points per direction", type=int, default=2)
  parser.add_argument('--density', help="write a density.xml to the given filename with density_res")
  parser.add_argument('--density_res', help="resolution for density",type=int, default=60)
  parser.add_argument('--vtk', help="write vtk file for given name (w/o extenstion)")
  parser.add_argument('--vtk_res', help="resolution for vtk export", type=int, default=200)
  parser.add_argument('--vtk_detailed', help="additional vtk output", action='store_true')
  parser.add_argument('--vtk_sens', help="additional sensitvity output via vtk", action='store_true')
  parser.add_argument('--cfseval', help="dry run cfs calls once, mainly for profiling", action='store_true')
  parser.add_argument('--lineplot', help="plots the distance value for the horizontal line crossing H1 in the given res", type=int)
  parser.add_argument('--noshow', help="don't show the image", action='store_true')  

  args = parser.parse_args()
  
  if not os.path.exists(args.input):
    print("error: cannot find '" + args.input + "'")
    os.sys.exit()
  
  shapes = read_xml(args.input, args.set, args.radius, args.cfseval)
  
  glob.shapes = shapes
  glob.rhomin = args.rhomin
  glob.transition = args.transition
  glob.boundary = args.boundary
  glob.combine = args.combine
  glob.order = args.order
  glob.n = [args.density_res, args.density_res, 1]

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

  if args.cfseval:
    dict = {
      "order": args.order,
      "silent": 1}
    design = ['density', 'rotAngle']
    cfs_init(args.rhomin, args.radius, args.boundary, args.transition, args.combine, args.density_res, args.density_res, 1, design, dict)
    des = cfs_map_to_design()
    #if len(design)>1:
      #des = np.reshape(des, (args.density_res*args.density_res,len(design)), 'C')
    ot.write_multi_design_file(args.input[0:-12] + '.eval.density.xml', des, ['density', 'rotAngle'])
    dummy_drho_vec = np.ones(args.density_res*args.density_res)
    drho = cfs_get_gradient(dummy_drho_vec, 'compliance')

  fig = plot_data(800,shapes,args.detail)
  if args.save:
    print("write '" + args.save + "'")
    fig.savefig(args.save)
  if not args.noshow:
    fig.show()
    input("Press Enter to terminate.")
  
else:
  # either a manual import shapghetti in python or a call from openCFS via SpaghettiDesign.cc
  # in the later case the cfs_init(), ... functions need to be provided
  #import optimization_tools as ot
  f = 'line.density.xml'
  X = [0.5,0.6]
  
  #glob.shapes = read_xml(f)
  #s = glob.shapes[0]

