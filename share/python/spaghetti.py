#!/usr/bin/env python
# this ia a tool for spaghetti visualization and calculation.
# code segments are based on LuFoVI/code/spaghetti.py by Daniel and Eloise 
import numpy as np
from numpy.linalg import norm 
import os
import cfs_utils as ut
from itertools import product
from operator import itemgetter
import copy 

# for automatic differentiation, replaced with normal numpy for testcase to work with cfs runners
#  import autograd
#  import autograd.numpy as agnp
import numpy as agnp

# for gradient check
import scipy.optimize as sciopt
from scipy.optimize import NonlinearConstraint

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
  
  plt.rcParams['text.usetex'] = True


# to conveniently access global values
class Global:
  # the default value are for easy debugging via import spaghetti
  def __init__(self):
    self.shapes = []         # array of Spaghetti 
    self.rhomin = 1e-30
    self.rhomax = 1
    self.radius = 0.25       # up to now constant for all spaghetti for arc
    self.boundary = 'poly'   # up to now only 'poly' and 'linear'
    self.transition = 0.05   # paramter for boundary: 2*h
    self.p = 3
    self.combine= 'max'      # up to now only 'max'
    self.orientation = 'rounded' # or 'straight'
    self.n = [10,10,1]       # [nx, ny, nz]
    self.opts = {} # there is an optional dictionary in xml in the python element
    
    # from opts
    self.silent = False # suppress command line output
    self.order = 6
    
    # these are arrays of size (nx+1,ny+1) with 
    self.idx_field = None # entries are int vectors of size shapes. -2 far outside, -1 far inside, >= 0 nearest part by shape 
    self.dist_field = None # entries are nodal closest distance by shape
    
    # this is an array of size (nx,ny,sum([len(s.optvar())) with the shape gradients at each element
    self.grad_field = None
    
    # this are fields of field elements of n-size with data to given to cfs via <result>.
    # this are gradient information but also arbitrary data can be added to this dictionary
    # the string key is reported to cfs via cfs_info_field_keys
    self.info_field = {}
    
    # this is the sparsity pattern for overlap constraints
    jac = None
  
  # total number of variables (not cached)
  def total(self):
    return sum([len(s.optvar()) for s in self.shapes])  
  
  # give optimzation variables of all spaghetti as defined array such that we can easily differentiate
  # @return p_x,p_y,q_x,q_y,p,a_1,...a_n-1,p_x,p_y,q_x,q_y,p,a_1,...a_n-1,...
  def var_all(self):
    vars = []
    for i in range(len(self.shapes)):
      vars.extend(self.shapes[i].optvar())
    return vars
    
glob = Global()    
    
## This functions are called from openCFS via SpaghettiDesign.cc

# called from SpaghettiDesign.cc constructor or via --cfseval
# @settings dict of key/string from openCFS or from command line
# @design tupel with design names as strings, usually only 'density'
# @dict dictionary transparently given from the xml file to python
def cfs_init(settings, design, dict):
  # non-zero value avoids divide by 0 in autograd. Seems to also work with 0 though
  glob.rhomin = float(settings['rhomin'])
  glob.rhomax = float(settings['rhomax']) if 'rhomax' in settings else 1.0
  if 'radius' in settings: # not set when called via --cfseval
    glob.radius = float(settings['radius'])
  glob.boundary = settings['boundary']
  glob.transition = float(settings['transition'])
  glob.combine =  settings['combine']
  if 'orientation' in settings:
    glob.orientation = settings['orientation']
  glob.n = np.array(eval(settings['n']),dtype=int)
  if 'dx' in settings:
    glob.dx = round(float(settings['dx']),8)
  else:
    glob.dx = 1./glob.n[0]
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
  glob.gradient_check = False
  if 'gradient_check' in dict:
    glob.gradient_check = dict['gradient_check'] == '1'
  
  if not glob.silent:
    print('cfs_init designs:',glob.design)
    print('cfs_init called: ', glob.n, dict) 
  
  
## set spaghetti. Initially create it, otherwise update
# @param id the 0-based index of the spaghetti, corresponds index within glob.shapes
# @param points tuple with two points where a point is list of two doubles 
def cfs_set_spaghetti(id, points, a_list, p):
  assert id < int(len(glob.shapes) + 1) # we may add another spaghetti
  assert len(points) == 2
  assert len(points[0]) == 2
  assert len(points[1]) == 2
  
  P = points[0]  
  Q = points[1]
  
  if id >= len(glob.shapes):
    base = sum([len(s.optvar()) for s in glob.shapes])
    # def __init__(self,id,base,radius,P,Q,a,p):
    glob.shapes.append(Spaghetti(id, base, glob.radius, P, Q, a_list, p))
    if not glob.silent:
      print('cfs_set_spaghetti: create ', glob.shapes[-1])
  else:
    glob.shapes[id].set(P, Q, a_list, p)
    if not glob.silent:
      print('cfs_set_spaghetti: update ', glob.shapes[id])
  
# @res [nx, ny] if given, otherwise glob.n is used
def cfs_map_to_design(res = None):
  if not glob.silent and hasattr(glob, 'design'):
    print('cfs_map_to_design: called for designs',glob.design)
  
  nx = res[0] if res else glob.n[0] 
  ny = res[1] if res else glob.n[1]
  
  if glob.idx_field is None or glob.dist_field is None:
    glob.idx_field, glob.dist_field, glob.idx_field_shapes_only = create_idx_field() 

  var = glob.var_all()

  rho = np.ones((nx, ny))
  angles = np.ones((nx, ny))
  for i in range(nx):
    for j in range(ny):
      tt = combine_designs(var, i, j, derivative=False)
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
  ny = n[1]
  # see density() for comment on layout
  assert len(ddes_vec) == np.prod(n)*len(glob.design)

  # prepare data, reset in Spaghetti.set()
  if glob.idx_field is None or glob.dist_field is None:
    glob.idx_field, glob.dist_field, glob.idx_field_shapes_only = create_idx_field()
  
  # quick check of basic derivatives
  if glob.gradient_check:
    eps = np.sqrt(np.finfo(float).eps)
    print('checking gradient with FD step: ' + str(eps))
    for shape in glob.shapes:
      optvar = shape.optvar()
      which =  ['u', 'nu', 'u0vert']
      for st in ['H', 't', 'tvert', 'alpha', 'C']:
        for ii in range(len(shape.a)):
          if ((st == 't' or st == 'tvert') and ii == len(shape.T)) or st == 'alpha' and (ii == 0 or ii == len(shape.a)-1):
            #print('skip ' + st + str(ii))
            continue
          which.append(st+str(ii))
      for qname in which:
        for xory in [0,1]:
          if qname == 'nu' or qname.startswith('alpha'):
            if xory == 1:
              continue
          fd_grad = sciopt.approx_fprime(optvar, shape.func, eps, [xory, qname])
          grad = shape.grad(optvar, [xory, qname])
          if norm(fd_grad-grad) > 1e-6:
            print(qname + ', variable: ' + str(xory))
            print("FD approximation:  " + str(fd_grad))
            print("Analytic gradient: " + str(grad))
            print("L2 error: " + str(norm(fd_grad-grad)))
    print('\n\n')

  var = glob.var_all()
  var = np.array(var)
  
  # prepare gradient, reset in Spaghetti.set()
  if glob.grad_field is None:
    glob.grad_field = np.empty((nx,ny,len(glob.design),glob.total()))
    # d_rho / d_shape
    #df =  autograd.jacobian(combine_designs,argnum=0)
    for i in range(nx):
      for j in range(ny):
         #tt = df(var,i, j,derivative='autograd') # autograd gives no nice vector
        tt = combine_designs(var,i, j,derivative=True)
        for k in range(len(glob.design)):
          if glob.design[k] == 'density':
            glob.grad_field[i,j,k] = tt[0:len(var)] # make a handleable vector
          elif glob.design[k] == 'rotAngle':
            glob.grad_field[i,j,k] = np.array(tt[len(var):]) # make a handleable vector
          else:
            print('design ' + str(glob.design[k]) + ' not permitted!')
            assert(False)
          # the components of t are nonzero only where the boundary is gray
          if glob.gradient_check:
            fd_grad = sciopt.approx_fprime(var, combine_designs_fd, eps, i, j, glob.design[k])
            if norm(fd_grad-glob.grad_field[i,j,k])>2e-3:# and norm(fd_grad-glob.grad_field[i,j,k])<1e6:
              print(glob.design[k])
              dx = glob.dx
              X = [(i+0.5) * dx, (j+0.5)*dx]
              print("X: " + str(X))
              combine_designs(var,i, j, derivative=False, verbose=True)
              ddes_k = np.reshape(ddes_vec[k*np.prod(glob.n):(k+1)*np.prod(glob.n)], (n[0],n[1]),order='F')
              print("CFS sensitivity: " + str(ddes_k[i,j]))
              print("FD approximation:  " + str(fd_grad))
              print("Analytic gradient: " + str(glob.grad_field[i,j,k]))
              print("combine_designs " + glob.design[k] + " gradient L2 error: " + str(norm(fd_grad-glob.grad_field[i,j,k])) + "\n")
  sens = np.zeros(len(var))   # return vector 
  # drho_vec represents a 2D matrix, glob.grad_field is three dimensional with the last dimension the variable
  for snum, s in enumerate(glob.shapes):
    var = glob.shapes[snum].optvar()
    optvar_names = glob.shapes[snum].optvar_names()   
    for i in range(len(var)):
      for k in range(len(glob.design)):
         ds = glob.grad_field[:,:,k,s.base+i].reshape(np.prod(glob.n), order='F') # the shape sensitivity as vector
         # sens_field = element wise the sensitivity of the function times the shape variables.
         # for optmization we need only the sum, we store it here for optional debug output via cfs
         sens_field = ddes_vec[k*np.prod(glob.n):(k+1)*np.prod(glob.n)] * ds   
         # '/' is no valid name for the cfs hdf5 writer but it works with gridToVTK?!
         if k == 0:
           glob.info_field['d_' + label + '_' + glob.design[k] + '_by_d_s' + str(snum) + '_' + optvar_names[i]] = np.concatenate((sens_field, np.zeros((np.prod(glob.n))))) if len(glob.design) == 2 else sens_field
         elif k == 1:
           glob.info_field['d_' + label + '_' + glob.design[k] + '_by_d_s' + str(snum) + '_' + optvar_names[i]] = np.concatenate((np.zeros((np.prod(glob.n))),sens_field))
         sens[s.base+i] += sens_field.sum() # equals for i in range(nx): for j in range(ny): sens +=  drho[i,j] * glob.grad_field[i,j]

  # make a field
  if glob.gradvtk:
    from pyevtk.hl import gridToVTK
    drho = ddes_vec.reshape((n[0],n[1]),order='F') # see density()
    total = glob.total()
    shapes = glob.shapes
    drho_ad = np.zeros((total,nx,ny,1))
    d_c_d_shapes_loc = np.zeros((total,nx,ny,1))
    drho_ar = np.expand_dims(drho,axis=2)
    pd={"d_c / d_rho": drho_ar}
    for i in range(nx):
      for j in range(ny):
        drho_ad[:,i,j,0] =  glob.grad_field[i,j]
        d_c_d_shapes_loc[:,i,j,0] =  drho[i,j] * glob.grad_field[i,j]
    for s in shapes:
      for e, n in enumerate(s.optvar_names()):
        pd["d_rho / d_s" + str(s.id) + '_' + n] = drho_ad[s.base+e]
        pd["d_c / d_s" + str(s.id) + '_' + n] = d_c_d_shapes_loc[s.base+e]
    x_ls = np.linspace(0,1,nx+1)
    y_ls = np.linspace(0,1,ny+1)
    name = str(glob.n[0]) + 'x' + str(glob.n[1]) + '-' + glob.boundary + '-tr_' + str(glob.transition) + '-rm_' + str(glob.rhomin) + '-order_' + str(glob.order)
    gridToVTK(name, x_ls,x_ls,np.zeros(1) , cellData=pd)
    if not glob.silent:
      print("# wrote '" + name + ".vtr'")

  if not glob.silent:
    print('cfs_get_gradient',label,'->',sens)

  return sens

# get constraint value to prevent arcs moving too close together and overlapping
def cfs_get_constraint_arc_overlap(constraint_num):
  shapenum, a_num, cfs_jac = glob.jac[constraint_num]
  shape = glob.shapes[shapenum]
  return shape.get_constraint_arc_overlap(a_num)

# get constraint gradient
def cfs_get_gradient_arc_overlap(constraint_num):
  shapenum, a_num, cfs_jac = glob.jac[constraint_num]
  shape = glob.shapes[shapenum]
  gradient_full = shape.get_gradient_arc_overlap(a_num)
  return gradient_full[cfs_jac-shape.base]

# get constraint sparsity pattern
# implemented exact sparsity pattern here, might be easier to always just add all 'a_i' variables even though some might be zero
def cfs_get_sparsity_arc_overlap(opt):
  glob.jac = [] # list of arrays of 0-based variable indices
  cfs_jac = []
  for snum, s in enumerate(glob.shapes):
    print(s.namedidx)
    if len(s.idx_a) == 0:
      continue
    base = s.base
    #sparse = np.arange(base,base+s.num_optvar)
    for i in range(1,s.n+1):
      sparse = [base+s.namedidx['px'], base+s.namedidx['py'], base+s.namedidx['qx'], base+s.namedidx['qy']]
      #for j in range(len(s.idx_a)):
      #  sparse.append(base+s.namedidx['a'+str(j+1)])
      #cfs_jac.append(np.array(sparse))
      #glob.jac.append((snum, i, cfs_jac[-1]))
      
      if i == 1 or i == 2: # first two segments, need just one, two or three a's instead of four
        sparse.append(base+s.namedidx['a1'])
        if len(s.idx_a) > 1:
          sparse.append(base+s.namedidx['a2'])
        if i == 2 and len(s.idx_a) > 2:
          sparse.append(base+s.namedidx['a3'])
        cfs_jac.append(np.array(sparse))
        glob.jac.append((snum, i, cfs_jac[-1]))
        if len(s.idx_a) == 1: # only two segments -> overlap symmetric on both sides and only needs one constraint
          break
      elif i == (s.n-1) or i == s.n: # last two segments, need just one, two or three a's instead of four
        if i == (s.n-1) and len(s.idx_a) > 2:
          sparse.append(base+s.namedidx['a'+str(s.n-3)])
        if len(s.idx_a) > 1:
          sparse.append(base+s.namedidx['a'+str(s.n-2)])
        sparse.append(base+s.namedidx['a'+str(s.n-1)])
        cfs_jac.append(np.array(sparse))
        glob.jac.append((snum, i, cfs_jac[-1]))
      else: # need to append four a's
        sparse.append(base+s.namedidx['a'+str(i-2)])
        sparse.append(base+s.namedidx['a'+str(i-1)])
        sparse.append(base+s.namedidx['a'+str(i)])
        sparse.append(base+s.namedidx['a'+str(i+1)])
        cfs_jac.append(np.array(sparse))
        glob.jac.append((snum, i, cfs_jac[-1]))
  print(cfs_jac)
  return cfs_jac

def get_vector_arc_overlap(var_all):
  const = []
  for snum, s in enumerate(glob.shapes):
    var = var_all[s.base:s.base+len(s.optvar())]
    cfs_set_spaghetti(s.id, [[var[0], var[1]], [var[2], var[3]]], var[5:], var[4])
    if len(s.idx_a) == 0:
      continue
    for i in range(1,s.n+1):
      const.append(s.get_constraint_arc_overlap(i))
  return np.array(const)

def get_jacobian_arc_overlap(var_opt):
  c_jac = np.zeros((len(get_vector_arc_overlap(var_opt)),len(var_opt)))
  cnum = 0
  for snum, s in enumerate(glob.shapes):
    if len(s.idx_a) == 0:
      continue
    for i in range(1,s.n+1):
      c_jac[cnum][s.base:s.base+s.num_optvar]=(s.get_gradient_arc_overlap(i))
  return c_jac

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

  dpi_x = res / ((maximal[0] - minimal[0]) * 100.0) 

  fig = matplotlib.pyplot.figure(dpi=100, figsize=(dpi_x*round(max(1,maximal[0])), dpi_x*round(max(1,maximal[1]))))
  ax = fig.add_subplot(111)
  ax.set_xlim(min(0,minimal[0]), max(1,maximal[0]))
  ax.set_ylim(min(0,minimal[1]), max(1,maximal[1]))
  return fig, ax

def dump_shapes(shapes):
  for s in shapes:
     print(s)   

colors = ['b', 'g', 'r', 'c', 'm', 'y', 'tab:orange', 'tab:brown', 'cornflowerblue', 'lime', 'tab:gray']
# transforms ids from 0 to 6 to color codes 'b' to 'k'. Only for matplotlib, not for vtk!
def matplotlib_color_coder(id):
  return colors[id % len(colors)]
    
class Spaghetti: 
  # @param id starting from 0
  # @param base sum of all len(optvar()) of all previous shapes. Starts with 0
  def __init__(self,id,base,radius,P,Q,a,p):
    self.id = id
    self.base = base
    self.color = matplotlib_color_coder(id)
    self.radius = radius
    self.namedidx = {'px': 0,'py': 1,'qx': 2,'qy': 3,'p': 4} # helper to iterate through readable names
    self.idx_a = {} # helper to get index of a in paper notation
    for i in range(len(a)):
      self.namedidx['a' + str(i+1)] = 5+i
      self.idx_a['a' + str(i+1)] = i+1
    self.set(P,Q,a,p)

    # derivative of fast_dist_ad by automatic differentiation, used for output, not for cfs
    #self.grad_fast_dist = autograd.grad(self.fast_dist_ad,argnum=0)

  
  # does not only set the variables but also computes the internal helpers (C, ...) and deletes glob index fields
  def set(self, P,Q,a,p, reset_fields=True):
    
    if reset_fields:
      glob.idx_field = None
      glob.idx_field_shapes_only = None
      glob.dist_field = None
      glob.grad_field = None
    
    self.num_optvar = 5+len(a)
    self.sens = {'grad_u': np.zeros((self.num_optvar,2)),
             'grad_norm_u': np.zeros((self.num_optvar)),
             'grad_u0vert': np.zeros((self.num_optvar,2)),
             'grad_H': np.zeros((len(a)+2, self.num_optvar,2)),
             'grad_t': np.zeros((len(a)+2, self.num_optvar,2)), # extended with empty first gradient to be consistent with paper notation
             'grad_tvert': np.zeros((len(a)+2, self.num_optvar,2)), # extended with empty first gradient to be consistent with paper notation
             'grad_alpha': np.zeros((len(a)+1,self.num_optvar)), # extended with empty first gradient to be consistent with paper notation
             'grad_C': np.zeros((len(a)+2,self.num_optvar,2))} # extended with empty first gradient to be consistent with paper notation
    
    assert len(P) == 2 and len(Q) == 2  
    self.P = np.array(P) # start coordinate
    self.Q = np.array(Q) # end coordinate
    self.a = (0.,)+tuple(a)+(0.,) # vector of normals with 0 added as artificial a for P=H_0 and Q=H_n. Might be just (0, 0)
    self.p = p # profile is full structure thickness
    self.w = p/2. # half structure
    
    # helper
    self.n = len(self.a) - 1 # number of segments
  
    self.U = self.Q - self.P # base line p -> q
    norm_U = norm(self.U)
    if norm_U < 1e-20:
      self.U[:] = [1e-20,0]
      norm_U = norm(self.U)
    self.U0 = self.U / norm_U
    self.sens['grad_u'][0:5,:] = [[-1, 0],
                [0, -1],
                [1, 0],
                [0, 1],
                [0, 0]]
    #for i in len(self.a): #zero initialized anyway
    #  self.sens['grad_u'][5+i]=[0, 0]
    self.sens['grad_norm_u'][0:5] = [-self.U[0]/norm_U, -self.U[1]/norm_U, self.U[0]/norm_U, self.U[1]/norm_U, 0]
    #for i in len(self.a): #zero initialized anyway
    #  self.sens['grad_norm_u'][5+i]=0
  
    self.V = np.array((-self.U[1],self.U[0])) # normal to u
    self.V0 = self.V / norm(self.V)
    nU3 = norm_U**3
    self.sens['grad_u0vert'][0:5,:] = [[-self.U[0]*self.U[1]/nU3, -self.U[1]**2/nU3],
                [self.U[0]**2/nU3, self.U[0]*self.U[1]/nU3],
                [self.U[0]*self.U[1]/nU3, self.U[1]**2/nU3],
                [-self.U[0]**2/nU3, -self.U[0]*self.U[1]/nU3],
                [0, 0]]
    #for i in len(self.a): #zero initialized anyway
    #  self.sens['grad_u0vert'][5+i]=[0, 0]
    
    self.E = [] # endpoints within p-q line but without p and q
    self.H_int = [] # real summits without p and q
    for i in range(1,self.n):
      self.E.append(self.P + i * self.U/self.n)
      self.H_int.append(self.E[-1] + self.a[i] * self.V0)
      self.sens['grad_H'][i][0:5,:] = [[1-i/self.n, 0]+self.a[i]*self.sens['grad_u0vert'][self.namedidx['px']],
                [0, 1-i/self.n]+self.a[i]*self.sens['grad_u0vert'][self.namedidx['py']],
                [i/self.n, 0]+self.a[i]*self.sens['grad_u0vert'][self.namedidx['qx']],
                [0, i/self.n]+self.a[i]*self.sens['grad_u0vert'][self.namedidx['qy']],
                [0, 0]]
      self.sens['grad_H'][i][self.namedidx['a'+str(i)],:] = self.V0 # other entries are zero
    self.sens['grad_H'][0][0:2,:] = [[1, 0], [0,1]] # add derivatives for H_0 = P and H_n = Q
    self.sens['grad_H'][self.n][2:4,:] = [[1, 0], [0,1]]
    assert len(self.E) == (len(self.a)-2)

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

    # list of segment directions as normalized vectors T
    # Careful: segments are 0-based here and 1-based in paper!
    self.T = []
    for i in range(1,self.n+1): # segment index from 1 to n as in paper 
      seg = self.segs[i-1] # get segment corresponding to paper notation!
      svec = seg[1]-seg[0]
      norm_svec = norm(svec)
      self.T.append(svec/norm(svec))
      for namedvar in self.optvar_names():
        idx = self.namedidx[namedvar]
        if (namedvar == 'px') or (namedvar == 'py') or (namedvar == 'qx') or (namedvar == 'qy') or (namedvar == 'p'):
          self.sens['grad_t'][i][idx,:] = ((1/self.n*self.sens['grad_u'][idx]
                                          +(self.a[i]-self.a[i-1])*self.sens['grad_u0vert'][idx])*norm_svec
                                          -svec*norm_U*self.sens['grad_norm_u'][idx]/(self.n**2*norm_svec))/norm_svec**2
        elif (self.idx_a[namedvar] == i-1): # a_{i-1}
          self.sens['grad_t'][i][idx,:] = (-self.V0*norm_svec+svec*(self.a[i]-self.a[i-1])/norm_svec)/norm_svec**2
        elif (self.idx_a[namedvar] == i): # a_i
          self.sens['grad_t'][i][idx,:] = (self.V0*norm_svec-svec*(self.a[i]-self.a[i-1])/norm_svec)/norm_svec**2
        self.sens['grad_tvert'][i][idx,0] = -self.sens['grad_t'][i][idx,1]
        self.sens['grad_tvert'][i][idx,1] = self.sens['grad_t'][i][idx,0]

    # list of normals M for each segment
    # Careful: normals are 0-based here and 1-based in paper!
    self.M = []
    for T in self.T:
      self.M.append(np.array([-T[1], T[0]]))

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
    self.alpha = [0] # artificial alpha_0, should not be used
    r = self.radius
    self.sens['grad_C'][0][0:2,:] = [[1, 0], [0,1]] # add derivatives for C_0 = P and C_n = Q
    self.sens['grad_C'][self.n][2:4,:] = [[1, 0], [0,1]]
    for i in range(1,len(self.T)):
      H = self.H[i] # summit point
      ti = self.T[i-1]     # normalized segment before H (0-based in code!!) 
      tiplus1 = self.T[i]    # normalized segment after H

      # get angle between two segments, standard form for angle between two vectors
      # cos(alpha) = -ti*tiplus1/ (||ti|| * ||tiplus1||)
      # If a=0 -> alpha = 0 and numerics goes crazy for the AD case -> see solution there
      cosa = np.dot(-ti,tiplus1)
      cosa = np.clip(cosa,-1,1)
      alpha = np.arccos(cosa)/2
      self.alpha.append(alpha)
      for namedvar in self.optvar_names():
          idx = self.namedidx[namedvar]
          self.sens['grad_C'][i][idx,:] = self.sens['grad_H'][i][idx,:] # first part of sensitivity independent of alpha
      if agnorm(-ti+tiplus1) > 1e-12:
        # set grad_alpha only here to prevent division by 0
        # grad_alpha only appears in the derivatives of arc segments, which vanish for ti = tiplus1
        # should there still be an integration point in the numerically small remaining area,
        # the spaghetti will be practically straight and dC^i/da^i = dH^i/da^i without the curve should be good enough as approximation
        for namedvar in self.optvar_names():
          idx = self.namedidx[namedvar]
          self.sens['grad_alpha'][i][idx] = (np.dot(ti,self.sens['grad_t'][i+1][idx,:])
                                             +np.dot(tiplus1,self.sens['grad_t'][i][idx,:]))/(2*np.sqrt(1-np.dot(ti,tiplus1)**2))
          # for ti != tiplus1 modify derivative of C^i to include arc information
          self.sens['grad_C'][i][idx,:] = self.sens['grad_C'][i][idx,:] + self.radius*(
            (self.sens['grad_t'][i+1][idx,:]-self.sens['grad_t'][i][idx,:])/(np.sin(alpha)*norm(tiplus1-ti))
            -np.cos(alpha)/(np.sin(alpha)**2)*self.sens['grad_alpha'][i][idx]*(tiplus1-ti)/norm(tiplus1-ti)
            -(tiplus1-ti)/(np.sin(alpha)*norm(tiplus1-ti)**3)*np.dot(tiplus1-ti,self.sens['grad_t'][i+1][idx,:]-self.sens['grad_t'][i][idx,:]))
        scaling = r/agnp.sin(alpha) # would be r/sin(0)
        B = -ti + tiplus1
        B0 = B/agnorm(B)
      else:
        B0 = agnp.array([ti[1], -ti[0]])
        scaling = r

      # from function of hypotenuse of a right-angled triangle
      C = H + scaling * B0
     # for namedvar in self.optvar_names():
     #   idx = self.namedidx[namedvar]
     #   self.sens['grad_C'][i][idx,:] = self.sens['grad_H'][i][idx,:]-

      # projection onto segment
      K1 = H-ti*np.dot(-ti,C-H)
      K2 = H+tiplus1*np.dot(tiplus1,C-H)
    
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
  def optvar(self):
    r = [self.P[0],self.P[1],self.Q[0],self.Q[1],self.p]
    for a in self.a[1:-1]:
      r.append(a)

    return r  

  # helper which gives the name of the variables returned by optvar()
  def optvar_names(self):
    r = ['px','py','qx','qy','p']
    for i in range(len(self.a[1:-1])):
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
        if np.matmul((X-C1), t) <= 0 and np.matmul((X-C0), t) >= 0:
          # (X-g) @ M is distance to segment. Positive on M side, negative on the other side. Instead of g, we can use H
          minimal = idx_min(minimal, (abs(np.matmul((X-H), M)) - w,i)) # don't forget abs!
    
    # arcs
    if not where or where == 'c':
      # n-1 arcs n ... 2n-2
      r = self.radius
      for i, C in enumerate(self.C[1:-1]): # arcs only are around interior C
        # check if we are in the cone (base is the center)
        v1 = self.T[i] 
        v2 = -self.T[i+1]
        XC = X - C
    
        if np.matmul(v1, XC) > 0 and np.matmul(v2, XC) > 0: # we are in the cone
          
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
      return (minimal[0], np.zeros(self.num_optvar)) # add fake gradient to be compatible with fast_dist
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
    # returning 0 gradient is technically wrong but saves time, as implemented boundary functions will lead to 0 gradient anyway 
    if idx == -1:
      return (-3 * glob.transition, np.zeros(self.num_optvar))
    if idx == -2:
      return (3 * glob.transition, np.zeros(self.num_optvar))
    
    # constants
    n = self.n
    w = self.w # p/2

    P = self.P
    Q = self.Q
    
    # endpoints are the most easy ones
    if idx == 2*n-1:
      gradp = (P-X)/norm(X-P) if norm(X-P) > 1e-15 else np.ones(2)
      return (agnorm(X-P) - w, np.concatenate((gradp,np.zeros(2),np.array([-0.5]),np.zeros(self.num_optvar-5))))
    if idx == 2*n:
      gradq = (Q-X)/norm(X-Q) if norm(X-Q) > 1e-15 else np.ones(2)
      return (norm(X-Q) - w, np.concatenate((np.zeros(2),gradq,np.array([-0.5]),np.zeros(self.num_optvar-5))))

    # the index of our segment or arc
    i = idx+1 if idx < n else idx - n + 1 # index for segments is still zero-based but not for gradient
    assert i >= 1 and i <= 2*n
    H = self.H[i]
    M = self.M[i-1] # zero based normals... gotta change this -_-

    # segments: abs((X-g) @ M) - w but instad of g we can use H
    if idx < n:
      grad = []
      for namedvar in self.optvar_names():
        if namedvar == 'p':
          grad.append(-0.5)
        else:
          gradidx = self.namedidx[namedvar]
          sig = 1 if np.sign(np.dot(X-H,M)) >=0 else -1
          grad.append(sig*(np.dot(X-H,self.sens['grad_tvert'][i][gradidx,:])
                                            -np.dot(self.sens['grad_H'][i][gradidx,:],M)))
      return (abs(np.matmul((X-H), M)) - w, np.array(grad))
    
    assert i >= 1 and i <= n-1
    r = self.radius
    C = self.C[i]
    grad = []
    for namedvar in self.optvar_names():
      if namedvar == 'p':
        grad.append(-0.5)
      else:
        gradidx = self.namedidx[namedvar]
        sig = 1 if np.sign(norm(X-C)-r) >= 0 else -1
        grad.append(sig*np.dot(C-X,self.sens['grad_C'][i][gradidx,:])/norm(X-C))
    return (abs(norm(X-C) - r) - w, np.array(grad))

  # for gradient check only
  def func(self, var,args):
    P = agnp.array(var[0:2])
    Q = agnp.array(var[2:4])
    p = var[4]
    a = var[5:]
    self.set(P, Q, a, p, reset_fields=False)
    if str(args[-1]) == 'u0vert':
      return self.V0[args[0]]
    elif str(args[-1]) == 'u':
      return self.U[args[0]]
    elif str(args[-1]) == 'nu':
      return norm(self.U)
    elif str(args[-1]).startswith('H'):
      return self.H[int(args[-1][-1])][args[0]]
    elif str(args[-1]).startswith('t'):
      return self.T[int(args[-1][-1])][args[0]]
    elif str(args[-1]).startswith('tvert'):
      return self.M[int(args[-1][-1])][args[0]]
    elif str(args[-1]).startswith('alpha'):
      return self.alpha[int(args[-1][-1])]
    elif str(args[-1]).startswith('C'):
      return self.C[int(args[-1][-1])][args[0]]
    res = boundary(self.fast_dist(args[0], args[1]), derivative = False)
    #res = self.fast_dist(args[0], args[1])
    #return res[0]
    return res

  # for gradient check only
  def grad(self, var,args):
    P = agnp.array(var[0:2])
    Q = agnp.array(var[2:4])
    p = var[4]
    a = var[5:]
    self.set(P, Q, a, p, reset_fields=False)
    if str(args[-1]) == 'u0vert':
      return self.sens['grad_u0vert'][:,args[0]]
    elif str(args[-1]) == 'u':
      return self.sens['grad_u'][:,args[0]]
    elif str(args[-1]) == 'nu':
      return self.sens['grad_norm_u']
    elif str(args[-1]).startswith('H'):
      return self.sens['grad_H'][int(args[-1][-1])][:,args[0]]
    elif str(args[-1]).startswith('t'):
      return self.sens['grad_t'][int(args[-1][-1])+1][:,args[0]]
    elif str(args[-1]).startswith('tvert'):
      return self.sens['grad_tvert'][int(args[-1][-1])+1][:,args[0]]
    elif str(args[-1]).startswith('alpha'):
      return self.sens['grad_alpha'][int(args[-1][-1])]
    elif str(args[-1]).startswith('C'):
      return self.sens['grad_C'][int(args[-1][-1])][:,args[0]]
    res = boundary(self.fast_dist(args[0], args[1]), derivative = True)
    #res = self.fast_dist(args[0], args[1])
    #print(res)
    return res[1]

  # gives the closest distance by part idx - see second return value from dist()!
  # does not check if the part is applicable to the coordinate (e.g. if we are for arc within the cone).
  # is meant as basis for automatic differentiation, not fast as it used no precompiled
  # param optvar the variables which are p_x,p_y,q_x,q_y,a_1,...a_n-1,w -> optvar()
  # param X the coordinate as array, 
  # param idx  0..n-1 = segment. n...2*n-2=arc, 2n-1=P,2n=Q -> see nice_idx()
  # only the distance itself and a fake 0 gradient to be compatible with fast_dist
  def fast_dist_ad(self,var,X,idx):
    # main changes for autograd
    # - use autograd.numpy (agnp)
    # - any array for computation needs to be agnp.array, e.g. V = [-U[1],U[0]] cannot be differentiated

    if idx == -1:
      return (-3 * glob.transition, np.zeros(self.num_optvar)) # return artificial gradient to conform with fast_dist()
    if idx == -2:
      return (3 * glob.transition, np.zeros(self.num_optvar))

    # design variables
    assert len(var) == 4 + 1 + len(self.a)-2
   
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
      return (agnorm(X-P) - w, np.zeros(self.num_optvar))
    if idx == 2*n:
      return (agnorm(X-Q) - w, np.zeros(self.num_optvar))

    # the index of our segment or arc
    i = idx if idx < n else idx - n
    assert i >= 0 and i < n 
    U = Q - P
    V0 = agnp.array([-U[1],U[0]]) / agnorm(U) # normal to U and normalized
    assert len(self.a) == n+1
 
    H_s = P if i == 0   else P + i/n * U + self.get_a(a, i) * V0    # summit of begin of segment
    H_e = Q if i >= n-1 else P + (i+1)/n * U + self.get_a(a, i+1) * V0  # summit of end of segment, which is Q for last segment
      
    T = H_e - H_s # segment as vector   
    M = agnp.array([-T[1], T[0]]) / agnorm(T) # normal

    # segments: abs((X-g) @ M) - w but instad of g we can use H_e
    if idx < n:
      return (agnp.abs(agnp.dot((X-H_e), M)) - w, np.zeros(self.num_optvar))
    
    # in the arc case, we need two segments and such three summits. we now use s(start=old s), c(center=old e), f(final, new)
    assert i >= 0 and i <= n-2
    r = self.radius
    H_c = H_e
    H_f = Q if i == n-2 else P + (i+2)/n * U + self.get_a(a, i+2) * V0
    
    v1 = H_s - H_c
    v2 = H_f - H_c
    v1 = v1 / agnorm(v1)
    v2 = v2 / agnorm(v2)

    # the scaling is based on the angle between the v, If a=0 -> alpha = 0 and numerics goes crazy
    if agnorm(v1+v2) > 1e-12:
      B = v1 + v2
      B0 = B/agnorm(B)
      cosa = agnp.dot(v1,v2)
      assert cosa >= -.99999999999999 and cosa <= .99999999999999
      #cosa = agnp.clip(cosa,-1,1) 
      alpha = agnp.arccos(cosa)
      scaling = r/agnp.sin(alpha/2)
    else:
      B0 = agnp.array([-v1[1], v1[0]])
      scaling = r # this is the case when the v and B are aligned

    C = H_c + scaling * B0
   
    # arcs: abs(norm(XC) - r) - w
    return (abs(agnorm(X-C) - r) - w, np.zeros(self.num_optvar))
  
  # helper to access a_0 to a_n without extending it to keep autograd's ArrayBoxes working
  def get_a(self,a,idx):
    if idx == 0:
      return 0
    elif idx == self.n:
      return 0
    else:
      return a[idx-1]

  # get constraint value to prevent arcs moving too close together and overlapping
  def get_constraint_arc_overlap(self, segnum):
    assert(len(self.C)>segnum)
    return np.dot(self.H[segnum]-self.H[segnum-1], self.C[segnum]-self.C[segnum-1])

  # get gradient value
  def get_gradient_arc_overlap(self, segnum):
    assert(len(self.a)>segnum)
    return np.dot(self.sens['grad_H'][segnum]-self.sens['grad_H'][segnum-1], self.C[segnum]-self.C[segnum-1]) + \
      np.dot(self.sens['grad_C'][segnum]-self.sens['grad_C'][segnum-1], self.H[segnum]-self.H[segnum-1])
   
  # give gradient via autograd for given X and idx for all parameters
  # @return array of size optvar
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
    return "id=" + str(self.id) + " P=" + str(self.P) + " Q=" + str(self.Q) + " a=" + str(self.a[1:-1]) \
         + " p=" + str(self.p) + " radius=" + str(self.radius) + " color=" + self.color  

  # shape info for a given index
  def to_string(self, idx):
    return "shape=" + str(self.id) + " color=" + str(self.color) + " idx=" + str(idx) + " a=" + str(self.a[idx+1])
   
   
   
# returns the boundary modelling. 
# @param dist positive inside, negative outside, boundary at 0
# @param transition is 2*h in feature mapping paper
# @return if outside |transition| -> rhomin or 1 
def boundary(dist, derivative=False):
  # in the non-derivative case the first (or in the vtk case, only) value used. Whatever dist[1] is otherwise?!
  phi = -dist if not derivative and type(dist) is float or type(dist) is np.float64 else  -dist[0] # positive inside, negative outside
  h = glob.transition/2.0
  #phi = phi-h
  rm = glob.rhomin
  rmx = glob.rhomax
  if phi <= -h:
    rho = rm
  elif phi >= h:
    rho = rmx
  elif glob.boundary == 'linear':
    rho = .5*((rmx-rm)*phi/h+rmx+rm)
  elif glob.boundary == 'poly':
    rho = 3.0/4.0*(rmx - rm) * (phi/h - phi**3/(3*h**3)) + .5 * (rmx+rm)
  else:
    print("Error: boundary type '" + glob.boundary + "' not implemented!")
    os.sys.exit()
  if derivative != True:
    return rho
  else:
    if phi <= -h or phi >= h:
      return (rho, 0*dist[1])
    elif glob.boundary == 'linear':
      return (rho, .5*((rmx-rm)/h) * dist[1])
    elif glob.boundary == 'poly':
      return (rho, -3.0/4.0*(rmx - rm)*(1/h - phi**2/(h**3)) *  dist[1])
   
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
  
  val, grad = s.fast_dist(X,idx)
  return boundary(val)


# differentiate fast_rho_ad w.r.t. to the optvar vector -> gives a vector
#grad_fast_rho = autograd.grad(fast_rho_ad,argnum=0) 
   
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
  Nx = discretization if discretization else glob.n[0]+1
  Ny = discretization if discretization else glob.n[1]+1
  
  
  h = glob.transition *.55 # is save to add a little buffer
  
  x_ls = np.linspace(0,glob.dx*glob.n[0],Nx)
  y_ls = np.linspace(0,glob.dx*glob.n[1],Ny)
  
  idx = np.ones((Nx,Ny,len(shapes)),dtype=int) 
  idx_shapes_only = np.ones((Nx,Ny,len(shapes)),dtype=int) 
  dist  = np.ones((Nx,Ny,len(shapes)))

  for i, x in enumerate(x_ls):
    for j, y in enumerate(y_ls):
      X = [x,y]
      for si, s in enumerate(shapes):
        d, k = s.dist(X)
        dist[i,j,si] = d
        idx_shapes_only[i,j,si] = k
        idx[i,j,si]  = k if d > -h and d < h else (-1 if d < h else -2)
        #idx[i,j,si]  = k if d > -2*h and d < 0.01 else (-1 if d <= -2*h else -2)
        # print('cif: i,j,X,d,k',i,j,X,d,k)   

  return idx, dist, idx_shapes_only
      

# integrate a rho for element with indices i and j (however the ordering is)
# uses glob.idx_field and glob.idx_field_shapes_only
# @param if ad fast_dist_ad evaluated, otherwise glob.dist_field is used 
def integrate_rho(var_all, shape, i, j, derivative = False):
  shape_num = shape.id
  optvar = var_all[shape.base:shape.base+len(shape.optvar())]
  idx_field = glob.idx_field
  order = glob.order
  
  # we take the indices  -2 far outside, -1 far inside
  idx1 = idx_field[i,j][shape_num]
  idx2 = idx_field[i+1,j][shape_num]
  idx3 = idx_field[i,j+1][shape_num]
  idx4 = idx_field[i+1,j+1][shape_num]
 
  # we quickly deal with elements inside of or far away from single shapes
  if idx1 == idx2 == idx3 == idx4 == -1:
    if derivative != True:
      return agnp.ones((order*order)) * glob.rhomax
    else:
      return (agnp.ones((order*order)), np.zeros((order*order,shape.num_optvar)))
  elif idx1 == idx2 == idx3 == idx4 == -2:
    if derivative != True:
      return glob.rhomin*agnp.ones((order*order))
    else:
      return (glob.rhomin*agnp.ones((order*order)), np.zeros((order*order,shape.num_optvar)))
  
  # we take non-cropped indices as we need them for higher order integration 
  idx_field = glob.idx_field_shapes_only
  idx1 = idx_field[i,j][shape_num]
  idx2 = idx_field[i+1,j][shape_num]
  idx3 = idx_field[i,j+1][shape_num]
  idx4 = idx_field[i+1,j+1][shape_num]
 
  dx = glob.dx
  deta = dx/(order-1) if order > 1 else None
  x = i * dx
  y = j * dx
  if order < 1: 
    raise('smallest order is 1 (midpoint')
  XX = [[x+k[0]*deta,y+k[1]*deta] for k in product(range(order),repeat=2)] if order >= 2 else [[x+.5*dx,y+.5*dx]]
  
  if idx1 == idx2 == idx3 == idx4:
    if derivative == False:
      return np.array([boundary(shape.fast_dist(X,idx1)) for X in XX])
    elif derivative == 'autograd' or derivative == 'grad_check': # for gradient check derivative == False cannot be used as it only uses cached values! 
      return agnp.array([boundary(shape.fast_dist_ad(optvar,X,idx1)) for X in XX])
    else: # analytical derivative
      if glob.gradient_check:
        for X in XX:
          err = sciopt.check_grad(shape.func, shape.grad, optvar, [X, idx1])
          if err > 1e-5:
            i = idx1+1 if idx1 < shape.n else idx1 - shape.n + 1 # index for segments is still zero-based but not for gradient
            eps = np.sqrt(np.finfo(float).eps)
            fd_grad = sciopt.approx_fprime(optvar, shape.func, eps, [X, idx1])
            grad = shape.grad(optvar, [X, idx1])
            print("index: " + shape.nice_idx(idx1))
            print("X: " + str(X))
            print("FD approximation:    " + str(fd_grad))
            print("Analytical gradient: " + str(grad))
            print("check_grad L2 error: " + str(err))
      rho, grad_ip = zip(*[boundary(shape.fast_dist(X,idx1),derivative) for X in XX])
      return rho, grad_ip
  else: # not all idx same
    #vi = [shape.dist(X) for X in XX]
    #print(vi)
    if derivative == False:
      return np.array([boundary(shape.dist(X, None, 'distance')) for X in XX])
    elif derivative == 'autograd' or derivative == 'grad_check':
      return agnp.array([boundary(shape.fast_dist_ad(optvar,X,shape.dist(X, None, 'index'))) for X in XX])
    else:
      if glob.gradient_check:
        for X in XX:
          d, idxx = shape.dist(X, None, 'all')
          err = sciopt.check_grad(shape.func, shape.grad, optvar, [X, idxx])
          if err > 1e-5:
            eps = np.sqrt(np.finfo(float).eps)
            fd_grad = sciopt.approx_fprime(optvar, shape.func, eps, [X, idxx])
            grad = shape.grad(optvar, [X, idxx])
            print("Non-unique element index: " + shape.nice_idx(idxx))
            print("X: " + str(X))
            print("FD approximation:    " + str(fd_grad))
            print("Analytical gradient: " + str(grad))
            print("check_grad L2 error: " + str(err))
      rho, grad_ip = zip(*[boundary(shape.fast_dist(X,shape.dist(X, None, 'index')),derivative) for X in XX])
      return rho, grad_ip
    
  return rho

# get material rotation angle for element with indices i and j (however the ordering is)
# uses glob.idx_field and glob.idx_field_shapes_only
# @param if ad fast_dist_ad evaluated, otherwise glob.dist_field is used 
def get_material_rotation(var_all, shape_num, i, j, derivative=False, verbose=False):
  shape = glob.shapes[shape_num]
  var = var_all[shape.base:shape.base+len(shape.optvar())]
  dx = glob.dx
  X = [(i+0.5) * dx, (j+0.5)*dx]
  
  idx1 = glob.idx_field_shapes_only[i,j][shape_num]
  idx2 = glob.idx_field_shapes_only[i+1,j][shape_num]
  idx3 = glob.idx_field_shapes_only[i,j+1][shape_num]
  idx4 = glob.idx_field_shapes_only[i+1,j+1][shape_num]
  if idx1 == idx2 == idx3 == idx4:
    idx = idx1
  else:
    d, idx = shape.dist(X)
  if verbose:
    print("index: " + shape.nice_idx(idx))
  
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
    if glob.orientation == 'straight':
      idx = 0 # use orientation of first segment
    else: # should be 'rounded'
      XP = P-X
      vec =  agnp.array([XP[1],-XP[0]])
      if derivative != True:
        return agnp.arctan2(vec[1],vec[0])
      else:
        if verbose:
          print("vector P-X: ", vec, ' with norm ', nvec2)
          print('derivative: ', datan2)
        nvec2 = vec[0]**2+vec[1]**2
        datan2 = np.array((-vec[1]/nvec2, vec[0]/nvec2)) if nvec2 > 1e-15 else np.zeros((2)) # nondifferentiable if X=Q, but should be inside of spaghetti anyways
        return (np.arctan2(vec[1],vec[0]), np.concatenate(([-datan2[1], datan2[0]],np.zeros(shape.num_optvar-2))))
  if idx == 2*n:
    if glob.orientation == 'straight':
      idx = n-1 # use orientation of last segment
    else: # glob.orientation == 'rounded'
      XQ = Q-X
      vec = agnp.array([XQ[1],-XQ[0]])
      if derivative != True:
        return agnp.arctan2(vec[1],vec[0])
      else:
        nvec2 = vec[0]**2+vec[1]**2
        datan2 = np.array((-vec[1]/nvec2, vec[0]/nvec2)) if nvec2 > 1e-15 else np.zeros((2)) # nondifferentiable if X=Q, but should be inside of spaghetti anyways
        if verbose:
          print("vector Q-X: ", vec, ' with norm ', nvec2)
          print('derivative: ', datan2)
        return (np.arctan2(vec[1],vec[0]), np.concatenate((np.zeros(2),[-datan2[1], datan2[0]],np.zeros(shape.num_optvar-4))))
  
  # the index of our segment or arc
  idx_seg_or_arc = idx+1 if idx < n else idx - n + 1 # index for segments is still zero-based but not for gradient
  assert idx_seg_or_arc >= 1 and idx_seg_or_arc <= n
  U = Q - P
  if norm(U) < 1e-20:
    U = [1e-20, 0]
  V0 = agnp.array([-U[1],U[0]]) / agnorm(U) # normal to U and normalized
  
  H_s = P if idx_seg_or_arc == 1   else P + (idx_seg_or_arc-1)/n * U + shape.get_a(a, idx_seg_or_arc-1) * V0    # summit of begin of segment
  H_e = Q if idx_seg_or_arc >= n else P + idx_seg_or_arc/n * U + shape.get_a(a, idx_seg_or_arc) * V0  # summit of end of segment, which is Q for last segment

  t = H_e - H_s # segment as vector
  
  # segments: parallel to segment
  if idx < n:
    vec = t/norm(t)
    if derivative != True:
      return agnp.arctan2(vec[1],vec[0])
    else:
      nvec2 = vec[0]**2+vec[1]**2
      datan2 = np.array((-vec[1]/nvec2, vec[0]/nvec2)) if nvec2 > 1e-15 else np.zeros((2)) # nondifferentiable if segment has zero length
      return (np.arctan2(vec[1],vec[0]), np.dot(shape.sens['grad_t'][idx_seg_or_arc],datan2))

  # arcs: perpendicular to X-C
  assert idx_seg_or_arc >= 1 and idx_seg_or_arc <= n-1
  M = agnp.array([-t[1], t[0]]) / agnorm(t) # normal
  r = shape.radius
  H_c = H_e
  H_f = Q if idx_seg_or_arc == n-1 else P + (idx_seg_or_arc+1)/n * U + shape.get_a(a, idx_seg_or_arc+1) * V0

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
    assert cosa >= -.99999999999999 and cosa <= .99999999999999
    #cosa = agnp.clip(cosa,-1,1)
    alpha = agnp.arccos(cosa)
    scaling = r/agnp.sin(alpha/2)

  C = H_c + nc*(scaling * B0)
  XC = C-X
  vec = agnp.array([XC[1],-XC[0]])
  if derivative != True:
    return agnp.arctan2(vec[1],vec[0])
  else:
    nvec2 = vec[0]**2+vec[1]**2
    if nvec2 < 1e-15:
      nvec2 = 1e-15 # only happens when X=C, which should lie in void
    return (agnp.arctan2(vec[1],vec[0]), (shape.sens['grad_C'][idx_seg_or_arc][:,1]*XC[0]-shape.sens['grad_C'][idx_seg_or_arc][:,0]*XC[1])/nvec2)

# compute average material rotation angle weighted by density field, if passed
def combine_angles(angles, density=None, derivative=False):
  if derivative == True:
    grad_a = angles[1]
    angles = angles[0]
  # use polar coordinates
  X = agnp.cos(angles)
  Y = agnp.sin(angles)
  # scale averaging importance using density
  if density is not None:
    if derivative == True:
      grad_dens = density[1]
      density = density[0]
    assert len(density) == len(angles)
    X = agnp.multiply(density, X)
    Y = agnp.multiply(density, Y)
  # get largest vector sum
  # 180° flips are allowed, as they don't change material properties
  # largest vector sum has least cancellations avoidable by 180° flips
  m = -1
  iter = product([1,-1],repeat=len(angles)-1)
  for xx in iter:
    flip = list(xx)
    flip.append((1))
    flip = agnp.array(flip)
    flipped_X = agnp.multiply(flip,X)
    flipped_Y = agnp.multiply(flip,Y)
    sumvecX = agnp.sum(flipped_X)
    sumvecY = agnp.sum(flipped_Y)
    # maximum of squared norm gives same vector as maximum of norm 
    n = sumvecX**2+sumvecY**2
    if n > m:
      # save maximum
      m = n
      flip_max = flip
      flipped_max_X = flipped_X
      flipped_max_Y = flipped_Y
      vec_max_sum_X = sumvecX
      vec_max_sum_Y = sumvecY
  if derivative != True:
    return agnp.arctan2(vec_max_sum_Y,vec_max_sum_X)
  else:
    # derivative of [vec_max_sum_X, vec_max_sum_Y] is [-flipped_max_X[s], flipped_max_Y[s]], of atan2 [-vec_max_sum_Y/norm(vec), vec_max_sum_X/norm(vec)]
    nvec2 = vec_max_sum_X**2+vec_max_sum_Y**2
    datan2 = [-vec_max_sum_Y/nvec2, vec_max_sum_X/nvec2]
    grad_max = []
    for s in range(len(angles)):
      grad_s = np.zeros(grad_a[s].shape)
      grad_s = np.multiply((datan2[0]*(-flipped_max_Y[s])+datan2[1]*flipped_max_X[s]),grad_a[s])
      if density is not None:
        grad_s = grad_s + np.multiply(((datan2[0]*flip_max[s]*np.cos(angles[s])+datan2[1]*flip_max[s]*np.sin(angles[s]))),grad_dens[s])
      grad_max.append(grad_s)
    return grad_max

# only changing argument order for gradient check of combine_angles w.r.t. density variable
def combine_angles_fd_density(density, angles):
  return combine_angles(angles, density, derivative=False)

# integrate rho for element with indices i and j (however the ordering is) for all shapes and take (smooth) maximum
# get rotation angles and take average angle between shapes
# uses glob.idx_field
# @param if ad fast_dist_ad evaluated
def combine_designs(var, i, j, derivative, verbose=False):
  num_shapes = len(glob.shapes)
  p = glob.p
  order = glob.order
  # optvar needs to be passed for autograd to work 
  if derivative != True:
    rho_shapes_ip = agnp.array([integrate_rho(var, s, i, j, derivative) for s in glob.shapes])
  else:
    rho_shapes_ip = []
    grad_shapes_ip = []
    for s in glob.shapes:
      rr, gg = integrate_rho(var, s, i, j, derivative)
      rho_shapes_ip.append(rr)
      grad_shapes_ip.append(gg)
    rho_shapes_ip = np.array(rho_shapes_ip)

  if glob.combine == 'p-norm':
    rho_shapes_ip_p = agnp.power(rho_shapes_ip,p)
    rho_shape_weights = agnp.sum(rho_shapes_ip_p,1)/(order*order)
    rho = agnp.sum(agnp.power(agnp.sum(rho_shapes_ip_p,0),(1/p)))/(order*order)
  elif glob.combine == 'softmax':
    exp = agnp.exp(p*rho_shapes_ip)
    sum_exp = agnp.sum(exp,0)
    rho_shape_weights = rho_shapes_ip*exp
    rho_ip = agnp.sum(rho_shape_weights,0)/sum_exp
    rho_shape_weights = agnp.sum(rho_shape_weights,1)/(order*order)
    rho = agnp.sum(rho_ip)/(order*order)
  elif glob.combine == 'KS':
    exp = agnp.exp(p*rho_shapes_ip)
    sum_exp = agnp.sum(exp,0)
    rho_shape_weights = agnp.sum(exp,1)/(order*order)-1
    rho = agnp.sum(agnp.log(sum_exp))/(p*order*order)
  else:
    rho = agnp.sum(agnp.max(rho_shapes_ip,0))/(order*order)
    #rho_ip = agnp.sum(rho_shapes)/(order*order)
  if derivative == True:
    grad_dens = []
    grad_rho_shape_weights = []
    if glob.combine == 'p-norm':
      sum_s = np.sum(rho_shapes_ip_p,0)
      for ip in range(len(sum_s)):
        sum_s[ip] = np.power(sum_s[ip],1/p-1) if sum_s[ip] > 1e-30 else 0 # avoid non-differentiability of norm as gradient in void is zero anyways
      rho_s_p = np.power(rho_shapes_ip, p-1)
      tmp = np.expand_dims(sum_s,axis=0)*rho_s_p
      for s in range(len(glob.shapes)):
        grad_s = np.expand_dims(tmp[s],axis=1)*grad_shapes_ip[s]
        grad_dens.append(np.sum(grad_s,axis=0)/(order*order))
        grad_rho_shape_weights.append(p*np.sum(np.expand_dims(rho_s_p[s],axis=1)*grad_shapes_ip[s],0)/(order*order))
    elif glob.combine == 'softmax':
      for s in range(len(glob.shapes)):
        grad_rho_shape_weights.append(np.sum(np.expand_dims(exp[s]+p*exp[s]*rho_shapes_ip[s],1)*grad_shapes_ip[s],0)/(order*order))
        grad_s = np.expand_dims(exp[s]/sum_exp*(1+p*(rho_shapes_ip[s]-rho_ip)),axis=1)*grad_shapes_ip[s]
        grad_dens.append(np.sum(grad_s,axis=0)/(order*order))
    elif glob.combine == 'KS':
      for s in range(len(glob.shapes)):
        grad_s = np.expand_dims(np.exp(p*rho_shapes_ip[s])/sum_exp,axis=1)*grad_shapes_ip[s]
        grad_dens.append(np.sum(grad_s,axis=0)/(order*order))
        grad_rho_shape_weights.append(np.sum(np.expand_dims(p*exp[s],axis=1)*grad_shapes_ip[s],0)/(order*order))
    elif glob.combine == 'max':
      rho_argmax = np.argmax(rho_shapes_ip,axis=0)
      for s in range(len(glob.shapes)):
        grad_s = np.zeros((len(grad_shapes_ip[s][0])))
        for ip in np.nonzero(rho_argmax==s)[0]:
          grad_s += grad_shapes_ip[s][ip]
        grad_dens.append(grad_s/(order*order))
    else:
      print("gradient for feature combination ", glob.combine, " not implemented!")
      assert(False)

  # angle computation
  angle = 0 # dummy result for void
  grad_a = np.zeros((len(var)))
  if len(glob.design) > 1:
    if glob.combine == 'max':
      rho_shape_weights = agnp.sum(rho_shapes_ip,1)/(order*order) # weights for angle average
    # for angle average use only non-zero values
    sidx = agnp.nonzero(rho_shape_weights > 1e-14)
    if len(sidx[0]) != 0:
      angles_shapes = agnp.array([get_material_rotation(var, idx, i, j, derivative) for idx in sidx[0]])
      if derivative == True:
        ang, grad_ang = zip(*angles_shapes)
        if verbose:
          print("gradient before combine_angles: ", grad_ang)
          [get_material_rotation(var, idx, i, j, derivative, verbose=True) for idx in sidx[0]]
        if glob.combine == 'max': # just take argmax close to center. Other "more precise" implementation possible as well
          rho_argmax = np.argmax(rho_shapes_ip[sidx[0]][0][round(0.5*(order-1)*(order+1))])
          shape = glob.shapes[sidx[0][rho_argmax]]
          grad_a[shape.base:shape.base+len(shape.optvar())] = grad_ang[rho_argmax]
          return np.concatenate((np.concatenate(grad_dens), grad_a))
        eps = np.sqrt(np.finfo(float).eps)
        if glob.gradient_check:
          for ind,s in enumerate(sidx[0]):
            fd_grad = sciopt.approx_fprime(var, get_material_rotation, eps, s, i, j)
            shape = glob.shapes[s]
            fd_grad = fd_grad[shape.base:shape.base+len(shape.optvar())]
            err = norm(fd_grad-grad_ang[ind])
            if err>3e-5:
              dx = glob.dx
              X = [(i+0.5) * dx, (j+0.5)*dx]
              d, idx = shape.dist(X)
              print("index: " + shape.nice_idx(idx))
              print("X: " + str(X))
              print("FD approximation:    " + str(fd_grad))
              print("Analytical gradient: " + str(grad_ang))
              print("get_material_rotation gradient L2 error: " + str(err))
        if glob.gradient_check:
          grad_max_ang = combine_angles((ang,np.ones(len(sidx[0]))), (rho_shape_weights[sidx[0]],np.zeros(len(sidx[0]))), derivative)
          fd_grad_ang = sciopt.approx_fprime(ang, combine_angles, eps, rho_shape_weights[sidx[0]])
          err = norm(fd_grad_ang-grad_max_ang)
          if err>3e-5:
            dx = glob.dx
            X = [(i+0.5) * dx, (j+0.5)*dx]
            print("X: " + str(X))
            print("FD approximation:    " + str(fd_grad_ang))
            print("Analytical gradient: " + str(grad_max_ang))
            print("combine_angles angle gradient L2 error: " + str(err))
          grad_max_dens = combine_angles((ang,np.zeros(len(sidx[0]))), (rho_shape_weights[sidx[0]],np.ones(len(sidx[0]))), derivative)
          fd_grad_dens = sciopt.approx_fprime(rho_shape_weights[sidx[0]], combine_angles_fd_density, eps, ang)
          err = norm(fd_grad_dens-grad_max_dens)
          if err>3e-5:
            dx = glob.dx
            X = [(i+0.5) * dx, (j+0.5)*dx]
            print("X: " + str(X))
            print("FD approximation:    " + str(fd_grad_dens))
            print("Analytical gradient: " + str(grad_max_dens))
            print("combine_angles density gradient L2 error: " + str(err))
        grad_dens_sidx = []
        for s in sidx[0]:
          grad_dens_sidx.append(grad_rho_shape_weights[s])
        grad_max = combine_angles((ang,grad_ang), (rho_shape_weights[sidx[0]],grad_dens_sidx), derivative)
        #grad_max = combine_angles((ang,grad_ang), None, derivative)
        if verbose:
          print("gradient after combine_angles: ", grad_max)
        for ind,s in enumerate(sidx[0]):
          shape = glob.shapes[s]
          grad_a[shape.base:shape.base+len(shape.optvar())] = grad_max[ind]
      else:
        if glob.combine == 'max': # just take argmax close to center. Other "more precise" implementation possible as well
          rho_argmax = np.argmax(rho_shapes_ip[sidx[0]][0][round(0.5*(order-1)*(order+1))])
          angle = angles_shapes[rho_argmax]
        else:
          angle = combine_angles(angles_shapes, rho_shape_weights[sidx[0]])
          #angle = combine_angles(angles_shapes, None)
          if verbose:
            print("active shapes " + str(sidx[0]) + " with densities " + str(rho_shape_weights[sidx[0]]))
            print("angles in degrees: " + str(180*angles_shapes/np.pi))
            print("weighted angles in degrees: " + str(180*(angles_shapes*rho_shape_weights[sidx[0]])/(np.max(rho_shape_weights[sidx[0]])*np.pi)) )
            print("averaged angle in degrees: " + str(180*angle/np.pi))
    # else: void, set above

  #X2_scaled = agnp.multiply(rho_shapes,agnp.cos(2*angles_shapes))
  #Y2_scaled = agnp.multiply(rho_shapes,agnp.sin(2*angles_shapes))
  #angle = 0.5*agnp.arctan2(agnp.sum(Y2_scaled,axis=0),agnp.sum(X2_scaled,axis=0))
  if derivative != True:
    return agnp.array([rho, angle])
  else:
    return np.concatenate((np.concatenate(grad_dens), grad_a))

# helper for gradient check
def combine_designs_fd(var,i, j, which = 'rotAngle'):
  res = combine_designs(var,i, j,derivative='grad_check')
  return res[1] if which == 'rotAngle' else res[0]
 
# generates a density map for a quare. 
# this is a trivial implementation, serving for reference which whall be deleted in near future     
def density(size):
  assert len(size) == 2 # [nx, ny]
  rho = cfs_map_to_design(size)
  
  return rho.reshape((nx,nx),order='F')
  return rho.reshape(size,order='F')

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
      
# reads 2D and returns list of Spaghetti and domain. Also sets some values to glob!
# @param radius if given overwrites the value from the xml header
# @return list of spaghettis and either [[min_x, min_y], [max_x, max_y]] if in density.xml or None
def read_xml(filename, set = None, radius = None, cfseval = False):
 
  xml = ot.open_xml(filename)
 
  shapes = []
  sq = 'last()' if set == None else '@id="' + str(set) + '"'

  if not radius:
    glob.radius = float(ot.xpath(xml, '//header/spaghetti/@radius')) 
    radius = glob.radius

  glob.n[0] = int(ot.xpath(xml, '//header/mesh/@x'))
  glob.n[1] = int(ot.xpath(xml, '//header/mesh/@y'))
  glob.n[2] = int(ot.xpath(xml, '//header/mesh/@z'))

  assert(glob.n[2] == 1)

  domain = None # this feature is only written by modern cfs (solar_heater, 12.2022) but we also want to read old files
  pn = xml.xpath('//header/coordinateSystems/domain')
  if len(pn) == 1: # we assume a single coordinate system
    dic = pn[0].attrib
    domain = [[float(dic['min_x']),float(dic['min_y'])],[float(dic['max_x']),float(dic['max_y'])]]

  glob.design = ['density']
  if len(xml.xpath('//element[@type="rotAngle"]')) > 0:
    glob.design.append('rotAngle')

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
    # width of noodle is 2*w -> don't confuse with P
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

    base = sum([len(s.optvar()) for s in shapes])
    noodle = Spaghetti(id=idx, base=base, radius=radius, P=(Px,Py), Q=(Qx,Qy), a=a, p=p)
    shapes.append(noodle)
    if cfseval:
      cfs_set_spaghetti(idx, [[Px, Py], [Qx, Qy]], a, p)
    print('# read noodle', noodle)
      
  return shapes, domain

# reads 2D and returns list of Spaghetti
# @param radius if given overwrites the value from the xml header
def write_xml(filename, shapes, remove_ghosts=False, padnormals=0):

  out = open(filename, "w")
  out.write('<?xml version="1.0"?>\n')
  out.write('<cfsErsatzMaterial>\n')
  out.write('  <header>\n')
  out.write('    <mesh x="' + str(glob.n[0]) + '" y="' + str(glob.n[1]) + '" z="' + str(glob.n[2]) + '"/>\n')
  out.write('    <spaghetti radius="' + str(glob.radius) + '"/>\n')
  out.write('  </header>\n')
  out.write('  <set id="spaghetti.py">\n')

  nr = 0
  shapeid = 0
  for shape in shapes:
    if remove_ghosts and shape.p < 0.0005:
      continue
    out.write('    <shapeParamElement nr="' + str(nr) + '" type="node" dof="x" tip="start" shape="' + str(shapeid) + '" design="' + str(shape.P[0]) + '"/>\n')
    out.write('    <shapeParamElement nr="' + str(nr+1) + '" type="node" dof="y" tip="start" shape="' + str(shapeid) + '" design="' + str(shape.P[1]) + '"/>\n')
    out.write('    <shapeParamElement nr="' + str(nr+2) + '" type="node" dof="x" tip="end" shape="' + str(shapeid) + '" design="' + str(shape.Q[0]) + '"/>\n')
    out.write('    <shapeParamElement nr="' + str(nr+3) + '" type="node" dof="y" tip="end" shape="' + str(shapeid) + '" design="' + str(shape.Q[1]) + '"/>\n')
    out.write('    <shapeParamElement nr="' + str(nr+4) + '" type="profile" shape="' + str(shapeid) + '" design="' + str(shape.p) + '"/>\n')
    nr += 5
    for idx, normal in enumerate(shape.a[1:-1]):
      out.write('    <shapeParamElement nr="' + str(nr) + '" type="normal" shape="' + str(shapeid) + '" design="' + str(shape.a[idx+1]) + '"/>\n')
      nr += 1
    for i in range(len(shape.a)-2, padnormals):
      out.write('    <shapeParamElement nr="' + str(nr) + '" type="normal" shape="' + str(shapeid) + '" design="0.000000001"/>\n')
      nr += 1
    shapeid += 1


  out.write('  </set>\n')
  out.write(' </cfsErsatzMaterial>\n')
  out.close()


    
# creates a matplotlib figure
# @domain is [[min_x, min_y],[max_x,max_y]] which is in modern density.xml     
# return fig
def plot_data(res, shapes, detail, domain = None):

  if domain:
    assert len(domain) == 2
    minimal = domain[0]
    maximal = domain[1]
  else:
    # could respect non-unit regions and out of bounds movement
    minimal = [0,0]
    min_dim = min((glob.n[0],glob.n[1]))
    maximal = [glob.n[0]/min_dim,glob.n[1]/min_dim] # normalize smaller dimension to 1, as there is no other element information in .density.xml

  lineopacity = 0.5 if args.gray else 1 # opacity value for plotting lines and points    

  fig, sub = create_figure(res, minimal, maximal)
  
  for s in shapes:
    
    if s.p < 1e-10: # omit zero-width shapes
      continue

    if detail > 0:
      # plot tangent lines with extended summits
      for num, seg in enumerate(s.segs):
        p1 = seg[0]
        p2 = seg[1]
        l = plt.Line2D((p1[0],p2[0]),(p1[1],p2[1]), alpha=lineopacity, color=s.color)
        sub.add_line(l)
        if detail > 2:
          p3 = 0.5*(p1+p2)
          v = p2-p1
          n = np.array((-v[0],v[1]))
          p4 = np.array((p3[0]+0.03*n[0], p3[1]-0.03*n[1]+0.01))
          angle = np.arctan2(v[1],v[0])*180/np.pi
          trans_angle = plt.gca().transData.transform_angles(np.array((angle,)),p4.reshape((1, 2)))[0]
          t = plt.text(p4[0], p4[1], r'$t^'+str(num+1)+'$', fontsize=16, rotation=angle, rotation_mode='anchor', alpha=lineopacity, color=s.color)

      # plot rectangles
      for num, L in enumerate(s.L):
        L1 = L[0]
        L2 = L[1]
        M = s.M[num]
        w = s.w
        
        verts = [L1-w*M, L2-w*M, L2+w*M, L1+w*M, L1-w*M,]
        codes = [Path.MOVETO, Path.LINETO, Path.LINETO, Path.LINETO,Path.CLOSEPOLY,]
        path = Path(verts, codes)
        col = s.color #if detail < 2 else [1,1,1]
        alph = .3 #if detail < 2 else 1
        filll = True #if detail < 2 else False
        patch = patches.PathPatch(path,facecolor=col, edgecolor=s.color, lw=1, alpha=alph, fill=filll)
        sub.add_patch(patch)
        patch2 = patches.PathPatch(path,facecolor=col, edgecolor=s.color, alpha=lineopacity, lw=1, fill=False)
        sub.add_patch(patch2)
        if detail > 1 and num == 0:
          p1 = 0.8*L1+0.2*L2
          sign = np.sign(s.a[1]) if (len(s.a) > 2 and s.a[1] != 0) else 1
          p2 = p1+sign*w*M
          sub.add_line(plt.Line2D((p1[0],p2[0]),(p1[1],p2[1]), alpha=lineopacity, color= 'red'))
          plt.annotate('p/2', 0.7*p2+0.3*p1, fontsize=20, xytext=(-5,6), textcoords='offset points', alpha=lineopacity, color = 'red')

      # start and endpoint
      fig.gca().add_artist(plt.Circle(s.P, 0.01, alpha=lineopacity, color = s.color))
      fig.gca().add_artist(plt.Circle(s.Q, 0.01, alpha=lineopacity, color = s.color))
      #sub.add_line(plt.Line2D((s.P[0],s.Q[0]),(s.P[1],s.Q[1]), alpha=lineopacity, color= 'gray'))
      if detail == 2:
        plt.annotate('$P^0$', s.P, fontsize=26, xytext=(3,3), textcoords='offset points', alpha=lineopacity, color = s.color)
        plt.annotate('$P^'+str(s.n)+'$', s.Q, fontsize=26, xytext=(3,3), textcoords='offset points', alpha=lineopacity, color = s.color)
      if detail > 2:
        plt.annotate('$P^0=C^0$', s.P, fontsize=26, xytext=(3,3), textcoords='offset points', alpha=lineopacity, color = s.color)
        plt.annotate('$P^'+str(s.n)+'=C^'+str(s.n)+'$', s.Q, fontsize=26, xytext=(3,3), textcoords='offset points', alpha=lineopacity, color = s.color)
        plt.annotate('$u$', 0.5*(s.P+s.Q), fontsize=16, xytext=(1,1), textcoords='offset points', alpha=lineopacity, color = 'gray')

      for num, E in enumerate(s.E):
        #fig.gca().add_artist(plt.Circle(E, 0.005, alpha=lineopacity, color = 'gray'))
        if detail > 2:
          H = s.H[num+1]
          sub.add_line(plt.Line2D((E[0],H[0]),(E[1],H[1]), alpha=lineopacity, color= 'green'))
          plt.annotate('$a_'+str(num+1)+'$', .5*(H+E), fontsize=16, xytext=(3,-6), textcoords='offset points', alpha=lineopacity, color = 'green')
          vec = s.P-s.Q
          gamma2 = 180/np.pi*np.arctan2(vec[1],vec[0])
          if s.a[num+1] > 0:
            gamma1 = gamma2-90
          else:
            gamma1 = gamma2
            gamma2 = gamma1+90
          sub.add_patch(patches.Arc(E, 0.05, 0.05, theta1=gamma1, theta2=gamma2, edgecolor='green', alpha=lineopacity, lw=1))
          angle = np.pi/180*.5*(gamma1+gamma2)
          center = E + 0.0125*np.array((np.cos(angle),np.sin(angle)))
          fig.gca().add_artist(plt.Circle(center, 0.002, alpha=lineopacity, color = 'green'))
        if detail > 2:
          plt.annotate('$E^'+str(num+1)+'$', E, fontsize=16, xytext=(1,1), textcoords='offset points', alpha=lineopacity, color = 'black')
          if num == 0:
            nvec = vec/norm(vec)
            p1 = E-0.1*nvec
            vec = s.Q-s.P
            nvec = vec/norm(vec)
            normal = np.array((-nvec[1],nvec[0]))
            plt.arrow(p1[0],p1[1], 0.1*nvec[0],0.1*nvec[1], head_width=0.01, alpha=lineopacity, color='black')
            plt.annotate(r'$u^0$', p1+.06*nvec, fontsize=16, xytext=(1,0), textcoords='offset points', alpha=lineopacity, color = 'black')
            plt.arrow(p1[0],p1[1], 0.1*normal[0],0.1*normal[1], head_width=0.01, alpha=lineopacity, color='black')
            plt.annotate(r'$u_0\bot$', p1+.04*normal, fontsize=16, xytext=(-16,1), textcoords='offset points', alpha=lineopacity, color = 'black')

      for num, H in enumerate(s.H_int): # the outer H which is P and Q is already in L
        fig.gca().add_artist(plt.Circle(H, 0.005, alpha=lineopacity, color = 'tab:gray' if args.gray else 'blue'))
        if detail > 1:
          plt.annotate('$P^'+str(num+1)+'$', H, fontsize=26, xytext=(1,1), textcoords='offset points', alpha=lineopacity, color = 'blue')
        
      for num, C in enumerate(s.C[1:-1]): # arcs are only around interior C
        if detail > 1:
          plt.annotate('$C^'+str(num+1)+'$', C, fontsize=26, xytext=(1,1), textcoords='offset points', alpha=lineopacity, color = 'gray')
        if detail > 1:
          fig.gca().add_artist(plt.Circle(C, 0.005, alpha=lineopacity, color = 'blue'))
          L=s.L[num+1][0]
          vec = C-L
          plt.annotate('$r^'+str(num+1)+'$', L+.5*vec, fontsize=20, xytext=(0,0), textcoords='offset points', alpha=lineopacity, color = 'blue')
          w = s.w
          L = L#+w*vec/norm(vec)
          sub.add_line(plt.Line2D((L[0],C[0]),(L[1],C[1]), alpha=lineopacity, color='gray', lw=1))
          L=s.L[num][1]
          vec = C-L
          w = s.w
          L = L#+w*vec/norm(vec)
          sub.add_line(plt.Line2D((L[0],C[0]),(L[1],C[1]), alpha=lineopacity, color='gray', lw=1))
        if detail > 3:
          L=s.L[num][1]
          vec = C-L
          sub.add_line(plt.Line2D((L[0],C[0]),(L[1],C[1]), alpha=lineopacity, color='dodgerblue', lw=1, linestyle='-'))
          H = s.H[num+1]
          sub.add_line(plt.Line2D((L[0],H[0]),(L[1],H[1]), alpha=lineopacity, color='dodgerblue', lw=1.5, linestyle='-'))
          sub.add_line(plt.Line2D((H[0],C[0]),(H[1],C[1]), alpha=lineopacity, color='dodgerblue', lw=1, linestyle='--'))
          v1 = C-H
          a1 = np.arctan2(v1[1], v1[0])*180/np.pi # Orientation of arc alpha/2
          v2 = s.H[num] - H
          a2 = np.arctan2(v2[1], v2[0])*180/np.pi # Orientation of arc alpha/2  
          alpha2 = np.min((np.abs(a1-a2),360.-np.abs(a1-a2)))

          sub.add_patch(patches.Arc(H, 0.1, 0.1, theta1=a1, theta2=a1+alpha2, alpha=lineopacity, color = 'dodgerblue', edgecolor='dodgerblue', facecolor='dodgerblue', lw=1))
          angle = a1+.5*(alpha2)
          center = H + 0.05*np.array((np.cos(angle),np.sin(angle)))
          leg = r'$\alpha$'+str(num+1)
          plt.annotate(leg, center, fontsize=16, xytext=(0,0), textcoords='offset points', alpha=lineopacity, color = 'dodgerblue')
  
      for L in s.L:  
        fig.gca().add_artist(plt.Circle(L[0], 0.005, alpha=lineopacity, color = 'gray'))
        fig.gca().add_artist(plt.Circle(L[1], 0.005, alpha=lineopacity, color = 'gray'))
      
      # plot normals
      if detail < 0:
        assert len(s.T) == len(s.M)
        for i, T in enumerate(s.T):
          M = s.M[i]
          p1 = .5*(s.H[i] + s.H[i+1])
          p2 = p1 + .1*M 
          
          sub.add_line(plt.Line2D((p1[0],p2[0]),(p1[1],p2[1]), alpha=lineopacity, color= 'red'))

    # plot half circles for start and end of noodle 
    M = s.M[0] # normal of the first segment tells us where to draw the radius
    angle = np.arctan2(M[1], M[0])*180/np.pi
    sub.add_patch(patches.Arc(s.P, 2*s.w, 2*s.w, theta1=angle, theta2=angle-180, edgecolor=s.color, alpha=lineopacity, lw=1))
    arc_patch(s.P, 0, s.w, theta1=angle, theta2=angle+180, facecolor=s.color, edgecolor=s.color, lw=1, alpha=0.3, fill=True)
    M = s.M[-1]
    angle = np.arctan2(M[1], M[0])*180/np.pi
    sub.add_patch(patches.Arc(s.Q, 2*s.w, 2*s.w, theta1=angle+180, theta2=angle, edgecolor=s.color, alpha=lineopacity, lw=1))
    arc_patch(s.Q, 0, s.w, theta1=angle, theta2=angle-180, facecolor=s.color, edgecolor=s.color, lw=1, alpha=0.3, fill=True)
  
    # plot arcs
    r = 2*s.radius
    for C,cosa, K1, K2 in s.arc:
      v1 = K1-C
      v2 = K2-C
      alpha = np.arctan2(v1[1], v1[0])*180/np.pi # Orientation of arc defined by C->K2
      beta = np.arctan2(v2[1], v2[0])*180/np.pi # Orientation of arc defined by C->K2

      gamma1 = min(alpha,beta)
      gamma2 = max(alpha,beta)
      if (gamma2-gamma1 > 180): # for edge case of arctan2 domain, angle difference should never be this large
        tmp = gamma1
        gamma1 = gamma2
        gamma2 = tmp+360

      if detail > 1:
        sub.add_patch(patches.Arc(C, r, r, theta1=0, theta2=360, edgecolor='gray', alpha=lineopacity, lw=1))
      sub.add_patch(patches.Arc(C, r-2*s.w, r-2*s.w, theta1=gamma1, theta2=gamma2, edgecolor=s.color, alpha=lineopacity, lw=1))
      arc_patch(C, s.radius-s.w, s.radius+s.w, theta1=gamma1, theta2=gamma2, facecolor=s.color, edgecolor=s.color, lw=1, alpha=0.3, fill=True)
      sub.add_patch(patches.Arc(C, r+2*s.w, r+2*s.w, theta1=gamma1, theta2=gamma2, edgecolor=s.color, alpha=lineopacity, lw=1))
      sub.add_patch(patches.Arc(C, r, r, theta1=gamma1, theta2=gamma2, edgecolor=s.color, alpha=lineopacity, lw=1.5))

  return fig

# discretized polygonal arc as patches.Arc can't be used with fill=True
def arc_patch(center, inner_radius, outer_radius, theta1, theta2, ax=None, resolution=50, **kwargs):
    # make sure ax is not empty
    if ax is None:
        ax = plt.gca()
    # generate the points
    theta = np.linspace(np.radians(theta1), np.radians(theta2), resolution)
    ateht = np.linspace(np.radians(theta2), np.radians(theta1), resolution)
    points = np.vstack((np.concatenate((inner_radius*np.cos(theta) + center[0], outer_radius*np.cos(ateht) + center[0])), 
                        np.concatenate((inner_radius*np.sin(theta) + center[1], outer_radius*np.sin(ateht) + center[1]))))
    # build the polygon and add it to the axes
    poly = patches.Polygon(points.T, closed=True, **kwargs)
    ax.add_patch(poly)
    return poly

# do L2-tracking of density
def track(shapes, design_track, saveall = False, track_iter = 30):
  rho = cfs_map_to_design()
  dist = np.sum((design_track.reshape(np.prod(glob.n), order='F')-rho)**2)
  var_all = glob.var_all()
  glob.iter = 0
  print("l2-distance: ", dist)
  #err = sciopt.check_grad(eval_l2,grad_l2,var_all,design_track)
  if glob.gradient_check:
    print("gradient check of l2-tracking: ", sciopt.approx_fprime(var_all, eval_l2, 1.5e-8, [design_track,False])-grad_l2(var_all,[design_track,False]))
  bound = sciopt.Bounds(lb=np.zeros(len(var_all)),ub=2*np.ones(len(var_all)))
  #cfs_get_sparsity_arc_overlap(None)
  n_const = len(get_vector_arc_overlap(var_all))
  if n_const > 0:
    nonlinear_constraints = NonlinearConstraint(get_vector_arc_overlap, np.zeros(n_const), np.inf*np.ones(n_const), jac=get_jacobian_arc_overlap)
    res = sciopt.minimize(eval_l2, var_all, args=[design_track, saveall], method="trust-constr", jac=grad_l2, constraints=nonlinear_constraints, bounds=bound, options={"disp": True, "verbose": 2, "maxiter": track_iter})
  else:
    res = sciopt.minimize(eval_l2, var_all, args=[design_track, saveall], method="trust-constr", jac=grad_l2, bounds=bound, options={"disp": True, "verbose": 2, "maxiter": track_iter})
  return res

# function value for L2 tracking optimization
def eval_l2(var_all, arg): # args = [density_track, save_figures]
  for shape in glob.shapes:
    shape_num = shape.id
    var = var_all[shape.base:shape.base+len(shape.optvar())]
    cfs_set_spaghetti(shape_num, [[var[0], var[1]], [var[2], var[3]]], var[5:], var[4])
  des = cfs_map_to_design()
  if arg[1]:
    ot.write_density_file('giffiles/dens' + str(glob.iter).zfill(4) + '.density.xml', np.reshape(des, (glob.n[0],glob.n[1]), 'F'))
  return np.sum((des - arg[0].reshape(np.prod(glob.n), order='F'))**2)

# gradient for L2 tracking optimization
def grad_l2(var_all, arg): # args = [density_track, save_figures]
  for shape in glob.shapes:
    shape_num = shape.id
    var = var_all[shape.base:shape.base+len(shape.optvar())]
    cfs_set_spaghetti(shape_num, [[var[0], var[1]], [var[2], var[3]]], var[5:], var[4])
  if arg[1]:
    try:
      fig = plot_data(800,shapes,args.detail)
      fig.savefig('giffiles/track' + str(glob.iter).zfill(4) + '.png')
      plt.close(fig)
      glob.iter += 1
    except:
      print("Could not save file in iteration ", glob.iter)
    # read the xml file to memory
    xml = ut.open_xml(args.input)
    for num_var, var in enumerate(var_all):
      # modify the xml file
      ut.replace(xml, '//shapeParamElement[@nr='+str(num_var)+']/@design', str(var))

    # write to a new name
    xml.write('giffiles/track' + '.density.xml')
  drho = 2*(cfs_map_to_design() - arg[0].reshape(np.prod(glob.n), order='F'))
  sens = cfs_get_gradient(drho, "L2-tracking")
  return sens

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
  
  total = sum([len(s.optvar()) for s in shapes])
  assert total == shapes[-1].base + len(shapes[-1].optvar()) 
  ddist         = np.zeros((total,N,N,1)) # ad for fast_dist_ad
  drho_ad       = np.zeros((total,N,N,1)) # ad for fast_rho_ad 
  
  for s in shapes:
    var = s.optvar()
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
          for e in range(len(s.optvar())):
            drho_ad[s.base+e,i,j,0] = g[e]

        if derivative and detailed:
          g = s.ddist(var,X,idx)
          for e in range(len(s.optvar())):
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
      for e, n in enumerate(s.optvar_names()):
        pd["d_rho / d_s" + str(s.id) + '_' + n] = drho_ad[s.base+e]

  if derivative and detailed:
    for s in shapes:
      for e, n in enumerate(s.optvar_names()):
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
  var = np.copy(s.optvar())
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
  dv = s.optvar()
  dap   = np.ones(res) # using dist
  dapad = np.ones(res) # using fast_dist_ad
  for i, x in enumerate(np.linspace(0,1,res)):
    X = [x,y]
    dap[i] = s.dist(X)[0]
    dapad[i] = s.fast_dist_ad(dv,X,idx[i])
  s.set(s.P, s.Q, org, s.p) # reset   
  
  print('# lineplot of distance for height ',y, ' optvar=',var)
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
  parser.add_argument('--saveall', help="save all sets as image with the given format. Might be png, pdf, eps, vtp", action='store_true')
  parser.add_argument('--detail', help="level of technical details for spaghetti plot", choices=[0, 1, 2, 3, 4], default=1, type=int)
  parser.add_argument('--rhomin', help="void outside the feature", type=float, default=1e-6)
  parser.add_argument('--rhomax', help="density of solid inside feature", type=float, default=1.0)
  parser.add_argument('--transition', help="size of the transition zone (2*h)", type=float, default=.1)
  parser.add_argument('--boundary', help="type of boundary modelling ('poly' or 'linear')", choices=['poly', 'linear'], default='poly')
  parser.add_argument('--combine', help="type of (smooth) maximum function for combination of shapes", choices=['softmax', 'KS', 'p-norm', 'max'], default='softmax')
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
  parser.add_argument('--gradient_check', help="check internal spaghetti derivatives", action='store_true')
  parser.add_argument('--track', help="track given density file to represent topology using spaghetti shapes")
  parser.add_argument('--track_iter', help="number of optimization iterations for tracking", type=int, default=30)
  parser.add_argument('--padnormals', help="pad normals with additional normals (zero-valued) so there is a total of x normals", type=int, default=0)
  parser.add_argument('--gray', help="plot grayscale image", action='store_true')
  parser.add_argument('--noticks', help="omit axis tick labels", action='store_true')

  args = parser.parse_args()
  
  if not os.path.exists(args.input):
    print("error: cannot find '" + args.input + "'")
    os.sys.exit()

  colors = ['tab:gray'] if args.gray else ['b', 'g', 'r', 'c', 'm', 'y', 'tab:orange', 'tab:brown', 'cornflowerblue', 'lime', 'tab:gray']
  
  # sets some values in glob!
  shapes, domain = read_xml(args.input, args.set, args.radius, args.cfseval)
  
  glob.shapes = shapes
  glob.rhomin = args.rhomin
  glob.rhomax = args.rhomax
  glob.transition = args.transition
  glob.boundary = args.boundary
  glob.combine = args.combine
  glob.order = args.order
  glob.gradient_check = args.gradient_check
  glob.design = ["density"]
  dict = {"order": args.order,
    "silent": '1',
    "p": 8}
  if args.density:
    glob.n = [100, 50, 1]

  if args.lineplot:
    if not args.noshow:
      print("error: use lineplot with --noshow")
      os.sys.exit()
    lineplot(args.lineplot)

  if args.track:
    to_track = ot.read_density(args.track, attribute="design", fill=1.0)
    #design = ['density', 'rotAngle']
    design = ["density"]
    cfs_init(args.rhomin, args.radius, args.boundary, args.transition, args.combine, glob.n[0], glob.n[1], glob.n[2], 1/min(glob.n[0:2]), 'straight', design, dict)
    res = track(shapes, to_track, args.saveall, args.track_iter)
    shapes = glob.shapes
    write_xml(args.track[0:-12] + '_tracked.density.xml', shapes, remove_ghosts=True, padnormals=args.padnormals)

  if args.vtk:
    write_vtk(args.vtk, args.vtk_res, args.vtk_detailed, args.vtk_sens)

  if args.density or args.cfseval:
    tmp = copy.deepcopy(vars(args)) # many settings are same from cfs as from command line
    settings = dict((k,v) for k, v in tmp.items() if v) # remove key with None value
    # we have already read nx and ny from spaghetti.py
    scale = args.density_res / glob.n[0] # rescale, eather fine for --density or corarse for debugging --cfseval
    glob.n = [int(scale * glob.n[0]), int(scale * glob.n[1])]
    glob.dx = 1.0/glob.n[0]
    
    dict = { "order": args.order, "silent": 1}

    cfs_init(settings, glob.design, dict)
    
    if args.density:
      rho = density(glob.n)
      ot.write_density_file(args.density,rho)
    else:
      des = cfs_map_to_design()
      if len(design)>1:
        des = np.reshape(des, (glob.n[0]*glob.n[1],2), 'F')
      #des = np.reshape(des, (args.density_res*args.density_res,len(design)), 'C')
      ot.write_multi_design_file(args.input[0:-12] + '.eval.density.xml', des, glob.design)
      dummy_drho_vec = np.ones(glob.n[0] * glob.n[1])
      drho = cfs_get_gradient(dummy_drho_vec, 'compliance')

  if args.saveall:
    xml = ut.open_xml(args.input)
    sets = []
    for set in xml.xpath('//set'):
      sets.append(int(set.attrib['id']))
      #print(etree.tostring(set.xpath('//[@id]')))
    print('read', len(sets),'sets from',args.input + ': ',end='') # only python3
    path = os.getcwd()
    dir = os.path.join(path, 'giffiles')
    if not os.path.exists(dir):
      os.mkdir(dir)
    for i in sets:
      print(i,' ',end='' if i < sets[-1] else '\n',flush=True)
      shapes, domain = read_xml(args.input, i, args.radius, args.cfseval)
      glob.shapes = shapes
      fig = plot_data(800,shapes,args.detail, domain)
      if args.noticks:
        plt.tick_params(
          axis='both',          # changes apply to x-axis and y-axis
          which='both',      # both major and minor ticks are affected
          bottom=True,      # ticks along the bottom edge are off
          top=False,         # ticks along the top edge are off
          labelleft=False,   # labels along the left edge are off
          labelbottom=False) # labels along the bottom edge are off
      fig.savefig('giffiles/' + str(i).zfill(4) + '.png')
      plt.close(fig)

  fig = plot_data(800,shapes,args.detail, domain)
  if args.noticks:
    plt.tick_params(
        axis='both',          # changes apply to x-axis and y-axis
        which='both',      # both major and minor ticks are affected
        bottom=True,      # ticks along the bottom edge are off
        top=False,         # ticks along the top edge are off
        labelleft=False,   # labels along the left edge are off
        labelbottom=False) # labels along the bottom edge are off
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

