import numpy as np

print('hallo from local kernel.py')

import cfs


# global variaable set by slope_sparsity()
jac = None

def slope_sparsity(opt):
  print('kernel.py: slope_sparsity')

  dim, nx, ny, nz = cfs.getDims()
  assert dim == 2
  
  global jac
  jac = [] # list of arrays of 0-based varible indices
  
  # we make no reserve slope constraints, hence add for both directions
  for y in range(ny):
    for x in range(nx):
      idx = y * nx + x
      if x < nx-1: # place for right neighbor
        jac.append(np.array([idx, idx+1])) # even index
        jac.append(np.array([idx, idx+1])) # odd
      if y < ny -1: # place for upper neighbor
        jac.append(np.array([idx,idx+nx])) # even
        jac.append(np.array([idx,idx+nx])) # odd   
  
  return jac
  
def slope_eval(loc_idx):
  #print('kernel.py: slope_eval')
  slope  = jac[loc_idx]
  # instead of requesting each value individually we could have a callback hook for opt_eval_func
  # and get the whole design to be served by the individual slope_eval() calls (see also python_function test)
  # but we are not for efficiency but for rapid prototyping
  first  = cfs.get_opt_design_value(slope[0])
  second = cfs.get_opt_design_value(slope[1])
  even   = (loc_idx % 2) == 0
  return first - second if even else second - first
  

def slope_grad(loc_idx):
  #print('kernel.py: slope_grad')
  even = (loc_idx % 2) == 0
  return np.array([+1.0, -1.0]) if even else np.array([-1.0, +1.0]) 
