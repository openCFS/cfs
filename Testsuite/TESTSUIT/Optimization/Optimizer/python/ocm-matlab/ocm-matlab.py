#!/usr/bin/env python

import matlab.engine
import os, cfs
import numpy as np

matlab_eng = None
# handle for a matlab ocm object
ocm_obj = None
# number of elements
n_elems = 0
# number of constraints
n_constr = 0
# number of max. iterations
max_iters = 0
# upper volume bound
vbound = 0

# optional for cfs. 
# n: number of design variables
# m: number of active constraints
# @return an integer with 0 for ok
def setup():
  global matlab_eng
  # start engine
  matlab_eng = matlab.engine.start_matlab()
  # set location of matlab file - currently same location as this script
  matlab_eng.addpath(os.path.dirname(__file__),nargout=0)
  
  print('ocm.py: setup() called')
  return 1

def init(n, m, iters, options):
  global n_elems, n_constr, max_iters, vbound
  n_elems = n
  n_constr = m
  max_iters = iters
  
  assert(n_elems > 0)
  assert(n_constr == 1) # necessary for ocm
  assert(max_iters > 0)
  
  # bounds
  xl = np.zeros(n_elems)
  xu = np.zeros(n_elems)
  gl = np.zeros(n_constr)
  gu = np.zeros(n_constr)
  
  cfs.bounds(xl,xu,gl,gu)
  
  # make sure we have same values for all elements -> store only one scalar
  all([e == xl[0] for e in xl])
  all([e == xu[0] for e in xu])
  
#   print('gbounds',gl, gu)
  assert(gl[0] <= gu[0])
  vbound = gu[0]
  assert(vbound > 0)
  
  # call matlab script that declares
  #matlab_eng.init_oc(nargout=0)
  
  global ocm_obj
  ocm_obj = matlab_eng.OCM(float(xl[0]),float(xu[0]),float(options["move"]),float(options["damp"]))
  
  print('ocm-matlab.py: init(n, m, options) called')
  
  # get properties from matlab ocm class
  # works only if properties are publix
  #matlab_eng.workspace["ocm_obj"] = ocm_obj
  #print("xl: ", matlab_eng.eval("ocm_obj.xl"), " xu:", matlab_eng.eval("ocm_obj.xu"), " move:", matlab_eng.eval("ocm_obj.move"), " damp:", matlab_eng.eval("ocm_obj.damp"))
  
  return 1

# expected and called by cfs
# @return an integer with 0 for ok
def solve():
  print('ocm.py: solve() called n=', n_elems, ' m=', n_constr, ' max_iter=', max_iters)

  # design
  x = np.zeros(n_elems)
  cfs.initialdesign(x)
  
  # volume constraint gradient only required once
  dv = np.zeros(n_elems) 
  cfs.evalgradconstrs(x, dv)
  
  dc = np.zeros(n_elems)
  cfs.evalgradobj(x,dc)
  iter = 0
  
  while iter < max_iters and not cfs.dostop():
    x, lmbd = oc(x, dc, dv)
    cfs.evalgradobj(x, dc)
    cfs.evalgradconstrs(x, dv)
    #print('ocm.py: next iter', np.average(x), lmbd)
    iter += 1    
    
  return 0

# calc volume
def calc_g(x):
  g = np.zeros(1)
  cfs.evalconstrs(x, g)
  return g[0]

# calc the next oc-step.  
# @param x the 1D numpy design array
# @param bound constraint value (volume constraint)
# @param dc objective function gradient (compliance gradient)
# @param dv constraint gradient (volume gradient). Assumed to be constant over iterations, homogeneous on regular mesh
# @return new design and lmbd
def oc(x, dc , dv):
  import sys
  global matlab_eng, vbound, ocm_obj
  assert(matlab_eng)
  
  l1 = 0
  l2 = 1e5

  xnew = x.copy()
  
  # convert to matlab data structure
  m_x = matlab.double(x.tolist())
  m_dc = matlab.double(dc.tolist())
  m_dv = matlab.double(dv.tolist())
  
  while (l2-l1) > 1e-6:
    lmbd = 0.5 * (l2+l1)
    
    xnew = matlab_eng.oc_update(ocm_obj, m_x, m_dc, m_dv, lmbd, nargout=1)
    
    # matlab returns matrix of size (1, nelems)
    xnew = np.asarray(xnew[0])
    vol = calc_g(xnew)
    cmp = vbound - vol
     
    if cmp > 0:
      l2=lmbd
    else:
      l1=lmbd
    
#     print("\nlower: %2.6f, upper: %2.6f, lambda: %2.6f, diff: %2.6f, vol: %2.6f" % (l1,l2,lmbd,cmp,vol))
#   
#   print("lower:", l1, " upper:", l2)    
  return xnew, lmbd