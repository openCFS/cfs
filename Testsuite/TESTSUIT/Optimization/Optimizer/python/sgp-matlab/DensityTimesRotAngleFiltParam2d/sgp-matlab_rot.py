#!/usr/bin/env python

import matlab.engine
import os, cfs
import numpy as np
from timeit import default_timer as timer

# define global variables
matlab_eng = None
# dictionary with all optimization parameters
params = None
# object handle for sub-problem solver
subsolver = None
# sgp log file
log = None
# filter matrix from matlab
F = None
# core stiffness tensor
E0 = None
# power-law exponent
p = 1
# sensitivity filtering?
filt_sens = 1

start_time = timer()

# optional for cfs. 
# n: number of design variables
# m: number of active constraints
# @return an integer with 0 for ok
def setup():
  global matlab_eng
  # start engine
  matlab_eng = matlab.engine.start_matlab()
  assert(matlab_eng)
  
  print('python: setup() called')
  return 1

def init(n_design, n_constr, max_iter, sim_name, options):
  global params, log, subsolver, F, filt_sens
  params = options
  ## make sure these parameters are provided by cfs
  # parameters for hierarichal grid in sub problems
  assert(params['levels'])
  assert(params['samples'])
  assert(params['pen_filter'])
 
  import os.path
  
  # set location of matlab file - currently same location as this script
  if not os.path.exists(params['matlabFile']):
    Exception("MATLAB file " + str(params['matlabFile']) + " does not exist!")
  
  matlab_eng.addpath(os.path.dirname(params['matlabFile']),nargout=0)
  
  # set location of mtimesx directory
  if not os.path.exists(params['mtimesxDir']):
    Exception("directory " + str(params['mtimesxDir']) + " does not exist!")
  
  # genpath() generates path name for recursive inclusion 
  matlab_eng.addpath(matlab_eng.genpath(params['mtimesxDir']),nargout=0)  
  
  log = open(sim_name+".sgp.log", "w")
  
  params['max_iter'] = max_iter
  params['dim'], nx, ny, nz = cfs.getDims()
  print(params['dim'], nx, ny, nz)
  n_elems = nx*ny*nz
  params['n_elems'] = n_elems
  
  print("python init(): dim=",params['dim'], " n_elems=", params['n_elems'], " n_design=", n_design, " max_iter=",params['max_iter'])
  
  assert(n_design > 0)
  params['n_design'] = n_design
  
  # filter radius
  assert('rad_filter' in params) 
  print("python init(): filter radius=", params['rad_filter'], " nelems=", params['n_elems'], " n_design=", n_design, " maxIter=", params['max_iter'])
  
  #### get bounds from cfs
  xl = np.zeros(n_design)
  xu = np.zeros(n_design)
  gl = np.zeros(n_constr)
  gu = np.zeros(n_constr)
  
  # we have designs: density \rho, angle \alpha
  # n_elems*2
  cfs.bounds(xl,xu,gl,gu)
  
  assert(len(xl) == n_elems * 2)
  assert(len(xu) == n_elems * 2)
  assert(gl[0] <= gu[0])
  assert(gu[0] > 0)
  
  # convert parameters to float for matlab
  params['levels'] = float(params['levels'])
  params['samples'] = float(params['samples'])
  params['pen_filter'] =float(params['pen_filter'])
  params['rad_filter'] = float(params['rad_filter'])
  params['xl'] = float(xl[0])
  params['xu'] = float(xu[0])
  params['al'] = float(xl[n_elems])
  params['au'] = float(xu[n_elems])
  params['v_bound'] = gu[0]
  
  params['max_bisec'] = int(params['max_bisec'])  if 'max_bisec' in params else 50
  
  print("optimization bounds: xl=",  params['xl'] , " xu=",  params['xu'], " al=",  params['al'], "au=",  params['au'], " vbound=", gu[0])
  
  if 'filter' not in params:
    params['filter'] = "parameter" if params['pen_filter'] > 0 else "sensitivity"
  
  print(params['filter'])
  assert(params['filter'] == "sensitivity" or params['filter'] == "parameter")
  filt_sens = 1 if params['filter'] == "sensitivity" else 0
  
  
  
  #### get core stiffness tensor cfs
  stiffness = np.zeros((6,6)) if params['dim'] == 3 else np.zeros((3,3))
  # for debugging
  #stiffness[0][0] = 3.0
  #stiffness[1][1] = 0.1
  #stiffness[2][2] = 1.0
  cfs.getOrgStiffness(stiffness)
  print("core stiffness tensor:\n", stiffness, "\n")
  global E0, p
  E0 = stiffness
  # check if stiffness is symmetric
  assert(np.allclose(stiffness, stiffness.T))
  
  #### get SIMP exponent from cfs
  p = cfs.getSimpExponent()
  
  # lower asymptote value
  L = 0.0
  
  print("p=",p)
  
  subsolver = matlab_eng.subproblem_solver(float(nx), float(ny), float(nz), params['xl'], params['xu'], params['al'], params['au'], params['levels'], params['samples'], L, filt_sens, params['pen_filter'], params['rad_filter'], matlab.double(stiffness.tolist()), float(p))
  # filter matrix
  F = np.asarray(matlab_eng.get_filter(subsolver, params['rad_filter'], nargout=1))  
  print('python: init optimization parameteres', params)
  return 1

# expected and called by cfs
# @return an integer with 0 for ok
def solve():
  global params, matlab_eng
  
  print('sgp-matlab.py: solve() called\n')
  #print('SGP parameters:', params)
  
  n_design = params['n_design']

  # design
  x = np.zeros(n_design)
  cfs.initialdesign(x)
  #print("initial design:",x)
  f_phys = cfs.evalobj(x)
  cfs.commitIteration()
  
  iter = 1
  max_iter = params['max_iter']
  log.write('{:<5} | {:<12} | {:<12} | {:<12} | {:<12} | {:<8}| {:<6} | {:<12} || {:<10} | {:<10} | {:<10} | {:<15}\n'.format('iter', 'f_merit', 'f_phys', 'f_reg', 'f_vol', 'vol', '#bisecs', 'lambda', ' t_iter', 't_vol', 't_filt', 't_avg_subprob'))
  log.write("-"*190 + "\n")
  #log.write('{:<100}\n'.format('='))
  log.write('{:<5d} | {:<12} | {:<12.6f} | {:<12} | {:<12} | {:<8} | {:<6} | {:<8} || {:<10} | {:<10} | {:<10} | {:<15}\n'.format(0, " ", f_phys, " ", " ", " ", " ", " ", " ", " ", " " ," "))
  
  # current iterate in terms of stiffness tensors
  rows = 3 if params['dim'] == 2 else 6
  
  global start_time
  t = timer()
  
  print("time until iteration loop:",t-start_time)
  f_old = 999999
  eps = 1e-7
  while iter <= max_iter and not cfs.dostop():
    dfdH = get_dfdH(x)
    x, f_new = sgp(iter, x, dfdH, )
    cfs.commitIteration()
    
    # finalize if change in objective is sufficiently small
    if abs(f_new-f_old) < eps:
      break
    
    f_old = f_new
    iter += 1
    
  
#   print("solution:",x)  
  
  matlab_eng.quit()
  return 1

def get_dfdH(x):
  global params

  n_elems = params['n_elems']
  rows = 3 if params['dim'] == 2 else 6
  dfdH = np.ones((n_elems,rows,rows))
  
  cfs.get_dfdH(dfdH)
  
#   export2matlab(dfdH)
#   dump(x,dfdH)
  
  return dfdH

def export2matlab(dfdH):
  from scipy.io import savemat
  # swapaxis
  # matlab expectes array of size 3 x 3 x nelems
  tmp = np.moveaxis(dfdH, 1, -1)
  savemat("dfdH.mat",{"dfdH":dfdH})

def dump(x,dfdH):
  from matviz_rot import get_rot_3x3
  global E0,p
  
  print("x:",x)
  nelems = dfdH.shape[0]
  alpha = x[nelems:]
  for e in range(nelems):
    Q = get_rot_3x3(-alpha[e])
    rot_tens = x[e]**p * np.dot(Q, np.dot(E0, Q.transpose()) )
    print("e=",e," rho[e]:",x[e], " alpha[e]=",alpha[e])
    print("rot tens:",rot_tens)  
  
  print("\ndfdH:")
  for e in range(nelems):
    print("e=",e,"\n")
    print(dfdH[e])
  sys.exit()
  

# x: current design values
# dfdH: gradient of objective 
def sgp(iter, x, dfdH):
  t_start = timer()
  
  global params, log, F
  
  n_levels = params['levels']
  n_samples = params['samples']
  n_elems = params['n_elems']
  n_design = params['n_design']
  p_filt = params['pen_filter']
  v_bound = params['v_bound']
  
  #print("size dfdH",dfdH.shape)
  
  # convert to matlab data structure
  dfdH = matlab.double(dfdH.tolist())
  x = matlab_eng.transpose(matlab.double(x.tolist()))
  
  x, lmbd, vol, f_model,  t_vol, n_bisecs, t_avg_sp = volume_bisection(x, dfdH)
  
#   print("x:",x)
  
  f_phys = cfs.evalobj(x)
  
  assert(len(x) == n_design) 
  assert(F is not None)
    
  t_filt_start = timer() 
  dens_filt = F.dot(x[:n_elems])
  alpha_filt = F.dot(x[n_elems:])
  # concatenate vertically
  x_filt = np.concatenate((dens_filt,alpha_filt),axis=0)
  t_filt_end = timer()
  
  f_reg = np.linalg.norm(x_filt - x)**2
  f_vol = vol - v_bound
  
  # merit function: compliance + regularization + volume constraint (+ proximal point term)
  f_merit = f_phys + 0.5* p_filt * f_reg + lmbd * f_vol
  
  t_end = timer()
  
  log.write('{:<5d} | {:<12.6f} | {:<12.6f} | {:<12.6f} | {:<12.6f} | {:<8.6f} | {:<6d} | {:<8.6f} || {:<10.6f} | {:<10.6f} | {:<10.6f} | {:<15.6f}\n'.format(iter, f_merit, f_phys, f_reg, f_vol, vol, n_bisecs, lmbd, t_end-t_start, t_vol, t_filt_end-t_filt_start,t_avg_sp))   
  return x, f_merit

# from previous iteration:
# x: design values
# dfdH: tensor derivatives
# material tensors H related to 'x' computed in matlab
def volume_bisection(x, dfdH):
  t_start = timer()
  
  global params, subsolver
  xl = params['xl']
  xu = params['xu']
#   n_bisections = params['max_bisect']
  n_bisections = params['max_bisec']
  assert(xl > 0 and xu < 1.01 and n_bisections > 0)
  
  diff_vol = 9999
  eps = 1e-7
  
  l1 = 0
  l2 = 100
  
  iter = 0
  # 1st column: density values
  # 2nd column: angle values
  x_new = x
  
  sum_s_time = 0
  # perform bisection as lang as volume residual is > eps
  while abs(diff_vol) > eps and iter < n_bisections:
    # volume penalty / Langrange multiplier
    lmbd = 0.5 * (l2+l1)
    
    t_sp_start = timer()
    x_new, f_min = matlab_eng.solve(subsolver,x,dfdH,lmbd,nargout=2)
    t_sp_end = timer()
    # averaged time spent in subsolver
    sum_s_time  = sum_s_time + (t_sp_end - t_sp_start)
#     print("x_new:",x_new)
    x_new = np.asarray(x_new[0])
    # TODO: check data structure of x_new
    vol = calc_g(x_new)
    diff_vol = params['v_bound'] - vol
    
    if diff_vol > 0:
      l2 = lmbd
    else:
      l1 = lmbd
    
    iter = iter + 1
    
  #print("xnew=",x_new, " vol=",vol, " lambda=",lmbd)
    
  #print("lambda=",lmbd, " diff_vol=",diff_vol)
  if diff_vol > 1e-3:
    print("WARNING: volume bisection diff_vol=",diff_vol)
    
#   print("volume bisection: iter=",iter)
  t_end = timer()
  
  return x_new, lmbd, vol, f_min, t_end - t_start, iter, sum_s_time/iter

# calc volume
def calc_g(x):
  g = np.zeros(1)
  cfs.evalconstrs(x, g)
  return g[0]
