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
  global params, log
  params = options
  params['max_bisec'] = int(params['max_bisec'])
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
  
  params['n_design'] = n_design  
  params['dim'], nx, ny, nz = cfs.getDims()
  n_elems = nx * ny * nz
  params['n_elems'] = n_elems
  
  print("nx,ny,nz,elems:",nx,ny,nz,n_elems)
  
  params['max_iter'] = max_iter
  # filter radius
  assert('rad_filter' in params) 
  print("python init(): filter radius=", params['rad_filter'], " nelems=", n_elems, " maxIter=", params['max_iter'], " maxBisec=", params['max_bisec'])
  
  #### get bounds from cfs
  xl = np.zeros(n_design)
  xu = np.zeros(n_design)
  gl = np.zeros(n_constr)
  gu = np.zeros(n_constr)
  
  cfs.bounds(xl,xu,gl,gu)
  
  print("init(): bounds on rho: lower=",xl[0], " upper=",xu[0])
  
  # convert parameters to float for matlab
  params['levels'] = float(params['levels'])
  params['samples'] = float(params['samples'])
  params['pen_filter'] =float(params['pen_filter'])
  params['rad_filter'] = float(params['rad_filter'])
  params['lower'] = float(xl[0]) # we're only interested in bounds for \rho
  params['upper'] = float(xu[0])
  params['v_bound'] = gu[0]
  
  #### get core stiffness tensor cfs
  stiffness = np.zeros((6,6)) if params['dim'] == 3 else np.zeros((3,3))
  cfs.getOrgStiffness(stiffness)
  stiffness = stiffness
  #print("core stiffness tensor:\n", stiffness, "\n")
  # check if stiffness is symmetric
  assert(np.allclose(stiffness, stiffness.T))
  
  #### get SIMP exponent from cfs
  params['p'] = cfs.getSimpExponent()
  #params['p'] = float(params['p'])
  assert(params['p'] > 0)
  
  # lower asymptote value
  L = 0.0
  
  global subsolver, F
  subsolver = matlab_eng.subproblem_solver(float(nx), float(ny), float(nz), params['lower'], params['upper'], params['levels'], params['samples'], L, params['pen_filter'], params['rad_filter'], matlab.double(stiffness.tolist()), float(params['p']))

  # filter matrix
  F = np.asarray(matlab_eng.calc_filter(subsolver, params['rad_filter'], nargout=1))  
  #print('python: init optimization parameteres', params)
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
  cfs.evalobj(x)
  cfs.commitIteration()
  #print("initial design:",x)
  #dfdH = get_dfdH(x)
  
  iter = 1
  max_iter = params['max_iter']
  log.write('{:<5} | {:<12} | {:<12} | {:<12} | {:<12} | {:<8}| {:<6} | {:<12} || {:<10} | {:<10} | {:<10} | {:<15}\n'.format('iter', 'f_merit', 'f_phys', 'f_reg', 'f_vol', 'vol', '#bisecs', ' lambda', ' t_iter', 't_vol', 't_filt', 't_avg_subprob'))
  log.write("-"*200 + "\n")
  #log.write('{:<100}\n'.format('='))
  
  # current iterate in terms of stiffness tensors
  rows = 3 if params['dim'] == 2 else 6
  
  global start_time
  t = timer()
  
  print("time until iteration loop:",t-start_time)
  
  while iter <= max_iter and not cfs.dostop():
    dfdH = get_dfdH(x)
    x = sgp(iter, x, dfdH)
    cfs.commitIteration()
    iter += 1
  
#   print("solution:",x)  
  
  matlab_eng.quit()
  return 1

def get_dfdH(x):
  global params
  
  n_design= params["n_design"]
  n_elems= params["n_elems"]
  
  rows = 3 if params['dim'] == 2 else 6
  dfdH = np.ones((n_elems,rows,rows))
  
  #cfs.evalobj(x)
  cfs.get_dfdH(dfdH)
  
#   print("\ndfdH:")
#   for e in range(n_elems):
#     print("e=",e,"\n")
#     print(dfdH[e])
    
  return dfdH

# x: current design values
# dfdH: gradient of objective 
def sgp(iter, x, dfdH):
  t_start = timer()
  
  global params, log, F
  assert(F is not None)
  
  n_levels = params['levels']
  n_samples = params['samples']
  n_design = params['n_design']
  n_elems = params['n_elems']
  p_filt = params['pen_filter']
  v_bound = params['v_bound']
  
  # convert to matlab data structure
  dfdH = matlab.double(dfdH.tolist())
  x = matlab_eng.transpose(matlab.double(x.tolist()))
  
  x, lmbd, vol, f_model,  t_vol, n_bisecs, t_avg_sp = volume_bisection(x, dfdH)
  
#   print("lambda:",lmbd)
  f_phys = cfs.evalobj(x)
  
  assert(len(x) == n_design) 
  
  rho = x[:n_elems]
#   print("x:",x)
#   print("rho:",rho)
  
  t_filt_start = timer() 
  rho_filt = F.dot(rho)
  t_filt_end = timer()
  
  f_reg = np.linalg.norm(rho_filt - rho)**2
  f_vol = vol - v_bound
  
  # merit function: compliance + regularization + volume constraint (+ proximal point term)
  f_merit = f_phys + 0.5* p_filt * f_reg + lmbd * f_vol
  
  t_end = timer()
  
  log.write('{:<5d} | {:<12.6f} | {:<12.6f} | {:<12.6f} | {:<12.6f} | {:<8.6f} | {:<6d} | {:<8.6f} || {:<10.6f} | {:<10.6f} | {:<10.6f} | {:<15.6f}\n'.format(iter, f_merit, f_phys, f_reg, f_vol, vol, n_bisecs, lmbd, t_end-t_start, t_vol, t_filt_end-t_filt_start,t_avg_sp) )
  
  return x

# from previous iteration:
# x: design values
# dfdH: tensor derivatives
# material tensors H related to 'x' computed in matlab
def volume_bisection(x, dfdH):
  t_start = timer()
  
  global params, subsolver
  xl = params['lower']
  xu = params['upper']
  n_elems = params['n_elems']
  n_bisections = params['max_bisec']
  #n_bisections = 50
  assert(xl > 0 and xu < 1.01 and n_bisections > 0)
  
  diff_vol = 9999
  eps = 1e-6
  
  l1 = 0
  l2 = 1e4
  
  iter = 0
  x_new = x
  
  sum_s_time = 0
  # perform bisection as lang as volume residual is > eps
  while abs(diff_vol) > eps and iter < n_bisections:
    # volume penalty / Langrange multiplier
    lmbd = 0.5 * (l2+l1)
    
    t_sp_start = timer()
    # subsolver returns only density values
    x_new, f_min = matlab_eng.solve(subsolver,x,dfdH,lmbd,nargout=2)
    t_sp_end = timer()
    sum_s_time  = sum_s_time + (t_sp_end - t_sp_start)
    
    # TODO: check data structure of x_new
    x_new = np.append(x_new[0],x[n_elems:])
    # append constant tensor entries for cfs
    #print("xnew:",x_new)
    
    vol = calc_g(x_new)
    diff_vol = params['v_bound'] - vol
    
    if diff_vol > 0:
      l2 = lmbd
    else:
      l1 = lmbd
    
    iter = iter + 1
    
  #print("x_new:",x_new[:n_elems],"diff_vol:",diff_vol, "vol:",vol, "\n")  
  t_end = timer()
  
  if diff_vol > eps:
    print("WARNING: volume constraint violated by",diff_vol)
  
  return x_new, lmbd, vol, f_min, t_end - t_start, iter, sum_s_time/iter

# calc volume
def calc_g(x):
  g = np.zeros(1)
  cfs.evalconstrs(x, g)
  return g[0]
