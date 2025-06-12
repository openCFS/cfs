# this implements Peter Dunning's auto projection filter update
# see Dunning, Achieving binary topology optimization solutions via automatic projection parameter increase
from pickle import NONE

print('hallo from auto_filter.py')

import cfs
import numpy as np

import auto_filter_tools as atf

# some global variables
# -----------------------
# which optimizer
ocm = None
mma = None
org_move_limit = None
python = None
snopt = None
# last objective for auto_filter
last_obj = None
last_last_obj = None
# beta updte speed
gamma = None
delta_beta_max_factor = None
# general
adapt_move_limits = None
divider = 6
tz_divider = None 
auto_divider = None
halt_on_osc = None
max_beta = None

# called as hook from cfs
def my_opt_post_init(opt):
  print('my_opt_post_init hook called')
  global gamma, delta_beta_max_factor, ocm, mma, snopt, adapt_move_limits, tz_divider, auto_divider, halt_on_osc, org_move_limit, max_beta

  gamma = float(opt['gamma']) if 'gamma' in opt else 1e-4
  delta_beta_max_factor = float(opt['delta_beta_max_factor']) if 'delta_beta_max_factor' in opt else 0.2
  
  print('check adapt_move_limits',opt['adapt_move_limits'], (eval(opt['adapt_move_limits'])),bool(eval(opt['adapt_move_limits'])))
  
  adapt_move_limits = bool(eval(opt['adapt_move_limits'])) if 'adapt_move_limits' in opt else True 
  tz_divider = float(opt['tz_divider']) if 'tz_divider' in opt else 6
  auto_divider = bool(eval(opt['auto_divider'])) if 'auto_divider' in opt else False
  halt_on_osc = bool(eval(opt['halt_on_osc'])) if 'halt_on_osc' in opt else False
  max_beta = float(eval(opt['max_beta'])) if 'max_beta' in opt else 500
  print('settings: gamma:',gamma,'delta_beta_max_factor',delta_beta_max_factor, 'adapt_move_limits',adapt_move_limits, 'tz_divider',tz_divider)

  # register logging beta
  filter = cfs.get_opt_filter_values()
  tf = int(filter['total_filters'])
  print('total_filters',tf)
  cfs.opt_register_log_property('beta',filter['beta'])
  cfs.opt_register_log_property('delta_beta',"-1")
  cfs.opt_register_log_property('f_change',"-1")
  if auto_divider:
    cfs.opt_register_log_property('divider',"6")

  opt = cfs.optimizer_get_properties()
  ocm = opt['optimizer'] == 'ocm'
  mma = opt['optimizer'] == 'mma'
  snopt = opt['optimizer'] == 'snopt'
  if not ocm and not snopt:
    print('Warning: unknown optimizer used', opt['optimizer'])
  if ocm or mma:
    cfs.opt_register_log_property('transition_zone',"0")
    org_move_limit = float(opt["move_limit"])

  if adapt_move_limits is True and not (ocm or mma):
    print("warning: no optimizer selected which can handle adapt_move_limits!",opt['optimizer'])
    adapt_move_limits = False
  if adapt_move_limits:
    cfs.opt_register_log_property('move_limit',"0")
  
    
  objs, cnstrs = cfs.get_opt_function_values()
  print('objective:', list(objs.keys()))
  print('constraints and observes:', list(cnstrs.keys()))
 
  #if not 'greyness_phy'

def opt_post_iter_dunning_auto_filter(opt):
  global ocm, mma, snopt, last_obj, last_last_obj, gamma, delta_beta_max_factor, adapt_move_limits, tz_divider, auto_divider, divider, halt_on_osc, max_beta
  
  # obtain current objective value  
  objs, cnstrs = cfs.get_opt_function_values()
  assert len(objs) == 1

  f_k = float(list(objs.values())[0])
  key = 'greyness' if 'greyness' in cnstrs.keys() else [k for k in cnstrs.keys() if 'greyness' in k][0] # for robust it is 'physical_greyness_f_2' 
  g = float(cnstrs[key])
  func = cfs.get_opt_function_properties(key)
  stop = cfs.get_opt_stopping_rules()
  if 'belowFunction_greyness' in stop:
    stop_val = float(stop['belowFunction_greyness'])
  elif 'belowFunction_greyness_f_2' in stop: 
    stop_val = float(stop['belowFunction_greyness_f_2'])
  else:
    print('WARNING: no belowFunction_greyness[_f_2] found! take .01')
    stop_val = .01
      
  #print('greyness',g,'stop',stop['belowFunction_greyness'])
  
  filter = cfs.get_opt_filter_values()
  beta_k = float(filter["beta"])

  # from Dunning "Achieving binary topology optimization solutions via automatic projection parameter increase"
  k = cfs.get_opt_iteraton()
  if k > 0: 
    assert last_obj is not None
    f_p = last_obj
    delta_beta = max((-gamma/2) * ((f_k + f_p) / (f_k - f_p)), 0) if f_k != f_p else 0
    if g > stop_val:
      beta_n = beta_k + min(delta_beta, delta_beta_max_factor*beta_k)
      beta_n = min(beta_n, max_beta)
    else:
      #print('grayness',g,'<',stop_val)
      beta_n = beta_k  
    #beta_n = 1.3**(int(k/20))
  else:
    beta_n = beta_k
    f_p = -1
    delta_beta = -1 
  
  tf = int(cfs.get_opt_filter_values()['total_filters'])
  for fi in range(tf): # set for all robust filters
    cfs.set_opt_filter_values(fi, beta_n)
  cfs.opt_set_log_property("beta", str(beta_n))
  # let cfs compute scaling and offset with new beta
  filter = cfs.get_opt_filter_values()
  cfs.opt_set_log_property("delta_beta", str(delta_beta))

  cfs.opt_set_log_property('f_change',str(((f_k + f_p) / (f_k - f_p)) if f_k != f_p else 0))
  #print('iter',cfs.get_opt_iteraton(),'f_k',f_k,'f_k-1',f_p,'d_f',f_p-f_k,'f_change',(f_k + f_p) / (f_k - f_p),'rdb',beta_n/beta_k,'delta_beta',delta_beta,'-> beta', beta_n)  

  if(k > 1 and adapt_move_limits and (ocm or mma)):
    assert last_last_obj is not None and f_p is not None and f_k is not None
    f_pp = last_last_obj
  
    eta = float(filter["eta"])
    tz = atf.transition_zone_grad(.1, filter['density'], 3, beta_n, eta)
    #if filter['density'] == 'tanh': 
    #  tz = atf.transition_zone_func(.1, filter['density'], beta_n, eta)
    #else:
    #  tz = atf.transition_zone_grad(.1, filter['density'], 3, beta_n, eta)  
    
    if auto_divider:
      osc = (f_p - f_k) * (f_pp - f_p) < 0
      if osc:
        divider = min(divider + .3, tz_divider)
      else:
        divider = max(divider - .2, 6)
    else:
      divider = tz_divider       
    
    ml_org = float(cfs.optimizer_get_properties()["move_limit"])
    
    ml = str(min(tz/divider,org_move_limit)) #  
    #print('tz',tz,'cand',tz/divider,'org',org_move_limit,'->',ml)
    cfs.optimizer_set_property("move_limit",ml)

    cfs.opt_set_log_property('move_limit',ml)
    cfs.opt_set_log_property('transition_zone',str(tz))
    if auto_divider:
      cfs.opt_set_log_property('divider',str(divider))
  if halt_on_osc and k > 2:
     assert last_last_obj is not None and f_p is not None and f_k is not None
     if (f_p - f_k) * (last_last_obj - f_p) < 0:
       cfs.opt_stop(False, 'objective oscillates')
  
  #if beta_n > 256:
  #  cfs.opt_stop(False, 'maximal beta reached: ' + str(beta_n))
  #  print('reached maximal beta',beta_n) 
  
  last_last_obj = last_obj  
  last_obj = f_k