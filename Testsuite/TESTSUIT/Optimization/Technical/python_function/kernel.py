import numpy as np

print('hallo from kernel.py')

import cfs

def vol_init(opt):
  print('kernel.py: vol_init', opt)

def vol_eval(opt):
  # we could also extract n once in init()
  n = cfs.get_opt_design_size()
  rho = np.zeros(n)
  cfs.get_opt_design_values(rho)
  #print(np.mean(rho), rho)
  return np.mean(rho)

def vol_grad(opt):
  #print('kernel.py: vol_grad')
  n = cfs.get_opt_design_size()
  grad = np.ones(n) / n
  return grad
