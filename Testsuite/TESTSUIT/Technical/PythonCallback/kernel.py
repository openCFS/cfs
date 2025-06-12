print('hallo from kernel.py')

import cfs

def my_post_grid(opt):
  print('my_post_grid with opt',opt)    
  print('cfs.single_pde_names()',cfs.single_pde_names())


def my_post_solve_problem(opt):
  print('my_post_solve_problem with opt', opt)
  print('cfs.single_pde_names()',cfs.single_pde_names())

def coef_init(opt):
  print('coef_init opt',opt)

def coef_eval(coord):
  print('coef_eval coord',coord)
  return [0.1,.2]