print('hallo from coef.py')

import cfs

def coef_init(opt):
  assert 'center' in opt
  cs = opt['center']
  global c 
  c = eval(cs)
  assert type(c) == list
  assert len(c) == 2
  

def coef_eval(coord):
  assert type(coord) == list
  assert len(coord) == 2
  # this is the same formula as in NodelForceExpression
  return [coord[0] - c[0], coord[1] - c[1]] 
