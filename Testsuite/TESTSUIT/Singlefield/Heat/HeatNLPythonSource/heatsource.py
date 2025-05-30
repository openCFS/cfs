print('hallo from heatsource.py')

import cfs
import numpy as np
from numpy import exp # to be compatible with gnuplot
try:
  import matplotlib.pyplot as plt
except ImportError:
  print('cannot import matplotlib.pyplot, the option "plot" with fail')

def coef_init(opt):
  global c, r, fac, beta, switch_temp, do_plot, do_center_output, sol_type
  fac = eval(opt['nominal']) # the nonimal heat source vlaue
  c = np.array(eval(opt['center'])) # center of the source
  assert len(c) == 2
  r = eval(opt['radius']) # radius fo the source
  switch_temp = eval(opt['switch_temp']) # temperature where the nonlinar source roughly decreases
  beta = eval(opt['beta']) # sharpness of the nonlinear function
  do_plot = 'plot' in opt
  do_center_output = 'center_output' in opt


  sps = cfs.single_pde_solutions(0)
  sol_type = sps[0][0]

  fe = cfs.fe_function_names()
  print('FE',fe, 'fac',fac)
 
num_eqns = None  

# this is an array of nodal rhs_return value, the damping value. By each assemle_rhs() the data is shifted
cur_rhs = None
prev_rhs = None
prev_prev_rhs = None


# called by a hook in Assemble::AssembleRHSLinForms()
def assemble_rhs(opt):
  
  global num_eqns, cur_rhs, prev_rhs, prev_prev_rhs, sol

  if prev_rhs is not None:
    prev_prev_rhs = np.copy(prev_rhs)
  if cur_rhs is not None:
    prev_rhs = np.copy(cur_rhs)

  if not num_eqns: # coef_init() is called too early
    neqns = cfs.fe_function_total_equations()
    assert len(neqns) == 1
    num_eqns = neqns[0]
    print('we have', num_eqns, 'equations')

    cur_rhs = np.zeros((num_eqns,2)) # heat and damping
    sol = np.zeros(num_eqns) # just for output

  assert num_eqns != None  
  
  old = np.copy(sol)
  cfs.fe_function_solution(0,True,sol)
  print('assemble_rhs: new sol min',np.mean(sol),'mean',np.mean(sol),'diff',np.linalg.norm(sol - old))  

# get temperature of current solution by node number. One could also implement a direct callable function and
# skip the solution and mapping handling
def get_temp_by_node(node):
  # the global solution sol and mapping is quite complex :(
  ns = cfs.fe_function_nodal_solution(node,0,sol_type)
  return ns[0] 
  
# nonlinar temp to heat map
def heat(x):
  global switch_temp
  global fac
  global beta
  #print(switch_temp,fac,beta)
  return fac * max((1-1/(1+exp(-beta * (x - switch_temp)))),0) 

hist = []

# give back the damped head value and its damping value based on the history
def damping(node, h_cur):
  global prev_rhs, prev_prev_rhs

  if prev_rhs is None:
    return h_cur, 0.0
    
  h_prev = prev_rhs[node-1,0]  
  if prev_prev_rhs is None:
    return (h_prev + h_cur) / 2.0, .5

  assert prev_prev_rhs is not None
  h_prev_prev = prev_prev_rhs[node-1,0]
  d_prev = prev_rhs[node-1,1]
  oscillate = (h_cur - h_prev) * (h_prev - h_prev_prev) < 0

  d = min(d_prev*1.2, .9) if oscillate else d_prev * .8

  return (1-d) * h_cur + d * h_prev, d  


def show_plot():
  max_t = max(np.max(hist),100)
  x = np.linspace(0, 1.2 * max_t, 100)
  y = [heat(v) for v in x]
  plt.plot(x,y)
  plt.plot(hist,[heat(x) for x in hist],'bo')
  plt.show()


center_hist = []

def coef_eval(number):
    
  # get coord by number. We could ask for the coord directly by leaving/setting the xml argument by_coord to true  
  coord = np.empty(2) 
  cfs.grid_node_coord(number,coord)

  temp = get_temp_by_node(number)
  h_cur = heat(temp)

  h, d = damping(number, h_cur)    
  cur_rhs[number-1,0] = h  
  cur_rhs[number-1,1] = d

  if do_center_output and list(coord) == [.5, .5]:
    prev_d = 0 if prev_rhs is None else prev_rhs[number-1][1]
    center_hist.append((h_cur-h, d - prev_d))
    if prev_rhs is None:
      print('curr',cur_rhs[number-1], 'org',h_cur,center_hist[-1])
    elif prev_prev_rhs is None:
      print('prev',prev_rhs[number-1],'curr',[h, d], 'org',h_cur,center_hist[-1])
    else:
      print('prev_prev',prev_prev_rhs[number-1], 'prev',prev_rhs[number-1],'curr',[h, d], 'org',h_cur,center_hist[-1])
    
    #print('node',number,'temp',temp,'heat',h,'fac',fac)
    global hist
    hist.append(temp)
    if do_plot:
      show_plot()

  loc = max(r - np.linalg.norm(c-coord),0)/r
  #res = h if loc > 0 else 0
  res = loc * h
  #print(r, c, coord,np.linalg.norm(c-coord), r - np.linalg.norm(c-coord),res)
  return res
