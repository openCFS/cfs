# this contains common functions for auto_filter.py, auto_filter_ocm.py and auto_filter_mma.py

import numpy as np

# this is the tanh transition zone without penalization
def transition_zone_func(eps, beta, eta):
  def f(x, beta, eta):
    return (np.tanh(beta*eta) + np.tanh(beta*(x-eta))) / (np.tanh(beta*eta) + np.tanh(beta*(1-eta)))
  l = 0
  u = 1
  while u-l > 1e-5:
    c = .5 * (u + l) 
    r_c  = f(c,beta,eta)-eps
    
    if r_c > 0:
      u = c
    else:
      l = c
  return 2*(eta-0.5 * (u + l))        

# the tanh implementation in cfs which requires non_lin_scale and non_lin_offet
def f(x,beta,eta):
  return 1-1/(np.exp(2*beta*(x-eta))+1)

# the tanh from Peter's paper which is 0 at 0 and 1 at 1
def t(x, beta, eta):
  return (np.tanh(beta*eta) + np.tanh(beta*(x-eta))) / (np.tanh(beta*eta) + np.tanh(beta*(1-eta)))

# the derivative of t
# wolframalpha
def dt(x, beta, eta):
    return beta/(np.cosh(beta*(x - eta))**2 *(np.sinh(beta*eta)/np.cosh(beta*eta) + np.sinh(beta - beta*eta)/np.cosh(beta - beta*eta)))

# wolframalpha
def inverse(y, beta, eta):
  return 1/beta * (np.arctanh((y-1)*np.tanh(beta*eta) + y * np.tanh(beta-beta*eta)) + beta*eta)

# solid heaviside as in cfs, would requice scale and offset
def sh(x, beta):
  return 1-np.exp(-beta*x) + x * np.exp(-beta)

# derivative of sh
def dsh(x,beta):
  return beta*np.exp(-beta*x) + np.exp(-beta)

# search where the derivative of tanh**3 or solid heaviside**3 is a given fraction of the maximum
# @param frac the fraction we search, e.g. 0.1 must be in (0;1]
# @param func 'tanh' or 'solid_heaviside'
# @return two values within (0;1) but one value for frac=1
def grad_val(frac, func, power, beta, eta):
  def ftanh(val):
    return power * t(val,beta,eta)**(power-1) * dt(val,beta,eta)
  
  def fsh(val):
    return power * sh(val, beta)**(power-1) * dsh(val, beta)
  
  def bisec(val, f, lower, upper):
    cnt = 0
    assert lower == 0 or upper == 1
    grow = lower == 0
    while upper-lower > 1e-5:
      c = .5 * (upper+lower)
      #print(val,cnt,'l',lower, 'u',upper,'delta',upper-lower, 'c',c, 'f',f(c), f(c) > val)
      # we need to know if we are montonous increasing or decreasing to step to the proper side
      if f(c) > val if grow else f(c) < val:
          upper = c
      else:
        lower = c
      cnt +=1
      if cnt > 100:
        break
    return .5 * (upper+lower)    
  
  assert func in ['tanh', 'solid_heaviside', 'expression']
  f = ftanh if func != 'solid_heaviside' else fsh

  # search maximum
  ls = np.linspace(0,1,300)
  samples = [f(v) for v in ls]
  pi = np.argmax(samples) # peak index
  pv = samples[pi] # peak value
  
  assert frac > 0 and frac <= 1
  if frac == 1:
    return ls[pi]
  
  #plt.plot(ls,samples)
  #plt.show()
  
  first = bisec(pv*frac, f, 0, ls[pi])
  second = bisec(pv*frac, f, ls[pi],1)
  return (first,second)  

def transition_zone_func(frac, func, beta, eta):
  assert func in ['tanh', 'expression']
  return 2*(eta-inverse(frac, beta, eta))

def transition_zone_grad(frac, func, power, beta, eta):
  first, second = grad_val(frac, func, power, beta, eta)
  return second - first

