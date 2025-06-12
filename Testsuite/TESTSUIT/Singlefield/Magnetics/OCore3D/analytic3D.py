import os;
#import sys;
#sys.path.append( os.path.realpath("..") );

import numpy as np;
import matplotlib.pyplot as plt;
from matplotlib.patches import ConnectionPatch;
from params import *;

wd = os.getcwd();

# for analytic solution
I = 1.0;
r_i =r_core_inner;
r_a =r_core_outer; 
mu_air= 4*pi*1.0e-7; 
mu_Fe = 4*pi*1.0e-7*5.0e3;
# functions of the analytic solution
def H( r ):
  """
  H_r(r) for a current carrying wire in z-direction.
  Relation obtained from Ampere's Law:
    
    \oint_{C(I)} H(r)*ds = \int\limits_{0}^{2*pi} H_r*r*dphi 
    = 2*H_r*pi*r = I -> H_r(r) = I/(2*pi*r) 
  """
  return I/(2*pi*r)

def mu( r ):
  """
  Function that gives the permeability as a function of the r-coordinate to obtain:
  B_r(r) = mu(r)*H_r(r)
  """
  if r >= r_i and r <= r_a:
    return mu_Fe;
  else: return mu_air;

b = wd + "/history/B.csv-1";
d = np.loadtxt( b, dtype=float, delimiter=',', skiprows=1, ); 

r = np.sqrt( d[:,1]*d[:,1] +  d[:,2]*d[:,2] + d[:,3]*d[:,3] );
B = np.sqrt( d[:,3]*d[:,3] +  d[:,4]*d[:,4] + d[:,5]*d[:,5] );

fig, ax = plt.subplots();
ax.plot(r, B, 'o',label='numeric' );

r = np.arange( 1e-4, r_domain, 1e-3 );
mu = np.array( [ mu(i) for i in r[1:] ] );
ax.plot(r[1:], mu*H(r[1:]), label='analytic' );

ax.legend();
ax.set(xlabel='r (m)', ylabel='B (T)');
plt.show();
fig.savefig("B_analytic.pdf");
