#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Calculate the magnetic field of cylindric coil along the center line
and compute the inductance as well.

For additional references, see

    http://de.wikipedia.org/wiki/Zylinderspule

In addition this script reads in the history results and compares them.
"""

from __future__ import division
from pylab import *
from os.path import *

# ===========
#  Variables
# ===========
t  = 0.1e-3          # thickness of coil
r_inner = 5e-3       # inner radius of coil
r = r_inner + t/2    # mean radius of coil
l = 40e-3            # length of coil
mu0 = 4 * pi * 1e-7  # permeability of air
I = 1                # current
N = 2                # number of turns
A_cros = t * l       # area of coil current

print "========================="
print " FLUX DENSITY COMPARISON "
print "========================="
# calculate axial flux density at arbitrary depth
def Bz(z):
    B = I*N/(2*l) * ( 
        (l/2-z) / (sqrt( (l/2 - z)**2 + r**2) )+
        (l/2+z) / (sqrt( (l/2 + z)**2 + r**2) ) )*mu0
    return B

# average flux density
B_avg=I*N/(l*2)*mu0
print "H_avg is ",B_avg    

# analytical solution
z = linspace(0,50e-3)
B_z = Bz(z)
plot(z, B_z*1e3,'r-')
grid()
xlabel(r"$z$-distance (m)",size=15)
ylabel(r"$|B_z|$ in mT",size=15)
xlim(0,5e-2)
title("Center flux of cylindric, air-filled coil",size=16)
labels = ["Analytical"]

# load history file (if present)
if exists("history/CylindricCoilEdge-magFluxDensity-hist-1.hist"):
    b_edge_hist = loadtxt("history/CylindricCoilEdge-magFluxDensity-hist-1.hist")
    plot(b_edge_hist[:,3],b_edge_hist[:,6]*1e3,'o')
    labels.append("FEM-Edge (elem midpoint)")
    edge_err = abs((b_edge_hist[0,6]-B_z[0])  / B_z[0])
    print "B_max error for edge formulation: ",edge_err*100,"%"
    
# load interpolated file (if present)
if exists("b-line.txt"):
    b_inter = loadtxt("b-line.txt")
    plot(b_inter[:,3],b_inter[:,6]*1e3)
    labels.append("FEM-Edge (interpolated)")

legend(labels)

print
print "=========================="
print " COMPARISON OF INDUCTANCE "
print "=========================="

# Formula for "short" cylindric coils (correction factor 1.1), taken from
# M. Albach, Grundlagen der Elektrotechnik 1, 
# 2nd Edition, p. 224
L = mu0 * N**2 * (pi* r_inner**2 ) / (l+1.1*r)
print "Analytical inductance is ",L,"H"

# define files with energy
wFiles = ["history/CylindricCoilEdge-magEnergy-region-air.hist",
           "history/CylindricCoilEdge-magEnergy-region-coil.hist",
           "history/CylindricCoilEdge-magEnergy-region-core.hist"]

# loop over all files and try to load energy per region
allExist = True
wTotal = 0.0
for f in wFiles:
    if not exists(f):
        allExist = False
        break
    else:
        wTotal += loadtxt(f)[1]

if allExist:
    # as we model only one eigth of the complete coil
    L_FEM = 2 * wTotal*8 / I**2
    print "Simulated inductance is ",L_FEM,"H"
    diff = abs(L_FEM-L)/L*100
    print "Difference: %2.2f"%diff,"%"

show()
