#!/usr/bin/python
# -*- coding: utf-8 -*-
"""

Calculate analytical capacitance of cylindric capacitor
with two dielectrica.
If exisiting, we plot also the simulated electric field
over the radius.

@author: ahauck
"""

# ensure to interpret every number as float
from __future__ import division
from math import *
from pylab import *
from os.path import *

# ------------------
#  Veriables 
# -------------------
r1=1e-3
r2=2e-3
r3=5e-3
l=1e-2
eps1=1.992E-11  # permittivity of PE
eps2=8.859E-12  # permittivity of air 
U=1             # <!-- potential on hot electrode


# Calculate capacitance
C = 2*pi*l / ( 1.0/eps1 * log(r2/r1) + 1.0/eps2 * log(r3/r2))
print "C is %g"%C
if exists("history/CylCap3D-elecEnergy-region-dielec-1.hist")\
and exists("history/CylCap3D-elecEnergy-region-dielec-2.hist"):
    wSim1=loadtxt("history/CylCap3D-elecEnergy-region-dielec-1.hist")[1]
    wSim2=loadtxt("history/CylCap3D-elecEnergy-region-dielec-2.hist")[1]
    # Attention: As we model only 1/4 of the whole setup, we have to 
    # multiply the energy / capacitance by a factor of 4.
    w_total = 4 * (wSim1 + wSim2)
    cSim = 2 * w_total / U**2
    diff = abs(cSim-C)/C*100
    print("C simulated is {0} ({1:g} %error)".format(cSim,diff))
    


# define function for calculatinf the analytical field
def E(r):
    # distinguish if we are in eps1 or eps2
    if (r < r2):
        return U/(eps1*r)*1/(1/eps1*log(r2/r1)+1/eps2*log(r3/r2))
    else:
        return U/(eps2*r)*1/(1/eps1*log(r2/r1)+1/eps2*log(r3/r2))    


r = linspace(r1,r3,50)
field = []
for actR in r:
    field.append(E(actR))

# plot analytical solution
plot(r*1e3,array(field),'-g')
labels = ["analytical"]

# if the field-interolation file is present, we can read this as well
# and compare it
if exists("field-line-1.txt"):
    eSim = loadtxt("field-line-1.txt")
    labels.append("simulated")
    # plot E(x) vs. x
    plot(sqrt(eSim[:,1]**2+eSim[:,2]**2)*1e3, sqrt(eSim[:,4]**2+eSim[:,5]**2),'ro-')

if exists("field-line-2.txt"):
    eSim = loadtxt("field-line-2.txt")
    labels.append("simulated")
    # plot E(x) vs. x
    plot(sqrt(eSim[:,1]**2+eSim[:,2]**2)*1e3, sqrt(eSim[:,4]**2+eSim[:,5]**2),'bd-')


axvline(r2*1e3,color="k")
xlabel("radius (mm)")
ylabel(r"$|E|$ (V/m)")
title("Radial electric field at half height")
legend(labels)
grid()
show()
