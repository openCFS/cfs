#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Script to calculate the Electric field and 
total energy of a L-shaped domain.

Here, a Poisson problem is solved with a pre-scribed
excitation of

    u(r,phi) = r^(2/3)* 2/3 * sin(2*phi/3 + pi/3)
    
We compute symbolically the electric field strength
(to be used as Neumann-type boundary condition),
as well as the total energy.

Created on Fri Sep 19 18:49:59 2014

@author: ahauck
"""
from __future__ import division, print_function
from sympy import *
from sympy.mpmath import quad

x = Symbol("x")
y = Symbol("y")


# define permittivity
eps = 1

# size of l-shape
l=1


print("===========================")
print(" Energy for L-shaped setup ")
print("===========================")

# define potential distribution function
r = sqrt(x**2 + y**2)
phi = atan2(y,x)

#u = r**Rational(2,3) * sin(2*phi/3+pi/3)
u = r**(2.0/3.0) * sin(2*phi/3+pi/3)
print("potential function: u(x,y) =", u, "V")

# define electric field
E = Matrix((diff(u,x), diff(u,y)))

print("E_x:", E[0])
print("E_y:", E[1])

# assemble expression for energyy density
w_dens = 0.5 *( E[0]**2 + E[1]**2) * eps
w_dens= w_dens.simplify()
print("energy density: w_dens(x,y) =", w_dens, "Ws/m^2")


print("performing integration ...")

# perform firast integration in x-direction only
w_x = integrate(w_dens,(x,0,l))
W_1 = integrate(w_x,(y,0,l)).evalf()
pprint(integrate(w_x,(y,0,l)))

w_x = integrate(w_dens,(x,-l,0))
W_2 = integrate(w_x,(y,0,l)).evalf()


W_total = 2*W_2+W_1
print("")
print("--------------------------------------------")
print("Total ElecEnergy (Ws): ", W_total , "Ws")
print("Total Energy (H1-norm):", sqrt(W_total) )
print("--------------------------------------------")

