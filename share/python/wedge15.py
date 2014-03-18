# -*- coding: utf-8 -*-
"""
Created on Mon Mar 25 19:05:11 2013

@author: ahauck
"""
import numpy as np
from sympy import *
import matplotlib.pylab as pl


x = Symbol("x")
y = Symbol("y")
z = Symbol("z")

# ===========================================================================
#  3D Wedge Functions   
# ===========================================================================
# define corner coordinates
coords = {
    1 : (0.0,  0.0, -1.0),
    2: (1.0,  0.0, -1.0),
    3: (0.0,  1.0, -1.0),
    4: (0.0,  0.0,  1.0),
    5: (1.0,  0.0,  1.0),
    6: (0.0,  1.0,  1.0),
    7: (0.5,  0.0, -1.0),
    8: (0.5,  0.5, -1.0),
    9: (0.0,  0.5, -1.0),
    10: (0.5,  0.0,  1.0),
    11: (0.5,  0.5,  1.0),
    12: (0.0,  0.5,  1.0),
    13: (0.0,  0.0,  0.0),
    14: (1.0,  0.0,  0.0),
    15: (0.0,  1.0,  0.0)
}

# define shape functions

# --- Original implementation ---
#funcs = {
# 1 :  0.5 * z * (1 - z) *  (1 - x - y)  * (2 * x + 2*y -1),
# 2 :  0.5 * z * (1 - z) *   x * (1 - 2 * x),
# 3 :  0.5 * z * (1 - z) *   y * (1 - 2 * y),
#
# 4 : -0.5 * z * (1 + z) * (1 - x - y) * (2*x + 2 * y -1),
# 5 : -0.5 * z * (1 + z) * x * (1 - 2 * x),
# 6 : -0.5 * z * (1 + z) * y * (1 - 2 * y),
#
# 7 : -2 * z * (1 - z) * x * (1 - x - y),
# 8 : -2 * z * (1 - z) * x * y,
# 9 : -2 * z * (1 - z) * y * (1 - x - y),
#
#10 :  2 * z * (1 + z) * x * (1 - x - y),
#11 :  2 * z * (1 + z) * x * y,
#12 :  2 * z * (1 + z) * y * (1 - x - y),
#
#13 : (1 - x - y) * (1 - z * z),
#14 : x * (1 - z * z),
#15 : y * (1 - z * z)   
#   }
   
# auxilliary variables
t = 1 -x -y
L1 =  1 -x -y
L2 = x
L3 = y

funcs = {
     1 : 0.5 * L1 * (2*L1 - 1) * (1 - z) - 0.5 * L1 * (1 - z*z),
     2 : 0.5 * L2 * (2*L2 - 1) * (1 - z) - 0.5 * L2 * (1 - z*z),
     3 : 0.5 * L3 * (2*L3 - 1) * (1 - z) - 0.5 * L3 * (1 - z*z),
    
     4 : 0.5 * L1 * (2*L1 - 1) * (1 + z) - 0.5 * L1 * (1 - z*z),
     5 : 0.5 * L2 * (2*L2 - 1) * (1 + z) - 0.5 * L2 * (1 - z*z),
     6 : 0.5 * L3 * (2*L3 - 1) * (1 + z) - 0.5 * L3 * (1 - z*z),
    
     7 : 2 * L1 * L2 * (1-z),
     8 : 2 * L2 * L3 * (1-z),
     9 : 2 * L1 * L3 * (1-z),
    
    10 : 2 * L1 * L2 * (1+z),
    11 : 2 * L2 * L3 * (1+z),
    12 : 2 * L1 * L3 * (1+z),
    
    13 : L1 * (1-z*z),
    14 : L2 * (1-z*z),
    15 : L3 * (1-z*z)
   }

# ===========================================================================
#  1D Line Functions   
# ===========================================================================

funcs1D = {
    1: 0.5 * x * (x-1),
    2: 1.0 - x*x,
    3: 0.5*x*(x+1)
 }
   
# define corner coordinates
coords1D = {
    1 : -1,
    2: 0,
    3: 1
}

def printFuncNumerical():
    for p in coords:
        print "point: ",p,":",coords[p]
        point = coords[p]
        
        for i in funcs:
            f = funcs[i]
            val = f.evalf(subs={x:point[0],y:point[1], z:point[2]},chop=True)
            print "val: ",val
            
        print
        
        
def printDerivNumerical():
    for p in coords:
        print "point: ",p,":",coords[p]
        point = coords[p]
    
        for i in funcs:
            f = funcs[i]
            grad = []
            grad.append(diff(f,x))
            grad.append(diff(f,y))
            grad.append(diff(f,z))
            
            # print number representation
            for g in grad:
                print "{:1.4e}".format(float(g.evalf(subs={x:point[0],y:point[1], z:point[2]},chop=True))),
            print
        print
    
    
def printFuncAnal():
        for i in funcs:
            f = funcs[i]
            grad = []
            grad.append(diff(f,x))
            grad.append(diff(f,y))
            grad.append(diff(f,z))
            
            # print number representation
            print "shape["+str(i-1)+"] = "+str(simplify(f))+";"
            
def printDerivAnal():
        for i in funcs:
            f = funcs[i]
            grad = []
            grad.append(diff(f,x))
            grad.append(diff(f,y))
            grad.append(diff(f,z))
            
            # print number representation
            d = 0
            for g in grad:
                print "deriv["+str(i-1)+"]["+str(d)+"] = "+str(simplify(grad[d]))+";"
                d +=1
            print
    
    
def plotVals1D():
    # print function 1 along zeta-direction
    colors = ["r","g","b"]
    funcIndex = [1,2,3]
    xi = np.linspace(-1,1,50)
    for i in range(0,3):
        v = np.zeros(xi.shape)
        vderiv = np.zeros(xi.shape)
        pos = 0
        for zet in xi:
            f = funcs1D[funcIndex[i]]
            grad = diff(f,x)
            v[pos] = float(f.evalf(subs={x:zet},chop=True))
            vderiv[pos] = float(grad.evalf(subs={x:zet},chop=True))
            pos += 1
        
        
        pl.plot(xi,v,color=colors[i])    
        pl.plot(xi,vderiv,'--',color=colors[i])
    
    pl.grid()
    pl.show()
    
    
def plotVals():
    # print function 1 along zeta-direction
    colors = ["r","g","b"]
    #funcIndex = [1,13,4]
    funcIndex=[1,13,4]
    zeta = np.linspace(-1,1,50)
    labels=[]
    for i in range(0,3):
        v = np.zeros(zeta.shape)
        vderiv = np.zeros(zeta.shape)
        pos = 0
        for zet in zeta:
            f = funcs[funcIndex[i]]
            grad = diff(f,z)
            point = coords[3]
            point = [0.25, 0.25, 0]
            v[pos] = float(f.evalf(subs={x:point[0],y:point[1], z:zet},chop=True))
            vderiv[pos] = float(grad.evalf(subs={x:point[0],y:point[1], z:zet},chop=True))
            pos += 1
        
        pl.plot(zeta,v,color=colors[i])    
        pl.plot(zeta,vderiv,'--',color=colors[i])
        labels.append("N"+str(funcIndex[i]))
        labels.append("N_deriv"+str(funcIndex[i]))
    
    pl.grid()
    pl.legend(labels)
    pl.show()
    
    
    
#plotVals1D()
#plotVals()
#printFuncNumerical()
#printDerivNumerical()
printFuncAnal()
print
print
printDerivAnal()

            
        
