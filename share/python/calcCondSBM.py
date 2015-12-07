#!/usr/bin/python

# Read in matrix in matrix market format 
# and plot sparsity pattern.
import sys
#from pylab import *
from scipy.sparse import *
from numpy import linalg as LA
from sbmtools import *

baseName = sys.argv[1]
n = int(sys.argv[2])
isSymmetric = True
A = sbmRead(baseName,n,n,isSymmetric)
a = A.todense()
kappa = LA.cond(a)
print "condition number is %e"% kappa
#ef = LA.eigvals(a)
#print "eigvals are",sort(ef)
savetxt("cond.txt",array([kappa]))



