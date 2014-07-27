#!/usr/bin/python

# Read in matrix in matrix market format 
# and plot sparsity pattern.
import sys
from pylab import *
from scipy.io import *
from numpy import linalg as LA

file = sys.argv[1]
print mminfo(file)
mat = mmread(file)

# Add first col/row for padding values with 0
[m,n] = mat.shape
#a = matrix(zeros((m+1,n+1)))
#a[1:m+1,1:n+1] = abs(mat.todense())
a = mat.todense()

#for i in arange(0,mat.shape[0]):
#	a[i,i] = 1.0
#	print "inverting row",i
#print a
kappa = LA.cond(a)
print "condition number is %e "%kappa
#ef = LA.eigvals(a)
#print "eigvals are",sort(ef)
savetxt("cond.txt",array([kappa]))



