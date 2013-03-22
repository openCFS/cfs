#!/usr/bin/python

# Read in matrix in matrix market format 
# and plot sparsity pattern.
from pylab import *
from scipy.io import *

file = sys.argv[1]
print mminfo(file)
mat = mmread(file)

# Add first col/row for padding values with 0
[m,n] = mat.shape
#a = matrix(zeros((m+1,n+1)))
#a[1:m+1,1:n+1] = abs(mat.todense())

# Reverse orientation of axis so we can visualize the array like a matrix

# set the locations of the xticks
spy(mat)
ax = gca()
ylim(-0.5,m-0.5)
xlim(n-0.5,-0.5)
#set_offset_position()
#yticks(0.5+arange(0,m,1),arange(1,m+1,1))
#xticks(-0.5+arange(0,n+1,1),arange(1,n+1,1))
yticks(arange(0,m,1),arange(m,0,-1))
xticks(arange(0,n+1,1),arange(n,0,-1))

print mat
grid(True,which="major")
show()

