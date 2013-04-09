#!/usr/bin/python

# Read in matrix in either matrix market format or 
# dense matrix format and plot colored entries.
from pylab import *
from scipy.io import *
from sbmtools import *

if len(sys.argv) != 3 and len(sys.argv) != 4:
    print "Print values of (dense / sparse) SBM matrix"
    print "Usage:\n\t ",sys.argv[0], " <matrix-file-basename> <numRows> <type=log|lin>"
    sys.exit()    
baseName = sys.argv[1]
numRows = int(sys.argv[2])
useLog = ""
if len(sys.argv) == 4:
    useLog = sys.argv[3]
    
mat = sbmRead(baseName, numRows, numRows, True)
[rowSize,colSize] = sbmInfo(baseName,numRows,numRows,True)
mat = mat.todense()

# Add first col/row for padding values with 0
[m,n] = mat.shape
if useLog == "log":
    # Calculate logarithm
    b=log10(abs(mat))
    # Remove entries with very small values (e.g. -inf due to 0 entries)
    #b[b<1e-14]=0
else:
    b = mat
    
# Reverse orientation of axis so we can visualize the array like a matrix
fig = matshow(b)

# show visually the different blocks
offset = 0
for i in arange(0,numRows-1):
    offset += rowSize[i]
    axhline(y=offset,color="k")

offset = 0    
for i in arange(0,numRows-1):
    offset += colSize[i]
    axvline(x=offset,color="k")
     

#spy(m)
#grid()
cb=colorbar()
if useLog == "log":
    cb.set_label("Exponent")
else:
    cb.set_label("Value")
show()


# ==== OLD VERSION WHICH DOES NOT COMPILE ANY MORE ===
#ax = gca()
#ax.invert_yaxis()
#title("Entries of '"+file+"' (log)")
#edgeColor="None"
##pcolormesh(b,e   dgecolors=edgeColor)

#x = arange(m)
#y = arange(n)
#X, Y = meshgrid(x,y)

#print "shape of b",b.shape
#print "shape of X",X.shape
#print "shape of Y",Y.shape
#pcolor(X,Y,b.T)
#pcolormesh(b)