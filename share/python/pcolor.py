#!/usr/bin/python

# Read in matrix in either matrix market format or 
# dense matrix format and plot colored entries.
from pylab import *
from scipy.io import *

if len(sys.argv) != 2 and len(sys.argv) != 3:
    print "Print values of (dense / sparse) matrix"
    print "Usage:\n\t ",sys.argv[0], " <matrix-file> <type=log|lin>"
    sys.exit()    
file = sys.argv[1]
useLog = ""
if len(sys.argv) == 3:
    useLog = sys.argv[2]
    
# First attempt: try to read matrix market format
try:
    print mminfo(file)
    mat = mmread(file)
    mat = mat.todense();
except ValueError:
    print "trying to read dense matrix .."
    try:
        mat = matrix(loadtxt(file))
    except ValueError:
        print "The file '"+\
        file+"' is neiter a matrix-market, nor a dense matrix file"
        sys.exit()
    

# Add first col/row for padding values with 0
[m,n] = mat.shape
#a = matrix(zeros((m+1,n+1)))
#a[1:m+1,1:n+1] = abs(mat)
if useLog == "log":
    # Calculate logarithm
    b=log10(abs(mat))
    # Remove entries with very small values (e.g. -inf due to 0 entries)
    #b[b<1e-14]=0
else:
    b = abs(mat)
    
# Reverse orientation of axis so we can visualize the array like a matrix

matshow(b)

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
