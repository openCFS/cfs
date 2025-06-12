from scipy.io import mmread, mmwrite
import scipy.sparse
import scipy.sparse.linalg
import numpy as np
import sys

# read in the arguments 
AFileName = sys.argv[1]
arg_sigma = complex(sys.argv[2])
arg_k = int(sys.argv[3])
arg_which = sys.argv[4]
arg_tol = sys.argv[5]

# import system matrices
A = mmread(AFileName)

# compute eigenvalues and corresponding eigenvectors
eigenvalues, eigenvectors = scipy.sparse.linalg.eigs(A,k=arg_k,sigma=arg_sigma,which=arg_which,tol=arg_tol)

# print result in console to test logging
print("The eigenvalues are: ", eigenvalues, "\n")
print("The corresponding eigenvectors are: ", eigenvectors, "\n")

# reshape the result arrays so they have only 1 column
eigenvalues = eigenvalues[np.newaxis].T
eigenvectors = eigenvectors.reshape((-1,1), order='F')

# export to matrix market array format
mmwrite(sys.argv[6], eigenvalues)
mmwrite(sys.argv[7], eigenvectors)

