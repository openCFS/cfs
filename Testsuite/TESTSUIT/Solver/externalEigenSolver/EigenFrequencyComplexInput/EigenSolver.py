from scipy.io import mmread, mmwrite
import scipy.sparse
import scipy.sparse.linalg
import numpy as np
import sys

# read in the arguments 
AFileName = sys.argv[1]
BFileName = sys.argv[2]
arg_sigma = complex(sys.argv[3])
arg_k = int(sys.argv[4])
arg_tol = sys.argv[5]

# import system matrices
A = mmread(AFileName)
B = mmread(BFileName)

# compute eigenvalues and corresponding eigenvectors
eigenvalues, eigenvectors = scipy.sparse.linalg.eigs(A,k=arg_k,M=B,sigma=arg_sigma,tol=arg_tol)

# reshape the result arrays so they have only 1 column, for eigenfrequencies real values are expected
eigenvalues = eigenvalues[np.newaxis].T
eigenvectors = eigenvectors.reshape((-1,1), order='F')

# export to matrix market coordinate format
eigenvalues = scipy.sparse.coo_matrix(eigenvalues)
eigenvectors = scipy.sparse.coo_matrix(eigenvectors)
mmwrite(sys.argv[6], eigenvalues)
mmwrite(sys.argv[7], eigenvectors)

