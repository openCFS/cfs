import numpy as np
from scipy.io import mmread, mmwrite
import scipy.sparse
import sys

if __name__ == "__main__":

    # get system
    A = mmread(sys.argv[2]).toarray()
    b = mmread(sys.argv[1]).toarray()

    # solve A^T x = b
    x = np.linalg.solve(A.T,b)

    # export
    x = scipy.sparse.coo_matrix(x)
    mmwrite(sys.argv[4],x)
    sys.exit(0)