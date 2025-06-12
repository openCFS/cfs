import numpy as np
from scipy.io import mmread, mmwrite
import scipy.sparse
import sys

if __name__ == "__main__":
    step = sys.argv[1]
    #print("Reading input values\n")
    matrix = mmread(sys.argv[2]).toarray()
    vector = mmread(sys.argv[3]).toarray()
    #print(f"Solving linear system:\n")
    solution = np.linalg.solve(matrix,vector)
    #print("Export solution\n")
    solution = scipy.sparse.coo_matrix(solution)
    mmwrite(sys.argv[4],solution)
    sys.exit(0)
