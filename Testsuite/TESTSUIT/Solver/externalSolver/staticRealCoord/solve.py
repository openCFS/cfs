import numpy as np
from scipy.io import mmread, mmwrite
import scipy.sparse
import sys

if __name__ == "__main__":
   print("Reading input values\n")
   matrix = mmread(sys.argv[1]).toarray()
   vector = mmread(sys.argv[2]).toarray()
   print(f"Solving linear system:\n")
   solution = np.linalg.solve(matrix,vector)
   solution = scipy.sparse.coo_matrix(solution)
   print("Export solution\n")
   mmwrite(sys.argv[3],solution)
   sys.exit(0)
