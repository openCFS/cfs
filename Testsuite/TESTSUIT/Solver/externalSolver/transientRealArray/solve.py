import numpy as np
from scipy.io import mmread, mmwrite
import sys

if __name__ == "__main__":
    print("Reading input values\n")
    matrix = mmread(sys.argv[2]).toarray()
    vector = mmread(sys.argv[3]).toarray()
    print(f"Solving linear system:\n")
    solution = np.linalg.solve(matrix,vector)
    print("Export solution\n")
    mmwrite(sys.argv[4],solution)
    sys.exit(0)
