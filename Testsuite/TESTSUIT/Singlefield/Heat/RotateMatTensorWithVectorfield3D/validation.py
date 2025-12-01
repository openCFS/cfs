import numpy as np
from numpy.linalg import norm
from scipy.spatial.transform import Rotation
import ast
import pandas as pd
from pathlib import Path

# ------------------
# Information regarding the coeffunctionRot.txt output!
# Each row accounts for one rotation.
# The columns contain the data as described in EXPECTED_COLS
# ------------------
EXPECTED_COLS = ["DBG_name","elemNum", "refVec", "targetVec", "T", "R", "Trot"]

def parse_vector_matrix(s: str) -> np.ndarray:
    """
    Parses a vector string formatted like: "[1, 0, 0]" or "[1.0, 2.0]".
    Parses a matrix string formatted like:
      "[[a, b], [c, d]]" (2x2) or
      "[[a, b, c], [d, e, f], [g, h, i]]" (3x3)

    :param s: string to be parsed
    """
    return np.array(ast.literal_eval(str(s).strip()))

def load_dbg_file(path):
    """
    Loads the debug .txt file with the above specified structure 
    and returns a dictionary containing the exptected cols (cols_to_parse) as keys! 
    
    :param path: path to the log file, str
    """
    cols_to_parse = ["elemNum", "refVec", "targetVec", "T", "R", "Trot"]
    # result dictionary contains all relevant information as a numpy arrays
    rotResult = {}
    # load the string into a pandas dataframe
    df = pd.read_csv(path, sep=';', header=None)
    for i, col in enumerate(cols_to_parse):     # iterate over each relevant column
        # read the vector or matrix to get a numpy array
        parsed_series = df.iloc[:, i+1].apply(parse_vector_matrix)
        rotResult[col] = np.stack(parsed_series.to_numpy())
    return rotResult

def calc_rotation(dbgResult: np.array) -> None:
    """
    Calculates each rotation using a scipy function ( scipy.spatial.transform.Rotation.align_vectors() ).

    :param dbgResult: parsed and structured results, dict
    """
    # get the initial material tensor
    T = dbgResult['T']
    num_rows = len(T)
    # load dimension of the problem
    dimension = np.size(T[0], 0)
    calcRotMatrix = np.zeros((num_rows, dimension, dimension))
    calcRotMatTensor = np.zeros_like(calcRotMatrix)
    for i in range(num_rows):   # iterate over all the rotations
        # load the vectors in 3d for scipy
        v1 = np.zeros(3)
        v2 = np.zeros(3)
        v1[0:dimension] = dbgResult['refVec'][i]
        v2[0:dimension] = dbgResult['targetVec'][i]
        # Rotation that aligns v1 → v2
        rot, _ = Rotation.align_vectors([v2], [v1])
        # transform to matrix shape and calculate the problem depending on the problem dimension (2D or 3D)
        calcRotMatrix[i,:,:] = rot.as_matrix()[:dimension, :dimension]
        calcRotMatTensor[i,:,:] = calcRotMatrix[i, :, :] @ T[i] @ calcRotMatrix[i].T    

    return calcRotMatrix, calcRotMatTensor

def validate_rotation(dbgResult: np.array, max_report: int = 20):
    """
    Compares the cfs results as to the calculated results using scipy.
    Outputs the relative errors, sorted by magnitude (largest first)
    
    :param dbgResult: parsed and structured results, dict
    :param max_report: defines the number printed relative errors, int
    """
    # calc scipy rotation
    calcRotMatrix, calcRotMatTensor = calc_rotation(dbgResult)
    # load cfs results
    cfsRotMatrix = dbgResult['R']
    cfsRotMatTensor = dbgResult['Trot']
    # get number of rows and preallocate the relative error arrays
    num_rows = len(cfsRotMatrix)
    rel_err_rotMatrix = np.zeros(num_rows)
    rel_err_rotMatTensor = np.zeros(num_rows)
    # calculation of the errors, taking the largest elementwise-calculated error of each matrix
    for k in range(num_rows):
        rel_err_rotMatrix[k] = np.max(np.abs(cfsRotMatrix[k,:,:] - calcRotMatrix[k,:,:]) / (np.abs(cfsRotMatrix[k,:,:])+1e-30))
        rel_err_rotMatTensor[k] = np.max(np.abs(cfsRotMatTensor[k,:,:] - calcRotMatTensor[k,:,:]) / (np.abs(calcRotMatTensor[k,:,:])+1e-30))
    # sorting the arrays
    rel_err_rotMatrix = np.sort(rel_err_rotMatrix)[::-1]
    rel_err_rotMatTensor = np.sort(rel_err_rotMatTensor)[::-1]

    print("Rotation validation summary")
    print(f"  Rows loaded:            {num_rows}")
    print('Errors by size (largest first):')
    for i in range(max_report):
        print(f"Rel error of the rotation Matrix number {i+1}: {rel_err_rotMatrix[i]}")
        print(f"Rel error of the rotated Material Tensor number {i+1}: {rel_err_rotMatTensor[i]}")

if __name__ == "__main__":
    
    # path to the logging output file
    PathToDebugFile = root_dir = Path(__file__).resolve().parent / r'coeffunctionRot.txt'
    dbgResult = load_dbg_file(PathToDebugFile)
    validate_rotation(dbgResult)
    

