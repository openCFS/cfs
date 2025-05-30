"""
  This script should test the interface between python and matlab - thus it's not performance optimized
  The matlab.engine must be installed via:
    - go to $MATLAB_ROOT_DIR/extern/engines/python
    - run 'python setup install'
    
  documentation: https://de.mathworks.com/help/matlab/matlab-engine-for-python.html
  
  >>> print(round(square(2.2),2))
  4.84
"""

def square(val):
  """ 
    squares a scalar or a matrix via matlab
    
    >>> str(square([[1,2,3], [4,5,6], [7,8,9]])).strip()
    '[[30.0,36.0,42.0],[66.0,81.0,96.0],[102.0,126.0,150.0]]'
  """
  
  import numpy as np
  import matlab.engine
  import sys, os
  
  if type(val) is list:
    # convert python list to matlab matrix
    val = matlab.double(val)
  
  # start engine
  eng = matlab.engine.start_matlab()
  
  # set location of matlab file - currently same location as this script
  eng.addpath(os.path.dirname(__file__),nargout=0)
  
  import io
  out = io.StringIO()
  err = io.StringIO()

  res = eng.square(val,nargout=1,stdout=out,stderr=err)
  
  eng.quit()
  
  return res

if __name__ == "__main__":
  import doctest
  from sys import exit
  
  result = doctest.testmod()
  # exit with number of failed tests as exit code
  exit(result.failed)
