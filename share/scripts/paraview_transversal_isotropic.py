# This is a Paraview programmable filter to visualize for a FMO result the stiffest and second stiffest (if any)
# direction and the corresponding magnitude.
#
# Apply a successive CellCenters Filter then the Glyph filter can be used. The Glyph filter for the second stiffest 
# vector needs to be scaled
#
# The current CFS-HDF5 Reader does not work with Paraview programmable filter. The result needs to exported (e.g. VTK)
# first.


import numpy
import numpy.linalg

from numpy import dot
from numpy import sin
from numpy import cos
from numpy import sqrt

# give the rotation matrix for a 2D elasticity tensor in Hill-Mandel notation
# @param angle: 0 ... PI
# @return:  a 3*3 matrix (numpy)
def rotation_matrix_mandel(angle):
   #assert(angle >= 0 and angle <= numpy.pi)
   result = numpy.zeros((3,3))
   
   result[0,0] = cos(angle)**2
   result[0,1] = sin(angle)**2
   result[0,2] = -1.0 * sqrt(2.0)/2.0 * sin(2.0 * angle)
   result[1,0] = sin(angle)**2
   result[1,1] = cos(angle)**2
   result[1,2] = sqrt(2.0)/2.0 * sin(2.0 * angle)
   result[2,0] = sqrt(2.0)/2.0 * sin(2.0 * angle)
   result[2,1] = -1.0 * sqrt(2.0)/2.0 * sin(2.0 * angle)
   result[2,2] = cos(2.0*angle)
   
   return result


# tests the strengs in x-direction of 2D Hill-Mandel tensor rotated by an angle
# @param tensor: a 2D Hill-Mandel tensor
# @param angle: 0 .. PI
# @return: the x-strength
def test_rotation_mandel(tensor, angle):
  assert(tensor.shape == (3,3))
  
  theta = rotation_matrix_mandel(angle)
  rot = dot(theta.transpose(),dot(tensor,theta))
  
  return rot[0,0]

# finds stiffest orientation for a 2D Hill-Mandel tensor
# @param steps: how many probes
# @return: idx of opt angle, array n*2 first angle, second value
def find_stiffest_rotation_mandel(tensor, steps, dump = None):
  m = -1e60
  opt = -1
  data = numpy.zeros((steps, 2))
  idx = 0
  for x in numpy.arange(0, numpy.pi, numpy.pi/steps):
    test = test_rotation_mandel(tensor, x)
    if not dump == None:
      print str(test) + "\t" + str(x) 
    if test > m:
      m = test
      opt = idx
    data[idx, 0] = x
    data[idx, 1] = test
    idx += 1
  return opt, data  

def to_tensor(input):
  assert(len(input) == 6)
  # "e11", "e22", "e33", "e23", "e13", "e12";
  #    0      1      2      3      4      5
  tensor = numpy.zeros((3,3))
  tensor[0,0] = input[0]
  tensor[0,1] = input[5]
  tensor[0,2] = input[4]
  tensor[1,0] = input[5]
  tensor[1,1] = input[1]
  tensor[1,2] = input[3]
  tensor[2,0] = input[4]
  tensor[2,1] = input[3]
  tensor[2,2] = input[2]
  
  return tensor
 
 
## give vector from angle
def to_vector(angle): 
  result = numpy.zeros((3))
  result[0] = cos(angle)
  result[1] = sin(angle)
  result[2] = 0
  return result
  
  

## for a paraview tensor in vector notation gives an vector with the stiffest orientation
# @param input: a vector of 6 components
# @param steps: the resolution for find the angle
# @return stiffest vector, stiffest value, secondary vector, secondary value
def fmo_stiffest_orientations(input, steps = 200):

  tensor = to_tensor(input)
  # idx is the stiffest direction
  idx, data = find_stiffest_rotation_mandel(tensor, steps)
  
  # find secondary direction, test lower
  inflection = find_inflection(data, idx, 0)
  if inflection > -1:
    scnd = find_max(data, inflection, 0)
  else:
    inflection = find_inflection(data, idx, len(data)-1)
    scnd = find_max(data, inflection, len(data)-1)
  
  second_vec = []
  second_val = 0
  
  # handle monotropic material
  if scnd > -1:
    second_vec = to_vector(data[scnd, 0])
    second_val = data[scnd, 1]
  else:
    second_vec = [0, 0, 0]
  
  return to_vector(data[idx, 0]), data[idx,1], second_vec, second_val

## find inflection from a given starting point in a given direction, if it exists!!
# @param data from find_stiffest_rotation_mandel
# @param start idx of the global maximum
# @param end either 0 or len(data)-1
# @return the index or -1 if not exits
def find_inflection(data, start, end):
  step = 1
  if end < start:
    step = -1
  for i in range(start, end, step):
    if data[i + step, 1] > data[i, 1]:
      return i + step
  return -1    

## searches for the second maximum
# @return the index within data
def find_max(data, start, end):
  m_v = -1e30
  m_i = -1
  
  step = 1
  if end < start:
    step = -1
    
  for i in range(start, end, step):
    test = data[i, 1]
    # check only if the next is smaller
    if test > data[i + step, 1]:
      if test > m_v:
        m_v = test
        m_i = i  

  return m_i


## reads the parametrized data
def readParams():
  dens = readParam('design_density')
  Ex = readParam('design_emodul-iso')
  Ey = readParam('design_emodul')
  Gxy = readParam('design_gmodul')
  rotAngle = readParam('design_rotAngle')
  
  return dens, Ex, Ey, Gxy, rotAngle


## reads the parameter param
# @return the filtered design, if the value exists. Otherwise the plain design
def readParam(param):
  name = param
  name += '_smart'
  read = inputs[0].CellData[name]
  if read is None:
    name = param
    name += '_plain'
    read = inputs[0].CellData[name]

  return read


## computes the tensor entries from the parameters
def computeTensors(penalty, theta, dens, Ex, Ey, Gxy, rotAngle):
  el = len(dens)
  tensordata = numpy.zeros((el, 6))
  for i in range(el):
    tensordata[i,0] = Ex[i]/(1.0-theta)
    tensordata[i,1] = Ey[i]/(1.0-theta)
    tensordata[i,2] = Gxy[i]
    tensordata[i,5] = sqrt(Ex[i]*Ey[i]*theta)/(1.0-theta)
    tensor = to_tensor(tensordata[i,:])
    Theta = rotation_matrix_mandel(rotAngle[i])
    rot = dot(Theta.transpose(),dot(tensor,Theta))    
    tensordata[i,0] = rot[0,0]
    tensordata[i,1] = rot[1,1]
    tensordata[i,2] = rot[2,2]
    tensordata[i,3] = rot[1,2]
    tensordata[i,4] = rot[0,2]
    tensordata[i,5] = rot[0,1]
  
  tensordata *= dens**penalty
  return tensordata
  

## computes the tensor entries from the parameters
def computeTensor(penalty, theta, dens, Ex, Ey, Gxy, rotAngle):
  tensor = numpy.zeros((3,3))
  tensor[0,0] = Ex/(1.0-theta)
  tensor[1,1] = Ey/(1.0-theta)
  tensor[2,2] = Gxy
  tensor[0,1] = sqrt(Ex*Ey*theta)/(1.0-theta)
  tensor[1,0] = tensor[0,1]
  Theta = rotation_matrix_mandel(rotAngle)
  rot = dot(Theta.transpose(),dot(tensor,Theta))    
  rot *= dens**penalty
  return rot
  
  
## takes takes the paraview tensor data and generates a set of vectors and the magnitudes
# @return array of normalized vector, vector of original length
def pv_tensor_directions(input, steps):
  el = len(input) 

  first_vec = zeros((el, 3))
  first_val = zeros((el))
  second_vec = zeros((el, 3))
  second_val = zeros((el))
  
  for i in range(el):
    t = input[i]
    first_vec[i], first_val[i], second_vec[i], second_val[i] = fmo_stiffest_orientations(t.transpose(), steps)

  return first_vec, first_val, second_vec, second_val


#***********************************************************************************

def main():
  data = inputs[0].CellData['mechTensor']
  if data is None:
    dens, Ex, Ey, Gxy, rotAngle = readParams()
    penalty = 3.0
    theta = 0.09
    data = computeTensors(penalty, theta, dens, Ex, Ey, Gxy, rotAngle)
  
  first_vec, first_val, second_vec, second_val = pv_tensor_directions(data, 90)
  output.CellData.append(first_vec,"stiffest_vec")
  output.CellData.append(first_val,"stiffest_mag")
  output.CellData.append(second_vec,"second_vec")
  output.CellData.append(second_val,"second_mag")
  
main()
