# This is a Paraview programmable filter to visualize for a FMO result the stiffest and second stiffest (if any)
# direction and the corresponding magnitude.
#
# Apply a successive CellCenters Filter then the Glyph filter can be used. The Glyph filter for the second stiffest 
# vector needs to be scaled
#
# The current CFS-HDF5 Reader does not work with Paraview programmable filter. The result needs to exported (e.g. VTK)
# first.


# visualize:
# import pylab
# pylab.plot(data[:,1])
# pylab.show()

import numpy
import numpy.linalg

from numpy import dot
from numpy import sin
from numpy import cos
from numpy import sqrt


e2d = numpy.zeros((2,2))
e2d[0,0] = 1.51
e2d[0,1] = 0.0
e2d[1,0] = 0.0
e2d[1,1] = 1.27

p0 = numpy.zeros((2,3))
p0[0,0] = 0.0
p0[0,1] = 0.0
p0[0,2] = 17.0
p0[1,0] = -6.5
p0[1,1] = 23.3
p0[1,2] = 0.0


## This rotates a 2*2 2D tensor via the third direction. As in Richter and CFS
def get_rot_2x2(angle):
  R = numpy.zeros((2,2))

  R[0][0] =  cos(angle)
  R[0][1] =  sin(angle)
  R[1][0] =  -sin(angle)
  R[1][1] =  cos(angle)
  
  return R

## This rotates a 3*3 2D tensor via the third direction. As in Richter and CFS
def get_rot_3x3(angle, R):
  
  Q = numpy.zeros((3,3))
  
  Q[0][0] = R[0][0]*R[0][0];
  Q[0][1] = R[0][1]*R[0][1];
  Q[0][2] = 2.0*R[0][0]*R[0][1];

  Q[1][0] = R[1][0]*R[1][0];
  Q[1][1] = R[1][1]*R[1][1];
  Q[1][2] = 2.0*R[1][0]*R[1][1];

  Q[2][0] = R[0][0]*R[1][0];
  Q[2][1] = R[0][1]*R[1][1];
  Q[2][2] = R[0][0]*R[1][1] + R[0][1]*R[1][0];
  
  return Q

## this performs the cfs rotation for 2D.
# works for 2*2 tensors (permitivity), 2*3 (piezo coupling) and 3*3 (elasticity in Voigt notation)
# @param tensor: we assume the cfs rotation alpha=-90 and gamma=-90 already to be done -> XML-Reference
# @return with the dimension of tensor  
def rotate_cfs(tensor, angle):

  out = numpy.zeros(tensor.shape)

  R = get_rot_2x2(angle)
  
  if tensor.shape == (2,2):
    return dot(R,dot(tensor,R.transpose()))
  
  Q = get_rot_3x3(angle, R)
  
  if tensor.shape == (2,3):
    return dot(R, dot(tensor, Q.transpose()))

  if tensor.shape == (3,3):
   return dot(Q, dot(tensor, Q.transpose()))

  assert(false) 


# performs a cfs-rotation study for elast (Voigt), elec and piezo tensors
# in the result the first two steps are repeated such the neighbor search is easier!
# @param steps: how many probes
# @return: array n * 3: first=angle, second=magnitude, third=norm for 0-positions in orthonormal 
def perform_cfs_rotation(tensor, steps, calc_norm=False):
  data = numpy.zeros((steps+2, 3))
  idx = 0
  for x in numpy.arange(0, numpy.pi, numpy.pi/steps):
    test = rotate_cfs(tensor, x)
    data[idx, 0] = x
    data[idx, 1] = test[0,0]
    if calc_norm:
      if tensor.shape == (2,2):
        data[idx, 2] = sqrt(2.0 * test[0,1]**2)
      if tensor.shape == (3,3):
        data[idx, 2] = sqrt(2.0 * (test[0,2]**2 + test[1,2]**2))
    idx += 1
    
  data[steps, 0] = data[0,0]  
  data[steps, 1] = data[0,1]
  data[steps, 2] = data[0,2]
  data[steps+1, 0] = data[1,0]  
  data[steps+1, 1] = data[1,1]
  data[steps+1, 2] = data[1,2]
  
  return data  

# finds the two largest maxima
# @param param: data vector where the last two entries shall be a double of the first ones
# @return the indices of the largest and second largest maxima. If there is no second then -1
def find_maxima(data):
  first = [-1, -1e60]
  second = [-1, -1e60]
  for i in range(1, len(data)-1):
    if data[i-1] <= data[i] and data[i] >= data[i+1]:
      if data[i] >= first[1]:
        second = first[:] # deep copy
        first[0] = i
        first[1] = data[i]
      elif data[i] >= second[1]:
        second[0] = i
        second[1] = data[i]        

  return first[0], second[0]

## @see find_maxima 
def find_minima(data):
  first = [-1, +1e60]
  second = [-1, +1e60]
  for i in range(1, len(data)-1):
    if data[i-1] >= data[i] and data[i] <= data[i+1]:
      if data[i] <= first[1]:
        second = first[:]
        first[0] = i
        first[1] = data[i]
      elif data[i] <= second[1]:
        second[0] = i
        second[1] = data[i]        

  return first[0], second[0]


# give the rotation matrix for a 2D elasticity tensor in Hill-Mandel notation
# @param angle: 0 ... PI
# @return:  a 3*3 matrix (numpy)
def rotation_matrix_mandel(angle):
   assert(angle >= 0 and angle <= numpy.pi)
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
# @return: the x-strength, orth error
def test_rotation_mandel(tensor, angle):
  assert(tensor.shape == (3,3))
  
  theta = rotation_matrix_mandel(angle)
  rot = dot(theta.transpose(),dot(tensor,theta))
  
  return rot[0,0], sqrt(rot[0,2]**2 + rot[1,2]**2)  

# finds stiffest orientation for a 2D Hill-Mandel tensor
# @param steps: how many probes
# @return: idx of opt angle, array n*2 first angle, second value
def find_stiffest_rotation_mandel(tensor, steps, dump = None):
  m = -1e60
  opt = -1
  data = numpy.zeros((steps, 2))
  err  = numpy.zeros((steps, 1))
  idx = 0
  for x in numpy.arange(0, numpy.pi, numpy.pi/steps):
    test, norm = test_rotation_mandel(tensor, x)
    if not dump == None:
      print str(test) + "\t" + str(x) + "\t" +str(norm) 
    if test > m:
      m = test
      opt = idx
    data[idx, 0] = x
    data[idx, 1] = test
    err[idx, 0] = norm
    idx += 1
  return opt, data  

## transforms Hill-Mandel to Voigt elasticity tensor 
def HillMandel2Voigt(tensor):
  ret = numpy.zeros((3,3))
  ret[0,0] = tensor[0,0]
  ret[0,1] = tensor[0,1]
  ret[0,2] = 1.0/sqrt(2.0) * tensor[0,2]
  ret[1,0] = tensor[1,0]
  ret[1,1] = tensor[1,1]
  ret[1,2] = 1.0/sqrt(2.0) * tensor[1,2]
  ret[2,0] = 1.0/sqrt(2.0) * tensor[2,0]
  ret[2,1] = 1.0/sqrt(2.0) * tensor[2,1]
  ret[2,2] = 0.5 * tensor[2,2]
  return ret

## transforms Voigt elasticity tensor to Hill-Mandel notation 
def Voigt2HillMandel(tensor):
  ret = numpy.zeros((3,3))
  ret[0,0] = tensor[0,0]
  ret[0,1] = tensor[0,1]
  ret[0,2] = sqrt(2.0) * tensor[0,2]
  ret[1,0] = tensor[1,0]
  ret[1,1] = tensor[1,1]
  ret[1,2] = sqrt(2.0) * tensor[1,2]
  ret[2,0] = sqrt(2.0) * tensor[2,0]
  ret[2,1] = sqrt(2.0) * tensor[2,1]
  ret[2,2] = 2.0 * tensor[2,2]
  return ret


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

# data = inputs[0].CellData['mechTensor'] 
# first_vec, first_val, second_vec, second_val = pv_tensor_directions(data, 90)
# output.CellData.append(first_vec,"stiffest_vec")
#output.CellData.append(first_val,"stiffest_mag")
#output.CellData.append(second_vec,"second_vec")
#output.CellData.append(second_val,"second_mag")