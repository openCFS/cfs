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


# creates a 2D elasticity tensor. To the HillMandel2Voigt conversion if necessary!
def to_mech_tensor(input):
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

# creates a piezoelectric coupling tensor
def to_piezo_tensor(input):
  assert(len(input) == 6)
  # "e11", "e12", "e13", "e21", "e22", "e23";
  #    0      1      2      3      4      5
  tensor = numpy.zeros((2,3))
  tensor[0,0] = input[0]
  tensor[0,1] = input[1]
  tensor[0,2] = input[2]
  tensor[1,0] = input[3]
  tensor[1,1] = input[4]
  tensor[1,2] = input[5]
  return tensor

# creates a permittivity tensor
def to_piezo_tensor(input):
  assert(len(input) == 3)
  # "e11", "e22", "e12"
  #    0      1      2 
  tensor = numpy.zeros((2,2))
  tensor[0,0] = input[0]
  tensor[0,1] = input[2]
  tensor[1,0] = input[2]
  tensor[1,1] = input[1]
  return tensor
 
## give vector from angle
def to_vector(angle): 
  result = numpy.zeros((3))
  result[0] = cos(angle)
  result[1] = sin(angle)
  result[2] = 0
  return result

## gives the two Poisson's ratios for a given tensor which needs to be in Voigt notation in elasticity notation
# @return v21, v12
def poissons_ratio(tensor):
  c = numpy.linalg.inv(tensor)
  return (-c[0,1] / c[1,1]), (-c[1,0] / c[0,0])

## find stiffest directions for a tensor
# @param tensor elast tensor in voigt notation, elec tensor or piezo tensor
# @return first, second the scaled vectors of stiffest directions. second is 0,0,0 if there is none, v21 and v12 if desired
def find_stiffest_orientation(tensor, steps, also_poissons_ratio=False):
  data = perform_cfs_rotation(tensor, steps, also_poissons_ratio)
  idx_first, idx_second = find_maxima(data[:,1])
  
  first = to_vector(data[idx_first,0]) * data[idx_first,1]
  
  second = [0,0,0]
  if idx_second > -1: 
    second = to_vector(data[idx_second,0]) * data[idx_second,1]

  if also_poissons_ratio:
    # find the one corresponding to the stronger direction
    # note, we might have zeros not in the strongest direction (-> PZT-5A!!)
    idx_first, idx_second = find_minima(data[:,2])
    idx = idx_first
    if data[idx_second,1] > data[idx_first, 1]:
      idx = idx_second
    v21, v12 = poissons_ratio(rotate_cfs(tensor, data[idx,0]))  
    return first, second, v21, v12
  else:
    return first, second


## takes takes the paraview tensor data and generates a set of vectors of orientations
# @return array of first and secondary stiffness tensors, v21 and v12
def pv_mech(input, steps):
  el = len(input) 

  first  = zeros((el, 3))
  second = zeros((el, 3))
  v21    = zeros((el, 3))
  v12    = zeros((el, 3))
  
  for i in range(el):
    t = input[i]
    first[i], second[i], v21[i], v12[i] = find_stiffest_orientation(HillMandel2Voigt(to_mech_tensor(t.transpose())), steps, True)

  return first, second, v21, v12

## process the paraview stuff
def paraview_show():
  data = inputs[0].CellData['mechTensor'] 
  first, second, v21, v12 = pv_mech(data, 90)
  output.CellData.append(first,"mech_stiffest")
  output.CellData.append(second,"mech_second")
  output.CellData.append(v21,"mech_v21")
  output.CellData.append(v12,"mech_v12")
  
  
# paraview_show()  
