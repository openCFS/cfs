# The file has been copied from paraview_fmo.py which is now depreciated.

# visualize:
# import pylab
# pylab.plot(data[:,1])
# pylab.show()

import sys
import numpy as np
import numpy.linalg

from numpy import dot
from numpy import sin
from numpy import cos
from numpy import sqrt

## print tensor nicely
def dump_tensor(tensor, toString=False):
  tensor = np.round(tensor,15)
  out = ""
  for y in range(tensor.shape[0]):
    msg = '{:2}: '.format(y+1)
    for x in range(tensor.shape[1]):
      msg += '{:10.4g} '.format(tensor[y][x])
    if toString:
      out += msg + "\n"
    else:
      print(msg)

  return out

## This rotates a 2x2 2D tensor via the third direction. As in Richter and CFS
## WARNING!!!! This rotates in clockwise direction and not counter-clockwise as in CFS!!!
def get_rot_2x2(angle):
  R = np.zeros((2,2))

  R[0][0] =  cos(angle)
  R[0][1] =  -sin(angle)
  R[1][0] =  sin(angle)
  R[1][1] =  cos(angle)

  return R

## This rotates a 3x3 2D material tensor in Voigt notation via the third direction. 
## Richter rotates clockwise and CFS counter-clockwise (default). 
## The rotation direction of this matrix depends on the direction of get_rot_2x2(angle)
def get_rot_3x3(angle):
  
  R = get_rot_2x2(angle)

  Q = np.zeros((3,3))

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


## This is a variant of DesignMaterial::SetOneAxisRotationMatrix
#  see Wikipedia: Kardan/Tait–Bryan angles, x-y-z (extrinsic rotation)
# @param axis is either 'x', 'y' or 'z'
# @param angle in rad ccw
def get_rot_3x3_3d_one_axis(axis, angle):
  if axis == 'x':
    return np.array([ [1,0,0], [0,cos(angle),-sin(angle)], [0,sin(angle),cos(angle)] ])
  if axis == 'y':
    return np.array([ [cos(angle),0,sin(angle)], [0,1,0], [-sin(angle),0,cos(angle)] ])
  if axis == 'z':
    return np.array([ [cos(angle),-sin(angle),0], [sin(angle),cos(angle),0], [0,0,1] ])
  
  assert False, 'inavild axis value'

  
## This rotates in 3D counterclockwise around three axes.
# @param axes 'xyz' rotates first by gamma around z and last by alpha around x and is default 
#  see Wikipedia: Kardan/Tait–Bryan angles, x-y-z (extrinsic rotation)
def get_rot_3x3_3d(alpha, beta, gamma, axes = 'xyz'):
  assert len(axes) == 3
  R1 = get_rot_3x3_3d_one_axis(axes[0], alpha) # x for axis = 'xyz' 
  R2 = get_rot_3x3_3d_one_axis(axes[1], beta)  # y for axis = 'xyz'
  R3 = get_rot_3x3_3d_one_axis(axes[2], gamma) # z for axis = 'xyz'

  return (R1.dot(R2)).dot(R3) # note (R1 @ R2) @ R3 == R1 @ (R2 @ R3)
  
  
## This rotates a 6x6 3D tensor counterclockwise around last axis by gamma,
#  then midde axis by beta and last first axis by alpha. By default 'xyz' this
#  is first gamma around z, then beta around y and last alpha around x
def get_rot_6x6(alpha, beta, gamma, axes = 'xyz'):

  R = get_rot_3x3_3d(alpha, beta, gamma, axes)

  Q = np.zeros((6,6))

  Q[0][0] = R[0][0]*R[0][0]
  Q[0][1] = R[0][1]*R[0][1]
  Q[0][2] = R[0][2]*R[0][2]
  Q[0][3] = 2.0*R[0][1]*R[0][2]
  Q[0][4] = 2.0*R[0][0]*R[0][2]
  Q[0][5] = 2.0*R[0][0]*R[0][1]

  Q[1][0] = R[1][0]*R[1][0]
  Q[1][1] = R[1][1]*R[1][1]
  Q[1][2] = R[1][2]*R[1][2]
  Q[1][3] = 2.0*R[1][1]*R[1][2]
  Q[1][4] = 2.0*R[1][0]*R[1][2]
  Q[1][5] = 2.0*R[1][0]*R[1][1]

  Q[2][0] = R[2][0]*R[2][0]
  Q[2][1] = R[2][1]*R[2][1]
  Q[2][2] = R[2][2]*R[2][2]
  Q[2][3] = 2.0*R[2][1]*R[2][2]
  Q[2][4] = 2.0*R[2][0]*R[2][2]
  Q[2][5] = 2.0*R[2][0]*R[2][1]

  Q[3][0] = R[1][0]*R[2][0]
  Q[3][1] = R[1][1]*R[2][1]
  Q[3][2] = R[1][2]*R[2][2]
  Q[3][3] = R[1][1]*R[2][2] + R[1][2]*R[2][1]
  Q[3][4] = R[1][0]*R[2][2] + R[1][2]*R[2][0]
  Q[3][5] = R[1][0]*R[2][1] + R[1][1]*R[2][0]

  Q[4][0] = R[0][0]*R[2][0]
  Q[4][1] = R[0][1]*R[2][1]
  Q[4][2] = R[0][2]*R[2][2]
  Q[4][3] = R[0][1]*R[2][2] + R[0][2]*R[2][1]
  Q[4][4] = R[0][0]*R[2][2] + R[0][2]*R[2][0]
  Q[4][5] = R[0][0]*R[2][1] + R[0][1]*R[2][0]

  Q[5][0] = R[0][0]*R[1][0]
  Q[5][1] = R[0][1]*R[1][1]
  Q[5][2] = R[0][2]*R[1][2]
  Q[5][3] = R[0][1]*R[1][2] + R[0][2]*R[1][1]
  Q[5][4] = R[0][0]*R[1][2] + R[0][2]*R[1][0]
  Q[5][5] = R[0][0]*R[1][1] + R[0][1]*R[1][0]

  return Q

## this rotates a tensor counterclockwise around z-axis by theta (and y-axis by phi in 3D)
# works for 2*2 tensors (permitivity), 2*3 (piezo coupling) and 3*3 (elasticity in Voigt notation)
# @param tensor: we assume the cfs rotation alpha=-90 and gamma=-90 already to be done -> XML-Reference
# @return with the dimension of tensor
def rotate_cfs(tensor, theta, phi = None):

  out = np.zeros(tensor.shape)

  # 3D only for pure elasticity
  if tensor.shape == (6,6):
    Q = get_rot_6x6(0.0, phi, theta)
    #print "theta=" + str(theta) + " phi=" + str(phi) + " -> " + str(dot(Q, dot(tensor, Q.transpose())))
    return dot(Q, dot(tensor, Q.transpose()))

  # 2D for elec, piezo and mech
  R = get_rot_2x2(theta)
  if tensor.shape == (2,2):
    return dot(R,dot(tensor,R.transpose()))

  Q = get_rot_3x3(theta)
  if tensor.shape == (2,3):
    return dot(R, dot(tensor, Q.transpose()))

  if tensor.shape == (3,3):
   return dot(Q, dot(tensor, Q.transpose()))

  assert(False)


## this performs a Hill-Mandel 2D elasticity tensor rotation.
# Hill-Mandel rotation is trace invariant
def rotate_hill_mandel(tensor, theta):

  assert(theta >= 0 and theta <= np.pi)

  out = np.zeros(tensor.shape)

  Q = np.zeros((3,3))
  Q[0,0] = cos(theta)**2
  Q[0,1] = sin(theta)**2
  Q[0,2] = -1.0 * sqrt(2.0)/2.0 * sin(2.0 * theta)
  Q[1,0] = sin(theta)**2
  Q[1,1] = cos(theta)**2
  Q[1,2] = sqrt(2.0)/2.0 * sin(2.0 * theta)
  Q[2,0] = sqrt(2.0)/2.0 * sin(2.0 * theta)
  Q[2,1] = -1.0 * sqrt(2.0)/2.0 * sin(2.0 * theta)
  Q[2,2] = cos(2.0*theta)

  return dot(Q.transpose(),dot(tensor,Q))

def test_rotation(tensor, steps, notation):

  res = []

  for x in np.arange(0, np.pi, np.pi/steps):
    test = 0
    if notation == "mandel":
      test = rotate_hill_mandel(tensor, x)
    else:
      test = rotate_cfs(tensor, x)
    res.append(test.trace())

  return res



# performs a cfs-rotation study for elast (Voigt), elec and piezo tensors
# in the result the pi is copied to the end!
# @param steps: how many probes
# @param aux: "default" = none, "ortho_norm", "mono_norm" for 3D elasticity, "e21_normed"
# @return: 1st: list of angles (angle or (theta, phi) in 3D), 2nd: vector of magnitudes, 3rd: aux vector of option
def perform_voigt_tensor_samping(tensor, steps, aux_data = "default"):
  angle = []
  data  = []
  aux   = []
  
  assert(np.ndim(tensor) == 2)

  if tensor.shape == (6,6):
    # 3D case for mech
    # we do not use the standard spherical coordinate system but in to_vector() stuff from our rotation matrix
    # therefore, as z = sin(theta) instead of cos(theta) we need to run from [pi/2, 3/2 * pi] to cover all z

    #Create points on Surface
    for theta in np.arange(0.5 * np.pi, 1.5001 * np.pi, np.pi / (steps-1)):
      for phi in np.arange(0, 2.0001 * np.pi, 2 * np.pi / (steps-1)):
        test = rotate_cfs(tensor, phi, theta)
        angle.append((phi, theta))
        data.append(test[0,0])

        if aux_data == "ortho_norm":
          aux.append(sqrt(2.0 * (test[0,3]**2 + test[1,3]**2 + test[2,3]**2 + test[0,4]**2 + test[1,4]**2 + test[2,4]**2 + test[3,4]**2 +
                                 test[0,5]**2 + test[1,5]**2 + test[2,5]**2 + test[3,5]**2 + test[4,5]**2)))
        if aux_data == "mono_norm":
          aux.append(sqrt(2.0 * (test[0,4]**2 + test[1,4]**2 + test[2,4]**2 + test[3,4]**2 +
                                 test[0,5]**2 + test[1,5]**2 + test[2,5]**2 + test[3,5]**2)))
  else:
    # 2D case for elec, piezo, mech
    # for mech it is enought to go from 0 to pi and duplicate, for piezo we need all. Do slow and easy
    for x in np.arange(0, 2.0*np.pi, 2.0*np.pi/steps):
      test = rotate_cfs(tensor, x)
      #print "angle=" + str(x) + " -> " + str(test)
      angle.append(x)
      data.append(test[0,0])
      #print str(test[0,2])

      if aux_data == "ortho_norm" or aux == "mono_norm":
        if tensor.shape == (2,2):
          aux.append(sqrt(2.0 * test[0,1]**2))
        if tensor.shape == (3,3):
          aux.append(sqrt(2.0 * (test[0,2]**2 + test[1,2]**2)))
        if tensor.shape == (2,3):
          aux.append(np.min((sqrt(test[0,0]**2 + test[0,1]**2 + test[1,2]**2), sqrt(test[1,0]**2 + test[1,1]**2 + test[0,2]**2))))
      if aux_data == "e21_normed":
        aux.append(test[2,2])
        # aux.append(np.abs(test[1,0]/test[1,1]))

  assert(len(data) == len(angle))
  assert(len(aux) == len(data) or len(aux) == 0)
  return angle, data, aux

# finds the two largest maxima
# @param param: data vector where the last two entries shall be a double of the first ones
# @return the indices of the largest and second largest maxima. If there is no second then -1
# THIS FUNCTION MIGHT NOT BE CORRECT!!!
def find_maxima(data):
  first = [-1, -1e60]
  second = [-1, -1e60]
  for i in range(1, (len(data)/2)+1):
    if data[i-1] <= data[i] and data[i] >= data[i+1]:
      if data[i] >= first[1]:
        second = first[:] # deep copy
        first[0] = i
        first[1] = data[i]
      elif data[i] >= second[1]:
        second[0] = i
        second[1] = data[i]

  return first[0], second[0]

## find all minima
# searches only in half+1 space
# @return unsorted list of (angle, value)
def find_minima(angle, data):
  result = []

  # 1D or 2D search?
  if type(angle[0]) == float:
    for i in range(1, (len(data)/2)+1):
      if data[i-1] >= data[i] and data[i] <= data[i+1]:
        result.append((angle[i], data[i]))
  else:
    samples = int(sqrt(len(angle)))
    half    = int(samples/2)
    assert(samples**2 == len(angle))

    for i in range(0, half):
      for j in range(0, samples):
        #print "i=" + str(i) + " j=" + str(j)
        idx = i * samples + j
        this  = data[idx]
        #print "this=" + str(i * samples + j) + " phi=" + str(angle[idx][0]) + " theta=" + str(angle[idx][1]) + " -> " + str(this)
        north = data[(i-1 if i > 0 else half) * samples + j]
        #print "north=" + str((i-1 if i > 0 else half) * samples + j) + " -> " + str(north)
        south = data[i+1 * samples + j] # range is samples/2
        #print "south=" + str(i+1 * samples + j) + " -> " + str(south)
        west  = data[i * samples + (j-1 if j > 0 else samples-1)]
        #print "west=" + str(i * samples + (j-1 if j > 0 else samples-1)) + " -> " + str(west)
        east  = data[i * samples + (j+1 if j < samples-1 else 0)]
        #print "east=" + str(i * samples + (j+1 if j < samples-1 else 0)) + " -> " + str(east)
        if this <= north and this <= south and this <= west and this <= east:
          #print "new minimum found"
          result.append((angle[idx], this))

  return result


## transforms Hill-Mandel to Voigt elasticity tensor
def HillMandel2Voigt(tensor):
  ret = tensor.copy()

  for i in range(len(ret)-1):
    ret[i, len(ret)-1] *= 1/sqrt(2.0)
    ret[len(ret)-1, i] *= 1/sqrt(2.0)

  ret[len(ret)-1,len(ret)-1] *= 0.5
  return ret

## transforms Voigt elasticity tensor to Hill-Mandel notation
def Voigt2HillMandel(tensor):
  ret = tensor.copy()

  for i in range(len(ret)-1):
    ret[i, len(ret)-1] *= sqrt(2.0)
    ret[len(ret)-1, i] *= sqrt(2.0)

  ret[len(ret)-1,len(ret)-1] *= 2.0
  return ret


# creates a 2D elasticity tensor. Do the HillMandel2Voigt conversion if necessary!
# we have a trace and a column style
# t2d = 11 22 33 23 13 12
# t3d = 11 22 33 44 55 66 23 13 12 34 24 14 45 35 25 15 56 46 36 26 16 
# c3d = 11 12 22 13 23 33 14 24 34 44 15 25 35 45 55 16 26 36 46 56 66
# @param style None, 'trace' or 'column'. None = default 2d is trace style and 3d is column style
# note that tensorviz.py has trace as default for style
def to_mech_tensor(input, style = None):
  assert style in [None, 'trace', 'column']
  assert len(input) == 6 or len(input) == 21

  if len(input) == 6:
    assert not style == 'column' # not implemented
    # "e11", "e22", "e33", "e23", "e13", "e12";
    #    0      1      2      3      4      5
    tensor = np.zeros((3,3))
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
  else:
    tensor = np.zeros((6,6))

    if style == 'trace':
      # 11 22 33 44 55 66 23 13 12 34 24 14 45 35 25 15 56 46 36 26 16
      #  0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20        
      tensor[0,0] = input[0]
      tensor[1,1] = input[1]
      tensor[2,2] = input[2]
      tensor[3,3] = input[3]
      tensor[4,4] = input[4]
      tensor[5,5] = input[5]
      tensor[1,2] = input[6]
      tensor[0,2] = input[7]
      tensor[0,1] = input[8]
      tensor[2,3] = input[9]
      tensor[1,3] = input[10]
      tensor[0,3] = input[11]
      tensor[3,4] = input[12]
      tensor[2,4] = input[13]
      tensor[1,4] = input[14]
      tensor[0,4] = input[15]
      tensor[4,5] = input[16]
      tensor[3,5] = input[17]
      tensor[2,5] = input[18]
      tensor[1,5] = input[19]
      tensor[0,5] = input[20]
    else:  
      # e11, e12, e22, e13, ..
      idx = 0
      for x in range(6):
        for y in range(6):
          if y > x:
             continue
          else:
            val = float(input[idx])
            idx += 1
            tensor[x][y] = val
            tensor[y][x] = val
    return tensor

## invert to_mech_tensor()
# 2D = "strain notation", 3D is by upper columns
# @param as_array as given by H5
def to_mech_vector(tensor, as_array=False):
  assert(tensor.shape == (3,3) or tensor.shape == (6,6))

  vec = []

  if len(tensor) == 3:
    vec.append(tensor[0,0])
    vec.append(tensor[1,1])
    vec.append(tensor[2,2])
    vec.append(tensor[1,2])
    vec.append(tensor[0,2])
    vec.append(tensor[0,1])
  else:
    for x in range(6):
      for y in range(6):
        if y > x:
           continue
        else:
          vec.append(tensor[y,x])

    assert(len(vec) == 21)

  if as_array:
    array = np.zeros((len(vec)))
    array[:] = vec
    return array
  else:
    return vec

# creates a piezoelectric coupling tensor
def to_piezo_tensor(input):
  # print "tpt: " + str(len(input)) + " -> " + str(input)
  assert(len(input) == 6)
  # "e11", "e12", "e13", "e21", "e22", "e23";
  #    0      1      2      3      4      5
  tensor = np.zeros((2,3))
  tensor[0,0] = input[0]
  tensor[0,1] = input[1]
  tensor[0,2] = input[2]
  tensor[1,0] = input[3]
  tensor[1,1] = input[4]
  tensor[1,2] = input[5]
  return tensor


## invert to_piezo_tensor()
def to_piezo_vector(tensor, as_array=False):
  assert(tensor.shape == (2,3))

  vec = []

  vec.append(tensor[0,0])
  vec.append(tensor[0,1])
  vec.append(tensor[0,2])
  vec.append(tensor[1,0])
  vec.append(tensor[1,1])
  vec.append(tensor[1,2])

  if as_array:
    array = np.zeros((len(vec)))
    array[:] = vec
    return array
  else:
    return vec


# creates a permittivity tensor
def to_elec_tensor(input):
  assert(len(input) == 3)
  # "e11", "e22", "e12"
  #    0      1      2
  tensor = np.zeros((2,2))
  tensor[0,0] = input[0]
  tensor[0,1] = input[2]
  tensor[1,0] = input[2]
  tensor[1,1] = input[1]


  return tensor

## give vector from angle
# angle a scalar or (phi,theta)
def to_vector(angle):

  result = np.zeros((3))
  if np.isscalar(angle):
    result[0] = cos(angle)
    result[1] = sin(angle)
    result[2] = 0
  else:
    assert(len(angle) == 2)
    # we use not the usual conversion for the spherical coordinate system but the first row of the rotation matrix!
    # x = radius * cos(phi) * sin(theta)
    # y = radius * sin(phi) * sin(theta)
    # z = radius * cos(theta)
    result[0] = cos(angle[1]) * cos(angle[0])
    result[1] = -1 * cos(angle[1]) * sin(angle[0])
    result[2] = sin(angle[1])

  return result

## gives the two Poisson's ratios for a given tensor which needs to be in Voigt notation in elasticity notation
# @return v21, v12
def poissons_ratio(tensor):
  c = numpy.linalg.inv(tensor)
  return (-c[0,1] / c[1,1]), (-c[1,0] / c[0,0])

## find stiffest directions for a tensor
# @param tensor elast tensor in voigt notation, elec tensor or piezo tensor
# @return first, second the scaled vectors of stiffest directions. second is 0,0,0 if there is none, v21 and v12 if desired
# THIS FUNCTION MIGHT NOT BE CORRECT!!!
def find_stiffest_orientation(tensor, steps, also_poissons_ratio=False):
  angle, data, aux = perform_voigt_tensor_samping(tensor, steps, also_poissons_ratio)
  idx_first, idx_second = find_maxima(data)

  first = to_vector(angle[idx_first]) * data[idx_first]

  second = [0,0,0]
  if idx_second > -1:
    second = to_vector(angle[idx_second]) * data[idx_second]

  if also_poissons_ratio:
    # find the one corresponding to the stronger direction
    # note, we might have zeros not in the strongest direction (-> PZT-5A!!)
    idx_first, idx_second = find_minima(aux)
    idx = idx_first
    if data[idx_second] > data[idx_first]:
      idx = idx_second
    v21, v12 = poissons_ratio(rotate_cfs(tensor, angle[idx]))
    return first, second, v21, v12
  else:
    return first, second


# e2d = np.zeros((2,2))
# e2d[0,0] = 1.51
# e2d[0,1] = 0.0
# e2d[1,0] = 0.0
# e2d[1,1] = 1.27
#
# p0 = np.zeros((2,3))
# p0[0,0] = 0.01
# p0[0,1] = 0.01
# p0[0,2] = 17.0
# p0[1,0] = -6.5
# p0[1,1] = 23.3
# p0[1,2] = 0.01
#
# e3d = np.zeros((6,6))
# e3d[0,0] = 1.0
#
# iso3d = to_mech_tensor(eval("[9.999406e-01,2.999666e-01,9.999406e-01,2.999666e-01,2.999666e-01,9.999406e-01,0.0,0.0,0.0,3.499870e-01,0,0,0,0,3.499870e-01,0,0,0,0,0,3.499870e-01]"))
#
# t3d = HillMandel2Voigt(to_mech_tensor(eval("[0.00401617093935,-2.56173585052e-07,0.00401617094245,2.20294623258e-07,2.2027694833e-07,0.00401566744551,3.38547287661e-07,3.38534101428-07,-5.93702112942e-07,0.00401568927622,-4.13754162037e-08,-6.18071141999e-08,6.94156633186e-08,2.81430319706e-07,0.00401628783183,-6.181071466e-08,-4.13602971057e-08,6.94018181029e-08,2.81408882188e-07,-1.66007908852e-07,0.00401628785427]")))
#
# triv3d = to_mech_tensor(eval("[1, 0, 0.5, 0, 0, 0.25, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ,0 ]"))
#
# orient_0_degrees = to_mech_tensor(eval("[2.022,0.615,0.0148,0.0,0.0,0.0949]"))
