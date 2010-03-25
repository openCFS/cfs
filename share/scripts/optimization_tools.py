# This collects some tool routines for optimization

import libxml2
import numpy
import math

from cfs_utils import *

## Reads a cubic (3D) density.xml file as 3D NDArray
# @param filename from which the last 'set' is used
def read_cubic_density(filename):
  doc = libxml2.parseFile(filename)
  xml = doc.xpathNewContext()
  # read the design attributes from the last set
  res = xml.xpathEval("//set[last()]/element/@design")
      
  # validation
  if len(res) == 0:
    raise RuntimeError("found no elements in the last 'set' of '" + filename + "'")
  raw_size = pow(len(res), (1.0/3.0)) 
  size = int(round(raw_size))

  if abs(size - raw_size) > 0.0001:
    raise RuntimeError("last set has " + str(len(res)) + " elements which appears to be not cubic")

  # parse data to array
  data = numpy.zeros((size, size, size))

  # copy data from linear list
  for i in range(size):
    for j in range(size):
      for k in range(size):
        idx = int(i * size * size + j * size + k)
        # if idx < 0 or idx >= len(res):
        #  raise RuntimeError("idx " + str(idx) + " is not within " + str(len(res)) + " stupid!")
        # print str(idx) + " -> " + str(res[idx].getContent())
        data[i, j, k]  = float(res[idx].getContent()) 

  return data  

## write the data to a density.xml file
# @param data_inp a ndata array of dim 2 or three or a list of data
# @param setname_inp the name of the set or a list of setnames
def write_density_file(filename, data_inp, setname_inp):
  out = open(filename, "w")
  out.write('<?xml version="1.0"?>\n')
  out.write('<cfsErsatzMaterial>\n')
  out.write('  <header>\n')
  out.write('    <design initial="0.5" lower="1e-3" name="density" region="mech" upper="1"/>\n')
  out.write('    <transferFunction application="mech" design="density" param="1" type="simp"/>\n')
  out.write('  </header>\n')
 
  # check if we deal with lists or not
  data_list = []
  setname_list = []
  if type(data_inp) is list:
    data_list = data_inp
    setname_list = setname_inp
  else:
    data_list.append(data_inp)
    setname_list.append(setname_inp)
 
  for i in range(len(data_list)):
    data    = data_list[i]
    setname = setname_list[i]
    out.write('  <set id="' + setname + '">\n')
    edge = data.shape[0] # be careful!
    counter = 1

    dim = data.ndim 

    for i in range(edge):
      for j in range(edge):
        for k in range(cond(dim == 2, 1, edge)):
           val = getNDArrayEntry(data, dim, i, j, k)
           out.write('    <element nr="' + str(counter) + '" type="density" design="' + str(val) + '"/>\n')
           counter = counter + 1       
         
    out.write('  </set>\n')

  out.write(' </cfsErsatzMaterial>\n')
  out.close()


## apply a treshold filter on a 2D/3D matrix
# @param data original 2D/3D data
# @param treshold all smaller threshold becomes min, otherwise max
# @return a new data array 
def threshold_filter(data, threshold, min, max):
  
  dim = data.ndim
  edge = data.shape[0]
  res = numpy.zeros((edge, edge, edge))
  
  for i in range(edge):
    for j in range(edge):
      for k in range(cond(dim == 2, 1, edge)):
        org = getNDArrayEntry(data, dim, i, j, k)
        val = cond(org < threshold, min, max)
        setNDArrayEntry(res, dim, i, j, k, val)        

  return res  



## rotate a 3d matrix
# exchange x and y
def rotate_matrix_x_y(data):
  
  dim = data.ndim
  edge = data.shape[0]
  res = numpy.zeros((edge, edge, edge))
  
  for i in range(edge):
    for j in range(edge):
      for k in range(cond(dim == 2, 1, edge)):
        val = getNDArrayEntry(data, dim, i, j, k)
        setNDArrayEntry(res, dim, j, i, k, val)        

  return res  

## rotate a 3d matrix
# exchange x and z
def rotate_matrix_x_z(data):
  
  dim = data.ndim
  edge = data.shape[0]
  res = numpy.zeros((edge, edge, edge))
  
  for i in range(edge):
    for j in range(edge):
      for k in range(cond(dim == 2, 1, edge)):
        val = getNDArrayEntry(data, dim, i, j, k)
        setNDArrayEntry(res, dim, k, j, i, val)        

  return res

## Enlarge a density matrix by multiples of the original data
# @param data a two or three dimensional ndata array
# @param x_times how often to repeat in x direction. 1 = no change
# @param z_times ignored when data is 2D
# @return a new ndata array (might be huge!!)
def enlarge_matrix(data, times_x, times_y = 1, times_z = 1):

  dim = data.ndim
  x_edge = data.shape[0]
  y_edge = data.shape[1]
  z_edge = cond(dim == 2, 0, data.shape[2])
  res = numpy.zeros((x_edge * times_x, y_edge * times_y, z_edge * times_z))
  
  for x in range(times_x):
    for y in range(times_y):
      for z in range(cond(dim == 2, 0, times_z)):
        # copy complete block
        x_base = x * x_edge
        y_base = y * y_edge
        z_base = z * z_edge
        
        for i in range(x_edge):
          for j in range(y_edge):
            for k in range(z_edge):
              val = getNDArrayEntry(data, dim, i, j, k)
              setNDArrayEntry(res, dim, x_base + i, y_base + j, z_base + k, val)

  return res                     
                     