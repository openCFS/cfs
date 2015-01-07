# -*- coding: utf-8 -*-
import libxml2
import math
import os
import string
import numpy

# replace a single xpath value -> must exsit once!
# xml is a xpathContext: doc = libxml2.parseFile("params.xml") -> xml = doc.xpathNewContext()
# save your doc later via doc.saveFile("param_tmp.xml")
def replace(xml, path, value):
  res = xml.xpathEval(path)
  if  len(res) == 0:
    raise RuntimeError(path + " not found")
  if len(res) > 1:
    raise RuntimeError(path + " has " + str(len(res)) + " hits")
  data = res[0]
  data.setContent(value)
  return    

## removes the defined xml entity.
# extend to attribute first by .hasProp('name2') == None stuff
def remove(xml, path):
  res = xml.xpathEval(path)
  if  len(res) == 0:
    raise RuntimeError(path + " not found")
  if len(res) > 1:
    str(res)
    raise RuntimeError(path + " has " + str(len(res)) + " hits")
  data = res[0]
  # TODO the node/attribute = prop stuff
  data.unlinkNode() 


# returns an xpath value
def xpath(xml, path):
  res = xml.xpathEval(path)
  if  len(res) == 0:
    raise RuntimeError(path + " not found")
  if len(res) > 1:
    str(res)
    raise RuntimeError(path + " has " + str(len(res)) + " hits")
  data = res[0]
  return data.getContent()    

# does at leas one element exist
def has(xml, path):
  res = xml.xpathEval(path)
  if  len(res) == 0:
   return false
  else:
   return true


  
# dump a xml node
def dump(xml, path):
  res = xml.xpathEval(path)
  if  len(res) == 0:
    raise RuntimeError(path + " empty")
  for i in res:
    print str(i) + ": " + i.getContent()
  
# mimic conditional operator
def cond(test, trueval, falseval):
  if test:
    return trueval
  else:
    return falseval

## apparently python has to string to bool casting
def toBool(string):
  c = string[0].upper()
  return c == "T" or c == "Y"

## transform the number to the significant digits. 
# done slowly by strings conversion as an /10.0 might lead to artefacts .:(
def digits(value, decimal_place):
  format = "%." + str(decimal_place) + "f"
  string = format % value
  return float(string) 

# covert a complex number "(a,b)" to two gnuplot compatible strings "a \t b" 
# -> remove brackets and replace the comma by a tab, nothing else
# return the string for gnuplot printing
def toGnuPlot(complex_string):
  ret = string.lstrip(complex_string, "(")
  ret = string.rstrip(ret, ")")
  return string.replace(ret, ",", "\t")

# execute cmd and rais error when not 0
def execute(cmd, output = False):
 if output:
   print cmd
 ret = os.system(cmd) <> 0
 if(ret <> 0):
   raise RuntimeError("execution of '" + cmd + "' -> " + str(ret))
 return     

# return the first line of a file
def first_line(file_name, append = ""):
   file = open(file_name, "r")
   line = file.readline()
   file.close() 
   if len(line) == 0:
     raise RuntimeError("file '" + file_name + "' is empty")
   return line.strip() + append + "\n"

# return the second line of a file, this skips the gnuplot header
def second_line(file_name, append = ""):
   file = open(file_name, "r")
   lines = file.readlines()
   file.close() 
   if len(lines) == 0:
     raise RuntimeError("file '" + file_name + "' is empty")
   return lines[1].strip() + append + "\n"

# return the last line of a file and append 'append'
def last_line(file_name, append = ""):
   file = open(file_name, "r")
   lines = file.readlines()
   file.close() 
   if len(lines) == 0:
     raise RuntimeError("file '" + file_name + "' is empty")
   string = lines[len(lines)-1]
   return string.rstrip() + append + "\n"

# there shall be a predefined class somewhere, I just didn't find it
class Coordinate:
  def __init__(self):
    self.x = 0.0
    self.y = 0.0 
    self.z = 0.0 

  def __init__(self, x, y, z):
    self.x = x
    self.y = y 
    self.z = z 
    
  # i, j, k ints from 0 to div-1  
  def toCoordinate(self, i, j, k, div):
    # shift the coordinate to the center of the elements
    self.x = i / float(div) + 0.5 / div
    self.y = j / float(div) + 0.5 / div 
    self.z = k / float(div) + 0.5 / div 

  # export to arry for numerical stuff
  def toArray(self):
    return [self.x, self.y, self.z]

  # other is also Coordinate
  def dist(self, other):
    return math.sqrt((self.x - other.x)**2 + (self.y - other.y)**2 + (self.z - other.z)**2)  
 
  def printline(self):
    print str(self.x) + ", " + str(self.y) + ", " + str(self.z)
    
  def toString(self):
    return str(self.x) + ", " + str(self.y) + ", " + str(self.z) 
       

# extracts an entry, if data is of lower dimension, the indices are ignored
def getNDArrayEntry(data, i, j, k, d = None):
  if data.ndim == 4:
    return data[i,j,k,d]
  if data.ndim == 3:
    return data[i,j,k]
  if data.ndim == 2:
    return data[i, j]
  if data.ndim == 1:
    return data[i]
  raise RuntimeError("cannot handle dimension")
    
# see getNDArrayEntry(data, dim, i, j, k)
# save_out_of_dim do nothing if infeasible dimensions are not 0
def setNDArrayEntry(data, i, j, k, value, save_out_of_dim = False):
  if data.ndim == 3:
    data[i,j,k] = value
    return
  if data.ndim == 2:
    if not save_out_of_dim or k == 0:
      data[i, j] = value
    return
  if data.ndim == 1:
    if not save_out_of_dim or (j == 0 and k == 0):
      data[i] = value
    return
  raise RuntimeError("cannot handle dimension")

## returns the x, y, and z dimension of a ndarray. z=1 for 2d 
# call x, y, z = getDim(data)
# call x, y, z, d = getDim(data, True)
def getDim(data, get4dims = False):
  x = data.shape[0]
  y = 1
  if data.ndim >= 2:
    y = data.shape[1]
  z = 1
  if data.ndim >= 3:
    z = data.shape[2]
  if not get4dims:
    return x, y, z
  else:
    d = None
    if data.ndim >= 4:
      d = data.shape[3]
    return x, y, z, d

## helps to clean an array with repeated entries as it happens hen nodes and elements are defined in cfs with a too small inc value
# @param data array which is a history file read by numpy.loadtxt()
def cleanOversampledArray(data):
  assert(len(data.shape) == 2)
  assert(data.shape[0] > 1)  
  
  # find unique indices
  unique = []
  # the first element is unique
  unique.append(0)
  line = data[:,0] 
  for i in range(1, data.shape[0]):
    if line[i] <> line[unique[len(unique)-1]]:
      unique.append(i)
  
  # copy unique data
  columns = data.shape[1]
  result = numpy.zeros((len(unique), columns))
  for i in range(len(unique)):
    for c in range(columns):
      val = data[unique[i], c]
      result[i][c] = val
  
  return result    


## convert a list to a numpy array 
def listToNDArray(data):
  ret = numpy.zeros((len(data)))
  
  for i in range(len(data)):
    ret[i] = data[i]
    
  return ret  


## finds a value in an ndarray
# @param silent if True -1,-1,-1 is returned, otherwise an error
# @return the coordinates x, y, z or an error, see silent
def findInNDArray(data, value, silent=False):
  x, y, z = getDim(data)
  for i in range(x):
    for j in range(y):
      for k in range(z):
        if data[i, j, k] == value:
          return i, j, k
 
  if not silent:
    raise RuntimeError(" value'" + str(value) + "' not found in data with x=" + str(x) + " y=" + str(y) + " z=" + str(z))
  else:
    return -1, -1, -1

## checks the status of a CFS problem run by the info.xml file
# @param problem string without '.info.xml' 
# return 'not_found', 'running', 'finished', 'aborted'. 'cannot_determine'
def check_cfs_status(problem):
  if os.path.exists(problem + ".info.xml"):
    try:
      doc = libxml2.parseFile(problem + ".info.xml")
      xml = doc.xpathNewContext()
      status = xpath(xml, "//cfsInfo/@status")
      return status
    except:
      return "cannot_determine"
  else:
    return "not_found"
     
