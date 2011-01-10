# -*- coding: utf-8 -*-
import libxml2
import math
import os
import string

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
def execute(cmd):
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
def getNDArrayEntry(data, i, j, k):
  if data.ndim == 3:
    return data[i,j,k]
  if data.ndim == 2:
    return data[i, j]
  if data.ndim == 1:
    return data[i]
  raise RuntimeError("cannot handle dimension")
    
# see getNDArrayEntry(data, dim, i, j, k)
def setNDArrayEntry(data, i, j, k, value):
  if data.ndim == 3:
    data[i,j,k] = value
    return
  if data.ndim == 2:
    data[i, j] = value
    return
  if data.ndim == 1:
    data[i] = value
    return
  raise RuntimeError("cannot handle dimension")

## returns the x, y, and z dimension of a ndarray. z=1 for 2d 
# call x, y, z = getDim(data)
def getDim(data):
  x = data.shape[0]
  y = 1
  if data.ndim >= 2:
    y = data.shape[1]
  z = 1
  if data.ndim >= 3:
    z = data.shape[2]
  return x, y, z

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
     
