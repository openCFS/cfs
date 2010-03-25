# -*- coding: utf-8 -*-
import libxml2
import math
import os

# replace a single xpath value -> must exsit once!
# xml is a xpathContext: doc = libxml2.parseFile("params.xml") -> xml = doc.xpathNewContext()
# save your doc later via doc.saveFile("param_tmp.xml")
def replace(xml, path, value):
  res = xml.xpathEval(path)
  if  len(res) == 0:
    raise RuntimeError(path + " not found")
  if len(res) > 1:
    raise RuntimeError(path + " has " + len(res) + " hits")
  data = res[0]
  data.setContent(value)
  return    

# returns an xpath value
def xpath(xml, path):
  res = xml.xpathEval(path)
  if  len(res) == 0:
    raise RuntimeError(path + " not found")
  if len(res) > 1:
    raise RuntimeError(path + " has " + len(res) + " hits")
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

  # other is also Coordinate
  def dist(self, other):
    return math.sqrt((self.x - other.x)**2 + (self.y - other.y)**2 + (self.z - other.z)**2)  
 
  def printline(self):
    print str(self.x) + ", " + str(self.y) + ", " + str(self.z)
    
  def toString(self):
    return str(self.x) + ", " + str(self.y) + ", " + str(self.z) 
       

# helper: return numpy.ndarray element with 2/3D tolerance
# if dim = 3 the k entry is ignored
def getNDArrayEntry(data, dim, i, j, k):
  if dim == 3:
    return data[i,j,k]
  if dim == 2:
    return data[i, j]
  raise " cannot handle dimension " + str(dim)
    
# see getNDArrayEntry(data, dim, i, j, k)
def setNDArrayEntry(data, dim, i, j, k, value):
  if dim == 3:
    data[i,j,k] = value
    return
  if dim == 2:
    data[i, j] = value
    return
  raise " cannot handle dimension " + str(dim)


