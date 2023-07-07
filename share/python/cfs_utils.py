# -*- coding: utf-8 -*-


# this is to prevent annoying  
# /usr/lib64/python3.6/site-packages/h5py/__init__.py:36: FutureWarning: Conversion of the second argument of issubdtype from `float` to `np.floating` is deprecated. In future, it will be treated as `np.float64 == np.dtype(float).type`.
#  from ._conv import register_converters as _register_converters
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

# libxml2 is not available for python3. lxml is available on python2 and python3 and is a successor of libxml2
# lxml is technically based on the libxml2 C-code bat hat the nicer python interface
import lxml
import lxml.etree
import six
import math
import os
import string
import numpy

# helper to write values and coordinates to a csv (comma separated values) file
# see element_to_node_2d in hdf5_tools
def write_cvs_file(filename, nodes, values):
  import pandas

  nodes_dim = len(nodes[0]) 
  val_dim = len(values[0])

  if nodes_dim == 2 and val_dim == 2:  
    df = pandas.DataFrame(list(zip(*[nodes[:, 0], nodes[:, 1], values[:,0], values[:,1]])))
    df.to_csv(filename, header=False, index=False)
  else:
    raise "implemented generic write_cvs_file"  


## trivial helper which only helps to avoid lxml.etree stuff
def open_xml(file):
  if not os.path.exists(file):
    raise RuntimeError("xml file '" + file + "' not found")
  
  xml = lxml.etree.parse(file, lxml.etree.XMLParser(remove_comments=True))
  return xml


# helper to get namespace for a lxml query. if 'cfs:' or 'mat:' is in the query a mapping is returned, else None
def namespace(query):
  if 'cfs:' in query:
    return {'cfs':'http://www.cfs++.org/simulation'}
  elif 'mat:' in query:
    return {'mat':'http://www.cfs++.org/material'}
  else:
    return None

# replace a single xpath value -> must exist once!
# the xpath shall contain a single result. e.g. '//cfs:materialData/@file
# internally we get via lxml the element by removing /@file from the expression
# if the original attribute value has '/nx' this will stay, even when value has not '/nx' set.
#@param unique if True there must be one match, if false there may be more than one - none is also an option
#@return the number of replaces attributes
def replace(xml, path, value, unique = True):
  res = xml.xpath(path, namespaces = namespace(path))
  if unique:  
    if len(res) == 0:
      raise RuntimeError(path + " not found")
    if len(res) > 1:
      raise RuntimeError(path + " not unique, has " + str(len(res)) + " hits")

  # in the attribute case we have to fake
  idx = path.rfind('/@')
  if idx > 0:
    # query element
    elem = xml.xpath(path[0:idx], namespaces = namespace(path))
    for e in elem:
      data = e.attrib[path[idx+2:]]
      by_nx = str(data).find('/nx') > 0 and str(value).find('/nx') == -1
      e.attrib[path[idx+2:]] = str(value) + ('/nx' if by_nx else '')
    return len(elem)
  else:
    for e in res:
      e.text = str(value)
    return len(res)

## removes the defined xml entity.
def remove(xml, path, unique = True):
  res = xml.xpath(path, namespaces = namespace(path)) 
  if unique: 
    if len(res) == 0:
      raise RuntimeError(path + " not found")
    if len(res) > 1:
      raise RuntimeError(path + " has " + str(len(res)) + " hits")
  for data in res:
    data.getparent().remove(data)

# returns an xpath value. Assumes we have lxml as xml tree, change your old libxml2 code.
# if 'cfs:' is in path we add a namespace mapping automatically
# example '//cfs:materialData/@file' for any materialData element;
# to get text inside a xml tag, e.g. material tensor from .info.xml:
# xpath(xml, "//iteration[last()]/homogenizedTensor/tensor/real/text()") gives a string
def xpath(xml, path):
  # assume lxml first
  try:
    res = xml.xpath(path, namespaces = namespace(path))  
    if  len(res) == 0:
      raise RuntimeError(path + " not found with ns='" + str(namespace(path)) + "'")
    if len(res) > 1:
      str(res)
      raise RuntimeError(path + " has " + str(len(res)) + " hits")
    data = res[0]
    return str(data)
  except AttributeError: # this happens when xml is from libxml2 and not lxml
    if six.PY2:
      raise RuntimeError('is your xml paramater from libxml2? You need to switch to lxml')
    else: 
      raise RuntimeError('parameter seems to be no lxml attribute ' + str(xml))

  
# does at least one element exist
def has(xml, path):
  try:
    res = xml.xpath(path, namespaces = namespace(path)) 
    if  len(res) == 0:
      return False
    else:
      return True
  except AttributeError: # this happens when xml is from libxml2 and not lxml
    if six.PY2:
      raise RuntimeError('is your xml paramater from libxml2? You need to switch to lxml')
    else: 
      raise RuntimeError('parameter seems to be no lxml attribute ' + str(xml))
  
# dump a xml node
def dump(xml, path):
  res = xml.xpath(path, namespaces = namespace(path))
  if  len(res) == 0:
    raise RuntimeError(path + " empty")
  for i in res:
    lxml.etree.dump(i)
  
# mimic conditional operator. depreciated, use python a = b if i > 2 else c
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

## the python function isdigit() checks only for 0...9, hence -3 and 3.2 are no digits
def isfloat(string):
  try:
    float(string)
    return True
  except ValueError:
    return False
  
  
# get the real part of a complex number string of type '(r,i)' as float
def getReal(complex_string):
  assert(complex_string[0] == '(')
  assert(complex_string.find(',') > 1)
  return float(complex_string[1:complex_string.find(',')])

# covert a complex number "(a,b)" to two gnuplot compatible strings "a \t b" 
# -> remove brackets and replace the comma by a tab, nothing else
# return the string for gnuplot printing
def toGnuPlot(complex_string):
  ret = string.lstrip(complex_string, "(")
  ret = string.rstrip(ret, ")")
  return string.replace(ret, ",", "\t")


## helps generate a shell script for submission via qsub on RRZE HPC systems
# It takes a template script, adds optinally a cd to current and the cmd line
#@param template use qsub_templtate.sh as base for your own template, it might be sufficient
#@param cmd the cfs call from run.py
#@param filename shall end with .sh
#@param silent if so suppress output#
#@return the qsub command 
def generate_qsub_script(template, cmd, filename, silent = False):
  if not os.path.exists(template):
    raise RuntimeError("qsub template not found '" + template + "'")
   
  with open(template) as f:
    lines = f.readlines()

  # do we have a cd?
  cd = [l for l in lines if l.startswith('cd')]
  if not cd:
    pwd = os.getcwd()
    if not silent:
      print(" add 'cd " + pwd + "'")
    lines.append('cd ' + pwd + '\n')
  
  # the job to be executed
  lines.append(cmd + '\n')
  
  # write the new qsub file
  if not silent:
    print("generate script '" + filename + "'")
  out = open(filename, "w")
  out.writelines(lines)
  out.close()
  return "qsub " + filename

# execute cmd and raise error when not 0 and not silent
# return error code, 0 for no problem
def execute(cmd, output = False, silent = False):
 if output:
   print(cmd)
 ret = os.system(cmd)
 if ret != 0 and not silent:
   raise RuntimeError("execution returns " + str(ret) + ": '" + cmd + "'") 
 return ret    

# return the first line of a file
def first_line(file_name, append = "", no_eol = False):
   file = open(file_name, "r")
   line = file.readline()
   file.close() 
   if len(line) == 0:
     raise RuntimeError("file '" + file_name + "' is empty")
   return line.strip() + append + ("\n" if not no_eol else "")

# return the second line of a file, this skips the gnuplot header
def second_line(file_name, append = "", no_eol = False):
   file = open(file_name, "r")
   lines = file.readlines()
   file.close() 
   if len(lines) == 0:
     raise RuntimeError("file '" + file_name + "' is empty")
   return lines[1].strip() + append + ("\n" if not no_eol else "")

# return the last line of a file and append 'append'
def last_line(file_name, append = "", no_eol = False):
   file = open(file_name, "r")
   lines = file.readlines()
   file.close() 
   if len(lines) == 0:
     raise RuntimeError("file '" + file_name + "' is empty")
   string = lines[len(lines)-1]
   return string.rstrip() + append + ("\n" if not no_eol else "")

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
    print(str(self.x) + ", " + str(self.y) + ", " + str(self.z))
    
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
    if line[i] != line[unique[len(unique)-1]]:
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
      xml = open_xml(problem + ".info.xml")
      status = xpath(xml, "//cfsInfo/@status")
      return status
    except:
      return "cannot_determine"
  else:
    return "not_found"
     
## takes a list of dicts assuming all dicts have the same structure and prints them as gnuplot table
def dicts_to_gnuplot(dicts):     
  header = "#"
  for idx, key in enumerate(dicts[0].keys()):
    header += '(' + str(idx+1) + ') ' + key + ' \t'
  print(header)

  for item in dicts:
    line = ""
    for idx, key in enumerate(dicts[0].keys()):
      line += str(item[key])  + ' \t'
    print(line)

# assumes a rectangualr regular mesh of size nx x ny and creates for a given displacement a vtk file.
# this allows to visualize displacement u with paraview without having to deal with hdf5
# u is an array of dim 2 with x as fast variable (left lower to right upper)
# u = read_displacement(hdf5_file)
def cfs_displacement_to_vtk(u, nx, ny, name):
  from pyevtk.hl import gridToVTK
  
  assert len(u)==(nx+1)*(ny+1)
  assert u.shape[1] == 2
  x=np.linspace(0,1,(nx+1))      
  y=np.linspace(0,1,(ny+1))
  
  ux = np.zeros((nx+1,ny+1,1))
  uy = np.zeros((nx+1,ny+1,1))
   # fake third component of displacement vector
  uz = np.zeros((nx+1,ny+1,1))
  
  for j in range(0,ny+1):
    for i in range(nx, -1, -1):
      idx = j * (ny+1) + i
      print(i,j,idx,u[idx][0])
      ux[i,j,0] = u[idx][0]
      uy[i,j,0] = u[idx][1]

  print('x', x)      
  print('y', y)
  print('ux', ux)
  print('uy', uy)
  # gridToVTK expects 3D data, thus we fake the third dimension
  gridToVTK(name, x,y,np.zeros(1),pointData={"displacement":(ux,uy,uz)})
  print("# wrote '" + name + ".vtr'") 

