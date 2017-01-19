# -*- coding: utf-8 -*-
# This collects some tools for cfs-grid handling.
# based on the grid-export by cfs -G (best together with -g)
# Who ever wants can add gid-mesh and hdf5 interfaces. 

import numpy
import numpy.linalg
import os
from lxml import etree

from cfs_utils import *

# this is an example application to process a structural optimization result on an acoustic grid.
#
# as the written density file does not contain the design and filter specification you need an
# optimization with e.g. evaluateInitialDesign
# 
# from optimization_tools import *
# from cfs_grid import *
# s = read_elements("piezo_grid.info.xml", "piezo")
# t = read_elements("acoustic_fine_grid.info.xml", "piezo")
# m = map_elements(s, t)
# o_d = read_density("struct.density.xml", False, 30, 30, 1)
# o_n = read_density("struct.density.xml", True, 30, 30, 1)
# r = apply_elmennr_mapping(o_n, m)
# write_density_file("shifted_struct.density.xml", o_d, "last", r)



## Element type
class Element:
  # x,y,z the barycenter
  def __init__(self, elemNum, x, y, z):
    self.elemNum = elemNum
    self.coord   = [x, y, z]
    
  def toString(self):
    return "elem=" + str(self.elemNum) + " coord=" + str(self.coord) 

## read the elements of an specified region
# @param info info.xml with exported grid (-G)
# @return list of Element 
def read_elements(info, region):
  if not os.path.exists(info):
    raise RuntimeError("file '" + info + "' doesn't exist")

  tree = etree.parse(info, etree.XMLParser(remove_comments=True))
  
  root = tree.getroot()
  set  = root.xpath("//grid/regionList/region[@name='" + region + "']")[0]
  
  result = [0]*len(set)
  counter = 0
  for elem in set:
    result[counter] = Element(int(elem.get("id")), float(elem.get("x")), float(elem.get("y")), float(elem.get("z")))
    counter += 1
    
  return result

## maps elements by their barycenter. E.g. when the same region has different element number in different meshes.
#  Handles also rigid shift of the regions!
# @returns a list of tuples with two element numbers
def map_elements(source, target):
  
  assert(len(source) == len(target))
  
  # identify shift
  s_lower, s_upper = find_corners(source)
  t_lower, t_upper = find_corners(target)
  
  # differences
  d_lower = numpy.subtract(t_lower, s_lower)
  d_upper = numpy.subtract(t_upper, s_upper)
  
  diff = numpy.linalg.norm(numpy.subtract(d_lower, d_upper))
  if diff > 1e-5:
    raise RuntimeError("the regions seem to have no common shift: " + str(d_lower) + " and " + str(d_upper) + " are " + str(diff) + " away") 
  
  shift = d_upper # rename to common shift
  result = [None]*len(source)

  # n^2 effort
  for s in range(len(source)):
    test = numpy.add(source[s].coord, shift)
    found = False 
    for t in range(len(target)):
      if numpy.linalg.norm(numpy.subtract(test, target[t].coord)) < 1e-5:
        result[s] = [source[s].elemNum, target[t].elemNum]
        found = True
        break
    if not found:
      raise RuntimeError("found no matching target element for " + source[s].toString() + " with shift " + str(shift))
  
  return result

## find lower and upper corner of the barycenters of an list of Element
# @return lower, upper array of size 3
def find_corners(elements):
  lower = [1e30, 1e30, 1e30]
  upper = [-1e30, -1e30, -1e30]  
  
  for e in range(len(elements)):
    test = elements[e].coord
    for c in range(3):
      if test[c] < lower[c]:
        lower[c] = test[c]
      if test[c] > upper[c]:
        upper[c] = test[c]
  
  return lower, upper       
    
        
