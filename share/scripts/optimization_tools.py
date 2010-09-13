# -*- coding: utf-8 -*-
# This collects some tool routines for optimization

import libxml2
import numpy
import math
import os
from lxml import etree

from cfs_utils import *
from distutils.command.build_scripts import first_line_re

## Reads a cubic (3D) density.xml file as 3D NDArray
# @param filename from which the last 'set' is used
def read_cubic_density(filename):
  vals = read_density_as_vector(filename)

  length = len(vals)
  size=int(math.pow(length+1, 1.0/3.0))
  print "length = " + str(length) + ', size = ' + str(size) + '^3'
  
  if math.pow(size, 3) == length:
    print "calculated length seems to be okay"
  else:
    raise RuntimeError("last set has " + str(length) + " elements which appears to be wrong")

  # parse data to array
  ret = numpy.zeros((size, size, size))
  
  # copy data from linear list
  for i in range(size):
    for j in range(size):
      for k in range(size):
        idx = int(i * size * size + j * size + k)
        ret[i, j, k] = vals[idx]
        
  print "have " + str(ret.shape) + " values"

  return ret
  
## Reads a density.xml file as vector
# @param filename from which the last 'set' is used
def read_density_as_vector(filename):
  if not os.path.exists(filename):
    raise RuntimeError("file '" + filename + "' doesn't exist")
  
  tree = etree.parse(filename, etree.XMLParser(remove_comments=True))
  
  root = tree.getroot()
  sett = root.xpath("//set[last()]")[0]
  print "reading set with id = " + sett.get("id")
  
  length = len(sett)
  
  counter = 0
  vals=[0]*length
  for element in sett:
    # traverse the elements and get the design
    vals[counter] = float(element.get("design"))
    counter = counter + 1
  
  print "found " + str(length) + " elements, read " + str(counter) + " elements"
  
  return vals
  
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

## evaluates the physical volume fraction
# @param data 2D/3D data
# @param penalty 1 for no penalty
def physical_volume(data, penalty):

  dim = data.ndim
  edge = data.shape[0]
  vol = 0.0
  
  for i in range(edge):
    for j in range(edge):
      for k in range(cond(dim == 2, 1, edge)):
        org = getNDArrayEntry(data, dim, i, j, k)
        vol = vol + pow(org, penalty)
        
  return vol / data.size              


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


## threshold towards a given final volume
# @param data original 2d/3d data
# @param min niminum final density in the result, ignores penalty, max = 1
# @param target volume fraction to be reached
# @param material_penalty penalty to interpret original data
# @return: the new data array as first value and the threshold as second value
def auto_threshold_filter(data, min, target, material_penalty):
  
  dim = data.ndim
  edge = data.shape[0]
  res = numpy.zeros((edge, edge, edge))

  lower = min
  upper = 1
  while upper - lower > 1e-14:
    mid = lower + 0.5 *(upper - lower)
    val = physical_volume(threshold_filter(data, mid, min, 1), material_penalty)
    print " lower=" + str(lower) + " mid=" + str(mid) + " upper=" + str(upper) + " val=" + str(val)
    if val > target:
      lower = mid
    else:
      upper = mid
  
  return threshold_filter(data, lower, min, 1), lower   


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
              
## asuming a regular 2d grid this function takes a density-file input
# and writes a finer version of it, where every element density is just
# spread to 4 elements
# element 1 density get written to elements 1,2,5,6 for a 2x2 infile
# @param infile the density file with the densities to be blown upper
# @param outfile name of the output file
def refine_density(infile, outfile):
  res = read_density_as_vector(infile)

  # validation
  if len(res) == 0:
    raise RuntimeError("found no elements in the last 'set' of '" + filename + "'")

  out = open(outfile, "w")
  out.write('<?xml version="1.0"?>\n')
  out.write('<cfsErsatzMaterial>\n')

  # get header from old density file
  header = extract_old_header(infile)
  out.write(header)
  
  out.write('<set id="refined">\n')

  new=[0]*4*len(res)
  print "length of new vector = " + str(len(new))

  num=int(math.sqrt(len(res)+1))
  print "copying " + str(num) + "^2 values"

  for r in range(num):
    for c in range(num):
      val =  res[r*num + c]

      start = 4*r*num + 2*c
      new[start + 0] = val
      new[start + 1] = val
      new[start + 0 + 2*num] = val
      new[start + 1 + 2*num] = val

  for i in range(len(new)):
    out.write('  <element nr="' + str(i+1)   + '" type="density" design="' + str(new[i]) + '"/>\n')

  out.write('</set>\n')
  out.write('</cfsErsatzMaterial>\n')
  out.close()

  return

## asuming a regular grid 3d this function takes a density-file input
# and writes a finer version of it, where every element density is just
# spread to 8 elements
# @see 2d version above
# @param infile the density file with the densities to be blown upper
# @param outfile name of the output file
def refine_density_3d(infile, outfile):
  res = read_density_as_vector(infile)

  # validation
  if len(res) == 0:
    raise RuntimeError("found no elements in the last 'set' of '" + filename + "'")

  out = open(outfile, "w")
  out.write('<?xml version="1.0"?>\n')
  out.write('<cfsErsatzMaterial>\n')

  # get header from old density file
  header = extract_old_header(infile)
  out.write(header)
  
  out.write('<set id="refined">\n')

  new=[0]*8*len(res)
  print "length of new vector = " + str(len(new))

  num=int(math.pow(len(res)+1, 1.0/3.0))
  num2=int(math.pow(num, 2))
  
  print "copying " + str(num) + "^3 values"
  print 'length of new vector = ' + str(len(new))

  for r in range(num):
    for c in range(num):
      for s in range(num):
        index = s*num2 + r*num + c
        val =  res[index]

        start = 8*s*num2 + 4*r*num + 2*c
        new[start + 0] = val
        new[start + 1] = val
        new[start + 0 + 2*num] = val
        new[start + 1 + 2*num] = val

        start = start + 4*num2
        new[start + 0] = val
        new[start + 1] = val
        new[start + 0 + 2*num] = val
        new[start + 1 + 2*num] = val

  for i in range(len(new)):
    out.write('  <element nr="' + str(i+1)   + '" type="density" design="' + str(new[i]) + '"/>\n')

  out.write('</set>\n')
  out.write('</cfsErsatzMaterial>\n')
  out.close()

  return

# Assuming there is a sequence of CFS runs by a paramerer study. This method finds the closes
# valid run
# @param filename: the filename with a special key for the running number. e.g. 'j_near_f_$_ah_60.info.xml'
# @param key: the key withing filename, here $
# @param start: the first number to search for the file (inclusive!)
# @param stop: the last number to check (inclusive). Smaller or larger start!
# @return: the number within start and end 
# @raise exception: if no info.xml with status 'finished' found    
def find_closes_info_xml(filename, key, start, end):
  checked = 0
  failed  = 0
  # we can search upwards or downwards
  direction = cond(end < start, -1, 1)
  for i in range(start, end + direction, direction):
    name = string.replace(filename, key, str(i))
    
    print  "check " + name + "\n"
    
    if not os.path.exists(name):
      continue
    
    doc = libxml2.parseFile(name)
    xml = doc.xpathNewContext()
    
    status = xpath(xml, "/cfsInfo/@status")
    print "file " + name + " has status " + status

    if status == "finished":
      return i
    else:
      failed = failed + 1
    
    checked = checked + 1
    
  raise RuntimeError("Could not find close info.xml to '" + filename + "': from " + str(checked) + " files " + str(failed) + " were not finished")       


# a function to restore missing element information in a density file 
# the element density information of the old file is read in and for missing 
# element numbers, the lower density from the parameter "lower" is restored
# @infile density file with missing element information
# @outfile file where the restored information is to be written 
# @targetlength how many elements does the target file have?
#    WARNING: this assumes that the first element number is 1!
# @lower the value of the lower density of the missing elements
def restore_lower_densities(infile, outfile, targetlength, lower):
  header = extract_old_header(infile)

  infi = open(infile, "r")
  counter = 0 

  out = open(outfile, "w")
  out.write('<?xml version="1.0"?>\n')
  out.write('<cfsErsatzMaterial>\n')
  out.write(header)
  out.write('\n<set id="restored">\n')

  for event, element in etree.iterparse(infi):
    if element.tag == "element":
      # traverse the elements
      # first we get the number and the design
      nr = int(element.get("nr")) 
      design = float(element.get("design"))
      counter = counter + 1 
      # if we didn't skip anything, just take the values and go for it
      if nr == counter:
        out.write('  <element nr="' + str(counter)   + '" type="density" design="' + str(design) + '"/>\n')
      else:
        # we missed some entries: restore them to lower bound
        missing = nr - counter
        for n in range(missing):
          out.write('  <element nr="' + str(counter)   + '" type="density" design="' + str(lower) + '"/>\n')
          counter = counter + 1 

        # the one design is missing, must insert here
        counter = nr
        out.write('  <element nr="' + str(counter)   + '" type="density" design="' + str(design) + '"/>\n')

  infi.close()

  # finished everything except one: if there are some elements missing 
  # at the end, we have to insert them too!
  if counter < targetlength:
    missing = targetlength - counter
    counter = counter + 1 
    for n in range(missing):
      out.write('  <element nr="' + str(counter)   + '" type="density" design="' + str(lower) + '"/>\n')
      counter = counter + 1 


  out.write('</set>\n')
  out.write('</cfsErsatzMaterial>\n')
  out.close()
  return

# sometimes when writing a new density file it is good to preserve the old header information 
# this function returns the header of a density file 
# @param infile the density file containing the header
def extract_old_header(infile):
  infi = open(infile, "r")
  for event, element in etree.iterparse(infi):
    if element.tag == "header":
      print "found header"
      # we use the header from the original file
      header = etree.tostring(element)
      break
  infi.close()
  return header


# this is a helper strucht for readLumpedMechDisplacement
class LumpedMechDisplacement:
  def __init__(self, elemNum, x, y, z):
    self.elemNum = elemNum
    self.x       = x
    self.y       = y
    self.z       = z
    
  def toString(self):
    return "elem=" + str(self.elemNum) + " x=" + str(self.x) + ", " + str(self.y) + ", " + str(self.z) 


# read the info.xml file and extracts the lumpedMechDisplacement from an eigenfrequency analysis
# @param info_xml the file name
# @param region the region name
# @return an list of tuples. the first item is the frequency, the next is a list of  LumpedMechDisplacement elements.
#        Multiple frequencies are not possible but the last one is used.
def readLumpedMechDisplacement(info_xml, region):
  doc = libxml2.parseFile(info_xml)
  xml = doc.xpathNewContext()
  # detect multiple eigenfrequencies
  last_freq = "-1.0"
  result = []
  for step in range(0,5000):
    xpath_base = "//result[@data='lumpedMechDisplacement'][@location='" + region + "']/item[@step='" + str(step) + "']"
    res = xml.xpathEval(xpath_base + "/@step_val")
    if len(res) == 0 or res[0].getContent() == last_freq:
      continue
    last_freq = res[0].getContent()
    item = []
    item.append(float(last_freq))
    # read data
    e = xml.xpathEval(xpath_base + "/@id")
    x = xml.xpathEval(xpath_base + "/@x")
    y = xml.xpathEval(xpath_base + "/@y")
    z = xml.xpathEval(xpath_base + "/@z")
    # assume len(res) = len(x) = len(y), eventually also len(z)
    data = []
    for i in range(len(e)):
      e_val = int(e[i].getContent())
      x_val = float(x[i].getContent())
      y_val = float(y[i].getContent())
      z_val = cond(len(z) == 0, 0.0, float(z[i].getContent()))
      data.append(LumpedMechDisplacement(e_val, x_val, y_val, z_val))
    item.append(data)
    result.append(item)
  return result  

# create an interpolated density file from a readLumpedMechDisplacement() result for positive z displacement.
# below the first frequency we assume full material, obove the highest frequency an exception is thrown
# @param data is the result from readLumpedMechDisplacement
# @param f the frequency to interpolate
# @param density_file the filename where the density file is written to
def interpolateLumpedMechDisplacementAsDensity(data, f, density_file):
  # find upper and lower frequency
  lower = 0.0
  upper = -1.0
  
  # the data array of interest
  lower_data = [] # stays empty if f is below first eigenfrequency, assume 1.0 values then
  upper_data = []
  
  for i in range(len(data)):
    test = data[i][0]
    if test < f:
      lower = test
      continue 
    else:
      upper  = test
      upper_data = data[i][1]
      if i > 0:
        lower_data = data[i-1][1]
      break
  if upper == -1.0:
    raise RuntimeError("Cannot find frequency " + str(f) + " within lumpedMechDisplacement data set from f=" + str(data[0][0]) + " to f=" + str(data[len(data)-1][0]))

  # determine interpolation parameter, normalized to 1
  beta = (f - lower) / (upper - lower)
  alpha = 1.0 - beta
  print "f=" + str(f) + " lower=" + str(lower) + " upper=" + str(upper) + " alpha=" + str(alpha) + " beta=" + str(beta) + "\n"

  # preliminary result, before normalizing to 1.0
  tmp = []
  # we don't know if the positive or negative real parts dominates. the firs mode has the same sign which might be neg
  max_v = -1e30   
  min_v = +1e30
  for i in range(len(upper_data)):
    lower_value = 1.0 # there is no read condition operator in python :( 
    if len(lower_data) > 0:
      lower_value = lower_data[i].z
    upper_value = upper_data[i].z
    v = alpha * lower_value + beta * upper_value
    # print "lv=" + str(lower_value) + " uv=" + str(upper_value) + " v=" + str(v) + "\n"
    tmp.append(v)
    max_v = cond(v > max_v, v, max_v)
    min_v = cond(v < min_v, v, min_v)
    
  # normalize and eliminate the right values
  for i in range(len(tmp)):
    if abs(max_v) >= abs(min_v) and max_v > 0:
      # print "a max=" + str(max_v) + " min=" + str(min_v) + " t=" + str(tmp[i]) + " -> " + str(cond(tmp[i] > 0, tmp[i] / max_v, 0.0)) + "\n" 
      tmp[i] = cond(tmp[i] > 0, tmp[i] / max_v, 0.0)
    else:
      # print "b max=" + str(max_v) + " min=" + str(min_v) + " t=" + str(tmp[i]) + " -> " + str(cond(tmp[i] < 0, tmp[i] / min_v, 0.0)) + "\n"      
      tmp[i] = cond(tmp[i] < 0, tmp[i] / min_v, 0.0)
    if tmp[i] < 0.0 or tmp[i] > 1.0:
      raise RuntimeError("invalid density " + str(tmp[i]))
    # prevent 0
    if tmp[i] < 1e-7:
      tmp[i] = 1e-7
      
  # write the stuff    
  out = open(density_file, "w")
  out.write('<?xml version="1.0"?>\n')
  out.write('<cfsErsatzMaterial>\n')
  out.write('  <header>\n')
  out.write('    <design initial="0.5" lower="1e-3" name="density" upper="1"/>\n')
  out.write('    <transferFunction application="mech" design="density" param="1" type="simp"/>\n')
  out.write('  </header>\n')
  out.write('  <set id="f_' + str(f) + '_' + str(alpha) + '*' + str(lower) + '+' + str(beta) + '*' + str(upper) + '_min_' + str(min_v) + '_max_' + str(max_v) + '">\n')
  for i in range(len(tmp)):
    out.write('    <element nr="' + str(upper_data[i].elemNum) + '" type="density" design="' + str(tmp[i]) + '"/>\n')
  
  out.write('  </set>\n')

  out.write(' </cfsErsatzMaterial>\n')
  out.close()


# helper for external use - create a 2D "window"
# strength is the size of the window in fraction. .5 is maximum
# return 2d array
def make2DWindow(divider, strength, lower):
  ret = numpy.zeros((divider, divider))
  
  if strength <= 0.01 or strength > 0.48:
    raise RuntimeError("invalid strength")
  
  #fraction
  f = strength
  # counter fraction
  cf = 1.0 - f
  
  for ii in range(divider):
    for jj in range(divider):
      i = float(ii)
      j = float(jj)
  
      ip = i / (divider-1)
      jp = j / (divider-1)
      
      #print str(ip) + " - " + str(jp) + " f " + str(f) + " cf " + str(cf)
      if ip > f and ip < cf and jp > f and jp < cf:
        ret[i][j] = lower
      else:
        ret[i][j] = 1.0

  return ret    
      
# extrudes an 2D array to the third dimenstion
def extrude(data_2d):
  edge = data_2d.shape[0] 
  if edge <> data_2d.shape[1]:
    raise RuntimeException("require quadratic input")
  
  ret = numpy.zeros((edge, edge, edge))
  
  for i in range(edge):
    for j in range(edge):
      val = data_2d[i][j]
      for k in range(edge):
        ret[i][k][j] = val
        
  return ret
