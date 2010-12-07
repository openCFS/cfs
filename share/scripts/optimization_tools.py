# -*- coding: utf-8 -*-
# This collects some tool routines for optimization

import libxml2
import numpy
import numpy.linalg
import math
import os
from lxml import etree

from cfs_utils import *
from distutils.command.build_scripts import first_line_re

## Read an arbitrary density file as NDArray
# Uses the <mesh x="30" y="20" z="1"/> element in the header of the density file
# if not the whole domain is design domain, the data is read as 1D array
# @param elemnr if False the design is read, otherwise nr
# @param x, y, z optional mesh size in case it is not given in the density file. Note, the smallest number is 1, not 0!!
def read_density(filename, elemnr=False, x = None, y = None, z = None):
  vals = read_density_as_vector(filename, elemnr)

  tree = etree.parse(filename, etree.XMLParser(remove_comments=True))
  root = tree.getroot()
  if x == None:
    x = int(root.xpath("//mesh/@x")[0])
  if y == None:  
    y = int(root.xpath("//mesh/@y")[0])
  if z == None:
    z = int(root.xpath("//mesh/@z")[0])

  assert(x > 0 and y > 0 and z > 0)  

  # density files where not the whole domain is design domain are read and re-written
  # as 1D arrays
  dim = cond(y > 1, cond(z > 1, 3, 2), 1)
  ret = 0
     
  if len(vals) > x*y*z:
    raise RuntimeError("density mesh information x=" + str(x) + " y=" + str(y) + " z=" + str(z)\
                        + " does not match " + str(len(vals)) + " elements in " + filename) 
  if len(vals) == x*y*z:
    # full array!
    if dim == 2:
      ret = numpy.zeros((x, y))
    else:
      ret = numpy.zeros((x, y, z))

    # copy data from linear list
    for k in range(z):
      for j in range(y):
        for i in range(x):
          idx = int(x*y*k + x*j + i)
          # print "i=" + str(i) + " j=" + str(j) + " k=" + str(k) + " idx=" + str(idx)
          setNDArrayEntry(ret, i, j, k, vals[idx])

  
  if len(vals) < x*y*z:
    # we need to be 1D
    print "read density file '" + filename + "' with " + str(len(vals)) + " element smaller x=" + str(x) \
         + " y=" + str(y) + " z=" + str(z) + " mesh as " + cond(elemnr, "elem nr", "design")
    x = len(vals)
    ret = numpy.zeros((x))
    for i in range(x):
      ret[i] = vals[i]    
        
  return ret


  
## Reads a density.xml file as vector
# @param filename from which the last 'set' is used
# @param elemnr if False the design is read on true nr
def read_density_as_vector(filename, elemnr=False):
  if not os.path.exists(filename):
    raise RuntimeError("file '" + filename + "' doesn't exist")
  
  tree = etree.parse(filename, etree.XMLParser(remove_comments=True))
  
  root = tree.getroot()
  sett = root.xpath("//set[last()]")[0]
  # print "reading set with id = " + sett.get("id")
  
  length = len(sett)
  
  attribute = cond(elemnr, "nr", "design")
  # print "check for attribute " + attribute
  counter = 0
  vals=[0]*length
  for element in sett:
    # traverse the elements and get the design
    vals[counter] = float(element.get(attribute))
    counter = counter + 1
  
  # print "found " + str(length) + " elements, read " + str(counter) + " elements"
  
  return vals
  
## write the data to a density.xml file
# @param data_inp a ndata array (1D, 2D or 3D) or a list of data
# @param setname_inp the name of the set or a list of setnames
# @param elemnr if set, the element number is taken from this elemnr ndarray.
#               The data can be obtained from read_density(...,elemnr=True)
def write_density_file(filename, data_inp, setname_inp, elemnr=None):
  # check if we deal with lists or not
  data_list = []
  setname_list = []
  if type(data_inp) is list:
    data_list = data_inp
    setname_list = setname_inp
  else:
    data_list.append(data_inp)
    setname_list.append(setname_inp)

  out = open(filename, "w")
  out.write('<?xml version="1.0"?>\n')
  out.write('<cfsErsatzMaterial>\n')
  out.write('  <header>\n')

  data    = data_list[0]
  x, y, z = getDim(data)
  out.write('    <mesh x="' + str(x) + '" y="' + str(y) + '" z="' + str(z) + '"/>\n')  
  out.write('    <design initial="0.5" lower="1e-3" name="density" region="mech" upper="1"/>\n')
  out.write('    <transferFunction application="mech" design="density" param="1" type="simp"/>\n')
  out.write('  </header>\n')

 
  for i in range(len(data_list)):
    data    = data_list[i]
    setname = setname_list[i]
    out.write('  <set id="' + setname + '">\n')
    dim = data.ndim
    x, y, z = getDim(data)
    nr = 1

    for k in range(z):
      for j in range(y):
        for i in range(x):    
           val = getNDArrayEntry(data, i, j, k)
           if elemnr <> None:
             nr = int(getNDArrayEntry(elemnr, i, j ,k))
           
           # print " i=" + str(i) + " j=" + str(j) + " k=" + str(k) + " idx=" + str(nr)
           out.write('    <element nr="' + str(nr) + '" type="density" design="' + str(val) + '"/>\n')
           nr = nr + 1       
         
    out.write('  </set>\n')

  out.write(' </cfsErsatzMaterial>\n')
  out.close()


## replaces the element numbers by new element numbers.
# @param org ndarray of element numbers from read_density(,elemnr=True)
# @param map element mapping as tupel list from cfs_grid.map_elements()
# @return again an ndarray which can be used as parameter for write_density_file()
def apply_elmennr_mapping(org, map):
  x, y, z = getDim(org)
  assert(x*y*z == len(map))

  result = numpy.zeros((x, y, z))

  for k in range(z):
    for j in range(y):
      for i in range(x):    
         val = getNDArrayEntry(org, i, j, k)
         # be save and search n^2 -> if too slow make it faster :)
         found = False
         for p in range(len(map)):
           if val == map[p][0]:
             setNDArrayEntry(result, i, j, k, int(map[p][1]))
             found = True
             break
         if not found:
           raise RuntimeException("could not find elemnr=" + str(val) + " from x=" + str(val) + " y=" + str(y) + " z=" + str(z))
            
  return result
  


## evaluates the physical volume fraction
# @param data 2D/3D data
# @param penalty 1 for no penalty
def physical_volume(data, penalty):

  dim = data.ndim
  x = data.shape[0]
  y = data.shape[1]
  z = 1
  if dim >= 3:
    z = data.shape[2]
  vol = 0.0
  
  for i in range(x):
    for j in range(y):
      for k in range(z):
        org = getNDArrayEntry(data, i, j, k)
        vol = vol + pow(org, penalty)
        
  return vol / data.size              


## apply a threshold filter on a 2D/3D matrix
# @param data original 2D/3D data
# @param treshold all smaller threshold becomes min, otherwise max
# @return a new data array 
def threshold_filter(data, threshold, min, max):
  
  dim = data.ndim
  x, y, z = getDim(data)
  res = numpy.copy(data)
    
  # handle the case that threshold is smaller than min
  barrier = cond(threshold < min, min, threshold)
    
  for i in range(x):
    for j in range(y):
      for k in range(z):
        org = getNDArrayEntry(data, i, j, k)
        val = cond(org < barrier, min, max)
        setNDArrayEntry(res, i, j, k, val)        

  return res  


## threshold towards a given final volume
# @param data original 2d/3d data
# @param min niminum final density in the result, ignores penalty, max = 1
# @param target volume fraction to be reached
# @param material_penalty penalty to interpret original data
# @return: the new data array as first value and the threshold as second value
def auto_threshold_filter(data, min, target, material_penalty):
  
  res = numpy.copy(data)

  lower = min
  upper = 1
  while upper - lower > 1e-14:
    mid = lower + 0.5 *(upper - lower)
    val = physical_volume(threshold_filter(data, mid, min, 1), material_penalty)
    print " lower=%15.15g mid=%15.15g upper=%15.15g val=%15.15g " % (lower , mid , upper , val)
    # print " lower=" + str(lower) + " mid=" + str(mid) + " upper=" + str(upper) + " val=" + str(val)
    if val > target:
      lower = mid
    else:
      upper = mid
  
  return threshold_filter(data, lower, min, 1), lower   

## interpolates a floating point coordinate within an numpy.array
# @param data an numpy.array
# @param position an arbitrary 2D/ 3D array from where we want to get the data
# @param quiet error if the data does not exist or evaluate periodic
# @param periodic if data does not exist it is extended periodically or default_val is taken
# @param default_value
def interpolate(data, position, quiet=False, periodic=True, default_value=0.001):

  #this is the corrected position pointer
  corrected = numpy.copy(position)
  
  for i in range(len(position)):
    x = int(position[i])
    # are we in the correct data range?
    if x > 0 and x < data.shape[i]:
      corrected[i] = int(x)
      continue
    if not quiet:
      raise RuntimeError("position " + str(b) + " not within data " + str(data.shape))
    if not periodic:
      return default_value
    # periodic case
    if x < 0:
      corrected[i] = data.shape[i] + x
    else:
      corrected[i] = x - data.shape[i]
  
  if len(corrected) == 2:
    val = data[corrected[0], corrected[1]]
  else:
    val = data[corrected[0], corrected[1], corrected[2]]    

  return val 

## rotates an ndarray by an arbitrary rotation matrix
# for all intermediate materials there will be data with no origin,
# these data is considered periodic
def rotate(data, rot):
  result = numpy.copy(data)
  x, y, z = getDim(data)
  
  center = [x * 0.5, y * 0.5, z * 0.5]
  
  for i in range(x):
    for j in range(y):
      for k in range(z):
        # algorithm for each target point rotate backwards to the original point
        target = [i-center[0], j - center[1], k - center[2]]
        brt = numpy.dot(rot.T, target) # back_rot_target
        # denormalize the center
        origin = brt + center
        # print "rotate " + str(target) + " back to " + str(brt) + " -> " + str(origin) 
        val = interpolate(data, origin, quiet=True, periodic=True)
        setNDArrayEntry(result, i, j, k, val) 
         
  return result       

## rotate a 3d matrix
# exchange x and y
def rotate_matrix_x_y(data):
  
  dim = data.ndim
  edge = data.shape[0]
  res = numpy.zeros((edge, edge, edge))
  
  for i in range(edge):
    for j in range(edge):
      for k in range(cond(dim == 2, 1, edge)):
        val = getNDArrayEntry(data, i, j, k)
        setNDArrayEntry(res, j, i, k, val)        

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
        val = getNDArrayEntry(data, i, j, k)
        setNDArrayEntry(res, k, j, i, val)        

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
              val = getNDArrayEntry(data, i, j, k)
              setNDArrayEntry(res, x_base + i, y_base + j, z_base + k, val)

  return res                     

## handles arbitrary 2d and 3d density files, doubles in each dimension 
# @param infile the density file with the densities to be blown upper
# @param outfile name of the output file
# @param returns also the new ndata array
def refine_density(infile, outfile):
  org = read_density(infile)

  x, y, z = getDim(org)
  dim = org.ndim
  # we cannot handle 3D density files with z is one layer as these are identified as 2D :(
  out = numpy.zeros((x*2, y*2, cond(z == 1, 1, z*2)))

  for i in range(x):
    for j in range(y):
      for k in range(z):
        val = getNDArrayEntry(org, i, j, k)
        # in 2D the z-component is ignored and we set the value twice
        setNDArrayEntry(out, i*2 + 0, j*2 + 0, k*2 + 0, val)
        setNDArrayEntry(out, i*2 + 0, j*2 + 0, k*2 + 1, val)
        setNDArrayEntry(out, i*2 + 0, j*2 + 1, k*2 + 0, val)
        setNDArrayEntry(out, i*2 + 0, j*2 + 1, k*2 + 1, val)
        setNDArrayEntry(out, i*2 + 1, j*2 + 0, k*2 + 0, val)
        setNDArrayEntry(out, i*2 + 1, j*2 + 0, k*2 + 1, val)
        setNDArrayEntry(out, i*2 + 1, j*2 + 1, k*2 + 0, val)
        setNDArrayEntry(out, i*2 + 1, j*2 + 1, k*2 + 1, val)
        
  write_density_file(outfile, out, "refined")
  return out 

# Assuming there is a sequence of CFS runs by a parameter study. This method finds the closes
# valid run
# @param filename: the filename with a special key for the running number. e.g. 'j_near_f_$_ah_60.info.xml'
# @param key: the key within filename, here $
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



# checks if file is a valid density file
# simple and stupid check, maybe better version in the future
# @param infile the density file to be checked
def is_valid_density_file(infile):
  ersatzfound = False
  headerfound = False
  setfound = False

  infi = open(infile, "r")
  for event, element in etree.iterparse(infi):
    if element.tag == "cfsErsatzMaterial":
      ersatzfound = True
    if element.tag == "header":
      headerfound = True
    if element.tag == "set":
      setfound = True

    # if all is found, then we break
    if ersatzfound and headerfound and setfound:
      break
  infi.close()

  return ersatzfound and headerfound and setfound


# do an ascii print of the density data
def ascii_print(data, threshold):
  x, y, z = getDim(data)
  assert()
