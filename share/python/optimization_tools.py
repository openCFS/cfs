# This collects some tool routines for optimization

import numpy
import numpy.linalg
import math
import os
import sys
from lxml import etree
from PIL import Image
from cfs_utils import *

# dump some information about a density file
def print_design_info(filename, attribute, set = None, fill = None):
  try: 
    dens = read_density(filename, attribute, set=set, fill=fill)
    print("for attribute '" + attribute + "' min=" + str(numpy.amin(dens)) + " max=" + str(numpy.amax(dens)))
  except:  
    print("failed to read '" + attribute + "'")
    #print  sys.exc_info()



# # Read an arbitrary density file as NDArray
# Uses the <mesh x="30" y="20" z="1"/> element in the header of the density file
# if not the whole domain is design domain, the data is read as 1D array
# @param attribute the scalar attribute: "design" (default), "physical", "nr"
# @param x, y, z optional mesh size in case it is not given in the density file. Note, the smallest number is 1, not 0!!
# @param set optional name of set, if not given use default
# @param fill if the density is a subset of the mesh return 1D or fill with fill_value if fill is set
def read_density(filename, attribute="design", x=None, y=None, z=None, set=None, fill = None):
  vals = read_density_as_vector(filename, attribute, set)
  mx, my, mz = read_mesh_info(filename, silent=True)

  if x == None:
    x = mx
  if y == None:  
    y = my
  if z == None:
    z = mz
  
  if x == None and y == None and z == None:
    x = len(vals)
    y = 1
    z = 1

  assert(x > 0 and y > 0 and z > 0)  
  
  # density files where not the whole domain is design domain are read and re-written
  # as 1D arrays
  dim = cond(y > 1, cond(z > 1, 3, 2), 1)
  ret = numpy.zeros((x, y)) if dim == 2 else numpy.zeros((x, y, z)) # might be overwritten to 1D 
       
  if len(vals) > x * y * z:
    raise RuntimeError("density mesh information x=" + str(x) + " y=" + str(y) + " z=" + str(z)\
                        + " does not match " + str(len(vals)) + " elements in " + filename) 
  if len(vals) == x * y * z:

    # copy data from linear list
    for k in range(z):
      for j in range(y):
        for i in range(x):
          idx = int(x * y * k + x * j + i)
          # print "i=" + str(i) + " j=" + str(j) + " k=" + str(k) + " idx=" + str(idx)
          setNDArrayEntry(ret, i, j, k, vals[idx])

  
  if len(vals) < x * y * z:
    if fill == None:
      # we need to be 1D
      print("read density file '" + filename + "' with " + str(len(vals)) + " element smaller x=" + str(x) \
           + " y=" + str(y) + " z=" + str(z) + " mesh as " + attribute)
      x = len(vals)
      ret = numpy.zeros((x)) # overwrite full array
      for i in range(x):
        ret[i] = vals[i]    
    else:
      print("fill with " + str(fill) + " values '" + attribute + "' from non complete file '" + filename \
           + "' with " + str(len(vals)) + " elements for x=" + str(x) + " y=" + str(y) + " z=" + str(z))

      ret += fill # was zeros, fill with default value
      num = read_density_as_vector(filename, 'nr', set)
      nt = x * y * z
      for c in range(len(num)):
        n = int(num[c] -1)
        coords = numpy.unravel_index(n,(z,y,x))
        setNDArrayEntry(ret, coords[2], coords[1], coords[0], vals[c]) # 'nr' is one base in general
        
  return ret

# # read arbitrary multi-design density file as numpy array
def read_multi_design(filename, design1, design2=None, design3=None, design4=None, design5 = None, design6 = None, matrix=False, attribute="design"):
  if not os.path.exists(filename):
    raise RuntimeError("file '" + filename + "' doesn't exist")
  tree = etree.parse(filename, etree.XMLParser(remove_comments=True))
  root = tree.getroot()
  if matrix:
    x, y, z = read_mesh_info(filename, silent=True)
  
    if x == None and y == None and z == None:
      x = 1
      y = 1
      z = 1
  
    assert(x > 0 and y > 0 and z > 0)  
  sett = root.xpath("//set[last()]")[0]
  #print(len(sett))
  
  designs = 1
  if design2:
    designs = 2
  if design3:
    designs = 3  
  if design4:
    designs = 4
  if design5:
    designs = 5
  if design6:
    designs = 6
    
  length = len(sett) / designs
  
  out = numpy.zeros((length, designs))
  for element in sett:
    nr = int(element.get("nr"))
    type = element.get("type")
    idx = -1
    if type == design1:
      idx = 0
    if design2 and type == design2:
      idx = 1
    if design3 and type == design3:
      idx = 2
    if design4 and type == design4:
      idx = 3
    if design5 and type == design5:
      idx = 4
    if design6 and type == design6:
      idx = 5
    if not idx == -1:
      tmp = element.get(attribute)
      if tmp is None:
        print("Could not read '" + attribute + "' for design " + type + "! Fallback to 'design'.")
        tmp = element.get("design")
      des = float(tmp)
      out[nr - 1, idx] = des
  if matrix:
    output = numpy.zeros((x, y, z, designs))
    for t in range(designs):
      count = 0
      for k in range(z):
        for j in range(y):
          for i in range(x):
            output[i][j][k][t] = out[count][t]
            count += 1
    out = output
  return out
  
## returns all set ids from a density xml
# return a list of string like stuff
def read_set_ids(filename):
  tree = etree.parse(filename)
  root = tree.getroot()
  sett = root.xpath("//set/@id")
  return sett
  
  
## tests for an attribute in the density.xml
# @param attribute you might want to check for physical which is not present when generated via create_density.py
# @return True if the attribute is given
def test_density_xml_attribute(filename, attribute, set = None):
  tree = etree.parse(filename, etree.XMLParser(remove_comments=True))      
  query = '//set[last()]' if set is None else '//set[@id="' + set + '"]'
  query += '/element/@' + attribute
  s = tree.getroot().xpath(query)
  return len(s) > 0
    
  
## parse the 'mesh' information of the header of a density.xml file
# @param silent if True and no mesh info is found return None, None, None otherwise raise exception
# return x, y, z as ints
def read_mesh_info(filename, silent):
  tree = etree.parse(filename)
  root = tree.getroot()
  mesh = root.xpath('//cfsErsatzMaterial/header/mesh')
  
  if len(mesh) == 0:
    if not silent:
      raise RuntimeError("file '" + filename + "' has no mesh information")
    return None, None, None
  
  else:
    assert(len(mesh) == 1)
    nx = int(mesh[0].get("x"))
    ny = int(mesh[0].get("y"))
    nz = int(mesh[0].get("z"))
    # return int(mesh[0].get("x")), int(mesh[0].get("y")), int(mesh[0].get("z"))   
    return nx, ny, nz    
       
  
# # Reads a density.xml file as vector
# @param filename from which the last 'set' is used
# @param attribute the scalar attribute: "design" (default), "physical", "nr"
# @param set optionally give a set name (string) when not given the last is used 
def read_density_as_vector(filename, attribute="design", set=None):
  if not os.path.exists(filename):
    raise RuntimeError("file '" + filename + "' doesn't exist")
  
  tree = etree.parse(filename, etree.XMLParser(remove_comments=True))
  root = tree.getroot()
  query = '//set[last()]/element' if set is None else '//set[@id="' + set + '"]/element'
  sett = root.xpath(query)
  
  # print "check for attribute " + attribute
  counter = 0
  vals = [0] * len(sett)
  for element in sett:
    # traverse the elements and get the design
    vals[counter] = float(element.get(attribute))
    counter = counter + 1
  
  # print "found " + str(length) + " elements, read " + str(counter) + " elements"
  
  return vals
  
# # write the data to a density.xml file
# @param data_inp a ndata array (1D, 2D or 3D) or a list of data
# @param setname_inp the name of the set or a list of setnames
# @param elemnr if set, the element number is taken from this elemnr ndarray.
#               The data can be obtained from read_density(...,elemnr=True)
def write_density_file(filename, data_inp, setname_inp="set", param=0, elemnr=None):
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

  data = data_list[0]
  x, y, z = getDim(data)
  out.write('    <mesh x="' + str(x) + '" y="' + str(y) + '" z="' + str(z) + '"/>\n')  
  out.write('    <design initial="0.5" lower="1e-3" name="density" region="mech" upper="1"/>\n')
  if param > 0:
    out.write('    <transferFunction application="mech" design="density" param="' + str(param) + '" type="simp"/>\n')
  else:
    out.write('    <transferFunction application="mech" design="density" param="1" type="simp"/>\n')
  out.write('  </header>\n')

 
  for i in range(len(data_list)):
    data = data_list[i]
    setname = setname_list[i]
    out.write('  <set id="' + setname + '">\n')
    dim = data.ndim
    x, y, z = getDim(data)
    nr = 1

    for k in range(z):
      for j in range(y):
        for i in range(x):    
           val = getNDArrayEntry(data, i, j, k)
           if elemnr is not None:
             nr = int(getNDArrayEntry(data, i, j , k))
            
           # print " i=" + str(i) + " j=" + str(j) + " k=" + str(k) + " idx=" + str(nr)
           if param > 0:
            out.write('    <element nr="' + str(nr) + '" type="density" design="' + str(val) + '" physical="' + str(val ** param) + '"/>\n')
           else:
            out.write('    <element nr="' + str(nr) + '" type="density" design="' + str(val) + '"/>\n')
           if elemnr is None:
             nr = nr + 1       
         
    out.write('  </set>\n')

  out.write(' </cfsErsatzMaterial>\n')
  out.close()


# # write simple multi-design density files
def write_multi_design_file(filename, data, designs, elemnr=None):

  # assert(data.shape[1] == len(designs))
  out = open(filename, "w")
  out.write('<?xml version="1.0"?>\n')
  out.write('<cfsErsatzMaterial>\n')
  out.write('  <header>\n')
  out.write('  </header>\n')
  out.write('  <set id="optimization_tools.py">\n')
    
  for d in range(len(designs)):
    for e in range(len(data)):
      enr = e + 1 if elemnr == None else int(elemnr[e])
      out.write('    <element nr="' + str(enr) + '" type="' + designs[d] + '" design="' + str(data[e, d]) + '"/>\n')
  out.write('  </set>\n')
  out.write(' </cfsErsatzMaterial>\n')
  out.close()

## reads partial domain within a full array. E.g. if the design domain is circular
# visualize via img = Image.fromarray(a, 'RGB'), img.show()
def read_density_as_full_array(filename, attribute='design', fill=0.0, set = None):
  msh  = read_mesh_info(filename, silent = False)
  vals = read_density_as_vector(filename, attribute, set)
  nrs  = read_density_as_vector(filename, 'nr', set)

  assert(msh[2] == 1) # implement 3D if you need it
  
  a = numpy.ones(msh[0:2]) * fill  
  
  nx = msh[0]
  ny = msh[1]
  
  assert(len(nrs) <= nx * ny)
  
  for i in range(len(nrs)):
    n = int(nrs[i]) -1 # the element number is 1 based!
    # assume a lexicographic order
    y = int(n / ny)
    x = n % ny
    #print 'i=' + str(i) + ' y=' + str(y) + ' x=' + str(x) 
    
    a[y, x] = vals[i]
      
  return a

# # replaces the element numbers by new element numbers.
# @param org ndarray of element numbers from read_density(,elemnr=True)
# @param map element mapping as tupel list from cfs_grid.map_elements()
# @return again an ndarray which can be used as parameter for write_density_file()
def apply_elmennr_mapping(org, map):
  x, y, z = getDim(org)
  assert(x * y * z == len(map))

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
           raise RuntimeError("could not find elemnr=" + str(val) + " from i=" + str(i) + " j=" + str(j) + " k=" + str(k))
            
  return result
  




# # evaluates the physical volume fraction
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


# # apply a threshold filter on a 2D/3D matrix
# @param data original 2D/3D data or vector
# @param treshold all smaller threshold becomes min, otherwise max
# @return new data of input type 
def threshold_filter(data, threshold, min, max):
  
  # handle the case that threshold is smaller than min
  barrier = numpy.max((threshold, min))
  res = None  
  
  if type(data) == list:
    res = [0.0] * len(data)
    for i in range(len(data)):
      res[i] = max if data[i] >= barrier else min
  else:
    
    dim = data.ndim
    x, y, z = getDim(data)
    res = numpy.copy(data)
  
    for i in range(x):
      for j in range(y):
        for k in range(z):
          org = getNDArrayEntry(data, i, j, k)
          val = cond(org < barrier, min, max)
          setNDArrayEntry(res, i, j, k, val)        

  return res  


# # threshold towards a given final volume
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
    mid = lower + 0.5 * (upper - lower)
    val = physical_volume(threshold_filter(data, mid, min, 1), material_penalty)
    print " lower=%15.15g mid=%15.15g upper=%15.15g val=%15.15g " % (lower , mid , upper , val)
    # print " lower=" + str(lower) + " mid=" + str(mid) + " upper=" + str(upper) + " val=" + str(val)
    if val > target:
      lower = mid
    else:
      upper = mid
  
  return threshold_filter(data, lower, min, 1), lower   

# # interpolates a floating point coordinate within an numpy.array
# @param data an numpy.array
# @param position an arbitrary 2D/ 3D array from where we want to get the data
# @param quiet error if the data does not exist or evaluate periodic
# @param periodic if data does not exist it is extended periodically or default_val is taken
# @param default_value
def interpolate(data, position, quiet=False, periodic=True, default_value=0.001):

  # this is the corrected position pointer
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

# # rotates an ndarray by an arbitrary rotation matrix
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
        target = [i - center[0], j - center[1], k - center[2]]
        brt = numpy.dot(rot.T, target)  # back_rot_target
        # denormalize the center
        origin = brt + center
        # print "rotate " + str(target) + " back to " + str(brt) + " -> " + str(origin) 
        val = interpolate(data, origin, quiet=True, periodic=True)
        setNDArrayEntry(result, i, j, k, val) 
         
  return result       

# # flip a 2d matrix
# exchange x and y
def flip_2d_matrix(data):
  
  assert(data.ndim == 2)
  
  x = data.shape[0]
  y = data.shape[1]
  res = numpy.zeros((y, x))
  
  for i in range(x):
    for j in range(y):
        val = getNDArrayEntry(data, i, j, 0)
        setNDArrayEntry(res, j, i, 0, val)        

  return res  

# # rotate a 3d matrix
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


# # rotate a 3d matrix
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

# # Enlarge a density matrix by multiples of the original data
# @param data a two or three dimensional ndata array
# @param x_times how often to repeat in x direction. 1 = no change
# @param z_times ignored when data is 2D
# @return a new ndata array (might be huge!!)
def enlarge_matrix(data, times_x, times_y=1, times_z=1):

  dim = data.ndim
  x_edge = data.shape[0]
  y_edge = data.shape[1]
  
  if dim == 3:
    z_edge = data.shape[2]
    res = numpy.zeros((x_edge * times_x, y_edge * times_y, z_edge * times_z))
    
    for x in range(times_x):
      for y in range(times_y):
        for z in range(times_z):
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
  else:
    res = numpy.zeros((x_edge * times_x, y_edge * times_y))
    
    for x in range(times_x):
      for y in range(times_y):
        # copy complete block
        x_base = x * x_edge
        y_base = y * y_edge
        
        for i in range(x_edge):
          for j in range(y_edge):
            val = getNDArrayEntry(data, i, j, 0)
            setNDArrayEntry(res, x_base + i, y_base + j, 0, val)
    
    return res     



# # handles arbitrary 2d and 3d density files, doubles in each dimension 
# @param infile the density file with the densities to be blown upper
# @param outfile name of the output file
# @param returns also the new ndata array
def refine_density(infile, outfile, design1=None, design2=None, design3=None, design4=None, design5 = None):
  if design1 is None:
    org = read_density(infile)
    x, y, z = getDim(org)
    ndes = 1
    dim = org.ndim
  else:
    org = read_multi_design(infile, design1, design2, design3, design4, design5, True)
    x, y, z, ndes = getDim(org, True)
    dim = 3
    if z == 1:
      dim = 2
    
  # we cannot handle 3D density files with z is one layer as these are identified as 2D :(
  
  out = numpy.zeros((x * 2, y * 2))

  # In 3D we need to overwrite with 3D array or else setNDArrayEntry will not work properly
  if(dim == 3):
    out = numpy.zeros((x * 2, y * 2, z * 2))
    
  output = numpy.zeros((out.size, ndes))

  # print "debug: x=" + str(x) + ", y=" + str(y) + ", z=" + str(z) + ", dim=" + str(dim)
  for d in range(ndes):
    for i in range(x):
      for j in range(y):
        for k in range(z):
          val = getNDArrayEntry(org, i, j, k, d)
          # in 2D the z-component is ignored and we set the value twice
          setNDArrayEntry(out, i * 2 + 0, j * 2 + 0, k * 2 + 0, val)
          setNDArrayEntry(out, i * 2 + 0, j * 2 + 0, k * 2 + 1, val)
          setNDArrayEntry(out, i * 2 + 0, j * 2 + 1, k * 2 + 0, val)
          setNDArrayEntry(out, i * 2 + 0, j * 2 + 1, k * 2 + 1, val)
          setNDArrayEntry(out, i * 2 + 1, j * 2 + 0, k * 2 + 0, val)
          setNDArrayEntry(out, i * 2 + 1, j * 2 + 0, k * 2 + 1, val)
          setNDArrayEntry(out, i * 2 + 1, j * 2 + 1, k * 2 + 0, val)
          setNDArrayEntry(out, i * 2 + 1, j * 2 + 1, k * 2 + 1, val)
    if design1 is None:
      write_density_file(outfile, out, "refined")
      return out
    else:
      output[:,d] = numpy.ravel(out,order='F')
  
  write_multi_design_file(outfile, output, (design1, design2, design3, design4, design5)[0:ndes])
  return output

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
  failed = 0
  # we can search upwards or downwards
  direction = cond(end < start, -1, 1)
  for i in range(start, end + direction, direction):
    name = string.replace(filename, key, str(i))
    
    print  "check " + name + "\n"
    
    if not os.path.exists(name):
      continue
    
    xml = open_xml(name)
    
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
        out.write('  <element nr="' + str(counter) + '" type="density" design="' + str(design) + '"/>\n')
      else:
        # we missed some entries: restore them to lower bound
        missing = nr - counter
        for n in range(missing):
          out.write('  <element nr="' + str(counter) + '" type="density" design="' + str(lower) + '"/>\n')
          counter = counter + 1 

        # the one design is missing, must insert here
        counter = nr
        out.write('  <element nr="' + str(counter) + '" type="density" design="' + str(design) + '"/>\n')

  infi.close()

  # finished everything except one: if there are some elements missing 
  # at the end, we have to insert them too!
  if counter < targetlength:
    missing = targetlength - counter
    counter = counter + 1 
    for n in range(missing):
      out.write('  <element nr="' + str(counter) + '" type="density" design="' + str(lower) + '"/>\n')
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


# this is a helper struct for readLumpedMechDisplacement
class LumpedMechDisplacement:
  def __init__(self, elemNum, x, y, z):
    self.elemNum = elemNum
    self.x = x
    self.y = y
    self.z = z
    
  def toString(self):
    return "elem=" + str(self.elemNum) + " x=" + str(self.x) + ", " + str(self.y) + ", " + str(self.z) 


# read the info.xml file and extracts the lumpedMechDisplacement from an eigenfrequency analysis
# @param info_xml the file name
# @param region the region name
# @return an list of tuples. the first item is the frequency, the next is a list of  LumpedMechDisplacement elements.
#        Multiple frequencies are not possible but the last one is used.
def readLumpedMechDisplacement(info_xml, region):
  xml = open_xml(info_xml)
  # detect multiple eigenfrequencies
  last_freq = "-1.0"
  result = []
  for step in range(0, 5000):
    xpath_base = "//result[@data='lumpedMechDisplacement'][@location='" + region + "']/item[@step='" + str(step) + "']"
    res = xml.xpath(xpath_base)
    if len(res) == 0 or not 'step_val' in res[0].attrib or res[0].attrib['step_val'] == last_freq:
      continue
    last_freq = res[0].attrib['step_val']
    item = []
    item.append(float(last_freq))
    # read data
    e = xml.xpath(xpath_base + "/@id")
    x = xml.xpath(xpath_base + "/@x")
    y = xml.xpath(xpath_base + "/@y")
    z = xml.xpath(xpath_base + "/@z")
    # assume len(res) = len(x) = len(y), eventually also len(z)
    data = []
    for i in range(len(e)):
      e_val = int(e[i])
      x_val = float(x[i])
      y_val = float(y[i])
      z_val = cond(len(z) == 0, 0.0, float(z[i]))
      data.append(LumpedMechDisplacement(e_val, x_val, y_val, z_val))
    item.append(data)
    result.append(item)
  return result  

# # create an interpolated density file from a readLumpedMechDisplacement() result for positive z displacement.
# below the first frequency we assume full material, above the highest frequency an exception is thrown
# @param data is the result from readLumpedMechDisplacement
# @param f the frequency to interpolate
# @param density_file the filename where the density file is written to
# @param mode 'min' or 'max' or 'auto'
# @return lower, upper, alpha, beta, min, max
def interpolateLumpedMechDisplacementAsDensity(data, f, density_file, mode="auto"):
  # find upper and lower frequency
  lower = 0.0
  upper = -1.0
  
  # the data array of interest
  lower_data = []  # stays empty if f is below first eigenfrequency, assume 1.0 values then
  upper_data = []
  
  for i in range(len(data)):
    test = data[i][0]
    if test < f:
      lower = test
      continue 
    else:
      upper = test
      upper_data = data[i][1]
      if i > 0:
        lower_data = data[i - 1][1]
      break
  if upper == -1.0:
    raise RuntimeError("Cannot find frequency " + str(f) + " within lumpedMechDisplacement data set from f=" + str(data[0][0]) + " to f=" + str(data[len(data) - 1][0]))

  # determine interpolation parameter, normalized to 1
  beta = (f - lower) / (upper - lower)
  alpha = 1.0 - beta
  print "f=" + str(f) + " lower=" + str(lower) + " upper=" + str(upper) + " alpha=" + str(alpha) + " beta=" + str(beta) + "\n"

  # preliminary result, before normalizing to 1.0
  tmp = []
  # we don't know if the positive or negative real parts dominates. the first mode has the same sign which might be neg
  max_v = -1e30   
  min_v = +1e30
  for i in range(len(upper_data)):
    lower_value = 1.0  # there is no real condition operator in python :( 
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
    if mode == "both":
      v = tmp[i]
      max_case = cond(v / max_v > 1e-7, v / max_v, 1e-7)  # max_case in [1e-7:1]
      min_case = cond(v / min_v > 1e-7, v / min_v, 1e-7)  # if min_v > 0: min_case >= 1, for min_v < 0 also [1e-7:+1]
      tmp[i] = cond(v >= 0, max_case, -1.0 * min_case)  # this is not solid plate for min_v > 0!
    if mode == "max" or (mode == "auto" and abs(max_v) >= abs(min_v) and max_v > 0):
      # print "a max=" + str(max_v) + " min=" + str(min_v) + " t=" + str(tmp[i]) + " -> " + str(cond(tmp[i] / max_v > 1e-7, tmp[i] / max_v, 1e-7)) + "\n"
      tmp[i] = cond(tmp[i] / max_v > 1e-7, tmp[i] / max_v, 1e-7)
    else:
      # print "b max=" + str(max_v) + " min=" + str(min_v) + " t=" + str(tmp[i]) + " -> " + str(cond(tmp[i] / min_v > 1e-7, tmp[i] / min_v, 1e-7)) + "\n"      
      tmp[i] = cond(tmp[i] / min_v > 1e-7, tmp[i] / min_v, 1e-7)
      tmp[i] = cond(tmp[i] > 1.0, 1.0, tmp[i])
    if tmp[i] < 1e-8 or tmp[i] > 1.0:
      raise RuntimeError("invalid density " + str(tmp[i]))
      
  # write the stuff    
  out = open(density_file, "w")
  out.write('<?xml version="1.0"?>\n')
  out.write('<cfsErsatzMaterial>\n')
  out.write('  <header>\n')
  out.write('    <design initial="0.5" lower="1e-3" name="density" upper="1" region="all" />\n')
  out.write('    <transferFunction application="mech" design="density" param="1" type="simp"/>\n')
  out.write('  </header>\n')
  out.write('  <set id="f_' + str(f) + '_' + str(alpha) + '*' + str(lower) + '+' + str(beta) + '*' + str(upper) + '_min_' + str(min_v) + '_max_' + str(max_v) + '">\n')
  for i in range(len(tmp)):
    out.write('    <element nr="' + str(upper_data[i].elemNum) + '" type="density" design="' + str(tmp[i]) + '"/>\n')
  
  out.write('  </set>\n')

  out.write(' </cfsErsatzMaterial>\n')
  out.close()
 
  return lower, upper, alpha, beta, min_v, max_v
 

# helper for external use - create a 2D "window"
# strength is the size of the window in fraction. .5 is maximum
# return 2d array
def make2DWindow(divider, strength, lower):
  ret = numpy.zeros((divider, divider))
  
  if strength <= 0.01 or strength > 0.48:
    raise RuntimeError("invalid strength")
  
  # fraction
  f = strength
  # counter fraction
  cf = 1.0 - f
  
  for ii in range(divider):
    for jj in range(divider):
      i = float(ii)
      j = float(jj)
  
      ip = i / (divider - 1)
      jp = j / (divider - 1)
      
      # print str(ip) + " - " + str(jp) + " f " + str(f) + " cf " + str(cf)
      if ip > f and ip < cf and jp > f and jp < cf:
        ret[i][j] = lower
      else:
        ret[i][j] = 1.0

  return ret    
      
# extrudes an 2D array to the third dimenstion
def extrude(data_2d):
  edge = data_2d.shape[0] 
  if edge <> data_2d.shape[1]:
    raise RuntimeError("require quadratic input")
  
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

  lc = 0
  for line in open(infile, 'r'):
    lc += 1
    ls = line.strip().split()
    if len(ls) < 1:
      continue

    ls0strip = ls[0].strip()

    # we assume the first line is <?xml version="1.0"?>
    if lc == 1:
      if ls0strip != '<?xml':
        return False
      continue

    if ls0strip == '<cfsErsatzMaterial' or ls0strip == '<cfsErsatzMaterial>':
      ersatzfound = True

    if ls0strip == '<header' or ls0strip == '<header>':
      headerfound = True

    # if we find a set, all requirements are met, so we break
    if ls0strip == '<set' or ls0strip == '<set>':
      setfound = True
      break

    # if all is found, then we break
    if (ersatzfound and headerfound and setfound):
      break

  return (ersatzfound and headerfound and setfound)

# do an ascii print of the density data
def ascii_print(data, threshold):
  x, y, z = getDim(data)
  assert()
  
# # convert textual tensordata to XML for rectangle homogenization
# @param tensordata: file with lines  a b 11 12 22 33 where a and b are percent
# @param xmlfile writes the xml file if given
# return a string with the text 
def convert_hom_rect_tensor_to_xml(tensordata, xmlfile=None):
  a = numpy.genfromtxt(tensordata)
  assert(a.shape[1] == 6)
  
  out = "<homRect>\n"
  for i in range(len(a)):
    data = a[i]
    # <data a="0.0" b="0.0" e11="" e12="" e22="" e33="" />
    out += '<data a="' + str(data[0] / 100.0) + '" b="' + str(data[1] / 100.0)
    out += '" e11="' + str(data[2]) + '" e12="' + str(data[3]) + '" e22="' + str(data[4]) + '" e33="' + str(data[5]) + '"/>\n'
  out += "</homRect>\n"    

  if xmlfile:
    file = open(xmlfile, "w")
    file.write(out)
    file.close()

  return out    

# # rotate if possible
# @param data 2d array 
# @param center list with sufficient size of indices
# @param angle counter-clock-wise in rad
# @return new array with org where no rotation was possible
def rotate(data, center, angle):
  ret = data.copy()
  
  cx, cy = center
  
  nx, ny = data.shape
  
  # the general algorithm for rotation around a given point is
  # - translate to center
  # - rotate around center
  # - translate back
  
  for y in range(ny):
    for x in range(nx):
       x2 = numpy.cos(angle) * (x - cx) - numpy.sin(angle) * (y - cy) + cx
       y2 = numpy.sin(angle) * (x - cx) + numpy.cos(angle) * (y - cy) + cx
       #print "x=" + str(x) + " y=" + str(y) + " v="  + str(data[y, x]) + " ->  x2=" + str(x2) + " y2=" + str(y2)
       h = 1.0/nx
       bx = 0.5*h + x*h
       by = 0.5*h + y*h
       bx2 = 0.5*h + x2*h
       by2 = 0.5*h + y2*h
       print "org=(" + str(bx) + ", " + str(by) + ") v="  + str(data[y, x]) + " ->  (" + str(bx2) + ", " + str(by2) + ")"
        
       if x2 >= 0 and x2 < nx and y2 >= 0 and y2 < ny:
         ret[y, x] = data[y2, x2] 
       
  return ret
  
  
## returns an image out of an numpy array.
# @param data numpy array
# @param scale if not given the original size is given 
# @return visualize with .show()
def get_image(data, scale = None):
  
  x, y = data.shape
    
  ret = numpy.zeros((y, x), dtype="uint8")
  
  # copy data from linear list
  for i in range(y):
    for j in range(x):
      ret[y-i-1][j] = 255 - int(255 * data[j][i])

  img = Image.fromarray(ret)
  if scale:
    img = img.resize((scale, y/x * scale), Image.NEAREST) # higher quality with ANTIALIAS
  return img


## calculate the normalized perimeter
# @param data array
# @param eps for continuation, 0.0 is ok
# @param order TV order, 2 = taxi cap or 4. See Petersson, Beckers and Duysinx, 1999
# @param normalze if not normalized we give the perimeter in m assuming a 1x1m domain  
def perimeter(data, eps = 0.0, order = 2, normalize = True):
  
  assert(order == 2 or order == 4)
 
  n2, n1 = data.shape
  per = 0.0

  h1 = 1.0/n1
  h2 = 1.0/n2

   
  for j in range(0,n2): 
    for i in range(1,n1):
      mine  = data[j,i]    # rho_ij
      other  = data[j,i-1] # rho_i-1,j
      scale = h2/(n1-1) if normalize else h2
      per += scale * (numpy.sqrt((mine-other)**2 + eps**2) - eps)    

  for j in range(1,n2):
    for i in range(0,n1):
      mine  = data[j,i]    # rho_ij
      other  = data[j-1,i] # rho_i,j-1
      scale = h1/(n2-1) if normalize else h1
      per += h1 * (numpy.sqrt((mine-other)**2 + eps**2) - eps)    

  if order == 4:
    # Petersson, Beckers and Dun2sinn1, "Almost Isotropic Perimeters in Topologn2 Optimization: Theoretical and Numerical Aspects", 1999  
    # -> (5)      
    h = numpy.sqrt(h1**2 + h2**2)

    for j in range(0,n2-1):
      for i in range(0,n1-1):
        mine  = data[j,i]      # rho_ij
        other  = data[j+1,i+1] # rho_i+1,j+1
        scale = 1./((n2-1)*(n1-1)) if normalize else (h1*h2/h) 
        per += scale * (numpy.sqrt((mine-other)**2 + eps**2) - eps)         
  
    for j in range(0,n2-1):
      for i in range(1,n1):
        mine  = data[j,i]      # rho_ij
        other  = data[j+1,i-1] # rho_i-1,j+1
        scale = 1./((n2-1)*(n1-1)) if normalize else (h1*h2/h)
        per += scale * (numpy.sqrt((mine-other)**2 + eps**2) - eps)         

  # normalization is done for TV_2, for TV_4 we need to correkt the diagonal elements -> c_4^-1 = 1 + 2 | cos(theta) | below (3)
  if not normalize and order ==4:
    per /= 1 + 2 * numpy.cos(numpy.pi / 4)
  
  return per
  
# add a save ghost cell around an array
# @param reproduce shall the outer boundary be copied 
def add_ghost_cells(data, reproduce = False):
  x, y = data.shape
  result = numpy.zeros((x + 2, y + 2))
  result[1:x+1, 1:y+1] = data

  # copy the boundary
  if reproduce:
    result[0,]    = result[1,]
    result[y+1,]  = result[y,]
    result[:,0]   = result[:,1]
    result[:,x+1] = result[:,x] 

  return result
  
  
    
# a = read_multi_design("fmomulti-40.density.xml", "stiff1", "stiff2", "rotAngle", "rotAngle2")
# a[:,0] *= 0.11
# a[:,1] *= 0.11
# b = a[:,0:3]
# write_multi_design_file("hom_rect_40.density.xml", b, "hom_rect_a", "hom_rect_b", "rotAngle")

#a = read_multi_design("hom_rect_ml_min_compl.density.xml", "stiff1", "stiff2", "rotAngle")
#numpy.savetxt("hom_rect_ml_min_compl.dat", a)
    
#a = read_multi_design("smooth-full-100.xml", "density", "stiff1, "stiff2", "rotAngle")
#b = numpy.zeros(a.shape[0],4);
#b[:,0] = a[:,0]*a[:,1]**.5;
#b[:,1] = a[:,0]*a[:,2]**.5;
#b[:,2] = a[:,3];
#write_multi_design_file("sqrt-smooth-full-100.xml", b, "tensor11", "tensor22", "tensor33",  "rotAngle")
#a[:,1] = a[:,1]**.5;
#a[:,2] = a[:,2]**.5;
#write_multi_design_file("sqrt-smooth-full-100.xml", a, "density", "tensor11", "tensor22", "tensor33", "rotAngle")
