#!/usr/bin/env python
# this ia a tool for shape map optimization 

import matplotlib
# necessary for remote execution, even when only saved: http://stackoverflow.com/questions/2801882/generating-a-png-with-matplotlib-when-display-is-undefined
#matplotlib.use('Agg')
matplotlib.use('tkagg')

import matplotlib.patches
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
from matplotlib.path import Path
from scipy import ndimage
import numpy as np
from optimization_tools import *
import argparse
import copy
from scipy.interpolate import interp1d

try:
  from matviz_vtk import *
except:
  print('could not import matviz_vtk, hope you do not need it: ' + sys.exc_info()[0])

def create_figure(res, minimal, maximal):

  dpi_x = res / 100.0 
  
  fig = matplotlib.pyplot.figure(dpi=100, figsize=(dpi_x, dpi_x))
  ax = fig.add_subplot(111)
  ax.set_xlim(min(0,minimal[0]), max(1,maximal[0]))
  ax.set_ylim(min(0,minimal[1]), max(1,maximal[1]))
  return fig, ax

def dump_shapes(shapes):
  for s in shapes:
     print(s)   

def find_shape(shapes, el):
  for s in shapes:
    if el >= s.el[0] and el <= s.el[-1]:
      return s
  print('none of the ' + str(len(shapes)) + ' shapes has an element with nr ' + str(el)) 
  return None     

def find_shape_by_dof(shapes, dof):
  res = []  
  for s in shapes:
    if s.dof == dof:
      res.append(s)  
  return res     

  
class Shape: 
  def __init__(self, id, dof):
    self.id = id
    self.dof = dof
    self.el = []
    self.val = []
    self.profile = []
    self.valid = [] # only necessary for image import
    if id == 0:
      self.color = 'b'  
    if id == 1:
      self.color = 'g'  
    if id == 2:
      self.color = 'r'  
    if id == 3:
      self.color = 'c'  
    if id == 4:
      self.color = 'm'  
    if id == 5:
      self.color = 'y'  
    if id == 6:
      self.color = 'k'
    assert(id <= 6)  
  
  # print the shape info
  def __str__(self):
    return "id=" + str(self.id) + " dof=" + str(self.dof) + " color=" + self.color + " #el=" + str(len(self.el)) + " el=" + str(self.el)

  # shape info for a given index
  def to_string(self, idx):
    return "shape=" + str(self.id) + " dof=" + str(self.dof) + " color=" + str(self.color) + " idx=" + str(idx) + " val=" + str(self.val[idx]) + " profile=" + str(self.profile[idx]) + " valid=" + str(self.valid[idx]) + ""
   
  #@return x, y for center point a
  def get_center(self, idx):   
    free = idx * float(1./(len(self.el)-1))
    val  = self.val[idx]
    # this is the a middle line
    x = free if self.dof == 1 else val
    y = val  if self.dof == 1 else free
    return x, y
    
  # averages the current profiles (only valid)
  def average_valid_profile(self):
    if len(self.profile) != len(self.valid):
      print(len(self.profile))
      print(len(self.valid))    
    assert(len(self.profile) == len(self.valid))   
    c = 0
    s = 0
    for i in range(len(self.profile)):
      if self.valid[i]:
         c += 1
         s += self.profile[i]
    
    assert(c > 0)
    return s / c

  # return the coordinates for both profile nodes
  #@return x1,y1,x2,y2
  def get_profiles(self, idx):
    x,y = self.get_center(idx)
    if self.dof == 0:
       return x -.5 * self.profile[idx], y, x + .5 * self.profile[idx], y
    else:    
      return x, y -.5 * self.profile[idx], x, y + .5 * self.profile[idx]   
      
  # return the line coordinates for the profile
  #@param left (True) or right (False)
  def get_profile(self, idx_1, idx_2, left):
    x_val = []
    y_val = []

    x,y = self.get_center(idx_1)
    if self.dof == 0:
      x += (-.5 if left else +.5) * self.profile[idx_1]
    else:    
      y += (-.5 if left else +.5) * self.profile[idx_1]  
    
    x_val.append(x)
    y_val.append(y)
    
    x,y = self.get_center(idx_2)
    if self.dof == 0:
      x += (-.5 if left else +.5) * self.profile[idx_2]
    else:    
      y += (-.5 if left else +.5) * self.profile[idx_2]  
    
    x_val.append(x)
    y_val.append(y)

    return x_val, y_val    
    
## symmetrize assues symmetry on x-axis (mirror)
def symmetrize(shapes):
  x = find_shape_by_dof(shapes, 0) 
  y = find_shape_by_dof(shapes, 1)
  assert(len(x) == 2)

  center = .5* (np.asarray(x[0].val) + np.asarray(x[1].val))
  offset = np.mean(center) - 0.5 # positive is shift to right
  shift = int(offset * len(y[0].val) + .5)

  print("center on x-axis is " + str(np.mean(center)) + " (shift " + str(shift) + " nodes) with var " + str(np.var(center))) 

  assert(np.var(center) < 0.001) 
  
  # avg and shift the left and right x-structures (from bottom to top)
  for i in range(len(x[0].val)):
    lv = x[0].val[i] - offset
    rv = x[1].val[i] - offset
    m = .5*(rv-lv)
    x[0].val[i] = .5 - m
    x[1].val[i] = .5 + m 
    
    lp = x[0].profile[i]
    rp = x[1].profile[i]
    mp = .5*(lp+rp)
    x[0].profile[i] = mp
    x[1].profile[i] = mp 

  # shift the horizontal lines and average left/right counterparts    
  for s in y:
    # we need the original data, otherwise rv might be replaced by shifted lv
    nodes = copy.deepcopy(s.val)
    profiles = copy.deepcopy(s.profile)

    for i in range(int(.5 * len(s.val))):  
      lv = nodes[i+shift]      
      rv = nodes[-(i+1-shift)]
      s.val[i] = s.val[-(i+1)] = .5*(lv+rv) 
      
      lp = profiles[i+shift]      
      rp = profiles[-(i+1-shift)]
      s.profile[i] = s.profile[-(i+1)] = .5*(lp+rp) 
     

# @param file withe grad.plot or density.xml
# @param profile use if not in file
def read_file(filename, profile):
  shapes = []

  if filename.endswith('.plot'):  
    all = open(filename).readlines()
    assert(len(all) > 2)
    comment = all[0].split('\t')
    # add first shape
    first = all[1].split('\t')
    shapes.append(Shape(id = int(first[2]), dof=int(first[3])))
    curr = shapes[0]
    for i in range(1,len(all)): # skip first comment
      line = all[i].split()
      id = int(line[2])
      if curr.id != id:
        shapes.append(Shape(id = id, dof=int(line[3])))
        curr = shapes[-1]  
      curr.el.append(int(line[0]))
      curr.val.append(float(line[4]))
      curr.profile.append(profile)
      curr.valid.append(True)    
      
  else:
    # might be too much when not the whole domain is design (e.g. lbm with boundary)
    nx, ny, nz = read_mesh_info(filename, silent=True)
    assert(nx == ny)  
      
    tree = etree.parse(filename, etree.XMLParser(remove_comments=True))
    root = tree.getroot()
    query = '//set[last()]/shapeParamElement' 
    sett = root.xpath(query)

    # check for full mesh and be tolerant against lbm meshes. 
    # we don't known if profiles are written
    nshapes = int(len(sett)/nx + .5)
    shape_elems = nx+1 # default case
    if nshapes > len(sett)/nx: 
      # apparently we have the lbm case
      assert(len(sett)/nshapes == int(len(sett)/nshapes))
      shape_elems = int(len(sett)/nshapes)
      print('assume not the full domain is design: ' + str(nshapes) + ' shapes with ' + str(shape_elems) + ' (' + str(nx+1) + ') elements')   
  
      
    curr = None
    for el in sett:
      nr = int(el.get('nr'))
      v =  float(el.get('design'))
      if nr == 0 or nr % (shape_elems) == 0:
        dof = 0 if el.get('dof') == 'x' else 1
        node = el.get('type') == 'node' or el.get('type') == None
        if node: 
          shapes.append(Shape(id = len(shapes), dof=dof))
          curr = shapes[-1]
        else:
          # we have a profile, hence search the corresponding shape
          if profile != None:
            print('error: profile data given in ' + filename + ' and concurrently via command line')
            sys.exit(-2)  
          curr = find_shape(shapes, nr - shapes[-1].el[-1])   
     
      if node:  
        curr.el.append(nr)
        curr.val.append(v)
      else:  
        curr.profile.append(v)  
        curr.valid.append(True)
   
    # check if there was no profile in the xml
    if len(shapes[-1].val) > 0 and len(shapes[-1].profile) == 0:
       if profile != None: # error message comes below
         for shape in shapes:
           shape.profile = [profile] * len(shape.el)
           shape.valid   = [True] * len(shape.el)       

  # for both input formats
  if len(shapes[-1].val) != len(shapes[-1].profile):
    print("error: no profiles in '" + filename + "' and not given via command line")
    sys.exit(-3)    
  return shapes  

# resamples given shapes. ignores valid
def resample(shapes, resample):
  res = []  
  org_space = np.linspace(0, 1.0, num=len(shapes[0].val), endpoint=True)
  new_space = np.linspace(0, 1.0, num=resample+1, endpoint=True)

  for o in shapes:
     s = Shape(o.id, o.dof)
     s.el = list(range(len(res) * (resample+1), (len(res)+1) * (resample+1)))
     v = interp1d(org_space, o.val, kind='cubic')
     s.val = v(new_space)
     p = interp1d(org_space, o.profile, kind='cubic')
     s.profile = p(new_space)
     s.valid = (resample+1) * [True]
     res.append(s)

  return res   

# searches for data (< 128) within a 1D line
# returns tupes with center and profile in meters (assume unit cube!)
# @param line vector with len == nx of content int 0 ... 255
def find_data_in_line(line):

  data = []

  # the token
  in_data = line[0] < 128
  start = 0 if in_data else -1
  end   = -1 # not found yet
  for i in range(1,len(line)):
    v = line[i]
    if v < 128 and not in_data: # data starts
      assert(start == -1)
      assert(end == -1)
      in_data = True
      start = i
      end   = -1
    elif v >= 128 and in_data: # data ends
      assert(start != -1)
      assert(end == -1)
      in_data = False
      end = i
      center = 1.0/len(line) * (start+end)/2.0
      profile = 1.0/len(line) * float(end-start)
      data.append((center, profile))   
      start = -1
      end = -1           
  # check if we ended with data
  if end != -1:
    center = 1.0/len(line) * (start+end)/2.0
    profile = 1.0/len(line) * float(end-start)
    data.append((center, profile))

  return data        
  
# fixes the shape data which is not valid by interpolating between the last valid
def repair_shapes(shapes):
  cnt = 0  
  for shape in shapes:
    start = -1
    end   = -1
    assert(len(shape.valid) == len(shape.val))
    # find invalid ranges
    in_valid = shape.val[0]
    assert(in_valid) # we shall start good
    for i in range(len(shape.valid)):
      if not shape.valid[i] and in_valid: # invalid starts
        assert(start == -1)        
        assert(end == -1)
        in_valid = False
        start = i
        end = -1
      elif shape.valid[i] and not in_valid: #invalid ends
        assert(start != -1)
        assert(end == -1)
        in_valid = True
        end = i-1
        #print 'last valid region: ' + shape.to_string(start-1)
        for j in range(end+1-start):
          #print start
          #print end
          #print end+2-start  
          #print 1.0/(end+2-start)
          bal = (j+1)*1.0/(end+2-start)
          val = shape.val[start-1] + bal * (shape.val[end+1] - shape.val[start-1])  
          #print 'ivalid: bal = ' + str(bal) + " val=" + str(val) + " -> " + shape.to_string(start+j)
          shape.val[start+j] = val
          shape.profile[start+j] = shape.average_valid_profile()
          cnt += 1
        #print 'invalid region: start ' + str(start) + " -> " + shape.to_string(start)
        #print 'invalid region: end ' + str(end) + " -> " + shape.to_string(end)
        #print 'first valid region: ' + shape.to_string(end+1)
        # reset
        start = -1
        end = -1      
  print("repaired " + str(cnt) + " times")        

# fills shape content for import_form_image based on the previous shapes data
def smart_append_data(n, data, shapes):
  assert(len(data) >= len(shapes))
  # easy case where there is no crossing from y to x or vice versa
  if len(data) == len(shapes):
    for i in range(len(data)):
       shape = shapes[i] 
       a, p = data[i] 
       # el is missing!
       shape.val.append(a)
       # validate profile, it might be too large
       if len(shape.profile) > 0 and p >= 1.3 * shape.average_valid_profile():
         shape.profile.append(shape.average_valid_profile())
         shape.valid.append(False)
       else:
         shape.profile.append(p)
         shape.valid.append(True)
  # we copy the last if we have not enough data
  if len(data) < len(shapes):
    for  s in shapes:
      s.val.append(s.val[-1])  
      s.profile.append(s.average_valid_profile())
      s.valid.append(False)
  # we search the best data if we have too much and throw it away if best is not good enough
  if len(data) > len(shapes):
    for shape in shapes:
      ref_a = shape.val[-1]  
      best = None  
      for dat in data:
        if best == None or abs(dat[0] - ref_a) < abs(best[0] - ref_a):
           best = dat  
      if abs(best[0] - ref_a) < 1.5/n:
        # use best found  
        # print 'best ' + str(best) + ' for ref ' + str(ref_a) + ' for line data ' + str(data)
        shape.val.append(best[0])  
        shape.profile.append(best[1])
        shape.valid.append(True)
      else:    
        #  print 'best ' + str(best) + ' is not good enough for ref ' + str(ref_a) + ' with err ' + str(abs(best[0] - ref_a)) + ' and criteria ' + str(1.5/n)
        # copy old stuff
        shape.val.append(shape.val[-1])  
        shape.profile.append(shape.average_valid_profile())
        shape.valid.append(False)

  
def import_from_image(filename, resample, repair):
  img = Image.open(args.input).transpose(Image.FLIP_TOP_BOTTOM)
  img = img.convert("L")
  dat = numpy.asarray(img, dtype="uint8")
  assert(len(dat.shape) == 2)
  nx, ny = dat.shape
  assert(nx == ny) # check what is all needed!
  if resample == None:
    resample = nx 
  shapes = []  
  # we assume that the first line determines the number of shapes
  # start with columns
  num = len(find_data_in_line(dat[:,0]))
  print("we assume " + str(num) + " horizontal shapes")
  for i in range(num):
    shapes.append(Shape(id = len(shapes), dof=1))
  for i in numpy.linspace(0,nx-1, resample+1): #resample +1 because for 10 elements we have 11 variables
     smart_append_data(resample, find_data_in_line(dat[:,int(i)]), shapes[-num:]) # the last num shapes just appended

  num = len(find_data_in_line(dat[0,:]))
  print("we assume " + str(num) + " vertical shapes")
  for i in range(num):
    shapes.append(Shape(id = len(shapes), dof=0))
  for i in numpy.linspace(0,nx-1, resample+1):
     smart_append_data(resample, find_data_in_line(dat[int(i),:]), shapes[-num:]) # the last num shapes just appended
     
  cnt = 0   

  # add the element numbers   
  for shape in shapes:
    shape.el = list(range(cnt, cnt + len(shape.val)))
    cnt += len(shape.val)
        
  if repair:       
    repair_shapes(shapes)       
  return shapes
    
# creates a matplotlib figure     
def plot_data(res, shapes, unit):
  # find extreme bounds to also visualize negative node positions
  minimal = [0.0]*2
  maximal = [1.0]*2
  if not unit:
    minimal = [1e9]*2
    maximal = [-1e9]*2
    for shape in shapes:
      minimal[shape.dof] = min(minimal[shape.dof], min(shape.val))
      maximal[shape.dof] = max(maximal[shape.dof], max(shape.val))
  
  fig, sub = create_figure(res, minimal, maximal)
  
  for shape in shapes:
    n = len(shape.el)

    # circles
    for i in range(n):
      x,y = shape.get_center(i)
      c=plt.Circle((x,y), 0.005, color = shape.color)
      fig.gca().add_artist(c)
      
    for i in range(0,n-1):
      x, y = shape.get_profile(i, i+1, True) # left          
      l = plt.Line2D(x,y, marker='.', color=shape.color)        
      sub.add_line(l)
      x, y = shape.get_profile(i, i+1, False) # right          
      l = plt.Line2D(x,y, marker='.', color=shape.color)                                    
      sub.add_line(l)
      
  return fig, sub

# create vtk polydata tesselation
def create_vtk(shapes):
  # create vtk cells and points
  points = vtk.vtkPoints()
  cells = vtk.vtkCellArray()
  
  for shape in shapes:
    # we initialize the loop with the -1 values
    cx,cy = shape.get_center(0)
    last_center = points.InsertNextPoint(cx, cy, 0.0)
  
    # last (-1) 'left' and 'right' profile nodes
    xl, yl, xr, yr = shape.get_profiles(0)
    last_left = points.InsertNextPoint(xl, yl, 0.0)
    last_right = points.InsertNextPoint(xr, yr, 0.0)
  
  
    for i in range(1,len(shape.el)):
      # this center node
      cx,cy = shape.get_center(i)
      this_center = points.InsertNextPoint(cx, cy, 0.0)
  
      # this 'left' and 'right' profile nodes
      xl, yl, xr, yr = shape.get_profiles(i)
      this_left = points.InsertNextPoint(xl, yl, 0.0)
      this_right = points.InsertNextPoint(xr, yr, 0.0)
    
      # two quadrialterals:
      # last_left--------this_left
      #    |                 |
      # last_center------this_center
      #    |                 |
      # last_right-------this_right
      # we divide the quadrilaterals by triangles
      tri = vtk.vtkTriangle()
      tri.GetPointIds().SetId(0, last_left)
      tri.GetPointIds().SetId(1, last_center)
      tri.GetPointIds().SetId(2, this_left)
      cells.InsertNextCell(tri)
  
      tri = vtk.vtkTriangle()
      tri.GetPointIds().SetId(0, last_center)
      tri.GetPointIds().SetId(1, this_center)
      tri.GetPointIds().SetId(2, this_left)
      cells.InsertNextCell(tri)
  
      tri = vtk.vtkTriangle()
      tri.GetPointIds().SetId(0, last_center)
      tri.GetPointIds().SetId(1, last_right)
      tri.GetPointIds().SetId(2, this_center)
      cells.InsertNextCell(tri)
  
      tri = vtk.vtkTriangle()
      tri.GetPointIds().SetId(0, last_right)
      tri.GetPointIds().SetId(1, this_right)
      tri.GetPointIds().SetId(2, this_center)
      cells.InsertNextCell(tri)
  
      last_center = this_center
      last_left = this_left
      last_right = this_right
  
    polydata = vtk.vtkPolyData()
    polydata.SetPoints(points)
    polydata.SetPolys(cells)

  return polydata

def export(shapes, filename, suppress_profile):
  out = open(filename, "w")
  out.write('<?xml version="1.0"?>\n')
  out.write('<cfsErsatzMaterial>\n')
  out.write('  <header>\n')
  nx = len(shapes[0].val) - 1
  out.write('    <mesh x="' + str(nx) + '" y="' + str(nx) + '" z="1"/>\n')  
  out.write('  </header>\n')
  out.write('  <set id="shape_map.py">\n')
  # <shapeParamElement nr="0" type="node" dof="x" design="0.3"/>
  for shape in shapes:
    for i in range(len(shape.val)):  
      out.write('    <shapeParamElement nr="' + str(shape.el[i]) + '" type="node" dof="' + ('x' if shape.dof == 0 else 'y') + '" design="' + str(shape.val[i]) + '"/>\n')
  if not suppress_profile:
    base = shapes[-1].el[-1]+1
    for shape in shapes:
      for i in range(len(shape.val)):  
        out.write('    <shapeParamElement nr="' + str(base + shape.el[i]) + '" type="profile" dof="' + ('x' if shape.dof == 0 else 'y') + '" design="' + str(shape.profile[i]) + '"/>\n')
        
           
  out.write('  </set>\n')
  out.write(' </cfsErsatzMaterial>\n')
  
# __name__ is 'shape_map' if imported or '__main__' if run as commandline tool
if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument("input", help="a .density.xml or .grad.plot file or a image")
  parser.add_argument("--profile", help="give the profile if it is not in input", type=float)
  parser.add_argument('--resample', help="resample to this resolution", type=int)
  parser.add_argument('--repair', help="interpolate unsure data (when parsing from image)", action='store_true')
  parser.add_argument('--symmetrize', help="mirror on x-axis", action='store_true')
  parser.add_argument('--export', help="write a density.xml file with shapeParam variables")
  parser.add_argument('--suppress_profile', help="do not export profile", action='store_true')
  parser.add_argument('--unit_bound', help="enforce the visualization on a unit square", action='store_true')
  parser.add_argument('--save', help="save the image to the given name with the given format. Might be png, pdf, eps or even vtp!")
  parser.add_argument('--noshow', help="don't show the image", action='store_true')
  args = parser.parse_args()
  shapes = []
  if args.input.endswith('.xml') or args.input.endswith('.plot'):
    shapes = read_file(args.input, args.profile)
    if args.resample:
      shapes = resample(shapes, args.resample)  
  else:
    shapes = import_from_image(args.input, args.resample, args.repair)
  
  print('average profile is ' + str(1.0/len(shapes) * sum([ s.average_valid_profile() for s in shapes])))
    
  if args.symmetrize:
    symmetrize(shapes)  
  
  if args.export:  
    export(shapes, args.export, args.suppress_profile)

  # vtp generation exclusivly triggerd by saving an vtp file
  if args.save and args.save.endswith('.vtp'):
    polydata = create_vtk(shapes)
    show_write_vtk(polydata,800,args.save)
    if not args.noshow:
      show_vtk(polydata,800,show_edges=True) # show edges
  else:
    fig, sub = plot_data(800, shapes, args.unit_bound)
    if args.save:
      fig.savefig(args.save)  
    if not args.noshow:
      fig.show()
      input("Press Enter to terminate.")
else:
  f = 'shape_map_mech_39.grad.plot'
  print(f)
  #shapes = read_file(f, profile=0.1)
  #dump_shapes(shapes)
  #fig, sub = plot_data(800, shapes)
  #fig.show()
