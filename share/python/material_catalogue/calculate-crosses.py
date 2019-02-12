#!/usr/bin/python
import platform
from PIL import Image, ImageDraw
import lxml
import sys, bz2, os
from scipy import ndimage
import numpy as np
import math
from optimization_tools import *
from cfs_utils import *
import os.path
import matplotlib.pyplot as plt
from mesh_tool import *
import argparse

# lena = misc.lena()
# blurred_lena = ndimage.gaussian_filter(lena, sigma=3)
# very_blurred = ndimage.gaussian_filter(lena, sigma=5)


def save_image_as_densfile(im, outfile):
  imlist = list(im.getdata())
  # out = open(outfile, "w")
  out = bz2.BZ2File(outfile, 'wb')
  out.write('<?xml version="1.0"?>\n <cfsErsatzMaterial infoWriteCounter="0" \
  infoRejectCounter="0">\n <header>\n <mesh x="' + str(minres) + '" y="' + str(minres) + '" \
  z="1"/>\n <design adapt_lower="true" bimaterial="janrescaledvoid" \
  constant="false" fixed="false" initial="0.5" lower="0" name="density" \
  region="mech" scale="false" upper="1"/>\n <transferFunction \
  application="mech" design="density" param="' + str(p) + '" type="simp"/>\n \
  </header>\n<set id="frompng">\n')

  el = 1
  for i in imlist:
    val = (255. - float(i)) / 255.0
    out.write('<element nr="' + str(el) + '" type="density" design="'\
      + str(val ** (1. / p)) + '" physical="' + str(val) + '"/>\n')
    el += 1


  out.write('</set>\n </cfsErsatzMaterial>')
  out.close()
def insert_modified_frame(array, minres, x, y, steps, void, number, steps_p = 0, modify=False, triangle=False):
  # creates a density file with one or multiple 2D frame structures, optional: frame with round corners (modify option)
  # array: density array; minres: resolution of comp. domain; x,y, steps: x/steps is the thickness of the bar in x-direction;
  # void: density of void material, eg. 1e-9
  
  # triangles or quads used for triangulation
  if triangle:
    array = np.ones((2 * minres, minres))
  else:
    array = np.ones((minres, minres))
  
  # creates density file
  eps = 1e-8
  if steps_p:
    offx = int((minres / (2.*number)) * ((float(x) * steps_p / steps)**3./steps_p) + 0.5 + eps)
    offy = int((minres / (2.*number)) * ((float(y) * steps_p / steps)**3./steps_p) + 0.5 + eps)
  else:
    offx = int((minres / (2.*number)) * (float(x) / steps) + 0.5 + eps)
    offy = int((minres / (2.*number)) * (float(y) / steps) + 0.5 + eps)
  for nx in range(0, number):
    for ny in range(0, number):
      rx = int(float(nx) / float(number) * minres + 0.5 + eps)
      ry = int(float(ny) / float(number) * minres + 0.5 + eps)
      endx = int(float(nx + 1) / float(number) * minres + 0.5 + eps) - offx
      endy = int(float(ny + 1) / float(number) * minres + 0.5 + eps) - offy
      for i in range(offx + rx, endx):
        for j in range(offy + ry, endy):
          if not triangle:
            array[i][j] = void
          else:
            array[2 * i][j] = void
            array[2 * i + 1][j] = void
      # TODO: does not work for triangles yet
      if modify:    
        # modify frame for stress minimization with smooth interior corners
        for i in range(offx, minres - offx):
          for j in range(offy, minres - offy):
            if math.ceil((minres - 2.*offx) / 3.) <= math.ceil((minres - 2.*offy) / 3.):
              r = math.ceil((minres - 2.*offx) / 3.)
            else:
              r = math.ceil((minres - 2.*offy) / 3.) 
            m = [offx + r, offy + r]
            if i - offx < r and j - offy < r and (i - m[0]) * (i - m[0]) + (j - m[1]) * (j - m[1]) >= r * r:
              array[i][j] = 1.
            m = [minres - offx - r - 1, offy + r]
            if i >= minres - offx - r and j - offy < r and (i - m[0]) * (i - m[0]) + (j - m[1]) * (j - m[1]) >= r * r:
              array[i][j] = 1.
            m = [minres - offx - r - 1, minres - offy - r - 1]
            if i >= minres - offx - r and j >= minres - offy - r and (i - m[0]) * (i - m[0]) + (j - m[1]) * (j - m[1]) >= r * r:
              array[i][j] = 1.
            m = [offx + r, minres - offy - r - 1]
            if i - offx < r and j >= minres - offy - r and (i - m[0]) * (i - m[0]) + (j - m[1]) * (j - m[1]) >= r * r:
              array[i][j] = 1.  
  return array  

def insert_rotated_bar(tx, ty, rot, array, minres, midx, midy):
  # rotation angle should be between [0,pi] 
  # note that thickness of the bars not scaled by rotation angle
  #tx: thicknesses; ty: thicknesses; rot: rotation angle; minres: resolution; midx,midy: midpoint 
  eps = 1e-8
  for angle in [rot, rot + math.pi / 2.]:
    if angle != rot:
      tx = ty
    for t in range(0, tx): 
      d = (math.tan(angle))
      if d <= 1. and d > eps:
        dtan = round(abs(1. / d))
        i = midx - 1 - t + tx / 2 + 1
        j = int(midy - 1)
      elif d > 1.:
        dtan = round(abs(d))
        j = midx - 1 - t + tx / 2 + 1
        i = int(midy - 1)
      elif d < -1.:
        dtan = round(abs(d))
        j = midx - 1 - t + tx / 2 + 1
        i = int(midy - 1)
      elif d < -eps and d >= -1.:
        dtan = round(abs(1. / d))
        i = midx - 1 - t + tx / 2 + 1
        j = int(midy - 1)
      else:
        dtan = round(abs(d)) 
        i = midx - 1 - t + tx / 2 + 1
        j = int(midy - 1)
      count = 1
      # first part of the line
      while i >= 0 and i < minres and j >= 0 and j < minres:
        array[i][j] = 1.
        if d <= 1. and d > eps:
          if count == dtan:
            count = 1
            i = i - 1
          else:
            count += 1
          j = j + 1
        elif d > 1.:   
          if count == dtan:
            count = 1
            j = j + 1
          else:
            count += 1
          i = i - 1
        elif d < -1.:
          if count == dtan:
            count = 1
            j = j - 1
          else:
            count += 1
          i = i - 1
        elif d < -eps and d >= -1.:
          if count == dtan:
            count = 1
            i = i - 1
          else:
            count += 1
          j = j - 1
        else:
          if count == dtan:
            count = 1
            i = i - 1
          else:
            count += 1
          j = j + 1
          
      if d <= 1. and d > eps:
        i = midx - 1 - t + tx / 2 + 1
        j = int(midy - 1)
      elif d > 1.:
        dtan = round(abs(d))
        j = midx - 1 - t + tx / 2 + 1
        i = int(midy - 1)
      elif d < -1.:
        dtan = round(abs(d))
        j = midx - 1 - t + tx / 2 + 1
        i = int(midy - 1)
      elif d < -eps and d >= -1.:
        dtan = round(abs(1. / d))
        i = midx - 1 - t + tx / 2 + 1
        j = int(midy - 1)
      else:
        dtan = round(abs(d)) 
        i = midx - 1 - t + tx / 2 + 1
        j = int(midy - 1) 
      count = 1
      # second part of the line
      while i >= 0 and i < minres and j >= 0 and j < minres:
        array[i][j] = 1.
        if d <= 1. and d > eps:
          if count == dtan:
            count = 1
            i = i + 1
          else:
            count += 1
          j = j - 1
        elif d > 1.:   
          if count == dtan:
            count = 1
            j = j - 1
          else:
            count += 1
          i = i + 1
        elif d < -1.:
          if count == dtan:
            count = 1
            j = j + 1
          else:
            count += 1
          i = i + 1
        elif d < -eps and d >= -1.:
          if count == dtan:
            count = 1
            i = i + 1
          else:
            count += 1
          j = j + 1
        else:
          if count == dtan:
            count = 1
            i = i + 1
          else:
            count += 1
          j = j - 1
  return array
 
parser = argparse.ArgumentParser()
parser.add_argument("stp", help="number of grid points in one direction", type=int)
parser.add_argument("dimension", help="Dimension of the problem", type=int, default=2)
parser.add_argument("res", help="grid size", type=int)
parser.add_argument("folder", help="specify the output folder")
parser.add_argument("--point", help="point list", type=np.array, default=np.zeros((1, 1)))
parser.add_argument("--msfem", help="name of msfem grid file on the fine scale")
parser.add_argument("--hom", help="name of grid file for homogenization cell problem")
parser.add_argument("--sparse_msfem", help="sparse msfem option true or false")
parser.add_argument("--shape", help="choose between frame or cross", choices=['frame', 'cross', 'frame_modified', 'frame_w_triangles', 'auxetic','ortho_bc'])
parser.add_argument("--triangle_msfem", help="true or false")
parser.add_argument("--filter", help="filtered densities on or off")
parser.add_argument("--void_material", help="set value for void material", type=float, default=1e-9)
parser.add_argument("--epsilon", help="number of frames/crosses in the cell problem", type=int, default=1)
parser.add_argument("--design", help="select single thicknesses s1,s2,s3 for debugging,e.g. 0.1,0.3,0.")
parser.add_argument("--oversampling", help="name of the mesh with size minres/epsilon including only one base cell")
parser.add_argument("--penalization", help="creates a penalized material catalogue in the interval [0, 1/steps_p], step_p has to be given",type=int)
parser.add_argument("--gmsh",help="folder with stp files of geometry")






args = parser.parse_args()
pl = parser.parse_args()
dim = args.dimension
minres = args.res
steps = args.stp
folder = args.folder
void = args.void_material

if args.shape == "auxetic" and not args.gmsh:
  print("Error: --gmsh option is necessary for auxetic structures!")

if dim == 2:
  #setup calculation directory
  if args.msfem:
    if args.triangle_msfem:
      outfile = str(folder) + '_/jobs' 
      print("computing with " + str(steps) + " steps!")
      if not os.path.exists(str(folder)):
        os.mkdir(str(folder))
      os.system('cp mat.xml ' + str(folder))
    else:
      outfile = str(folder) + '/jobs' 
      print("computing with " + str(steps) + " steps!")
      if not os.path.exists(str(folder) + "/"):
        os.mkdir(str(folder))
      os.system('cp mat.xml ' + str(folder) + '/')
  else:
    outfile = str(folder) + '/jobs'
    print("computing with " + str(steps) + " steps!")
    if not os.path.exists(str(folder) + ''):
      os.mkdir(str(folder) + "/")
    os.system('cp mat.xml ' + str(folder) + '/')
if dim == 3:
  outfile = str(folder) + '/jobs'
  print("computing with " + str(steps) + " steps!")
  if not os.path.exists(str(folder) + ''):
    os.mkdir(str(folder) + '')
  os.system('cp mat.xml ' + str(folder) + '/')
  
print(outfile)
jobfile = open(outfile, "w")
if dim == 2:
  if args.msfem:
    array = void * np.ones((minres, minres))
    midx = minres / 2
    midy = minres / 2
    rot = 0.
    x = 0
    # loops over all different frame/cross thicknesses
    while x < steps + 1:
      y = 0
      while y < steps + 1:
        if args.penalization:
          x_tmp = x
          y_tmp = y
          x = float(1./ x) if x > 0 else 0
          y = float(1./ y) if y > 0 else 0
          if args.design:
            tmp = args.design.split(',')
            x = float(steps * float(tmp[0]))
            y = float(steps * float(tmp[1]))  
        elif args.design:
          tmp = args.design.split(',')
          x = int(steps * float(tmp[0]))
          y = int(steps * float(tmp[1]))
        # create the design array of the chosen shape for MSFEM 2D
        if args.shape == 'frame_modified':
          array = insert_modified_frame(array, minres, y, x, steps, void, args.epsilon,0, modify = True)
        elif args.shape == 'cross':
          array = insert_rotated_bar(x, y, 0., array, minres, midx, midy)
        elif args.shape == 'frame':
          array = insert_modified_frame(array, minres, y, x, steps, void, args.epsilon,0, modify=False)
        elif args.shape == 'frame_w_triangles':
          array = insert_modified_frame(array, minres, y, x, steps, void, args.epsilon,0, modify=False, triangle = True)
        else:
          print('Warning: base cell type is undefined, set --shape [frame or frame_modified or cross]')
        
        densfilename = str(x) + "-" + str(y) + "_msfem.dens.xml"
        
        # MSFEM oversampling in 2D
        if args.oversampling:
          # Warning: density file created here is only used for calculation of global stiffness matrix Kglob
          overarray = void * np.ones((int(minres/args.epsilon+0.5+1e-6), int(minres/args.epsilon+0.5+1e-6)))
          if args.shape == 'frame':
            overarray = insert_modified_frame(overarray, int(minres/args.epsilon+0.5+1e-6), y, x, steps, void, 1, 0, False)
          else:
            print('Warning: other shapes not implemented yet for oversampling.')
          overdensfilename = str(x) + "-" + str(y) + "_msfem_oversample.dens.xml"
          if args.filter == 'on':
            #filtering of the data due to the theory of homogenization
            overarray_filter = ndimage.uniform_filter(overarray, size=6)
            plt.gray()
            plt.imshow(overarray_filter)
            plt.show()
          else:
            overarray_filter = overarray
          # write density file
          if not args.sparse_msfem:
            write_density_file(str(folder) + "/" + overdensfilename, overarray_filter, "set")
          else:
            print('Warning: Sparse meshes not implemented for oversampling yet.')
          overarray = void * np.ones((int(minres/args.epsilon+0.5+1e-6), int(minres/args.epsilon+0.5+1e-6)))   
        
        # filter for smoothing the frame/cross, necessary for homogenization theory  
        if args.filter == 'on':
          # filtering of the data
          array_filter = ndimage.uniform_filter(array, size=6)
          plt.gray()
          plt.imshow(array_filter)
          plt.show()
        else:
          # no filter is applied
          array_filter = array
        if not args.sparse_msfem:
          write_density_file(str(folder) + "/" + densfilename, array_filter, "set")
        array = void * np.ones((minres, minres))
        
        # WARNING: not tested any more
        # creating cfs .xml files for msfem cell problems using triangle 2D elements
        if args.triangle_msfem:
          # create xml file for cfs
          doc = lxml.etree.parse("triangle_msfem.xml", lxml.etree.XMLParser(remove_comments=True))
          #err  = doc.xpath("//orthotropy/@err")[0] 
          val = []
          for i in range(6):
            val.append("0.")
          func = []
          for i in range(3):
            func.append("0.")
          func[0] = "1.-x-y"
          func[1] = "x"
          func[2] = "y"
          index1 = 0
          index2 = 1
          index = 0
          for i in range(6):
            # set MSFEM boundary conditions for different cell problems
            val[index1] = func[index]
            val[index2] = func[index]
            if i % 2 == 0:
              # degree of freedom x
              replace(doc, '//cfs:dirichletInhom[@name="left"][@dof="x"]/@value', val[0])
              replace(doc, '//cfs:dirichletInhom[@name="bottom"][@dof="x"]/@value', val[1])
              replace(doc, '//cfs:dirichletInhom[@name="right"][@dof="x"]/@value', val[2])
              replace(doc, '//cfs:dirichletInhom[@name="left"][@dof="y"]/@value', str(0.))
              replace(doc, '//cfs:dirichletInhom[@name="bottom"][@dof="y"]/@value', str(0.))
              replace(doc, '//cfs:dirichletInhom[@name="right"][@dof="y"]/@value', str(0.))
            else:
              # degree of fredom y
              replace(doc, '//cfs:dirichletInhom[@name="left"][@dof="y"]/@value', val[0])
              replace(doc, '//cfs:dirichletInhom[@name="bottom"][@dof="y"]/@value', val[1])
              replace(doc, '//cfs:dirichletInhom[@name="right"][@dof="y"]/@value', val[2])
              replace(doc, '//cfs:dirichletInhom[@name="left"][@dof="x"]/@value', str(0.))
              replace(doc, '//cfs:dirichletInhom[@name="bottom"][@dof="x"]/@value', str(0.))
              replace(doc, '//cfs:dirichletInhom[@name="right"][@dof="x"]/@value', str(0.))
            val[index1] = "0."
            val[index2] = "0."
            if i % 2 != 0:
              if index1 == 2:
                index1 = -1
              if index2 == 2:
                index2 = -1
              index1 += 1
              index2 += 1         
            if i % 2 == 0:         
              doc.write(str(folder) + '/' + str(x) + "-" + str(y) + '_msfem' + str(index) + '_x.xml',pretty_print=True)
            else:
              doc.write(str(folder) + '/' + str(x) + "-" + str(y) + '_msfem' + str(index) + '_y.xml',pretty_print=True)
            if i % 2 == 0:
              # add new job to jobfile
              jobfile.write('cfs.rel -m ' + str(args.msfem) + ' -x ' + densfilename + ' ' + str(x) + "-" + str(y) + '_msfem' + str(index) + '_x \n')
              os.system('cp triangle_msfem.xml ' + str(folder) + '/' + str(x) + "-" + str(y) + '_msfem' + str(index) + '_x.xml')
            else:
              # add new job to jobfile
              jobfile.write('cfs.rel -m ' + str(args.msfem) + ' -x ' + densfilename + ' ' + str(x) + "-" + str(y) + '_msfem' + str(index) + '_y \n')
              os.system('cp triangle_msfem.xml ' + str(folder) + '/' + str(x) + "-" + str(y) + '_msfem' + str(index) + '_y.xml')
              index += 1
          print(str(x) + ' ' + str(y) + ' is done')
        else:
          # MSFEM for 2D quadrileteral elements 
          # create xml file for cfs solving MSFEM cell problems
          doc = lxml.etree.parse("compliance_plain.xml", lxml.etree.XMLParser(remove_comments=True))
          val = []
          for i in range(8):
            val.append("0.")
          func = []
          for i in range(4):
            func.append("0.")
          func[0] = "0.25*(1.-x)*(1.-y)"
          func[1] = "0.25*(1.+x)*(1.-y)"
          func[2] = "0.25*(1.+x)*(1.+y)"
          func[3] = "0.25*(1.-x)*(1.+y)"
          index1 = 0
          index2 = 1
          index = 0
          for i in range(8):
            # set MSFEM boundary conditions for different cell problems
            val[index1] = func[index]
            val[index2] = func[index]
            if i % 2 == 0:
              # degree of freedom x
              replace(doc, '//cfs:dirichletInhom[@name="left"][@dof="x"]/@value', val[0])
              replace(doc, '//cfs:dirichletInhom[@name="bottom"][@dof="x"]/@value', val[1])
              replace(doc, '//cfs:dirichletInhom[@name="right"][@dof="x"]/@value', val[2])
              replace(doc, '//cfs:dirichletInhom[@name="top"][@dof="x"]/@value', val[3])
              replace(doc, '//cfs:dirichletInhom[@name="left"][@dof="y"]/@value', str(0.))
              replace(doc, '//cfs:dirichletInhom[@name="bottom"][@dof="y"]/@value', str(0.))
              replace(doc, '//cfs:dirichletInhom[@name="right"][@dof="y"]/@value', str(0.))
              replace(doc, '//cfs:dirichletInhom[@name="top"][@dof="y"]/@value', str(0.))
            else:
              # degree of fredom y
              replace(doc, '//cfs:dirichletInhom[@name="left"][@dof="y"]/@value', val[0])
              replace(doc, '//cfs:dirichletInhom[@name="bottom"][@dof="y"]/@value', val[1])
              replace(doc, '//cfs:dirichletInhom[@name="right"][@dof="y"]/@value', val[2])
              replace(doc, '//cfs:dirichletInhom[@name="top"][@dof="y"]/@value', val[3])
              replace(doc, '//cfs:dirichletInhom[@name="left"][@dof="x"]/@value', str(0.))
              replace(doc, '//cfs:dirichletInhom[@name="bottom"][@dof="x"]/@value', str(0.))
              replace(doc, '//cfs:dirichletInhom[@name="right"][@dof="x"]/@value', str(0.))
              replace(doc, '//cfs:dirichletInhom[@name="top"][@dof="x"]/@value', str(0.))
            val[index1] = "0."
            val[index2] = "0."
            if i % 2 != 0:
              if index1 == 3:
                index1 = -1
              if index2 == 3:
                index2 = -1
              index1 += 1
              index2 += 1
            if i % 2 == 0:         
              doc.write(str(folder) + '/' + str(x) + "-" + str(y) + '_msfem' + str(index) + '_x.xml',pretty_print=True)
            else:
              doc.write(str(folder) + '/' + str(x) + "-" + str(y) + '_msfem' + str(index) + '_y.xml',pretty_print=True)
            # sparse is not working currently
            if args.sparse_msfem:
              mesh = str(x) + "-" + str(y) + '.mesh'
            else:
              mesh = args.msfem
            if i % 2 == 0:
              # add new job to jobfile
              jobfile.write('cfs.rel -m ' + mesh + ' -x ' + densfilename + ' ' + str(x) + "-" + str(y) + '_msfem' + str(index) + '_x \n')
              #os.system('cp compliance_plain.xml ' + str(folder) + '/' + str(x) + "-" + str(y) + '_msfem' + str(index) + '_x.xml')
            else:
              # add new job to jobfile
              jobfile.write('cfs.rel -m ' + mesh + ' -x ' + densfilename + ' ' + str(x) + "-" + str(y) + '_msfem' + str(index) + '_y \n')
              #os.system('cp compliance_plain.xml ' + str(folder) + '/' + str(x) + "-" + str(y) + '_msfem' + str(index) + '_y.xml')
              index += 1
          os.system('cp cellproblem_fine.xml ' + str(folder) + '/' + str(x) + "-" + str(y) + '_msfem.xml')
          if args.oversampling:
            jobfile.write('cfs.rel -m ' + str(args.oversampling) + ' -x ' + overdensfilename + ' ' + str(x) + "-" + str(y) + '_msfem \n')
          else:
            jobfile.write('cfs.rel -m ' + str(args.msfem) + ' -x ' + densfilename + ' ' + str(x) + "-" + str(y) + '_msfem \n')
          print(str(x) + ' ' + str(y) + ' is done')
        if args.design:
            # stop calculations if only one point is needed (debug)
            x = steps + 1
            y = steps + 1
        else:
            if args.penalization:
              x = x_tmp
              y = y_tmp
            y += 1
      if not args.design:
        x +=1    
  elif args.shape == "frame" or args.shape == "frame_modified" or args.shape == "auxetic":
    # Homogenization for 2D frame structures
    array = void * np.ones((minres, minres))
    joblist = ()
    steps_p = 0
    x = 0
    while x < steps + 1:
      y= 0
      while y < steps + 1:
        if args.penalization:
          steps_p = args.penalization
          x_tmp = x
          y_tmp = y
          x = float(x) / steps_p if steps_p > 0 else 0
          y = float(y) / steps_p if steps_p > 0 else 0
          if steps_p <= 0:
            print("ERROR: steps_p is 0 or smaller")
          if args.design:
            tmp = args.design.split(',')
            x = float(steps * float(tmp[0]))
            y = float(steps * float(tmp[1]))  
        elif args.design:
          tmp = args.design.split(',')
          x = int(steps * float(tmp[0]))
          y = int(steps * float(tmp[1]))
        if args.penalization:     
          densfilename = str(x_tmp) + "-" + str(y_tmp) + "_" + str(1./steps_p) + ".dens.xml"
          problem = str(x_tmp)+ "-" +str(y_tmp)  + "_" + str(1./steps_p)
        else:
          densfilename = str(x) + "-" + str(y) + ".dens.xml"
          problem = str(x)+ "-" +str(y) 
                    
        if args.shape == "frame_modified":
          array = insert_modified_frame(array, minres, y, x, steps, void, args.epsilon, steps_p, True)
        elif args.shape == "auxetic":
          #os.system("/Applications/Gmsh.app/Contents/MacOS/gmsh "+ str(args.gmsh)+ "/"+ problem +".stp -3 -optimize_netgen ")
          #os.system("/Applications/Gmsh.app/Contents/MacOS/gmsh "+ str(args.gmsh)+"/"+ problem +".msh -refine ")
          os.system("gmsh "+ str(args.gmsh)+ "/"+ problem +".stp -3 -optimize_netgen ")
          os.system("gmsh "+ str(args.gmsh)+"/"+ problem +".msh -refine ")
          os.system("gmsh "+ str(args.gmsh)+"/"+ problem +".msh -refine ")
          create_mesh_from_gmsh(str(args.gmsh)+"/" + problem,"aux_cells")
          os.system('cp ' +str(args.gmsh)+"/" + problem +".mesh " + str(folder) + '/' + problem + '.mesh')
        else:
          array = insert_modified_frame(array, minres, y, x, steps, void, args.epsilon, steps_p, False)
        # filtering of the data
        if args.filter == 'on':   
          array_filter = ndimage.uniform_filter(array, size=6)
        else:
          array_filter = array
        array = np.ones((minres, minres))
        if not args.gmsh:
          write_density_file(str(folder) + "/" + densfilename, array_filter, "set")
          # add new job to jobfile
          jobfile.write('cfs.rel -m ' + str(args.hom) + ' -x ' + densfilename + ' ' + problem + ' \n')
          # create xml file for cfs
          os.system('cp inv_tensor.xml ' + str(folder) + '/' + problem + '.xml')
        else:
          os.system('cp inv_tensor_3D_gmsh.xml ' + str(folder) + '/' + problem + '.xml')
          jobfile.write('cfs.rel ' + problem + ' \n')  
        print(problem  + ' is done')
        if args.design:
          # stop calculations if only one point is needed (debug)
          x = steps + 1
          y = steps + 1
        else:
          if args.penalization:
              x = x_tmp
              y = y_tmp
          y += 1
      if not args.design:
        x += 1 
  elif args.shape == "cross":
    # Homogenization for 2D cross structure
    array = void * np.ones((minres, minres))
    x = 0
    joblist = ()
    while x < steps + 1:
      y = 0
      while y < steps + 1:
        if args.penalization:
          steps_p = args.penalization
          x_tmp = x
          y_tmp = y
          x = float(x) / steps_p if steps_p > 0 else 0
          y = float(y) / steps_p if steps_p > 0 else 0
          if steps_p <= 0:
            print("ERROR: steps_p is 0 or smaller")
          if args.design:
            tmp = args.design.split(',')
            x = float(steps * float(tmp[0]))
            y = float(steps * float(tmp[1]))  
        elif args.design:
          tmp = args.design.split(',')
          x = int(steps * float(tmp[0]))
          y = int(steps * float(tmp[1]))
        if args.penalization:     
          densfilename = str(x_tmp) + "-" + str(y_tmp) + "_" + str(1./steps_p) + ".dens.xml"
          problem = str(x_tmp)+ "-" +str(y_tmp)  + "_" + str(1./steps_p)
        else:
          densfilename = str(x) + "-" + str(y) + ".dens.xml"
          problem = str(x)+ "-" +str(y) 
        if args.penalization:  
          offx = int((minres / 2.) * (1 - ((float(x) * steps_p / steps)**3./steps_p)) + 0.5)
          offy = int((minres / 2.) * (1 - ((float(y) * steps_p / steps)**3./steps_p)) + 0.5)
        else:
          offx = int((minres / 2.) * (1 - float(x) / (steps)) + 0.5)
          offy = int((minres / 2.) * (1 - float(y) / (steps)) + 0.5)
        for i in range(offx, minres - offx):
          for j in range(0, minres):
            array[j][i] = 1.
        for i in range(offy, minres - offy):
          for j in range(0, minres):
            array[i][j] = 1.
        # filtering of the data
        if args.filter == 'on':   
          array_filter = ndimage.uniform_filter(array, size=6)
        else:
          array_filter = array        
        write_density_file(str(folder) + "/" + densfilename, array_filter, "set")
        array = void * np.ones((minres, minres))
        # add new job to jobfile
        jobfile.write('cfs.rel -m ' + str(args.hom) + ' -x ' + densfilename + ' ' + problem + ' \n')
        # create xml file for cfs
        os.system('cp inv_tensor.xml ' + str(folder) + '/' + problem + '.xml')
        print(problem + ' is done')
        if args.design:
          # stop calculations if only one point is needed (debug)
          x = steps + 1
          y = steps + 1
        else:
          if args.penalization:
              x = x_tmp
              y = y_tmp
          y += 1
      if not args.design:
        x+=1
  else:
    print('option not defined')
elif dim == 3:
  # Homogenization for 3D cross structures
  joblist = ()
  x = 0
  while x < steps + 1:
    y = 0
    while y < steps + 1:
      z = 0
      while z < steps + 1:
        if args.penalization:
          steps_p = args.penalization
          x_tmp = x
          y_tmp = y
          z_tmp = z
          x = float(x) / steps_p if steps_p > 0 else 0
          y = float(y) / steps_p if steps_p > 0 else 0
          z = float(z) / steps_p if steps_p > 0 else 0
          if steps_p <= 0:
            print("ERROR: steps_p is 0 or smaller")
          if args.design:
            tmp = args.design.split(',')
            x = float(steps * float(tmp[0]))
            y = float(steps * float(tmp[1]))
            z = float(steps * float(tmp[2]))    
        elif args.design:
          tmp = args.design.split(',')
          x = int(steps * float(tmp[0]))
          y = int(steps * float(tmp[1]))
          z = int(steps * float(tmp[2]))
        if args.penalization:     
          densfilename = str(x_tmp) + "-" + str(y_tmp) + "-" + str(z_tmp) + "_" + str(1./steps_p) + ".dens.xml"
          problem = str(x_tmp)+ "-" +str(y_tmp)  + "-" + str(z_tmp) + "_" + str(1./steps_p)
        else:
          densfilename = str(x) + "-" + str(y) +"-" + str(z) + ".dens.xml"
          problem = str(x)+ "-" +str(y) + "-" + str(z) 
         
        if not args.gmsh and not args.shape == "ortho_bc":  
          array = void * np.ones((minres, minres, minres))
          if args.penalization:
            offx = int((minres / 2.) * (1. - ((float(x) * steps_p / steps)**3./steps_p)) + 0.5)
            offy = int((minres / 2.) * (1. - ((float(y) * steps_p / steps)**3./steps_p)) + 0.5)
            offz = int((minres / 2.) * (1. - ((float(z) * steps_p / steps)**3./steps_p)) + 0.5)
          else:
            offx = int((minres / 2.) * (1. - float(x) / (steps)) + 0.5)
            offy = int((minres / 2.) * (1. - float(y) / (steps)) + 0.5)
            offz = int((minres / 2.) * (1. - float(z) / (steps)) + 0.5)
          print("test: offx " + str(offx) + " test: offy " + str(offy) + " test: offz " + str(offz))
          for i in range(0, minres):
            for j in range(offx, minres - offx):
              for k in range(offx, minres - offx):
                array[i][j][k] = 1.
          for i in range(offy, minres - offy):
            for j in range(0, minres):
              for k in range(offy, minres - offy):
                array[i][j][k] = 1.
          for i in range(offz, minres - offz):
            for j in range(offz, minres - offz):
              for k in range(0, minres):
                array[i][j][k] = 1.
          if args.filter == 'on':   
            array_filter = ndimage.uniform_filter(array, size=2)
          else:
            array_filter = array        
          write_density_file(str(folder) + "/" + densfilename, array_filter, "set")
          array = np.ones((minres, minres, minres))                    
          # add new job to jobfile
          jobfile.write('cfs.rel -m ' + str(args.hom) + ' ' + problem + ' \n')
          # create xml file for cfs
          os.system('cp inv_tensor_3D.xml ' + str(folder) + '/' + problem + '.xml')
          file = str(folder) + '/' + problem + '.xml'
          if os.path.isfile(file):      
            doc = lxml.etree.parse(file, lxml.etree.XMLParser(remove_comments=True))
            replace(doc, "//cfs:loadErsatzMaterial/@file", problem + '.dens.xml')
            doc.write(file,pretty_print=True)
        else:
          if args.gmsh:             
            #os.system("/Applications/Gmsh.app/Contents/MacOS/gmsh "+ str(args.gmsh)+ "/"+ problem +".stp -3 -optimize_netgen ")
            #os.system("/Applications/Gmsh.app/Contents/MacOS/gmsh "+ str(args.gmsh)+"/"+ problem +".msh -refine ")
            os.system("gmsh "+ str(args.gmsh)+ "/"+ problem +".stp -3 -optimize_netgen ")
            os.system("gmsh "+ str(args.gmsh)+"/"+ problem +".msh -refine ")
            os.system('cp ' +str(args.gmsh)+"/" + problem +".msh " + str(folder) + '/' + problem + '.mesh')
            os.system('cp inv_tensor_3D_gmsh.xml ' + str(folder) + '/' + problem + '.xml')
            jobfile.write('cfs.rel ' + problem + ' \n')
          elif args.shape == "ortho_bc":
            os.system('cp inv_tensor_3D_ortho_bc.xml ' + str(folder) + '/' + problem + '.xml')
            # basecell script can only handle values 0< x < 1, therefore change x,y,z appropriately
            x1 = x
            y1 = y
            z1 = z
            x1 = 1e-5*steps if x == 0 else x1
            y1 = 1e-5*steps if y == 0 else y1
            z1 = 1e-5*steps if z == 0 else z1
            x1 = 0.9999*steps if x == steps else x1
            y1 = 0.9999*steps if y == steps else y1
            z1 = 0.9999*steps if z == steps else z1 
            #print("x1 = "+str(x1) + ", x1_cell = "+ str(float(x1)/float(steps)))
            os.system("../basecell.py --res "+str(args.res)+" --x1 "+str(float(x1)/float(steps))+" --y1 " +str(float(y1)/float(steps))+ " --z1 "+str(float(z1)/float(steps))+" --target volume_mesh --beta 7 --eta 0.6  --interpolation heaviside --bend 0.5 --save "+str(folder) + "/"+problem)
            jobfile.write('cfs.rel ' + problem + ' \n')
        print(problem + ' is done')
        if args.design:
          # stop calculations if only one point is needed (debug)
          x = steps + 1
          y = steps + 1
          z = steps + 1
        else:
          if args.penalization:
              x = x_tmp
              y = y_tmp
              z = z_tmp
          z += 1
      if not args.design:
        y += 1
    if not args.design:
      x += 1     
jobfile.close()
sys.exit()








