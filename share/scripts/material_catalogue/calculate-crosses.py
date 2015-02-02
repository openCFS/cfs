#!/usr/bin/python

import Image, ImageDraw, sys, bz2, os
from scipy import ndimage
import numpy as np
import math
from optimization_tools import *
import libxml2
from cfs_utils import *
import os.path
import matplotlib.pyplot as plt

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
def insert_modified_frame(array, minres, x, y, steps, void, number, modify=True, triangle=False):
  # 2D frame structure
  if triangle:
    array = np.ones((2 * minres, minres))
  else:
    array = np.ones((minres, minres))
  eps = 1e-8
  offx = int((minres / (2.*number)) * (float(x) / (steps)) + 0.5 + eps)
  offy = int((minres / (2.*number)) * (float(y) / (steps)) + 0.5 + eps)
  for nx in range(0, number):
    for ny in range(0,number):
      rx = int(float(nx)/float(number)* minres+0.5+eps)
      ry = int(float(ny)/float(number)* minres+0.5+eps)
      endx = int(float(nx+1)/float(number)*minres+ 0.5 +eps)-offx
      endy = int(float(ny+1)/float(number)*minres+ 0.5 +eps)-offy
      for i in range(offx+rx, endx):
        for j in range(offy+ry, endy):
          if not triangle:
            array[i][j] = void
          else:
            array[2 * i][j] = void
            array[2 * i + 1][j] = void
      # TODO: does not work for triangles yet
      if modify:    
        # modify frame for stress minimization
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
  # rotationg angle should be between [0,pi] 
  # note that thickness of the bars not scaled by rotation angle
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
parser.add_argument("--shape", help="choose between frame or cross", choices=['frame', 'cross', 'frame_modified', 'frame_w_triangles'])
parser.add_argument("--triangle_msfem", help="true or false")
parser.add_argument("--filter", help="filtered densities on or off")
parser.add_argument("--void_material", help="set value for void material", type=float, default=1e-9)
parser.add_argument("--epsilon", help="number of frames/crosses in the cell problem", type=int, default=1)
parser.add_argument("--design", help="select single thicknesses s1,s2,s3 for debugging,e.g. 0.1,0.3,0.")




args = parser.parse_args()
pl = parser.parse_args()
dim = args.dimension
minres = args.res
steps = args.stp
folder = args.folder
void = args.void_material

if dim == 2:
  if args.msfem:
    if args.triangle_msfem:
      outfile = str(folder) + '_/jobs' 
      print "computing with " + str(steps) + " steps!"
      if not os.path.exists(str(folder)):
        os.mkdir(str(folder))
      os.system('cp mat.xml ' + str(folder))
    else:
      outfile = str(folder) + '/jobs' 
      print "computing with " + str(steps) + " steps!"
      if not os.path.exists(str(folder) + "/"):
        os.mkdir(str(folder))
      os.system('cp mat.xml ' + str(folder) + '/')
  else:
    outfile = str(folder) + '/jobs'
    print "computing with " + str(steps) + " steps!"
    os.mkdir(str(folder) + "/")
    os.system('cp mat.xml ' + str(folder) + '/')
if dim == 3:
  outfile = str(folder) + '/jobs'
  print "computing with " + str(steps) + " steps!"
  os.mkdir(str(folder) + '')
  os.system('cp mat.xml ' + str(folder) + '/')
  
print outfile
jobfile = open(outfile, "w")
if dim == 2:
  if args.msfem:
    array = void * np.ones((minres, minres))
    # steps_rot = 4
    # drot = math.pi/float(steps_rot)
    midx = minres / 2
    midy = minres / 2
    # for rot in range(-1.57,drot,1.57):
    rot = 0.
    x = 0
    while x < steps + 1:
      y = 0
      while y < steps + 1:
        if args.design:
          tmp = args.design.split(',')
          x = int(steps * float(tmp[0]))
          y = int(steps * float(tmp[1]))
        if args.shape == 'frame_modified':
          array = insert_modified_frame(array, minres, y, x, steps, void, args.epsilon, True)
        elif args.shape == 'cross':
          array = insert_rotated_bar(tx, ty, 0., array, minres, midx, midy)
        elif args.shape == 'frame':
          array = insert_modified_frame(array, minres, y, x, steps, void, args.epsilon, False)
        elif args.shape == 'frame_w_triangles':
          array = insert_modified_frame(array, minres, y, x, steps, void, args.epsilon, False, True)
        else:
          print 'Warning: base cell type is undefined, set --shape [frame or frame_modified or cross]'

        densfilename = str(x) + "-" + str(y) + "_msfem.dens.xml"
        if args.filter == 'on':
          # filtering of the data
          array_filter = ndimage.uniform_filter(array, size=6)
          plt.gray()
          plt.imshow(array_filter)
          plt.show()
        else:
          array_filter = array
        if not args.sparse_msfem:
          write_density_file(str(folder) + "/" + densfilename, array_filter, "set")
        array = void * np.ones((minres, minres))
        # WARNING: not tested any more
        if args.triangle_msfem:
          # create xml file for cfs
          doc = libxml2.parseFile("triangle_msfem.xml")
          xml = doc.xpathNewContext()
          xml.xpathRegisterNs('cfs', 'http://www.cfs++.org')
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
              replace(xml, '//cfs:dirichletInhom[@name="left"][@dof="x"]/@value', val[0])
              replace(xml, '//cfs:dirichletInhom[@name="bottom"][@dof="x"]/@value', val[1])
              replace(xml, '//cfs:dirichletInhom[@name="right"][@dof="x"]/@value', val[2])
              replace(xml, '//cfs:dirichletInhom[@name="left"][@dof="y"]/@value', str(0.))
              replace(xml, '//cfs:dirichletInhom[@name="bottom"][@dof="y"]/@value', str(0.))
              replace(xml, '//cfs:dirichletInhom[@name="right"][@dof="y"]/@value', str(0.))
            else:
              # degree of fredom y
              replace(xml, '//cfs:dirichletInhom[@name="left"][@dof="y"]/@value', val[0])
              replace(xml, '//cfs:dirichletInhom[@name="bottom"][@dof="y"]/@value', val[1])
              replace(xml, '//cfs:dirichletInhom[@name="right"][@dof="y"]/@value', val[2])
              replace(xml, '//cfs:dirichletInhom[@name="left"][@dof="x"]/@value', str(0.))
              replace(xml, '//cfs:dirichletInhom[@name="bottom"][@dof="x"]/@value', str(0.))
              replace(xml, '//cfs:dirichletInhom[@name="right"][@dof="x"]/@value', str(0.))
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
              doc.saveFile(str(folder) + '/' + str(x) + "-" + str(y) + '_msfem' + str(index) + '_x.xml')
            else:
	      doc.saveFile(str(folder) + '/' + str(x) + "-" + str(y) + '_msfem' + str(index) + '_y.xml')
            if i % 2 == 0:
              # add new job to jobfile
              jobfile.write('cfs.rel -m ~/meshes/' + str(args.msfem) + ' -x ' + densfilename + ' ' + str(x) + "-" + str(y) + '_msfem' + str(index) + '_x \n')
              os.system('cp triangle_msfem.xml ' + str(folder) + '/' + str(x) + "-" + str(y) + '_msfem' + str(index) + '_x.xml')
            else:
              # add new job to jobfile
              jobfile.write('cfs.rel -m ~/meshes/' + str(args.msfem) + ' -x ' + densfilename + ' ' + str(x) + "-" + str(y) + '_msfem' + str(index) + '_y \n')
              os.system('cp triangle_msfem.xml ' + str(folder) + '/' + str(x) + "-" + str(y) + '_msfem' + str(index) + '_y.xml')
              index += 1
          print str(x) + ' ' + str(y) + ' is done'
        else:
          # MSFEM for quadrileteral elements 
          # create xml file for cfs
          doc = libxml2.parseFile("compliance_plain.xml")
          xml = doc.xpathNewContext()
          xml.xpathRegisterNs('cfs', 'http://www.cfs++.org')
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
              replace(xml, '//cfs:dirichletInhom[@name="left"][@dof="x"]/@value', val[0])
              replace(xml, '//cfs:dirichletInhom[@name="bottom"][@dof="x"]/@value', val[1])
              replace(xml, '//cfs:dirichletInhom[@name="right"][@dof="x"]/@value', val[2])
              replace(xml, '//cfs:dirichletInhom[@name="top"][@dof="x"]/@value', val[3])
              replace(xml, '//cfs:dirichletInhom[@name="left"][@dof="y"]/@value', str(0.))
              replace(xml, '//cfs:dirichletInhom[@name="bottom"][@dof="y"]/@value', str(0.))
              replace(xml, '//cfs:dirichletInhom[@name="right"][@dof="y"]/@value', str(0.))
              replace(xml, '//cfs:dirichletInhom[@name="top"][@dof="y"]/@value', str(0.))
            else:
              # degree of fredom y
              replace(xml, '//cfs:dirichletInhom[@name="left"][@dof="y"]/@value', val[0])
              replace(xml, '//cfs:dirichletInhom[@name="bottom"][@dof="y"]/@value', val[1])
              replace(xml, '//cfs:dirichletInhom[@name="right"][@dof="y"]/@value', val[2])
              replace(xml, '//cfs:dirichletInhom[@name="top"][@dof="y"]/@value', val[3])
              replace(xml, '//cfs:dirichletInhom[@name="left"][@dof="x"]/@value', str(0.))
              replace(xml, '//cfs:dirichletInhom[@name="bottom"][@dof="x"]/@value', str(0.))
              replace(xml, '//cfs:dirichletInhom[@name="right"][@dof="x"]/@value', str(0.))
              replace(xml, '//cfs:dirichletInhom[@name="top"][@dof="x"]/@value', str(0.))
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
              doc.saveFile(str(folder) + '/' + str(x) + "-" + str(y) + '_msfem' + str(index) + '_x.xml')
            else:
	      doc.saveFile(str(folder) + '/' + str(x) + "-" + str(y) + '_msfem' + str(index) + '_y.xml')
            # sparse is not working currently
            if args.sparse_msfem:
              mesh = str(x) + "-" + str(y) + '.mesh'
            else:
              mesh = '~/meshes/' + args.msfem
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
          jobfile.write('cfs.rel -m ~/meshes/' + str(args.msfem) + ' -x ' + densfilename + ' ' + str(x) + "-" + str(y) + '_msfem \n')
          print str(x) + ' ' + str(y) + ' is done'
        if args.design:
            # stop calculations if only one point is needed (debug)
            x = steps + 1
            y = steps + 1
        else:
            y += 1
      if not args.design:
        x +=1    
  elif args.shape == "frame" or args.shape == "frame_modified":
    # 2D frame structure
    array = void * np.ones((minres, minres))
    x = 0
    while x < steps + 1:
      y= 0
      while y < steps + 1:
        if args.design:
          tmp = args.design.split(',')
          x = int(steps * float(tmp[0]))
          y = int(steps * float(tmp[1]))
        if args.shape == "frame_modified":
          array = insert_modified_frame(array, minres, y, x, steps, void, args.epsilon, True) 
        else:
          array = insert_modified_frame(array, minres, y, x, steps, void, args.epsilon, False)     
        densfilename = str(x) + "-" + str(y) + ".dens.xml"
        # filtering of the data
        if args.filter == 'on':   
          array_filter = ndimage.uniform_filter(array, size=6)
        else:
          array_filter = array
        write_density_file(str(folder) + "/" + densfilename, array_filter, "set")
        array = np.ones((minres, minres))
        # add new job to jobfile
        jobfile.write('cfs.rel -m ~/meshes/' + str(args.hom) + ' -x ' + densfilename + ' ' + str(x) + "-" + str(y) + ' \n')
        # create xml file for cfs
        os.system('cp inv_tensor.xml ' + str(folder) + '/' + str(x) + "-" + str(y) + '.xml')  
        print str(x) + ' ' + str(y) + ' is done'
        if args.design:
          # stop calculations if only one point is needed (debug)
          x = steps + 1
          y = steps + 1
        else:
          y += 1
      if not args.design:
        x += 1 
  elif args.shape == "cross":
    # 2D cross structure
    array = void * np.ones((minres, minres))
    x = 0
    while x < steps + 1:
      y = 0
      while y < steps + 1:
        if args.design:
          tmp = args.design.split(',')
          x = int(steps * float(tmp[0]))
          y = int(steps * float(tmp[1]))
        offx = int((minres / 2.) * (1 - float(x) / (steps)) + 0.5)
        offy = int((minres / 2.) * (1 - float(y) / (steps)) + 0.5)
        for i in range(offx, minres - offx):
          for j in range(0, minres):
            array[j][i] = 1.
        for i in range(offy, minres - offy):
          for j in range(0, minres):
            array[i][j] = 1.
        densfilename = str(x) + "-" + str(y) + ".dens.xml"
        # filtering of the data
        if args.filter == 'on':   
          array_filter = ndimage.uniform_filter(array, size=6)
        else:
          array_filter = array        
        write_density_file(str(folder) + "/" + densfilename, array_filter, "set")
        array = void * np.ones((minres, minres))
        # add new job to jobfile
        jobfile.write('cfs.rel -m ~/meshes/' + str(args.hom) + ' -x ' + densfilename + ' ' + str(x) + "-" + str(y) + ' \n')
        # create xml file for cfs
        os.system('cp inv_tensor.xml ' + str(folder) + '/' + str(x) + "-" + str(y) + '.xml')
        print str(x) + ' ' + str(y) + ' is done'
        if args.design:
          # stop calculations if only one point is needed (debug)
          x = steps + 1
          y = steps + 1
        else:
          y += 1
      if not args.design:
        x+=1
  else:
    print 'option not defined' 
elif dim == 3:
  # 3D cross structure
  x = 0
  while x < steps + 1:
    y = 0
    while y < steps + 1:
      z = 0
      while z < steps + 1:
        if args.design:
          tmp = args.design.split(',')
          x = steps * int(tmp[0])
          y = steps * int(tmp[1])
          z = steps * int(tmp[2])
        array = void * np.ones((minres, minres, minres))
        offx = int((minres / 2.) * (1 - float(x) / (steps)) + 0.5)
        offy = int((minres / 2.) * (1 - float(y) / (steps)) + 0.5)
        offz = int((minres / 2.) * (1 - float(z) / (steps)) + 0.5)
        print "test: offx " + str(offx) + " test: offy " + str(offy) + " test: offz " + str(offz)
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
        densfilename = str(x) + "-" + str(y) + "-" + str(z) + ".dens.xml"
        if args.filter == 'on':   
          array_filter = ndimage.uniform_filter(array, size=2)
        else:
          array_filter = array        
        write_density_file(str(folder) + "/" + densfilename, array_filter, "set")
        array = np.ones((minres, minres, minres))
        # add new job to jobfile
        jobfile.write('cfs.rel -m ~/meshes/' + str(args.hom) + ' ' + str(x) + '-' + str(y) + '-' + str(z) + ' \n')
        # create xml file for cfs
        os.system('cp inv_tensor_3D.xml ' + str(folder) + '/' + str(x) + '-' + str(y) + '-' + str(z) + '.xml')
        file = str(folder) + '/' + str(x) + '-' + str(y) + '-' + str(z) + '.xml'
        if os.path.isfile(file):      
          doc = libxml2.parseFile(file)
          xml = doc.xpathNewContext()
          xml.xpathRegisterNs('cfs', 'http://www.cfs++.org')
          replace(xml, "//cfs:loadErsatzMaterial/@file", str(x) + '-' + str(y) + '-' + str(z) + '.dens.xml')
          doc.saveFile(file)
        print str(x) + ' ' + str(y) + ' ' + str(z) + ' is done'
        if args.design:
          # stop calculations if only one point is needed (debug)
          x = steps + 1
          y = steps + 1
          z = steps + 1
        else:
          z += 1
      if not args.design:
        y += 1
    if not args.design:
      x += 1      
jobfile.close()
sys.exit()








