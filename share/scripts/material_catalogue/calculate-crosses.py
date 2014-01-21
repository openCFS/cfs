#!/usr/bin/python

import Image, ImageDraw, sys, bz2, os
from scipy import ndimage
import numpy as np
import math
from optimization_tools import *
import libxml2
from cfs_utils import *

import argparse

#lena = misc.lena()
# blurred_lena = ndimage.gaussian_filter(lena, sigma=3)
# very_blurred = ndimage.gaussian_filter(lena, sigma=5)


def save_image_as_densfile(im, outfile):
  imlist = list(im.getdata())
  #out = open(outfile, "w")
  out = bz2.BZ2File(outfile, 'wb')
  out.write('<?xml version="1.0"?>\n <cfsErsatzMaterial infoWriteCounter="0" \
  infoRejectCounter="0">\n <header>\n <mesh x="'+str(minres)+'" y="'+str(minres)+'" \
  z="1"/>\n <design adapt_lower="true" bimaterial="janrescaledvoid" \
  constant="false" fixed="false" initial="0.5" lower="0" name="density" \
  region="mech" scale="false" upper="1"/>\n <transferFunction \
  application="mech" design="density" param="'+str(p)+'" type="simp"/>\n \
  </header>\n<set id="frompng">\n')

  el = 1
  for i in imlist:
    val = (255.-float(i))/255.0
    out.write('<element nr="' + str(el) + '" type="density" design="'\
      + str(val**(1./p)) + '" physical="' + str(val) + '"/>\n')
    el += 1


  out.write('</set>\n </cfsErsatzMaterial>')
  out.close()
 
parser = argparse.ArgumentParser()
parser.add_argument("stp", help="number of grid points in one direction",type=int)
parser.add_argument("dimension", help="Dimension of the problem",type=int,default=2)
parser.add_argument("--res", help="grid size", type=int, default=100)
parser.add_argument("--point", help="point list", type =np.array,default=np.zeros((1,1)))

args = parser.parse_args()
pl = parser.parse_args()
dim = args.dimension
minres = args.res
steps = args.stp

if minres ==100:
  if dim==3:
    minres=20
p = 1.0
outfile = str(steps) + '/jobs'
if dim==2:
  outfile = str(steps) + '/jobs'
  if os.path.isdir(str(steps)):
    print "already computed with " + str(steps) + " steps, just saying..."
  else:
    print "computing with " + str(steps) + " steps!"
    os.mkdir(str(steps))
    os.system('cp mat.xml ' + str(steps) + '/')
elif dim ==3:
  outfile = str(steps) + '_3D/jobs'
  if os.path.isdir(str(steps)+'_3D'):
    print "already computed with " + str(steps) + " steps, just saying..."
  else:
    print "computing with " + str(steps) + " steps!"
    os.mkdir(str(steps)+'_3D')
    os.system('cp mat.xml ' + str(steps) + '_3D/')
  

jobfile = open(outfile, "w")
if dim == 4:
  # OLD Implementation
  for x in range(steps):
    for y in range(steps):
      if x + y == 0:
	       continue
      print str(x) + " " + str(y)
      steps2 = 2*(steps)
      im = Image.new('L', (steps2, steps2), 0) # create a blank image
      draw = ImageDraw.Draw(im) # create draw object
      #draw.polygon([(x,y),(x,im.size[1]-y),(im.size[0]-x,im.size[1]-y),(im.size[0]-x,y)], fill='white')
      draw.rectangle([(steps2-x-1,y),(x,steps2-y-1)],fill='white')
      # now draw the lines into the picture
      #for xx in range(x):
      #  draw.line((0, xx) + (steps, xx), fill=0)
      #for yy in range(y):
      #  draw.line((yy, 0) + (yy, steps), fill=0)

      del draw 

      # resize image to proper resolution
      if steps < minres:
	         im = im.resize((minres, minres))
      # create a smooth image by a gaussian_filter	
      im = ndimage.gaussian_filter(im, 5)
      # save image
      imagename = str(steps) + "/" + str(x) + "-" + str(y) + ".png"
      #im.save(imagename)
      #imagename = str(steps) + "/" + str(x) + "-" + str(y) + "2.png"




      #Rescale to 0-255 and convert to uint8
      rescaled = (255.0 / im.max() * (im - im.min())).astype(np.uint8)

      im = Image.fromarray(rescaled)
      im.save(imagename)

      densfilename = str(x) + "-" + str(y) + ".dens.xml.bz2"
      save_image_as_densfile(im, str(steps) + "/" + densfilename)

      # add new job to jobfile
      jobfile.write('cfs.rel -m ~/meshes/' + str(minres) + '.mesh -x ' + densfilename + ' ' + str(x) + "-" + str(y) + ' \n')

      # create xml file for cfs
      os.system('cp inv_tensor.xml ' + str(steps) + '/' + str(x) + "-" + str(y) + '.xml')
elif dim == 2:
  #cross structure
  array = 1e-6*np.ones((minres,minres))
  for x in range(steps+1):
    for y in range(steps+1):
      offx = int((minres/2.)*(1-float(x)/(steps))+0.5)
      offy = int((minres/2.)*(1-float(y)/(steps))+0.5)
      for i in range(offx,minres-offx):
        for j in range(0,minres):
          array[i][j] = 1.
      for i in range(offy,minres-offy):
        for j in range(0,minres):
          array[j][i] = 1.
      densfilename = str(x) + "-" + str(y) +".dens.xml"
      print str(numpy.sum(array)/(minres*minres)) + " vs " + str(float(x)/steps + float(y)/steps - (float(x)/steps)*(float(y)/steps))
      # filtering of the data
      array_filter = ndimage.uniform_filter(array, size = 6)
      write_density_file(str(steps) + "/" + densfilename,array_filter,"set")
      array = 1e-6*np.ones((minres,minres))
      # add new job to jobfile
      jobfile.write('cfs.rel -m ~/meshes/' + str(minres) + '.mesh -x ' + densfilename + ' ' + str(x) + "-" + str(y)+' \n')
      # create xml file for cfs
      os.system('cp inv_tensor.xml ' + str(steps) + '/' + str(x) + "-" + str(y) + '.xml')    
elif dim == 10:    
  #frame structure
  array = np.ones((minres,minres))
  for x in range(steps+1):
    for y in range(steps+1):
      offx = int((minres/2.)*(float(x)/(steps))+0.5)
      offy = int((minres/2.)*(float(y)/(steps))+0.5)
      print "test: offx " + str(offx) + " test: offy " + str(offy)
      print "test: x " + str(x) + " test: y " + str(y) + "-> " + str(float(x)/steps) + ", "+ str(float(y)/steps)
      for i in range(offx,minres-offx):
        for j in range(offy,minres-offy):
          array[i][j] = 1e-6
      densfilename = str(x) + "-" + str(y) +".dens.xml"
      print str(numpy.sum(array)/(minres*minres)) + " vs " + str(float(x)/steps + float(y)/steps - (float(x)/steps)*(float(y)/steps))
      # filtering of the data
   
      array_filter = ndimage.uniform_filter(array, size = 6)
      write_density_file(str(steps) + "/" + densfilename,array_filter,"set")
      array = np.ones((minres,minres))
      # add new job to jobfile
      jobfile.write('cfs.rel -m ~/meshes/' + str(minres) + '.mesh -x ' + densfilename + ' ' + str(x) + "-" + str(y)+' \n')
      # create xml file for cfs
      os.system('cp inv_tensor.xml ' + str(steps) + '/' + str(x) + "-" + str(y) + '.xml')           
elif dim == 3:
  # 3D cross structure
  for x in range(steps+1):
    for y in range(steps+1):
      for z in range(steps+1):
        array = 1e-6*np.ones((minres,minres,minres))
        offx = int((minres/2.)*(1-float(x)/(steps))+0.5)
        offy = int((minres/2.)*(1-float(y)/(steps))+0.5)
        offz = int((minres/2.)*(1-float(z)/(steps))+0.5)
        print "test: offx " + str(offx) + " test: offy " + str(offy) + " test: offz " + str(offz)
        for i in range(0,minres):
          for j in range(offx,minres-offx):
            for k in range(offx,minres-offx):
              array[i][j][k] = 1.
        for i in range(offy,minres-offy):
          for j in range(0,minres):
            for k in range(offy,minres-offy):
              array[i][j][k] = 1.
        for i in range(offz,minres-offz):
          for j in range(offz,minres-offz):
            for k in range(0,minres):
              array[i][j][k] = 1.
        densfilename = str(x) + "-" + str(y) + "-" + str(z) +".dens.xml"
        array_filter = ndimage.uniform_filter(array, size = 2)
        write_density_file(str(steps) + "_3D/" + densfilename,array_filter,"set")
        array = np.ones((minres,minres,minres))
        # add new job to jobfile
        jobfile.write('cfs.rel -m ~/meshes/'+'20x20x20'+ '.mesh ' + str(x) + '-' + str(y) + '-' + str(z) + ' \n')
        # create xml file for cfs
        os.system('cp inv_tensor_3D.xml ' + str(steps) + '_3D/' + str(x) + '-' + str(y) + '-' + str(z) + '.xml')
        file = str(steps) + '_3D/' + str(x) + '-' + str(y) + '-' + str(z) + '.xml'
        if os.path.isfile(file):      
          doc = libxml2.parseFile(file)
          xml = doc.xpathNewContext()
          xml.xpathRegisterNs('cfs', 'http://www.cfs++.org')
          replace(xml, "//cfs:loadErsatzMaterial/@file", str(x) + '-' + str(y) + '-' + str(z) + '.dens.xml')
          doc.saveFile(file)
elif dim ==11:
  for x in range(steps+1):
    for y in range(steps+1):
      for z in range(steps+1):
        offx = int((minres/2.)*(float(x)/(steps))+0.5)
        offy = int((minres/2.)*(float(y)/(steps))+0.5)
        offz = int((minres/2.)*(float(z)/(steps))+0.5)
	tmp  = -1.
	if pl.size > 3:
         tmp = 1.
        print "test: offx "+str(offx)+" test: offy "+str(offy)+" test: offz "+str(offz)
        print "est: x " + str(x) + " test: y " + str(y) + " test: z " + str(z)
	#if tmp > 0:
	#  for i in range(0,pl
        for i in range(offx,minres-offx):
          for j in range(offy,minres-offy):
            for k in range(offz,minres-offz):
              array[i][j][k] = 1e-6
        densfilename = str(x) + "-" + str(y) + "-" + str(z) +".dens.xml"
        array_filter = ndimage.uniform_filter(array, size = 6)
        write_density_file(str(steps) + "_3D/" + densfilename,array_filter,"set")
        array = np.ones((minres,minres,minres))
        # add new job to jobfile
        jobfile.write('cfs.rel -m ~/meshes/'+'20x20x20'+ '.mesh ' + str(x) + '-' + str(y) + '-' + str(z) + ' \n')
        # create xml file for cfs
        os.system('cp inv_tensor_3D.xml ' + str(steps) + '_3D/' + str(x) + '-' + str(y) + '-' + str(z) + '.xml')
        file = str(steps) + '_3D/' + str(x) + '-' + str(y) + '-' + str(z) + '.xml'
        if os.path.isfile(file):      
          doc = libxml2.parseFile(file)
          xml = doc.xpathNewContext()
          xml.xpathRegisterNs('cfs', 'http://www.cfs++.org')
          replace(xml, "//cfs:loadErsatzMaterial/@file", str(x) + '-' + str(y) + '-' + str(z) + '.dens.xml')
          doc.saveFile(file)
jobfile.close()
sys.exit()








