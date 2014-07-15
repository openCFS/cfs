#!/usr/bin/env python

import sys, os, string
import numpy
from numpy import dot
from numpy import sin
from numpy import cos
from numpy import sqrt
from lxml import etree
import argparse

from optimization_tools import *





## This rotates a 3*3 2D tensor via the third direction. As in Richter and CFS
def get_rot_3x3(angle): 
  R = numpy.zeros((2,2))

  R[0][0] =  cos(angle)
  R[0][1] =  sin(angle)
  R[1][0] =  -sin(angle)
  R[1][1] =  cos(angle)
  
  Q = numpy.zeros((3,3))
  
  Q[0][0] = R[0][0]*R[0][0];
  Q[0][1] = R[0][1]*R[0][1];
  Q[0][2] = 2.0*R[0][0]*R[0][1];

  Q[1][0] = R[1][0]*R[1][0];
  Q[1][1] = R[1][1]*R[1][1];
  Q[1][2] = 2.0*R[1][0]*R[1][1];

  Q[2][0] = R[0][0]*R[1][0];
  Q[2][1] = R[0][1]*R[1][1];
  Q[2][2] = R[0][0]*R[1][1] + R[0][1]*R[1][0];
  
  return Q
parser = argparse.ArgumentParser()
parser.add_argument("stp", help="number of grid points in one direction",type=int)
parser.add_argument("dimension", help="Dimension of the problem",type=int,default=2)
args = parser.parse_args()

steps = args.stp
dim = args.dimension
if dim == 2:
  filename = "detailed_stats_"+str(steps)
  out = open(filename, "w")
  out.write('  '+str(steps)+'   '+str(steps)+'   0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 \n')
elif dim ==3:
  filename = "detailed_stats_"+str(steps)+'_3D'
  out = open(filename, "w")
  out.write('  '+str(steps)+'   '+str(steps)+ '  '+str(steps)+'   0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00    0.000000e+00 0.000000e+00 0.000000e+00\n')
  

  #out.write('# x   y 11           12           22           33\n')
if dim== 2:
  for x in range(steps+1):
    for y in range(x+1):
      infoxml = str(steps) + "/" + str(x) + "-" + str(y)+ ".info.xml"
      #print infoxml
      if os.path.isfile(infoxml):      
        doc = libxml2.parseFile(infoxml)
        xml = doc.xpathNewContext()
        # complex values!  
        print infoxml + ' -> '
        tensortext  = xpath(xml, "//homogenizedTensor/tensor/real")	
        print infoxml + ' -> ' + tensortext	
        ts = tensortext.split()
        #print ts[0]
        #print ts[1]
        #print ts[4]
        #print ts[8]
  
        out.write(str(x).rjust(3) + ' ' + str(y).rjust(3) + ' ' + ts[0] + ' ' + ts[1] + ' ' + ts[4] + ' ' + ts[8] + '\n')
        if x != y:
          out.write(str(y).rjust(3) + ' ' + str(x).rjust(3) + ' ' + ts[4] + ' ' + ts[1] + ' ' + ts[0] + ' ' + ts[8] + '\n')
  
        if False:
          tensor = numpy.zeros((3,3))
          tensor[0,0] = ts[0]
          tensor[0,1] = ts[1]
          tensor[0,2] = ts[2]
          tensor[1,0] = ts[3]
          tensor[1,1] = ts[4]
          tensor[1,2] = ts[5]
          tensor[2,0] = ts[6]
          tensor[2,1] = ts[7]
          tensor[2,2] = ts[8]
          Q = get_rot_3x3(numpy.pi/2.0)
          rott = dot(Q, dot(tensor, Q.transpose()))
          out.write(str(y).rjust(3) + ' ' + str(x).rjust(3) + ' {:.6e} {:.6e} {:.6e} {:.6e}\n'.format(rott[0,0],rott[0,1],rott[1,1],rott[2,2]))
      else:
        print 'file ' + infoxml + ' not found'
elif dim == 3:
  for x in range(steps+1):
    for y in range(steps+1):
      for z in range(steps+1):
        infoxml = str(steps) + "_3D/" + str(x) + "-" + str(y)+ "-"+ str(z)+ ".info.xml"
        #print infoxml
        if os.path.isfile(infoxml):      
          doc = libxml2.parseFile(infoxml)
          xml = doc.xpathNewContext()
          # complex values!  
          print infoxml + ' -> '
          tensortext  = xpath(xml, "//homogenizedTensor/tensor/real")  
          print infoxml + ' -> ' + tensortext  
          ts = tensortext.split()
          #print ts[0]
          #print ts[1]
          #print ts[4]
          #print ts[8]
          out.write(str(x).rjust(3) + ' ' + str(y).rjust(3) + ' '+ str(z).rjust(3) +' ' + ts[0] + ' ' + ts[1] + ' ' +ts[2] +' '+ ts[7] + ' '+ ts[8] + ' '+ts[14] + ' '+ts[21] + ' '+ ts[28]+ ' '+ts[35] + '\n')
          #if x != y:
          #  out.write(str(y).rjust(3) + ' ' + str(x).rjust(3) + ' ' + ts[4] + ' ' + ts[1] + ' ' + ts[0] + ' ' + ts[8] + '\n')
        else:
          print 'file ' + infoxml + ' not found'

#print "number of tensors computed = " + str(count) + " / " + str(countall) + " (= {0:2.2f}%)".format(100.0*float(count)/float(countall))
#out.write("number of tensors computed = " + str(count) + "\n")

#out.write(str(steps) + ' ' + str(steps) + ' 5.698006e+00 1.994302e+00 5.698006e+00 1.851852e+00')
out.close()
