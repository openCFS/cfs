#!/usr/bin/env python

import sys, os, string
import numpy as np
import scipy.io
from numpy import dot, empty, size
from numpy import sin
from numpy import cos
from numpy import sqrt
from lxml import etree
import argparse
import scipy, scipy.sparse, scipy.sparse.linalg 
from hdf5_tools import *
from optimization_tools import *
from mesh_tool import *
from decimal import *
import os.path
import os
from cfs_utils import *
from subprocess import PIPE
import subprocess
import time
from io import StringIO 

tensors = numpy.zeros((11,9))
volumes = numpy.zeros((11,1))
scale = 116000.
tensors[0,:] = ['1.346154e-06', '5.769231e-07', '5.769231e-07', '1.346154e-06', '5.769231e-07', '1.346154e-06', '3.846154e-07', '3.846154e-07', '3.846154e-07']
tensors[10,:] = ['1.430976e+00', '6.734007e-01', '6.734007e-01', '1.430976e+00', '6.734007e-01','1.430976e+00', '3.787879e-01', '3.787879e-01', '3.787879e-01']
volumes[10] = 1.
name = "R.B.Fuller" 
upper = 8
lower = 1
name = "V7" 
upper = 5
lower = 0 

for i in range(1,10):
  infoxml = name +"_" + str(i*10) + "%.info.xml"
  if os.path.isfile(infoxml) and i > lower and i < upper:      
    doc = lxml.etree.parse(infoxml, lxml.etree.XMLParser(remove_comments=True))
    matrix = xpath(docgg, "//homogenizedTensor/tensor/real/text()")    
    res = list(map(float, matrix.split())) # convert list with string elements to list with float elements    
    ts = np.asarray(res)
    print(infoxml + ' -> '+ str(ts))
    tensors[i,:] = [ts[0],ts[1],ts[2],ts[7],ts[8],ts[14],ts[21],ts[28],ts[35]]
    vol = xpath(doc, "//domain/@structure_volume")
    hull_vol = xpath(doc, "//domain/@hull_volume")
    if vol != "cannot_determine":
      volumes[i] = float(vol)/float(hull_vol)
  else:
    print('file ' + infoxml + ' not found')
if name == "R.B.Fuller":
  tensors[1,:] = tensors[0,:] + (tensors[2,:]- tensors[0,:]) / 2.
  volumes[1,0] = volumes[0,0] + (volumes[2,0] - volumes[0,0])/ 2.
  tensors[8,:] = tensors[7,:] + (tensors[10,:] - tensors[7,:]) / 3.
  volumes[8,0] = volumes[7,0] + (volumes[10,0] - volumes[7,0])/ 3.
  tensors[9,:] = tensors[8,:] + (tensors[10,:] - tensors[8,:]) / 2.
  volumes[9,0] = volumes[8,0] + (volumes[10,0] - volumes[8,0])/ 2.
if name == "V7":
  tensors[5,:] = tensors[4,:] + (tensors[10,:] - tensors[4,:]) / 6.
  volumes[5,0] = volumes[4,0] + (volumes[10,0] - volumes[4,0])/ 6.
  tensors[6,:] = tensors[5,:] + (tensors[10,:] - tensors[5,:]) / 5.
  volumes[6,0] = volumes[5,0] + (volumes[10,0] - volumes[5,0])/ 5.
  tensors[7,:] = tensors[6,:] + (tensors[10,:] - tensors[6,:]) / 4.
  volumes[7,0] = volumes[6,0] + (volumes[10,0] - volumes[6,0])/ 4.
  tensors[8,:] = tensors[7,:] + (tensors[10,:] - tensors[7,:]) / 3.
  volumes[8,0] = volumes[7,0] + (volumes[10,0] - volumes[7,0])/ 3.
  tensors[9,:] = tensors[8,:] + (tensors[10,:] - tensors[8,:]) / 2.
  volumes[9,0] = volumes[8,0] + (volumes[10,0] - volumes[8,0])/ 2.

tensors = scale * tensors
print(tensors)
print(volumes)
# Read homogenized material tensors from cell problems in 3D and create detailed_stats table
steps = 10
filename = "detailed_stats_" + name
out = open(filename, "w")
out.write('  ' + str(steps) + '   ' + str(steps) + '  ' + str(steps) + '   0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00    0.000000e+00 0.000000e+00 0.000000e+00\n')
filename = "detailed_stats_vol_" + name
out_vol = open(filename, "w")
out_vol.write('  ' + str(steps) + '   ' + str(steps) + '  ' + str(steps) + '   0.000000e+00\n')
x = 0
while x < steps + 1:
  y = 0
  while y < steps + 1:
    z= 0
    while z < steps + 1:
      out.write(str(x).rjust(3) + ' ' + str(y).rjust(3) + ' ' + str(z).rjust(3) + ' ' + str(tensors[x,0]) + ' ' + str(tensors[x,1]) + ' ' + str(tensors[x,2]) + ' ' + str(tensors[x,3]) + ' ' + str(tensors[x,4]) + ' ' + str(tensors[x,5]) + ' ' + str(tensors[x,6]) + ' ' + str(tensors[x,7]) + ' ' + str(tensors[x,8]) + '\n')
      out_vol.write(str(x).rjust(3) + ' ' + str(y).rjust(3) + ' ' + str(z).rjust(3) + ' ' + str(volumes[x,0]) + '\n')
      z += 1
    y += 1
  x += 1
out.close()
out_vol.close()