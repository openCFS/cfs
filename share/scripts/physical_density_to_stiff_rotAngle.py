#!/usr/bin/env python
from optimization_tools import *
import numpy
import argparse
import math

parser = argparse.ArgumentParser()
parser.add_argument("input", help="density.xml file where pyhsical is read")
parser.add_argument("output", help="density.xml file where pyhsical is written as stiff1 and stiff2 and optional rotAngle")
parser.add_argument("--threshold", help="if within (0,1) input data is thresholded", type=float, default=-1.0)
parser.add_argument("--rot", help="rotangle from different density.xml")
parser.add_argument("--refangle", help="reference angle given, swaps stiff1 and stiff2 if benefitial for smooth angle", type=float)
parser.add_argument("--langle", help="lower bound for angle", type=float, default=-1.57)
parser.add_argument("--uangle", help="upper bound for angle", type=float, default=1.57)
parser.add_argument("--angle", help="set rotAngle between -pi/2 and pi/2", type=float)


args = parser.parse_args()
d = read_density_as_vector(args.input, "design")
if args.rot:
  d2 = read_multi_design(args.rot, "stiff1", "stiff2", "rotAngle")
if args.threshold > 0:
  t = threshold_filter(d, args.threshold, 1e-7, 1.0)
  print 'perform threshold at ' + str(args.threshold) + ": avg density " + str(numpy.sum(d) / len(d)) + " -> " + str(numpy.sum(t) / len(d))  
  d = t
if args.refangle:
    d2 = read_multi_design(args.input, "stiff1", "stiff2", "rotAngle")
    for i in range(len(d2)):
      swap = 0.
      drot = abs(args.refangle - d2[i, 2])
      drot1 = abs(args.refangle - (d2[i, 2] + math.pi / 2.))
      drot2 = abs(args.refangle - (d2[i, 2] - math.pi / 2.))
      if drot > drot1:
        swap = 1.
        if drot1 > drot2:
          swap = 2.
      elif drot > drot2:
        swap = 2.
      if swap == 1.:
        d2[i, 2] = d2[i, 2] + math.pi / 2.
        tmp = d2[i, 1]
        d2[i, 1] = d2[i, 2]
        d2[i, 2] = tmp
      elif swap == 2.:
        d2[i, 2] = d2[i, 2] - math.pi / 2.
        tmp = d2[i, 1]
        d2[i, 1] = d2[i, 2]
        d2[i, 2] = tmp
    write_multi_design_file(args.output, d2, ["stiff1", "stiff2", "rotAngle"])
else:
  if args.rot:
    data = numpy.zeros((len(d), 3))
    for i in range(len(d2)):
      if args.angle:
          data[i, 0] = d2[i, 1]
          data[i, 1] = d2[i, 2]
          data[i, 2] = args.angle    
      else:
          data[i, 0] = d[i]
          data[i, 1] = d[i]
          data[i, 2] = d2[i, 2]         
    write_multi_design_file(args.output, data, ["stiff1", "stiff2", "rotAngle"])
  else:
    data = numpy.zeros((len(d), 3))
    for i in range(len(d)): 
        # data[i,0] = 1.-numpy.sqrt(1.-d[i])
        # data[i,1] = 1.-numpy.sqrt(1.-d[i])
        data[i, 0] = d[i]
        data[i, 1] = d[i]
        data[i, 2] = 0.
        if args.rot:
          data[i, 2] = d2[i, 2]
        if args.angle:
          data[i, 2] = args.angle      
    write_multi_design_file(args.output, data, ["stiff1", "stiff2", "rotAngle"])
