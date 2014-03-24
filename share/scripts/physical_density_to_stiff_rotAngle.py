#!/usr/bin/env python
from optimization_tools import *
import numpy
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("input", help="density.xml file where pyhsical is read")
parser.add_argument("output", help="density.xml file where pyhsical is written as stiff1 and stiff2")
parser.add_argument("--threshold", help="if within (0,1) input data is thresholded", type=float, default=-1.0)
parser.add_argument("--rot", help="rotangle from different density.xml")

args = parser.parse_args()
d = read_density_as_vector(args.input, "design")
if args.rot:
  d2 = read_multi_design(args.rot, "stiff1", "stiff2", "rotAngle")
if args.threshold > 0:
  t = threshold_filter(d, args.threshold, 1e-7, 1.0)
  print 'perform threshold at ' + str(args.threshold) + ": avg density " + str(numpy.sum(d) / len(d)) + " -> " + str(numpy.sum(t) / len(d))  
  d = t
data = numpy.zeros((len(d), 3))
for i in range(len(d)): 
    # data[i,0] = 1.-numpy.sqrt(1.-d[i])
    # data[i,1] = 1.-numpy.sqrt(1.-d[i])
    data[i, 0] = d[i]
    data[i, 1] = d[i]
    data[i, 2] = 0.
    if args.rot:
      data[i, 2] = d2[i, 2]      
write_multi_design_file(args.output, data, ["stiff1", "stiff2", "rotAngle"])
