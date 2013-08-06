#!/usr/bin/env python
from optimization_tools import *
import argparse
import numpy

parser = argparse.ArgumentParser()
parser.add_argument("input", help="density.xml file where pyhsical is read")
parser.add_argument("output", help="density.xml file where pyhsical is written as stiff1 and stiff2")
parser.add_argument("--threshold", help="if within (0,1) input data is thresholded", type=float, default=-1.0)
parser.add_argument("--angle",help="if within 0 and 2*pi a rotAngle is added",type=float,default=-1.)
parser.add_argument("--attribute",help="read design or physical",default='physical')
parser.add_argument("--vts",help="vts_start conversion",default=False)




args = parser.parse_args()

d  = read_density_as_vector(args.input, args.attribute)
nr = read_density_as_vector(args.input, "nr")
if args.threshold > 0.:
  t = threshold_filter(d, args.threshold, 1e-6, 1.0)
  print 'perform threshold at ' + str(args.threshold) + ": avg density " + str(numpy.sum(d)/len(d)) + " -> " + str(numpy.sum(t)/len(d))  
  d = t
if args.angle >= 0.:
    data = numpy.zeros((len(d), 3))
    for i in range(len(d)):
      if args.vts:
        data[i,0] = 1.-numpy.sqrt(1.-d[i])
        data[i,1] = 1.-numpy.sqrt(1.-d[i])
        data[i,2] = args.angle
      else:
        data[i,0] = d[i]
        data[i,1] = d[i]
        data[i,2] = args.angle
    write_multi_design_file(args.output, data, ["stiff1", "stiff2","rotAngle"], nr)
else:
    data = numpy.zeros((len(d), 2))
    if args.vts:
      for i in range(len(d)):
        data[i,0] = 1.-numpy.sqrt(1.-d[i])
        data[i,1] = 1.-numpy.sqrt(1.-d[i])
    else:
      data[:,0] = d
      data[:,1] = d
    write_multi_design_file(args.output, data, ["stiff1", "stiff2"], nr)