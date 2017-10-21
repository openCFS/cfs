#!/usr/bin/env python
from optimization_tools import *
import argparse
import numpy

parser = argparse.ArgumentParser()
parser.add_argument("input", help="density.xml file where pyhsical is read")
parser.add_argument("output", help="density.xml file where pyhsical is written as stiff1 and stiff2 and optional rotAngle")
parser.add_argument("--threshold", help="if within (0,1) input data is thresholded", type=float, default=-1.0)
parser.add_argument("--angle", help="if within 0 and 2*pi a rotAngle is added", type=float, default=-1.)
parser.add_argument("--attribute", help="read design or physical", default='physical')
parser.add_argument("--vts", help="vts_start conversion", default=False)
parser.add_argument("--hom_stiff", help="hom_start with stiff1 and stiff2", default=False)
parser.add_argument("--stiff1", help="dat file with values for stiff1")
parser.add_argument("--stiff2", help="dat file with values for stiff2")
parser.add_argument("--rotAngle", help="dat file with values for rotAngle")
parser.add_argument("--density", help="normal density file is written")


args = parser.parse_args()
if args.stiff1 and args.stiff2 and args.rotAngle:
  stiff1 = numpy.fromfile(args.stiff1, dtype=float, count=-1, sep='\n')
  stiff2 = numpy.fromfile(args.stiff2, dtype=float, count=-1, sep='\n')
  rotAngle = numpy.fromfile(args.rotAngle, dtype=float, count=-1, sep='\n')
  data = numpy.zeros((len(stiff1), 3))
  nr = numpy.zeros((len(stiff1), 1))
  for i in range(len(stiff1)):
    data[i, 0] = stiff1[i]
    data[i, 1] = stiff2[i]
    data[i, 2] = rotAngle[i]
    nr[i, 0] = i + 1
    write_multi_design_file(args.output, data, ["stiff1", "stiff2", "rotAngle"], nr)	
elif args.stiff1 and args.stiff2:
  stiff1 = numpy.fromfile(args.stiff1, dtype=float, count=-1, sep='\n')
  stiff2 = numpy.fromfile(args.stiff2, dtype=float, count=-1, sep='\n')
  data = numpy.zeros((len(stiff1), 2))
  nr = numpy.zeros((len(stiff1), 1))
  for i in range(len(stiff1)):
    data[i, 0] = stiff1[i]
    data[i, 1] = stiff2[i]
    nr[i, 0] = i + 1
  write_multi_design_file(args.output, data, ["stiff1", "stiff2"], nr)
else:
  if args.hom_stiff:
    d = read_multi_design(args.input, "stiff1", "stiff2")
    nr = numpy.zeros((len(d), 1))
  else:
    d = read_density_as_vector(args.input, args.attribute)
    nr = read_density_as_vector(args.input, "nr")
  if args.threshold > 0.:
    t = threshold_filter(d, args.threshold, 1e-6, 1.0)
    print 'perform threshold at ' + str(args.threshold) + ": avg density " + str(numpy.sum(d) / len(d)) + " -> " + str(numpy.sum(t) / len(d))  
    d = t
  if args.angle >= 0.:
    data = numpy.zeros((len(d), 3))
    count = 1
    for i in range(len(d)):
      if args.vts:
        data[i, 0] = 1. - numpy.sqrt(1. - d[i])
        data[i, 1] = 1. - numpy.sqrt(1. - d[i])
        data[i, 2] = args.angle
      elif args.hom_stiff:
        nr[i, 0] = count
        count = count + 1
        data[i, 0] = d[i, 0]
        data[i, 1] = d[i, 1]
        data[i, 2] = args.angle
      else:
        data[i, 0] = d[i]
        data[i, 1] = d[i]
        data[i, 2] = args.angle
    write_multi_design_file(args.output, data, ["stiff1", "stiff2", "rotAngle"], nr)
  else:
    if args.density:
      data = numpy.zeros((len(d), 1))
      ndata = numpy.zeros((len(nr), 1))
    else:
      data = numpy.zeros((len(d), 2))
    if args.vts:
      for i in range(len(d)):
        data[i, 0] = 1. - numpy.sqrt(1. - d[i])
        data[i, 1] = 1. - numpy.sqrt(1. - d[i])
    if args.density:
      for i in range(len(d)):
        data[i,0] = d[i]
        ndata[i,0] = nr[i]
      write_density_file(args.output, data,elemnr=ndata)
    else:
      write_multi_design_file(args.output, data, ["stiff1", "stiff2"], nr)
