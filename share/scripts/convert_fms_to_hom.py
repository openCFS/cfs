#!/usr/bin/env python

from optimization_tools import *
import argparse
import sys


def check(data, idx, min_v, max_v, scale):
  for i in range(len(data)):
    data[i,idx] *= scale
    data[i,idx] = max(data[i,idx], min_v)
    data[i,idx] = min(data[i,idx], max_v)

parser = argparse.ArgumentParser()
parser.add_argument("out", help="output density file w/o rotAngle2")
parser.add_argument("--density_input", help="a density file created from FMS with angle2 which will be removed")
parser.add_argument("--fms_stiff1", help="file with stiff1 data from FMS (provided my Michael)")
parser.add_argument("--fms_stiff2", help="file with stiff2 data from FMS (provided my Michael)")
parser.add_argument("--fms_angle", help="file with angle data from FMS (provided my Michael)")
parser.add_argument("--min", help="lower bound for stiff1 and stiff2 (default '0.01')", default=0.01)
parser.add_argument("--max", help="upper bound for stiff1 and stiff2 (default '1.0')", default=1.0)
parser.add_argument("--scale", help="scale stiff1 and stiff2", default=1.0)

args = parser.parse_args()

b = None

if args.density_input <> None:
  a = read_multi_design(args.density_input, "stiff1", "stiff2", "rotAngle", "rotAngle2")
  b = a[:,0:3]
else:
  assert(args.fms_stiff1 <> None and args.fms_stiff2 <> None and args.fms_angle <> None )
  # the files are plain data lists (vectors)  
  stiff1 = numpy.genfromtxt(args.fms_stiff1)
  stiff2 = numpy.genfromtxt(args.fms_stiff2)
  angle = numpy.genfromtxt(args.fms_angle)
  b = numpy.zeros((len(stiff1), 3))
  b[:,0] = stiff1
  b[:,1] = stiff2
  b[:,2] = angle
  
check(b, 0, float(args.min), float(args.max), float(args.scale))
check(b, 1, float(args.min), float(args.max), float(args.scale))

write_multi_design_file(args.out, b, "stiff1", "stiff2", "rotAngle")
