#!/usr/bin/env python

from optimization_tools import *
from hdf5_tools import *
import argparse
import sys


def check(data, idx, min_v, max_v, scale):
  for i in range(len(data)):
    data[i,idx] *= scale
    data[i,idx] = max(data[i,idx], min_v)
    data[i,idx] = min(data[i,idx], max_v)

parser = argparse.ArgumentParser()
parser.add_argument("input", help="a density file created from FMS")
parser.add_argument("out", help="output density file w/o rotAngle2")
parser.add_argument("--min", help="lower bound for stiff1 and stiff2 (default '0.01')", default=0.01)
parser.add_argument("--max", help="upper bound for stiff1 and stiff2 (default '1.0')", default=1.0)
parser.add_argument("--scale", help="scale stiff1 and stiff2", default=1.0)

args = parser.parse_args()

a = read_multi_design(args.input, "stiff1", "stiff2", "rotAngle", "rotAngle2")
b = a[:,0:3]
check(b, 0, float(args.min), float(args.max), float(args.scale))
check(b, 1, float(args.min), float(args.max), float(args.scale))

write_multi_design_file(args.out, b, "stiff1", "stiff2", "rotAngle")
