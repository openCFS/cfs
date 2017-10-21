#!/usr/bin/env python

# the eigenfrequencies are not necessary sorted by frequency. This script does the sorting for the bloch.dat file
import numpy
import os
import sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--dim", help="2 or three dimensions", type=int, required=True, choices=[2,3])
parser.add_argument('input', help="unsorted bloc.dat file")
args = parser.parse_args()

if not os.path.exists(args.input):
  print("cannot open '" + args.input + "'")
  sys.exit() 

org = numpy.loadtxt(args.input)
 
dim = args.dim 

offset = 3 if dim == 2 else 4 # step, k_x, k_y (,k_z)
    
# sort the stuff, where we do not sort number and x and y (and z)
for i in range(len(org)):
  org[i][offset:].sort()

# create the output of the new file
print("#this is the sorted file '" + args.input + "'")
# copy the comment lines form the old file
f = open(args.input)
raw = f.readlines()
f.close()

i = 0
while raw[i].startswith('#'):
  print(raw[i], end=' ')
  i += 1

# now the content
for r in range(org.shape[0]):
  for c in range(org.shape[1]):
    print(str(org[r][c]) + ' \t', end=' ')
  print('\n', end=' ')  
    