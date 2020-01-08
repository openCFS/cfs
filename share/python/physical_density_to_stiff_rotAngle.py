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
parser.add_argument("--rot_ordering", help="option for ordering of the rotation angles")
parser.add_argument("--out_dat", help="file name for table of design variables")
parser.add_argument("--dim", help="dimension of output file, default: 2 (2D)", type=int, default=2)
parser.add_argument("--input_dat", help="input text file name with matrix of design variables")
parser.add_argument("--upper_bound", help="upper_bound for lattice cells recogniztion, output: how many lattice,void,solid cells",type=float, default=0.99999999)
parser.add_argument("--lower_bound", help="lower_bound for lattice cells recogniztion, output: how many lattice,void,solid cells",type=float, default=1e-6)





args = parser.parse_args()
if args.dim ==2 and not args.input_dat and (not args.lower_bound or not args.upper_bound):
  d = read_density_as_vector(args.input, "design")
elif args.input_dat:
  d = numpy.loadtxt(args.input_dat)
  d[:,0] *= -1. 
  write_multi_design_file(args.output, d, ["rotAngle","stiff1", "stiff2"])
elif args.lower_bound or args.upper_bound:
  d = read_multi_design(args.input, "stiff1", "stiff2", "stiff3")
  lattice_count = 0
  void_count = 0
  solid_count = 0
  for i in range(len(d)):
    if d[i,0] < args.lower_bound and d[i,1] < args.lower_bound and d[i,2] < args.lower_bound:
      void_count += 1
    elif d[i,0] <= args.upper_bound and d[i,1] <= args.upper_bound and d[i,2] <= args.upper_bound: 
      lattice_count += 1
    elif d[i,0] > args.upper_bound or d[i,1] > args.upper_bound or d[i,2] > args.upper_bound: 
      solid_count +=1
    else:
      print('not counted: ' +str(d[i,0]) + ', ' +str(d[i,1]) + ', ' +str(d[i,2]) + '\n')
  print('Percent of lattice: ' + str(100.*float(lattice_count)/len(d)) + ', Percent of void: '+ str(100.*float(void_count)/len(d)) + ' Percent of solid: '+str(100.*float(solid_count)/len(d)) + '\n')
  print(' total number of elements: '+str(len(d)) + ', total number counted: '+str(lattice_count + void_count + solid_count) + '\n')  
else:
  d = read_density_as_vector(args.input, "density","physical")
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
        tmp = d2[i, 0]
        d2[i, 0] = d2[i, 1]
        d2[i, 1] = tmp
      elif swap == 2.:
        d2[i, 2] = d2[i, 2] - math.pi / 2.
        tmp = d2[i, 0]
        d2[i, 0] = d2[i, 1]
        d2[i, 1] = tmp
    write_multi_design_file(args.output, d2, ["stiff1", "stiff2", "rotAngle"])
elif args.rot_ordering:
  # TODO: NOT DONE YET
  d2 = read_multi_design(args.input, "stiff1", "stiff2", "rotAngle")
  prev = d2[0, 2]
  for i in range(len(d2)):
    rot = d2[i, 2]
    dist1 = abs(rot - prev)
    dist2 = abs(rot + math.pi / 2. - prev)
    dist3 = 9999.  # abs(rot - math.pi / 2. - prev)
    # dist4 = abs(rot + math.pi - prev)
    # dist5 = abs(rot - math.pi - prev)
    # if dist2 < dist1 and dist2 < dist3:  # and dist2 < dist4 and dist2 < dist5:
    #  tmp = d2[i, 0]
    #  d2[i, 0] = d2[i, 1]
    #  d2[i, 1] = tmp
    #  d2[i, 2] = rot + math.pi / 2.
    # elif dist3 < dist1 and dist3 < dist2:  # and dist3 < dist4 and dist3 < dist5:
    if abs(d2[i, 2] + math.pi / 2.) < 0.3:
      tmp = d2[i, 0]
      d2[i, 0] = d2[i, 1]
      d2[i, 1] = tmp
      d2[i, 2] = d2[i, 2] + math.pi / 2.
    # elif dist4 < dist1 and dist4 < dist2 and dist4 < dist3 and dist4 < dist5:
   #   d2[i, 2] = rot + math.pi
   # elif dist5 < dist1 and dist5 < dist2 and dist5 < dist3 and dist5 < dist4:
   #   d2[i, 2] = rot - math.pi
    prev = d2[i, 2]
  # step = 0
  # for i in range(len(d2)):
  #  if step + 80 < len(d2):
  #    step = step + 80
  #  else:
  #    step = (i+2) % 40
  # for i in range(len(d2)):
  #  if abs(abs(d2[i, 2]) - math.pi) < 0.01:
  #    d2[i, 2] = 0.
  write_multi_design_file(args.output, d2, ["stiff1", "stiff2", "rotAngle"])
else:
  if args.out_dat:
    d2 = read_multi_design(args.input, "stiff1", "stiff2", "rotAngle",attribute="physical")
    file = open(args.out_dat, "w")
    for i in range(len(d2)):
      file.write(str(d2[i, 0]) + '           ' + str(d2[i, 1]) + '          ' + str(d2[i, 2]) + '\n')
    file.close()
  if args.rot:
    data = numpy.zeros((len(d2), 3))
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
    if not args.input_dat and not args.lower_bound and not args.upper_bound:
      data = numpy.zeros((len(d), 3)) 
      for i in range(len(d)): 
        # data[i,0] = 1.-numpy.sqrt(1.-d[i])
        # data[i,1] = 1.-numpy.sqrt(1.-d[i])
        data[i, 0] = d[i]
        data[i, 1] = d[i]
        data[i, 2] = 0. if args.dim == 2 else d[i]
        if args.rot and args.dim != 2:
          data[i, 2] = d2[i, 2]
        if args.angle and args.dim !=2:
          data[i, 2] = args.angle     
      if args.dim == 2: 
        write_multi_design_file(args.output, data, ["stiff1", "stiff2", "rotAngle"])
      else:
        write_multi_design_file(args.output, data, ["stiff1", "stiff2", "stiff3"])
