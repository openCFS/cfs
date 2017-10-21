#!/usr/bin/env python

import numpy
import os
import sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('input', help=".plot.dat with non-monotonous values from snopt line-search")
parser.add_argument('--column', help="column with the objective value, 0-based (default=1)", default=1, type=int)
parser.add_argument('--change_iterations', metavar='True', help="number iterations consecutively starting with 0", default=True)
parser.add_argument('--percent_away', metavar='percentage', help="find iteration number 'percentage' away from last iteration value", type=float)
parser.add_argument('--iteration_fraction', metavar='fraction', help="return percentage away from optimum after 'fraction' part of iterations", type=float)
parser.add_argument('--design_change_repetitions', nargs='+', metavar='n m', help="find lowest iteration number for which the design change was at max n for m sequential iterations")
parser.add_argument('--relative_change_repetitions', metavar='n m', nargs='+', help="find lowest iteration number for which the relative function value changed at max n for m sequential iterations")
parser.add_argument('--silent', metavar='False', help="print file or not", default=False, type=bool)
args = parser.parse_args()

if not os.path.exists(args.input):
  print "cannot open '" + args.input + "'"
  sys.exit() 

org = numpy.loadtxt(args.input)
 
# create the output of the new file
if not args.silent: 
  print "#this is the file '" + args.input + "'" + "with monotously decreasing objective values"
# copy the comment lines form the old file
f = open(args.input)
raw = f.readlines()
f.close()

i = 0
while raw[i].startswith('#'):
  if not args.silent:
    print raw[i],
  i += 1

# now the content
cval = org[0][args.column]
count = 0
for r in range(org.shape[0]):
  for c in range(org.shape[1]):
    if org[r][args.column] <= cval:
      cval = org[r][args.column]
      if c == 0:
	if args.change_iterations:
	  if not args.silent:
	    print str(int(count)) + ' \t',
	  count += 1
	else:
	  if not args.silent:
	    print str(int(org[r][c])) + ' \t',
      else:
	if not args.silent:
	  print str(org[r][c]) + ' \t',
      if c == (org.shape[1]-1):
	if not args.silent:
	  print '\n',
	  
if args.percent_away is not None:
  last = org[-1][args.column]
  for r in range(org.shape[0]):
    if org[r][args.column]/last <= (1.0 + .01*args.percent_away):
      print "Iteration " + str(r+1) + " of " + str(org.shape[0]) + " (" + str((float(r+1))/org.shape[0]) + ") is " + str(args.percent_away) + "% from optimum!"
      break
    
if args.iteration_fraction is not None:
  it = int(round((org.shape[0]-1)*args.iteration_fraction))
  print "After " + str(100*args.iteration_fraction) + "% of iterations, iteration " + str(it+1) + " is " + str((org[it,args.column]/org[-1,args.column]-1.00)*100) + "% from optimum!"

      
if args.design_change_repetitions is not None:
  for r in range(org.shape[0]):
    if r >= int(args.design_change_repetitions[1]):
      mm = max(org[r-int(args.design_change_repetitions[1]):r+1,args.column+1])
      if mm < float(args.design_change_repetitions[0]):
        print "Iteration " + str(r+1) + " of " + str(org.shape[0]) + " had maximum " + str(mm) + " design change for " + args.design_change_repetitions[1] + " consecutive iterations!"
        break
        
if args.relative_change_repetitions is not None:
  bb = False
  for r in range(org.shape[0]):
    if r > int(args.relative_change_repetitions[1]):
      for rr in range(int(args.relative_change_repetitions[1])):
        if org[r-rr-1][args.column]/org[r-rr][args.column]-1.0 > float(args.relative_change_repetitions[0]):
	  break
	else:
	  if rr == int(args.relative_change_repetitions[1])-1:
	    bb = True
    if bb:
      print "Iteration " + str(r+1) + " of " + str(org.shape[0]) + " had maximum " + args.relative_change_repetitions[0] + " relative change for " + args.relative_change_repetitions[1] + " consecutive iterations!"
      break
