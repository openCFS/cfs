#!/usr/bin/env python

# the eigenfrequencies are not necessary sorted by frequency. This script does the sorting for the bloch.dat file
import numpy
import os
import sys
import argparse
from random import choice

gap_count = 1
offset = None


def check_gap(data, test_col, range_start, range_end, eps, gnuplot):
  global gap_count
  
  mi = min(data[range_start:range_end+1,test_col])
  ma = max(data[range_start:range_end+1,test_col-1])
  rel = (mi - ma)/((ma+mi)/2.0)
  
  #print 'check_gap tm=' + str(test_col) + ' rs=' + str(range_start) + ' re=' + str(range_end) + ' -> mi=' + str(mi) + ' ma=' + str(ma) + ' rel=' + str(rel) 
    
  if (mi - ma) > eps: # rel negative if mi < ma
    if gnuplot:
      print 'set object ' + str(gap_count) + ' rect from ' + str(range_start) + ',' + str(ma) + ' to ' + str(range_end) + ',' + str(mi)
      gap_count += 1       
    else:
      type = 'full' if range_end - range_start > data.shape[0] else 'partial'
      print type + ' band gap between ' + str(ma) + ' and ' + str(mi) + ' within ' + str(range_start) + ' -> ' + str(range_end) + ' between modes ' + str(test_col-offset) + ' and ' + str(test_col-offset+1),
      print ' size: ' + str(mi - ma) + ' rel.size: ' + str(rel)

parser = argparse.ArgumentParser()
parser.add_argument("bloch", help="a sorted bloch.dat file. You might use sort_bloch.py first")
parser.add_argument("--dim", help="2 or three dimensions", type=int, required=True, choices=[2,3])
parser.add_argument('--mingap', help="minimal absolute (partial) band gap size (default 0.0 = all gaps)", default=0.0)
parser.add_argument('--nopartial', action='store_true', help='handle only full band gaps')
parser.add_argument('--maxmode', help="maximal mode number to be considered (default 9999)", default=9999, type=int)
parser.add_argument('--info', action='store_true', help='show range for all modes')
parser.add_argument('--gnuplot', action='store_true', help='create gnuplot output')
parser.add_argument('--nolines', action='store_true', help='gnuplot: do not concatenate points by lines')
parser.add_argument('--commonsymbol', action='store_true', help='gnuplot: use the same line symbol for all lines')
parser.add_argument('--nicelabel', action='store_true', help='gnuplot: use nice labels')
args = parser.parse_args()

org = numpy.loadtxt(args.bloch)
assert(org[-1][0] == 0.0) # the last step shall be a repetition of the first step
 
dim = args.dim 

offset = 3 if dim == 2 else 4 # step, k_x, k_y (,k_z)

# we search also for partial band gaps. The points are X=0, M, G, (R for 3D), Y=X=0 (again)
M = (org.shape[0]-1) / (dim+1) # -1 because the real data is %3 and %4 but the first line is repeated
G = 2 * M
R = 3 * M # ignored in 2D
Y = org.shape[0] + 1 # one over the top

eps = float(args.mingap) 
max_mode = min((args.maxmode + offset + 1,org.shape[1]))
max_freq =  numpy.amax(org[:,offset:])
cnt = 1

if args.info and not args.gnuplot: 
  print "\nmode range:"
  print "-----------"
  for i in range(offset, max_mode):
    mi = min(org[:,i])
    ma = max(org[:,i])
    print 'mode ' + str(i-offset+1) + ' min:' + str(mi) + ' max:' + str(ma), # print step 1-based
    print ' size:' + str(ma - mi) + ' rel.size: ' + str((ma - mi)/((ma+mi)/2.0)) 


if args.gnuplot:
  if args.commonsymbol:
    print 'set yrange [0:*]'
  else:
    print 'set yrange [0:' + str(max_freq * 1.2) + ']' # leave space for the labels
  for d in range(1, dim+1):
    print 'set arrow ' + str(d) + '  from ' + str(d * org.shape[0]/offset) + ',0 to ' + str(d * org.shape[0]/offset) + ',' + str(max_freq) + ' nohead lt rgb "gray" lw 2'  
  
  if args.nicelabel:
     print 'set ylabel "eigenfrequency in Hz"'
     print 'set xlabel "wave vector (IBZ)"'
     print 'set xtics ("G" 0, "X" ' + str(org.shape[0]/offset),
     if dim == 2:
       print ', "M" ' + str(2 * org.shape[0]/offset) + ', "G" ' + str(org.shape[0]-1) + ')'
     else:
       print ', "M" ' + str(2 * org.shape[0]/offset) + ', "R" ' + str(3 * org.shape[0]/offset) + ', "G" ' + str(org.shape[0]-1) + ')'
  else:
    print 'unset ylabel' 
    print 'unset xlabel'     

  
  print 'set style rectangle back fc rgb "gray"'
  wl =   '' if args.nolines else ' with linespoints '
  lc = ' lc 7 lt 1 ' if args.commonsymbol else ''
  for i in range(offset,  max_mode): # 1-based
    title = ' notitle ' if args.commonsymbol else ' t "' + str(i-offset+1) + '. mode" ' 
    print ('plot' if i <= offset else '    ') + '"' + args.bloch + '" u ' + str(i+1) + title + wl + lc + (' ,\\' if i < max_mode -1  else '')
  print 'replot' # necessary to show the boxes   
  
# search for partial band gaps
if not args.nopartial:
  if not args.gnuplot:
    print "\npartial band gaps >= rel.size " + str(eps)
    print "---------------------------------------"
  for i in range(offset+1, max_mode):
    check_gap(org, i, 0, M, eps, args.gnuplot)
    check_gap(org, i, M, G, eps, args.gnuplot)
    if dim == 2:
      check_gap(org, i, G, Y, eps, args.gnuplot)
    else:
      check_gap(org, i, G, R, eps, args.gnuplot)
      check_gap(org, i, R, Y, eps, args.gnuplot)
  
# search for full band gaps
if not args.gnuplot:
  print "\nfull band gaps >= rel.size " + str(eps)
  print "---------------------------------------"
for i in range(offset+1,  max_mode): # we start comparing the second to the first
  check_gap(org, i, 0, Y, eps, args.gnuplot)
 
      
