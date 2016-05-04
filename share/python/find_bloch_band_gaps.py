#!/usr/bin/env python

# the eigenfrequencies are not necessary sorted by frequency. This script does the sorting for the bloch.dat file
import numpy
import os
import sys
import argparse
from random import choice
from lxml import etree

gap_count = 0
offset = None

def check_gap(data, test_col, range_start, range_end, eps, gnuplot, xml):
  global gap_count
  #print "check_gap(rs=" + str(range_start) + ", re=" + str(range_end) + ")"
  
  mi = min(data[range_start:range_end+1,test_col])
  ma = max(data[range_start:range_end+1,test_col-1])
  rel = (mi - ma)/((ma+mi)/2.0)
  
  #print 'check_gap tm=' + str(test_col) + ' rs=' + str(range_start) + ' re=' + str(range_end) + ' -> mi=' + str(mi) + ' ma=' + str(ma) + ' rel=' + str(rel) 
    
  if (mi - ma) > eps: # rel negative if mi < ma
    gap_count += 1       
    if gnuplot:
      print 'set object ' + str(gap_count) + ' rect from ' + str(range_start) + ',' + str(ma) + ' to ' + str(range_end) + ',' + str(mi)
    else:
      type = 'full' if range_end - range_start > data.shape[0] else 'partial'
      if xml is not None:
        node = etree.SubElement(xml, type)
        node.append(etree.Element("wave_vector", start=str(range_start), end=str(range_end)))
        node.append(etree.Element("frequency", start=str(ma), end=str(mi)))
        node.append(etree.Element("mode", start=str(test_col-offset), end=str(test_col-offset+1)))
        node.attrib["size"] = str(mi - ma)
        node.attrib["rel_size"] = str(rel)
        node.attrib["count"] = str(gap_count)
      else:
        print type + ' band gap between ' + str(ma) + ' and ' + str(mi) + ' within ' + str(range_start) + ' -> ' + str(range_end) + ' between modes ' + str(test_col-offset) + ' and ' + str(test_col-offset+1),
        print ' size: ' + str(mi - ma) + ' rel.size: ' + str(rel)

# tries to determine dimension. searches for k_z in comment if there is a comment. 
# if dim is given at arguement and there is a mismatch prints a warning 
def get_dim(args):
  
  header = ""
  with open(args.bloch) as f:
    for line in f:
       if line.startswith('#'):
         header += line
       else:
         break
  
  has_header = header.find("k_y") > -1  
  header_dim = 3 if header.find("k_z") > -1 else 2 # check has_header!

  if args.dim is None:
    if has_header:
      if not args.gnuplot:
        print "detect dim " + str(header_dim) + " from comment in file " + args.bloch
      return header_dim  
    else:  
      raise "no --dim given and no valid comment in '" + args.bloch + "' to derive the dimension"
  else:
    if has_header:
      if args.dim <> header_dim and not args.gnuplot:
        print "--dim " + str(args.dim) + " given but header in '" + args.bloch + "' suggests dimension " + str(header_dim)
    return args.dim                                                                                                     
  
parser = argparse.ArgumentParser()
parser.add_argument("bloch", help="a sorted bloch.dat file. Sort by sort_bloch.py")
parser.add_argument("--dim", help="2 or three dimensions", type=int, required=False, choices=[2,3])
parser.add_argument('--mingap', help="minimal absolute (partial) band gap size (default 0.0 = all gaps)", default=0.0)
parser.add_argument('--nopartial', action='store_true', help='handle only full band gaps')
parser.add_argument('--maxmode', help="maximal mode number to be considered (default 9999)", default=9999, type=int)
parser.add_argument('--info', action='store_true', help='show range for all modes')
parser.add_argument('--xml', help='export info to a xml file')
parser.add_argument('--gnuplot', help='create gnuplot output, specify the type', choices = ['eps', 'png', 'console'])
#parser.add_argument('--show',  help='use gnuplot to create a png and shows it', action='store_true',)
parser.add_argument('--nolines', action='store_true', help='gnuplot: do not concatenate points by lines')
parser.add_argument('--commonsymbol', action='store_true', help='gnuplot: use the same line symbol for all lines')
parser.add_argument('--nicelabel', action='store_true', help='gnuplot: use nice labels')
parser.add_argument('--horizontal', action='store_true', help='display only the horizontal part of the wavevector (G->X)')
args = parser.parse_args()

org = numpy.loadtxt(args.bloch)
 
dim = get_dim(args)

offset = 3 if dim == 2 else 4 # step, k_x, k_y (,k_z)

# check for unsorted
# we may not simpy sort here as the plot is done via the original data which needs also the sorted and saved to another file
for w in range(len(org)):
  wv = org[w]
  for ev in range(offset+1, len(wv)):
    if wv[ev] < wv[ev-1]:
      print '#unsorted eigenfrequency: wave_vector ' + str(w) + ' index ' + (str(ev)) + ' value ' + str(wv[ev]) + ' < ' + str(wv[ev-1])                 


segment = org.shape[0]/offset # number of k for G->X = X->M = M -> G


# we search also for partial band gaps. The points are X=0, M, G, (R for 3D), Y=X=0 (again)
M = (org.shape[0]-1) / (dim+1) # -1 because the real data is %3 and %4 but the first line is repeated
if args.horizontal:
  # if only horizontal wave vectors were calculated, we need to correct segment
  non_hor = numpy.where(org[:,2] > 0.0)[0].shape[0] # number of rows where k_y is > 0.0
  if non_hor == 0:
    segment = org.shape[0] # don't divide   
    M = segment

G = 2 * M
R = 3 * M # ignored in 2D
Y = org.shape[0] + 1 # one over the top

# print "M=" + str(M) + " G=" + str(G) + " R=" + str(G) + " Y=" + str(Y) + " shape=" + str(org.shape[0])

eps = float(args.mingap) 
max_mode = min((args.maxmode + offset + 1,org.shape[1]))
max_freq =  numpy.amax(org[:,offset:])
cnt = 1

# the info.xml root node
root = None if not args.xml else etree.Element('bloch', file=args.bloch)

if args.info and not args.gnuplot: 
  print "\nmode range:"
  print "-----------"
  for i in range(offset, max_mode):
    mi = min(org[:,i])
    ma = max(org[:,i])
    print 'mode ' + str(i-offset+1) + ' min:' + str(mi) + ' max:' + str(ma), # print step 1-based
    print ' size:' + str(ma - mi) + ' rel.size: ' + str((ma - mi)/((ma+mi)/2.0)) 

if args.xml:
  modes = etree.SubElement(root, "modes")
  for i in range(offset, max_mode):
    mi = min(org[:,i])
    ma = max(org[:,i])
    mode = etree.SubElement(modes, "mode")
    mode.attrib["nr"] = str(i-offset+1)
    mode.attrib["min"] = str(mi)
    mode.attrib["max"] = str(ma)
    mode.attrib["size"] = str(ma - mi)
    mode.attrib["rel_size"] = str((ma - mi)/((ma+mi)/2.0))

gaps = None if args.xml is None else etree.SubElement(root, "gaps")  
  
if args.gnuplot:
  if args.gnuplot == "eps":
    print 'set size ratio 1.0'
    print 'set terminal postscript eps enhanced "Helvetica, 20" color' #  change to monochrome for papers
    print 'set output "' + args.bloch[:-len(".bloch.dat")] + '.eps"'
  if args.gnuplot == "png":
    print 'set size ratio 1.0'
    print 'set terminal png font Helvetica size 1000,1000'
    print 'set output "' + args.bloch[:-len(".dat")] + '.png"' # leave it as bloch.png
        
    # print 'set output "tmp.eps"'
  # unset eventually older boxes
  for i in range(60):
     print 'unset object ' + str(i)
  print 'set style rectangle back fc rgb "gray"'         

  
# search for partial band gaps
if not args.nopartial:
  if not args.gnuplot and not args.xml:
    print "\npartial band gaps >= rel.size " + str(eps)
    print "---------------------------------------"
  for i in range(offset+1, max_mode):
    check_gap(org, i, 0, M, eps, args.gnuplot, gaps)
    if not args.horizontal:
      check_gap(org, i, M, G, eps, args.gnuplot, gaps)
      if dim == 2:
        check_gap(org, i, G, Y, eps, args.gnuplot, gaps)
      else:
        check_gap(org, i, G, R, eps, args.gnuplot, gaps)
        check_gap(org, i, R, Y, eps, args.gnuplot, gaps)
  
# search for full band gaps
if not args.horizontal:
  if not args.gnuplot and not args.xml:
    print "\nfull band gaps >= rel.size " + str(eps)
    print "---------------------------------------"
  for i in range(offset+1,  max_mode): # we start comparing the second to the first
    check_gap(org, i, 0, Y, eps, args.gnuplot, gaps)

if args.gnuplot:
  if args.commonsymbol:
    print 'set yrange [0:*]'
  else:
    print 'set yrange [0:' + str(max_freq * 1.2) + ']' # leave space for the labels
  if args.horizontal:
    print 'set xrange [0:' + str(segment) + ']'    
  else:
    for d in range(1, dim+1):
      print 'set arrow ' + str(d) + '  from ' + str(d * segment) + ',0 to ' + str(d * segment) + ',' + str(max_freq) + ' nohead lt rgb "gray" lw 2'  
  
  if args.nicelabel:
     print 'set ylabel "eigenfrequency in Hz"'
     print 'set xlabel "wave vector (' + ('horizontal ' if args.horizontal else '') + 'IBZ)"'
     print 'set xtics ("G" 0, "X" ' + str(segment),
     if dim == 2:
       print ', "M" ' + str(2 * segment) + ', "G" ' + str(org.shape[0]-1) + ')'
     else:
       print ', "M" ' + str(2 * segment) + ', "R" ' + str(3 * segment) + ', "G" ' + str(org.shape[0]-1) + ')'
  else:
    print 'unset ylabel' 
    print 'unset xlabel'     

  

  wl =   '' if args.nolines else ' with linespoints '
  lc = ' lc 7 lt 1 ' if args.commonsymbol else ''
  for i in range(offset,  max_mode): # 1-based
    title = ' notitle ' if args.commonsymbol else ' t "' + str(i-offset+1) + '. mode" ' 
    print ('plot' if i <= offset else '    ') + '"' + args.bloch + '" u ' + str(i+1) + title + wl + lc + (' ,\\' if i < max_mode -1  else '')
 
 
if args.xml:
  root.find("gaps").attrib["count"]=str(gap_count)
  file = open(args.xml, "w")
  file.write(etree.tostring(root, pretty_print=True))
  file.close()
          
