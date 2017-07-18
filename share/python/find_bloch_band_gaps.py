#!/usr/bin/env python

# the eigenfrequencies are not necessary sorted by frequency. This script does the sorting for the bloch.dat file
import numpy
import os
import sys
import argparse
from lxml import etree

gap_count = 0
offset = None

def check_gap(data, test_col, range_start, range_end, eps, gnuplot, xml):
  global gap_count
  #print "check_gap(rs=" + str(range_start) + ", re=" + str(range_end) + ")"
  range_start = int(range_start)
  range_end = int(range_end)

  mi = min(data[range_start:range_end+1,test_col])
  ma = max(data[range_start:range_end+1,test_col-1])
  rel = (mi - ma)/((ma+mi)/2.0)
  
  #print 'check_gap tm=' + str(test_col) + ' rs=' + str(range_start) + ' re=' + str(range_end) + ' -> mi=' + str(mi) + ' ma=' + str(ma) + ' rel=' + str(rel) 
    
  if (mi - ma) > eps: # rel negative if mi < ma
    gap_count += 1       
    if gnuplot:
      print('set object ' + str(gap_count) + ' rect from ' + str(range_start) + ',' + str(ma) + ' to ' + str(range_end) + ',' + str(mi))
    else:
      mytype = 'full' if range_end - range_start > data.shape[0] else 'partial'
      if xml is not None:
        node = etree.SubElement(xml, mytype)
        node.append(etree.Element("wave_vector", start=str(range_start), end=str(range_end)))
        node.append(etree.Element("frequency", start=str(ma), end=str(mi)))
        node.append(etree.Element("mode", start=str(test_col-offset), end=str(test_col-offset+1)))
        node.attrib["size"] = str(mi - ma)
        node.attrib["rel_size"] = str(rel)
        node.attrib["count"] = str(gap_count)
      else:
        print(mytype + ' band gap between ' + str(ma) + ' and ' + str(mi) + ' within ' + str(range_start) + ' -> ' + str(range_end) + ' between modes ' + str(test_col-offset) + ' and ' + str(test_col-offset+1) + ' size: ' + str(mi - ma) + ' rel.size: ' + str(rel))


# return the feader of bloch.dat or "" if none
def get_header(filename):
  header = ""
  with open(args.bloch) as f:
    for line in f:
       if line.startswith('#'):
         header += line
       else:
         break
    f.close()

  return header  
     
# tries to determine dimension. searches for k_z in comment if there is a comment. 
# if dim is given at arguement and there is a mismatch prints a warning 
def get_dim(args):
  
  header = get_header(args.bloch)
    
  has_header = header.find("k_y") > -1  
  header_dim = 3 if header.find("k_z") > -1 else 2 # check has_header!

  if args.dim is None:
    if has_header:
      if not args.gnuplot:
        print("detect dim " + str(header_dim) + " from comment in file " + args.bloch)
      return header_dim  
    else:  
      raise "no --dim given and no valid comment in '" + args.bloch + "' to derive the dimension"
  else:
    if has_header:
      if args.dim != header_dim and not args.gnuplot:
        print("--dim " + str(args.dim) + " given but header in '" + args.bloch + "' suggests dimension " + str(header_dim))
    return args.dim                                                                                                     

## determines the number of segments and segment size
# for horizrontal num = 1, for 2D it is 3 for symmetric and 4 for quadrant and in 3D one more
# uses --horizontal and the header comment (# <ibz dim="2" edges="4"/>)
#@return array of of size segments with the last wave number
def get_segments(args, data, dim): 
  num = dim+1 # default 
  if args.horizontal:
    num = 1
    if not args.gnuplot:
       print('assume horizontal segment')  
  else:
    header = get_header(args.bloch)    
    if 'edges="' in header:
      start = header.find('edges="') + 7 # add the string
      end   = header[start:].find('"')  
      num = int(header[start:start+end])
      if not args.gnuplot:
        print('read ' + str(num) + ' segments from header for ' + str(len(data)) + '-1 data lines')   
  
  if (len(data)-1) % num != 0:
    if not args.gnuplot:
      print('--numer of data lines (' + str(len(data)) + ' -1) incompatible with ' + str(num) + ' segments') 

  if num == 1:
   return [len(data)]

  size = (len(data)-1) / num
  result = [size]
  for i in range(1, num):
    result.append(result[-1]+size)  
      
  # segment = data.shape[0]/offset # number of k for G->X = X->M = M -> G
  # Y = org.shape[0] + 1 # one over the top

  return result    
  
parser = argparse.ArgumentParser()
parser.add_argument("bloch", help="a sorted bloch.dat file. Sort by sort_bloch.py")
parser.add_argument("--dim", help="2 or 3 dimensions, for old .bloch.dat files without #<ibz/>", type=int, required=False, choices=[2,3])
parser.add_argument('--mingap', help="minimal absolute (partial) band gap size (default 0.0 = all gaps)", default=0.0)
parser.add_argument('--nopartial', action='store_true', help='handle only full band gaps')
parser.add_argument('--maxmode', help="maximal mode number to be considered", default=9999, type=int)
parser.add_argument('--maxfrequency', help="maximal frequency", type=float)
parser.add_argument('--info', action='store_true', help='show range for all modes')
parser.add_argument('--xml', help='export info to a xml file')
parser.add_argument('--gnuplot', help='create gnuplot output, specify the type', choices = ['eps', 'png', 'console'])
parser.add_argument('--nolines', action='store_true', help='gnuplot: do not concatenate points by lines')
parser.add_argument('--commonsymbol', action='store_true', help='gnuplot: use the same line symbol for all lines')
parser.add_argument('--paper', action='store_true', help="tune for paper publishing (e.g. gray)")
parser.add_argument('--fontsize', type=int, help="gnuplot font size")
parser.add_argument('--title', action='store_true', help="gnuplot: add title, off by default")
parser.add_argument('--nicelabel', action='store_true', help='gnuplot: use nice labels')
parser.add_argument('--horizontal', action='store_true', help='display only the horizontal part of the wavevector (G->X)')
args = parser.parse_args()

if not os.path.exists(args.bloch):
  print("file not found '" + args.bloch + "'")
  sys.exit(1)
org = numpy.loadtxt(args.bloch)
 
dim = get_dim(args)

offset = 3 if dim == 2 else 4 # step, k_x, k_y (,k_z)

# check for unsorted
# we may not simpy sort here as the plot is done via the original data which needs also the sorted and saved to another file
for w in range(len(org)):
  wv = org[w]
  for ev in range(offset+1, len(wv)):
    if wv[ev] < wv[ev-1]:
      print('#unsorted eigenfrequency: wave_vector ' + str(w) + ' index ' + (str(ev)) + ' value ' + str(wv[ev]) + ' < ' + str(wv[ev-1]))                 

segments = get_segments(args, org, dim)

eps = float(args.mingap) 
max_mode = min((args.maxmode + offset + 1,org.shape[1]))
max_freq =  numpy.amax(org[:,offset:])
cnt = 1

# the info.xml root node
root = None if not args.xml else etree.Element('bloch', file=args.bloch)

if args.info and not args.gnuplot: 
  print("\nmode range:")
  print("-----------")
  for i in range(offset, max_mode):
    mi = min(org[:,i])
    ma = max(org[:,i])
    print('mode ' + str(i-offset+1) + ' min:' + str(mi) + ' max:' + str(ma) + ' size:' + str(ma - mi) + ' rel.size: ' + str((ma - mi)/((ma+mi)/2.0))) 

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
  fontsize = str(args.fontsize) if args.fontsize else "18"
  
  if args.gnuplot == "eps":
    print('set size ratio 1.0')
    print('set terminal postscript eps enhanced "Helvetica, ' + fontsize + '" color') #  change to monochrome for papers
    print('set output "' + args.bloch[:-len(".dat")] + '.eps"')
  if args.gnuplot == "png":
    print('set size ratio 1.0')
    print('set terminal png size 1000,1000 font "Helvetica, ' + fontsize + '"')
    print('set output "' + args.bloch[:-len(".dat")] + '.png"') # leave it as bloch.png
        
    # print 'set output "tmp.eps"'
  # unset eventually older boxes
  for i in range(60):
     print('unset object ' + str(i))
  print('set style rectangle back fc rgb "gray"')         

  
# search for partial band gaps
if not args.nopartial:
  if not args.gnuplot and not args.xml:
    print("\npartial band gaps >= rel.size " + str(eps))
    print("---------------------------------------")
  for s in range(0,len(segments)):
    for i in range(offset+1, max_mode):  
      check_gap(org, i, 0 if s == 0 else segments[s-1], segments[s], eps, args.gnuplot, gaps)
  
# search for full band gaps
if not args.horizontal:
  if not args.gnuplot and not args.xml:
    print("\nfull band gaps >= rel.size " + str(eps))
    print("---------------------------------------")
  for i in range(offset+1,  max_mode): # we start comparing the second to the first
    check_gap(org, i, 0, segments[-1], eps, args.gnuplot, gaps)

if args.gnuplot:
  if args.maxfrequency:
    print('set yrange [0:' + str(args.maxfrequency) + ']')
  else:  
    if args.commonsymbol:
      print('set yrange [0:*]')
    else:
      print('set yrange [0:' + str(max_freq * 1.2) + ']') # leave space for the labels
  if args.horizontal:
    print('set xrange [0:' + str(segments[-1]) + ']')    
  else:
    for s in range(0,len(segments)-1):
      print('set arrow ' + str(s+1) + '  from ' + str(segments[s]) + ',0 to ' + str(segments[s]) + ',' + str(max_freq) + ' nohead lt rgb "gray" lw 2')  
  
  if args.nicelabel:
     print('set ylabel "eigenfrequency in Hz" offset 1')
     print('set xlabel "wave vector (' + ('horizontal ' if args.horizontal else '') + 'IBZ)"')
     symbols = ['"{/Symbol G}"', '"X"', '"{/Symbol M}"', '"R"']
     xtics = 'set xtics (' + symbols[0] + ' 0'
     for i in range(len(segments)):
        xtics += ', ' + symbols[i+1 if i < len(segments)-1 else 0] + ' ' + str(segments[i])
     print(xtics + ')')    
  else:
    print('unset ylabel') 
    print('unset xlabel')     

  wl =   '' if args.nolines else ' with linespoints lw 2 ' if not args.commonsymbol else ' with lines lw 2' 
  lc = ' lc black ' if args.paper or args.commonsymbol else ''
  for i in range(offset,  max_mode): # 1-based
    title = ' notitle ' if args.commonsymbol or not args.title else ' t "' + str(i-offset+1) + '. mode" ' 
    print(('plot' if i <= offset else '    ') + '"' + args.bloch + '" u ' + str(i+1) + title + wl + lc + (' ,\\' if i < max_mode -1  else '')) 
 
 
if args.xml:
  root.find("gaps").attrib["count"]=str(gap_count)
  file = open(args.xml, "w")
  file.write(etree.tostring(root, pretty_print=True))
  file.close()
          
