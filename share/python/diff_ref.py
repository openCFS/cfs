#!/usr/bin/env python
import matplotlib
#matplotlib.use('tkagg')
from matplotlib import pyplot as plt
import argparse
import os
import sys
import numpy as np

# this tool compares the objective value from a .plot.dat file with the ref.plot dat file for a test case
parser = argparse.ArgumentParser(description = 'Compare columns of a .plot.dat file with another (reference) file')
parser.add_argument("first", help=".plot.dat file, default is <directory name>.plot.dat", nargs='?')
parser.add_argument("second", help="second .plot.dat file, default is <first> but .ref.plot.dat", nargs='?')
parser.add_argument('--col', help="the column to display (1-based)", type=int, default=2)
parser.add_argument('--marker', help="style for current data line", choices=['o'], default='o')
args = parser.parse_args()

first = args.first
if not first:
  pwd = os.getcwd()
  first = os.path.basename(pwd) + '.plot.dat'
else:
  if not ".plot.dat" in first:
    print("the first parameter needs to be .plot.dat filename")
    os.sys.exit()  
if not os.path.exists(first):
  print("cannot open '" + first + "'")
  os.sys.exit()

second = args.second
if not second:
  second = first.replace(".plot.dat", ".ref.plot.dat")
if not os.path.exists(second):
  print("cannot open '" + second + "'")
  os.sys.exit()

fd = np.loadtxt(first)
sd = np.loadtxt(second)

if fd.shape[1] != sd.shape[1]:
  print("warning, the number of columns in '" + first  + "' (" + str(fd.shape[1]) + ") differs from '" + second + "' (" + str(sd.shape[1]) + ")") 

col = args.col-1

if col > fd.shape[1] and col > sd.shape[1]:
  print("--col is chosen too large")
  os.sys.exit()

plt.plot(sd[:,0], sd[:,col], marker=args.marker, label=second)
plt.plot(fd[:,0], fd[:,col], label=first)

loc = 'upper right' if fd[-1,col] < sd[0,col] else 'lower right' 

plt.legend(shadow=False,prop={'size':12})
print('close window to exit')
plt.show() 
