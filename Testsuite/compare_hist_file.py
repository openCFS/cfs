#!/usr/bin/env python
import numpy as np
import os.path
import argparse 

#=======================================================
# Script to compare text files, such as cfs .hist files
#=======================================================
# define arguments
parser = argparse.ArgumentParser(description = "Compares text files (such as cfs .hist files).")
parser.add_argument("test", help="the file to compare against the reference")
parser.add_argument("reference", help="the reference hist / text file")
parser.add_argument("--eps", help="the tolerance of each entry", type=float, default=1e-10)
args = parser.parse_args()

# print info
print("compare_hist_files: \nref='" + args.reference + "'\ntest='" + args.test + "'\n")

# check if files exist
if not os.path.isfile(args.reference):
  print("error: reference file '" + args.reference + "' does not exist")
  os.sys.exit(1)
if not os.path.isfile(args.test):
  print("error: test file '" + args.test + "' does not exist")
  os.sys.exit(1)

# read in files, using '#' as comment symbol and ',' as a column separator.
ref = np.loadtxt(args.reference, comments='#', delimiter=',')
print(args.reference)
print("--------------------------------------------------------------------------")
print(ref)
print("--------------------------------------------------------------------------")
tst = np.loadtxt(args.test, comments='#', delimiter=',')
print(args.test)
print("--------------------------------------------------------------------------")
print(tst)
print("--------------------------------------------------------------------------")

# compute relative errors of all entries
cmp = np.abs(np.asarray(ref).astype(float) - np.asarray(tst).astype(float))

# no entry is allowed to exceed the tolerance
if np.all(cmp <= args.eps):
  print('++ all differences are below tol ' + str(args.eps) + ' ++')
  os.sys.exit(0)
else:
  print('** error: at least one difference exceeds the tolerance ' + str(args.eps) + ' **')
  os.sys.exit(1) 