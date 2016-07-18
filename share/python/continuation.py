#!/usr/bin/env python
import sys
import libxml2
import argparse

from cfs_utils import *

## performs heaviside continuation by doubling filter/density/beta, starting from 1

def continuation(initial, old_beta, beta, mesh, short_problem, executable, show, failsafe = False):
    
  beta_problem = short_problem + "-beta_" + str(beta)  
    
  # first without iniital  
  start = "" if old_beta < 1 else "-x " + short_problem + "-beta_" + str(old_beta) + ".density.xml"
  # now check for initial
  if old_beta < 1 and initial <> None:
    assert(start == "")
    start = "-x " + initial    
  
  doc = libxml2.parseFile(short_problem + ".xml")
  xml = doc.xpathNewContext()
  xml.xpathRegisterNs('cfs', 'http://www.cfs++.org')

  # we assume one hit or three for robust
  res = xml.xpathEval("//cfs:filter/cfs:density/@beta")
  if len(res) == 0:
    res = xml.xpathEval("//cfs:shapeMap/@beta")
    if len(res) == 0:  
      raise RuntimeError("beta not found for filter/density/@beta and shapeMap/@beta")
  for data in res:
    data.setContent(str(beta))
  
  doc.saveFile(beta_problem + ".xml")
  
  cmd = executable + " " + start + " -m " + mesh + " " + beta_problem
  execute(cmd, output=True, silent = failsafe)
      
  if show:
    execute("show_density.py " + beta_problem + ".density.xml --save " + beta_problem + ".png")

  return


parser = argparse.ArgumentParser(description='Make heaviside continuation. In the simple case just replace cfs by continuation.py')
parser.add_argument('problem', help="the problem xml without extension where '-beta_x' will be added")
parser.add_argument('-m', "--mesh", help="the mesh file with extension", required=True)
parser.add_argument('-x', '--initial', help="optional density.xml for initial beta (with extension)")
parser.add_argument('--start', help="initial beta", type=int, default=1)
parser.add_argument('--max', help="maxmum beta which will be calculated", type=int, default=64)
parser.add_argument('--inc', help="beta increment b += inc*b. inc=1 doubles", type=float, default=1.0)
parser.add_argument('--executable', help="what to call for cfs", default='cfs_rel')
parser.add_argument('--noshow', help="suppress calling show_density.py, e.g. for 3d!", action='store_true')
parser.add_argument('--failsafe', help="ignore cfs exiting with error", action='store_true')

args = parser.parse_args()

old = -1  
beta = args.start

while beta <= args.max:
 continuation(args.initial, old, digits(beta, 1), args.mesh, args.problem, args.executable, not args.noshow, args.failsafe)
 old = digits(beta, 1)
 beta += beta * args.inc
