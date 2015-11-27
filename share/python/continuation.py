#!/usr/bin/env python
import sys
import libxml2
import argparse

from cfs_utils import *

## performs heaviside continuation by doubling filter/density/beta, starting from 1

def continuation(initial, old_beta, beta, mesh, short_problem, executable, show):
    
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

  replace(xml, "//cfs:filter/cfs:density/@beta", str(beta))
  
  doc.saveFile(beta_problem + ".xml")
  
  cmd = executable + " " + start + " -m " + mesh + " " + beta_problem
  execute(cmd, output=True)
  if show:
    execute("show_density.py " + beta_problem + ".density.xml --save " + beta_problem + ".png")

  return


parser = argparse.ArgumentParser(description='Make heaviside continuation. In the simple case just replace cfs by continuation.py')
parser.add_argument('problem', help="the problem xml without extension where '-beta_x' will be added")
parser.add_argument('-m', "--mesh", help="the mesh file with extension", required=True)
parser.add_argument('-x', '--initial', help="optional density.xml for initial beta (with extension)")
parser.add_argument('--beta', help="maxmum beta which will be calculated", type=int, default=64)
parser.add_argument('--executable', help="what to call for cfs", default='cfs_rel')
parser.add_argument('--non_show', help="suppress calling show_density.py, e.g. for 3d!", action='store_true')

args = parser.parse_args()

old = -1  
beta = 1

while beta <= args.beta:
 continuation(args.initial, old, beta, args.mesh, args.problem, args.executable, not args.non_show)
 old = beta
 beta *= 2
