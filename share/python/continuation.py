#!/usr/bin/env python
import sys
import libxml2
import argparse

from cfs_utils import *

## performs heaviside continuation by doubling filter/density/beta, starting from 1

def continuation(initial, old_beta, beta, mesh, problem, executable):
    
  # first without iniital  
  start = "" if old_beta == 1 else "-x " + problem + "-beta_" + str(old_beta) + ".density.xml"
  # now check for initial
  if old_beta == 1 and initial <> None:
    assert(start == "")
    start = "-x " + initial    
  
  doc = libxml2.parseFile(problem + ".xml")
  xml = doc.xpathNewContext()
  xml.xpathRegisterNs('cfs', 'http://www.cfs++.org')

  replace(xml, "//cfs:filter/cfs:density/@beta", str(beta))
  
  doc.saveFile(problem + "-beta_" + str(beta) + ".xml")
  
  cmd = executable + " " + start + " -m " + mesh + " " + problem + "-beta_" + str(beta)
  execute(cmd, output=True)
  return


parser = argparse.ArgumentParser()
parser.add_argument("--mesh", help="the mesh file with extension", required=True)
parser.add_argument('--problem', help="the problem xml without extension where '-beta_x' will be added", required=True)
parser.add_argument('--initial', help="optional density.xml for initial beta (with extension)")
parser.add_argument('--max_beta', help="maxmum beta which will be calculated", type=int, default=64)
parser.add_argument('--executable', help="what to call for cfs", default='cfs_rel')

args = parser.parse_args()

old = 1
beta = 1

while beta <= args.max_beta:
 continuation(args.initial, old, beta, args.mesh, args.problem, args.executable)
 old = beta
 beta *= 2
