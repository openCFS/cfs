#!/usr/bin/env python
import sys
import libxml2

from cfs_utils import *

## performs heaviside continuation by doubling filter/density/beta, starting from 1

def continuation(old_beta, beta, mesh, problem):
  start = cond(old_beta == 1, "", "-x " + problem + "-beta_" + str(old_beta))
  
  doc = libxml2.parseFile(problem + ".xml")
  xml = doc.xpathNewContext()
  xml.xpathRegisterNs('cfs', 'http://www.cfs++.org')

  replace(xml, "//cfs:filter/cfs:density/@beta", str(beta))
  replace(xml, "//cfs:export/@file", problem + "-beta_" + str(beta))
  
  doc.saveFile(problem + "-beta_" + str(beta) + ".xml")
  
  execute("cfs_trunk " + start + " -m " + mesh + " " + problem + "-beta_" + str(beta))
  return

def main():
  if len(sys.argv) != 3:
    print "usage: continuation <mesh file> <problem>"
    print "  e.g. continuation cantilever2d_5.mesh mech2d"
    print ""
    print "Loops for beta=1 to beta=512 by doubling beta"
    print "Sets filter/density/@beta and export/@file"
    print 'You need export/@file to be present, e.g. file="[problem]"' 
  else:
    mesh = sys.argv[1]
    problem = sys.argv[2]
    old = 1
    beta = 1
    while beta <= 512:
      continuation(old, beta, mesh, problem)
      old = beta
      beta *= 2
  
  return    

main()