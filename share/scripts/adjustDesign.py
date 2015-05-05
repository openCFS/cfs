#!/usr/bin/python

import sys
from optimization_tools import *

dens_lower = 0.01

if len(sys.argv) < 2:
  print "Usage: " + sys.argv[0] + " <mode> <filename> <maxTrace> <theta>"
  print "<theta> is only used in transversal-isotropic formulations, theta = nu_13*nu_31"
  print "<mode> can be:"
  print " - trans-iso:   parametrization using engineering constants and a transversal-isotropic material model"
  print " - dxtrans-iso: simultaneous topology and material optimization using engineering constants and a transversal-isotropic material model"
  print " - ortho:	 free orthotropic material optimization"
  print " - dxortho:     simultaneous topology and free orthotropic material optimization"
  print " - laminates:   homogenized layered material"
  print " - dxlaminates: simultaneous topology and homogenized layered material optimization"
  sys.exit(1)

mode = sys.argv[1]
filename = sys.argv[2]
maxTrace = float(sys.argv[3])

#tmp = read_multi_design(filename, "stiff1", "stiff2", "rotAngle", "density")#, None, False, "physical")
#tmp[:,0] *= 20
#tmp[:,1] *= 20
des = numpy.zeros((tmp.shape[0], 4));
tmp = read_multi_design(filename, "emodul-iso", "emodul", "gmodul", "rotAngle")#,"density")
#des = numpy.zeros((tmp.shape[0], 5));
des[:,0:tmp.shape[1]] = tmp




#des[:,0]*=(des[:,3]**3)
#des[:,1]*=(des[:,3]**3)

if mode == 'dxtrans-iso' or mode == 'dxortho' or mode == 'dxlaminates':
  dens = (des[:,0]+des[:,1]+2*des[:,2])/maxTrace
  des[:,4] = numpy.maximum(dens, dens_lower)
  des[:,0] /= des[:,4]
  des[:,1] /= des[:,4]
  des[:,2] /= des[:,4]

if mode == 'trans-iso' or mode == 'dxtrans-iso':
  theta = float(sys.argv[4])
  des[:,0] *= (1-theta)
  des[:,1] *= (1-theta)
 # des[:,2] *= .5
elif mode == 'ortho' or mode == 'dxortho':
  des[:,0] = des[:,0]**.5
  des[:,1] = des [:,1]**.5
elif mode == 'laminates' or mode == 'dxlaminates':
  des[:,0] = .5*(des[:,0]+.5*des[:,3])/maxTrace
  des[:,1] = .5*(des[:,1]+.5*des[:,3])/maxTrace
else:
  print "<mode> not defined!"

if mode == 'ortho':
  write_multi_design_file(filename[:-9] + "-" + mode + ".dens.xml", des[:,0:4], ["tensor11", "tensor22", "tensor33", "rotAngle"])
elif mode == 'dxortho':
  write_multi_design_file(filename[:-9] + "-" + mode + ".dens.xml", des, ["tensor11", "tensor22", "tensor33", "rotAngle", "density"])
elif mode == 'laminates':
  write_multi_design_file(filename[:-9] + "-" + mode + ".dens.xml", des[:,0:3], ["stiff1", "stiff2",  "rotAngle"])
elif mode == 'dxlaminates':
  write_multi_design_file(filename[:-9] + "-" + mode + ".dens.xml", des, ["stiff1", "stiff2", "rotAngle", "density"])
elif mode == 'trans-iso':
  write_multi_design_file(filename[:-9] + "-" + mode + ".dens.xml", des[:,0:4], ["emodul-iso", "emodul", "gmodul", "rotAngle"])
elif mode == 'dxtrans-iso':
  write_multi_design_file(filename[:-9] + "-" + mode + ".dens.xml", des, ["emodul-iso", "emodul", "gmodul", "rotAngle", "density"])