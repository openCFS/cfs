#!/usr/bin/python

import sys
from optimization_tools import *

dens_lower = 0.01

mode = sys.argv[1]
filename = sys.argv[2]

if mode == 'ortho' or mode == 'dxortho':
  tmp = read_multi_design(filename, "tensor11", "tensor22", "tensor33", "rotAngle", "tensor12")#,"density")
  des = numpy.zeros((len(tmp[:,0]), 6))
  des[:,0:5] = tmp
else:
  tmp = read_multi_design(filename, "emodul-iso", "emodul", "gmodul", "rotAngle")#,"density")
  des = numpy.zeros((len(tmp[:,0]), 5))

pi = 3.1415265358
ny = int(round((des.shape[0]/2)**.5))
nx = 2*ny
for i in range(nx):
  for j in range(ny):
    if (i > .65*nx and j < .5*ny) or (i < .35*nx and j >= .5*ny):
      num = j*nx+i
      if des[num,3] < 1.9:
        des[num,3] = des[num,3] + .5*pi
        des[num,0] = tmp[num,1]
        des[num,1] = tmp[num,0]

maxjump = 2.5
for i in range(nx):
  if des[i,3]-des[i+nx,3] > maxjump:
    des[i,3] = des[i,3] - pi
  elif des[i,3]-des[i+nx,3] < -maxjump:
    des[i,3] = des[i,3] + pi
  if des[(ny-1)*nx+i,3]-des[(ny-2)*nx+i,3] > maxjump:
    des[(ny-1)*nx+i,3] = des[(ny-1)*nx+i,3] - pi
  elif des[(ny-1)*nx+i,3]-des[(ny-2)*nx+i,3] < -maxjump:
    des[(ny-1)*nx+i,3] = des[(ny-1)*nx+i,3] + pi
      

if mode == 'dxtrans-iso' or mode == 'dxortho':# or mode == 'dxlaminates':
  dens = (des[:,0]+des[:,1]+des[:,2])/maxTrace
  des[:,4] = numpy.maximum(dens, dens_lower)
  des[:,0] /= des[:,4]
  des[:,1] /= des[:,4]
  des[:,2] /= des[:,4]

if mode == 'trans-iso' or mode == 'dxtrans-iso':
  write_multi_design_file(filename[:-9] + "-tmp.dens.xml", des, "emodul-iso", "emodul", "gmodul", "rotAngle", "density")
elif mode == 'ortho':
  write_multi_design_file(filename[:-9] + "-tmp.dens.xml", des[:,0:5], "tensor11", "tensor22", "tensor33", "rotAngle", "tensor12")
elif mode == 'dxortho':
  write_multi_design_file(filename[:-9] + "-tmp.dens.xml", des, "tensor11", "tensor22", "tensor33", "rotAngle", "tensor12", "density")
elif mode == 'laminates':
  write_multi_design_file(filename[:-9] + "-tmp.dens.xml", des[:,0:4], "stiff1", "stiff2", "gmodul", "rotAngle")
elif mode == 'dxlaminates':
  write_multi_design_file(filename[:-9] + "-tmp.dens.xml", des, "stiff1", "stiff2", "gmodul", "rotAngle", "density")

