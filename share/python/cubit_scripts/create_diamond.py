#!python

import cubit_scripts.matviz_cubit as mc
import importlib
importlib.reload(mc)

import numpy as np


repetitions = 10
reps = np.arange(1,11,1)
reps = [1]

# if design is given, we draw a parallelogram for homogenization
# if design is a filename, it will be read from this h5 file
designs = [0.5] #np.linspace(0.1,0.95,18)

# radius of rounded corners
radius = 0.6

meshsz = np.linspace(0.001,0.01,10)
meshsz = [0.0005]


for meshsize in meshsz:
  for repetitions in reps:
    for param in designs:
      shape, mids = mc.show_triangle_grad(None, None, None, None, None, True, radius, None, param, repetitions)
      meshfilename = 'cubit_diamond_v{:.3f}_b{:.3f}_{:.6f}_{:d}'.format(param, radius, meshsize, repetitions)

#      mc.mesh_shape(shape, meshsize, meshfilename)

      mc.mesh_shape_with_triangle(shape, 0.001, meshfilename, mids)
