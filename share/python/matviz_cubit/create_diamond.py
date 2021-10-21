#!python

import numpy as np

import matviz_cubit.matviz_cubit as mc
import importlib
importlib.reload(mc)

repetitions = 10
reps = np.arange(1,11,1)
reps = [2]

# if design is given, we draw a parallelogram for homogenization
# if design is a filename, it will be read from this h5 file
designs = [0.05] #np.linspace(0.1,0.95,18)

# radius of rounded corners
radius = 0.6

meshsz = np.linspace(0.001,0.01,10)
meshsz = [0.0005]


for meshsize in meshsz:
  for repetitions in reps:
    for param in designs:
      shape = mc.show_triangle_grad(None, None, None, None, None, True, radius, None, param, repetitions)
      meshfilename = 'cubit_diamond_v{:.3f}_b{:.3f}_{:.6f}_{:d}'.format(param, radius, meshsize, repetitions)

      mc.mesh_shape(shape, meshsize, meshfilename)
