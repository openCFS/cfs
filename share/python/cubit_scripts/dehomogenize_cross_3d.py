#!python

# This script has to be run from Cubit's Python command line
# or the Journal Editor.

import os
import numpy
import h5py
import hdf5_tools

# this reloads matviz_cubit each time, when this script runs.
# necessary, if matviz_cubit has changed
import cubit_scripts.matviz_cubit as mc
import importlib
importlib.reload(mc)


input = '/home/daniel/test.cfs'

#meshsize = 0.02

samples = [5,2,2]  #   492 surfaces -> 5s
samples = [10,4,4] #  3136 surfaces -> 1m 26s
samples = [15,6,6] #  9812 surfaces -> 13m 58s
samples = [20,8,8] # 20616 surfaces -> 56m 9s

grad = 'linear'

##################

# create some filenames
inputpath, inputfile = os.path.split(input)
filename, _ = os.path.splitext(inputfile)
#meshfilename = os.path.join(inputpath, '{}_b{:.3f}_{:.6f}_{:d}'.format(filename, radius, meshsize, samples))
#savefilename = os.path.join(inputpath,'{}_{:d}'.format(filename, samples))
savefilename = os.path.join(inputpath,'{}_{:d}'.format(filename, 1))

# open file
fid = h5py.File(input, 'r')

# get the FE centers for all 2d regions
# if the mesh was created with Cubit, we might have 1d regions (e.g. to apply pressure loads)
centers = [[None, None, None]]
min_bb = [numpy.Inf, numpy.Inf, numpy.Inf]
max_bb = [-numpy.Inf, -numpy.Inf, -numpy.Inf]
for region in fid['/Mesh/Regions']:
  if fid['/Mesh/Regions/{}'.format(region)].attrs['Dimension'] != 3:
    continue
  reg_centers, reg_min_bb, reg_max_bb, elem_dim, _, _, _, _, _, _  = hdf5_tools.centered_elements(fid, region)
  centers = numpy.concatenate((centers, reg_centers))
  min_bb = numpy.min([min_bb, reg_min_bb], 0);
  max_bb = numpy.max([max_bb, reg_max_bb], 0);
centers = centers[1:,:]

coords = (centers, min_bb, max_bb, elem_dim)

# For buckling analysis, the design is stored in one step
# and the buckling modes are stored in following steps
# Thus, the design is not stored in the last step, but somewhere before it.
step = min((99999, hdf5_tools.last_h5_step(fid)))
while not hdf5_tools.has_element(fid, "design_stiff1_smart", step):
  step -= 1
design = hdf5_tools.get_element(fid, "design_stiff1_smart", "mech", step)

if samples is not None:
  samples = samples if isinstance(samples, (list, tuple)) else [int(samples), int(samples), int(samples)]

design = [[a] for a in numpy.linspace(0.0185,0.146,20)]
shape = mc.show_cross_3d(coords, design, grad, samples, thres=None, savefile=None)
