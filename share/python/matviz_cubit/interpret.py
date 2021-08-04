#!python

import os
import numpy
import h5py
import hdf5_tools

import matviz_cubit.matviz_cubit as mc
import importlib
importlib.reload(mc)


input = '/home/daniel/buckling/optimization/pillar/globBuck_locBuck.warmstart/results_hdf5/0.480_snopt.cfs'
#input = '/home/daniel/buckling/microbuckling/setup1a_new/results_hdf5/macro.cfs'

# radius of rounded corners
radius = 0.6

meshsize = 0.02

samples = 1

grad = 'linear'

##################

inputpath, inputfile = os.path.split(input)
filename, _ = os.path.splitext(inputfile)
meshfilename = os.path.join(inputpath, '{}_b{:.3f}_{:.6f}_{:d}'.format(filename, radius, meshsize, samples))
savefilename = os.path.join(inputpath,'{}_b{:.3f}_{:d}'.format(filename, radius, samples))

fid = h5py.File(input, 'r')

centers = [[None, None, None]]
min_bb = [numpy.Inf, numpy.Inf, numpy.Inf]
max_bb = [-numpy.Inf, -numpy.Inf, -numpy.Inf]
for region in fid['/Mesh/Regions']:
  reg_centers, reg_min_bb, reg_max_bb, elem_dim, _, _, _, _, _, _  = hdf5_tools.centered_elements(fid, region)
  centers = numpy.concatenate((centers, reg_centers))
  min_bb = numpy.min([min_bb, reg_min_bb], 0);
  max_bb = numpy.max([max_bb, reg_max_bb], 0);
centers = centers[1:,:]

coords = (centers, min_bb, max_bb, elem_dim)

dim_2D = min_bb[2] == max_bb[2]

step = min((99999, hdf5_tools.last_h5_step(fid)))
while not hdf5_tools.has_element(fid, "design_stiff1_smart", step):
  step -= 1
design = hdf5_tools.get_element(fid, "design_stiff1_smart", "mech", step)

if samples is not None:
  samples = samples if isinstance(samples, (list, tuple)) else [int(samples), int(samples)]

shape = mc.show_triangle_grad(coords, design, grad, samples, None, equilateral=True, radius=radius, savefile=savefilename)

print('Relative Surface Area: {}'.format(cubit.get_surface_area(shape) / (max_bb[0]-min_bb[0]) / (max_bb[1]-min_bb[1])))

mc.mesh_shape(shape, meshsize, meshfilename)
