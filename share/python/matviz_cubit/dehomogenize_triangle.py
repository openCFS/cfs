#!/usr/bin/env python3

# This script has to be run from Cubit's Python command line
# or the Journal Editor.

import os
import numpy
import h5py
import hdf5_tools

# this reloads matviz_cubit each time, when this script runs.
# necessary, if matviz_cubit has changed
import matviz_cubit as mc
import importlib
importlib.reload(mc)


input = '/home/daniel/buckling/optimization/pillar/globBuck_locBuck_C1.warmstart/results_hdf5/0.480_snopt.cfs'
#input = '/home/daniel/buckling/microbuckling/setup1a_new/results_hdf5/macro.cfs'
input = '/home/daniel/buckling/two_scale_study_BE/results_hdf5/opt.cfs'

# radius of rounded corners
radius = 0.6

meshsize = 0.00125

samples = 8
samp = [14]

grad = 'linear'

boundary_enforcement = False
mid_enforcement = True

##################
for samples in samp:

  # create some filenames
  inputpath, inputfile = os.path.split(input)
  filename, _ = os.path.splitext(inputfile)
  meshfilename = os.path.join(inputpath, '{}_b{:.3f}_{:.6f}_{:d}'.format(filename, radius, meshsize, samples))
  savefilename = os.path.join(inputpath,'{}_b{:.3f}_{:d}'.format(filename, radius, samples))

  # open file
  fid = h5py.File(input, 'r')

  # get the FE centers for all 2d regions
  # if the mesh was created with Cubit, we might have 1d regions (e.g. to apply pressure loads)
  centers = [[None, None, None]]
  min_bb = [numpy.Inf, numpy.Inf, numpy.Inf]
  max_bb = [-numpy.Inf, -numpy.Inf, -numpy.Inf]
  for region in fid['/Mesh/Regions']:
    if fid['/Mesh/Regions/{}'.format(region)].attrs['Dimension'] != 2:
      continue
    reg_centers, reg_min_bb, reg_max_bb, elem_dim, _, _, _, _, _, _  = hdf5_tools.centered_elements(fid, region)
    centers = numpy.concatenate((centers, reg_centers))
    min_bb = numpy.min([min_bb, reg_min_bb], 0);
    max_bb = numpy.max([max_bb, reg_max_bb], 0);
  centers = centers[1:,:]

  coords = (centers, min_bb, max_bb, elem_dim)

  dim_2D = min_bb[2] == max_bb[2]

  # For buckling analysis, the design is stored in one step
  # and the buckling modes are stored in following steps
  # Thus, the design is not stored in the last step, but somewhere before it.
  step = min((99999, hdf5_tools.last_h5_step(fid)))
  while not hdf5_tools.has_element(fid, "design_stiff1_smart", step):
    step -= 1
  design = hdf5_tools.get_element(fid, "design_stiff1_smart", "mech", step)

  # if samples is only a number, create a list (one entry for each dimension)
  if samples is not None:
    samples = samples if isinstance(samples, (list, tuple)) else [int(samples), int(samples)]

  # create the geometry
  shape = mc.show_triangle_grad(coords, design, grad, samples, thres=None, equilateral=True, radius=radius, savefile=savefilename)

  print('Relative Surface Area: {}'.format(cubit.get_surface_area(shape) / (max_bb[0]-min_bb[0]) / (max_bb[1]-min_bb[1])))

  if boundary_enforcement:
    meshfilename += '_be'
    cubit.silent_cmd('create surface rectangle width 0.05 height 5.2 zplane ')
    r1 = cubit.get_last_id("surface")
    cubit.silent_cmd('move Surface {} x -0.025 y 2.6 include_merged '.format(r1))
    cubit.silent_cmd('create surface rectangle width 0.05 height 5.2 zplane ')
    r2 = cubit.get_last_id("surface")
    cubit.silent_cmd('move Surface {} x 1.025 y 2.6 include_merged '.format(r2))
    cubit.silent_cmd('unite surface {} {} {} '.format(shape,  r1, r2))
    cubit.silent_cmd('move Surface {} x 0.025 y 0 include_merged '.format(shape))

  if mid_enforcement:
    meshfilename += '_me'
    cubit.silent_cmd('webcut body 1 with plane xplane offset 0.5 ')
    shape = cubit.get_last_id("surface")
    surfaces = cubit.get_entities("surface")
    right_shape = [a for a in surfaces if a != shape][0]
    cubit.silent_cmd('move Surface {}  x 0.1 include_merged '.format(right_shape))
    cubit.silent_cmd('create surface rectangle width 0.1 height 5.2 zplane ')
    r = cubit.get_last_id("surface")
    cubit.silent_cmd('move Surface {} x 0.55 y 2.6 include_merged '.format(r))
    cubit.silent_cmd('unite body all ')
    shape = cubit.get_entities("surface")[0]

  mc.name_regions_and_nodes(shape, None)

  # mesh the geometry
  mc.mesh_shape(shape, meshsize, meshfilename)