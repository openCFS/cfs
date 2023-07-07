#!/usr/bin/env python3

# this reloads matviz_cubit each time, when this script runs.
# necessary, if matviz_cubit has changed
import cubit_scripts.cubit_tools as mc
import importlib
importlib.reload(mc)

import time

import os
import numpy
import h5py
import hdf5_tools
from cubit_scripts.cubit_tools import cubit

input = '/home/daniel/manni/BE_study/homogenDesign_v0.4.cfs'
input = '/home/daniel/manni/optimization/pillar/worst_case/globBuck_locBuck/validation2/0.700_snopt.cfs'

# radius of rounded corners
radius = 0.6

meshsize = 0.1

samples = 2
samp = [8,12,16,20]
samp=[12]

grad = 'linear'

boundary_enforcement = False
mid_enforcement = False

##################
for samples in samp:

  # create some filenames
  inputpath, inputfile = os.path.split(input)
  filename, _ = os.path.splitext(inputfile)
  filename = filename + '_rot'
  meshfilename = os.path.join(inputpath, '{}_b{:.3f}_{:.6f}_{:d}'.format(filename, radius, meshsize, samples))
  savefilename = os.path.join(inputpath,'{}_b{:.3f}_{:d}'.format(filename, radius, samples))

  if os.path.isfile(savefilename + ".cub5"):
    cubit.cmd('open "{}.cub5"'.format(savefilename))
    shape = cubit.get_last_id("surface")
    mids = numpy.loadtxt(savefilename + '.hole')

  else:
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
    shape, mids = mc.show_triangle_grad(coords, design, grad, samples, thres=None, equilateral=True, radius=radius, savefile=savefilename)
  
    print('Relative Surface Area: {}'.format(cubit.get_surface_area(shape) / (max_bb[0]-min_bb[0]) / (max_bb[1]-min_bb[1])))

  assert not (boundary_enforcement & mid_enforcement)

  if boundary_enforcement:
    meshfilename += '_be'
    bb = cubit.get_bounding_box("surface", shape)
    cubit.silent_cmd('create surface rectangle width 0.05 height {} zplane '.format(bb[5]))
    r1 = cubit.get_last_id("surface")
    cubit.silent_cmd('move Surface {} x -0.025 y {} include_merged '.format(r1, bb[5]/2))
    cubit.silent_cmd('create surface rectangle width 0.05 height {} zplane '.format(bb[5]))
    r2 = cubit.get_last_id("surface")
    cubit.silent_cmd('move Surface {} x {} y {} include_merged '.format(r2, bb[2]+0.025, bb[5]/2))
    cubit.silent_cmd('create surface rectangle width {} height 0.1 zplane '.format(bb[2]+0.1))
    r3 = cubit.get_last_id("surface")
    cubit.silent_cmd('move Surface {} x {} y {} include_merged '.format(r3, bb[2]/2, bb[5]+0.05))
    cubit.silent_cmd('unite body all ')
    shape = cubit.get_entities("surface")[0]
    cubit.silent_cmd('move Surface {} x 0.05 y 0 include_merged '.format(shape))

  if mid_enforcement:
    meshfilename += '_me'
    bb = cubit.get_bounding_box("surface", shape)
    cubit.silent_cmd('webcut body 1 with plane xplane offset {} '.format(bb[2]/2))
    shape = cubit.get_last_id("surface")
    surfaces = cubit.get_entities("surface")
    right_shape = [a for a in surfaces if a != shape][0]
    cubit.silent_cmd('move Surface {}  x 0.1 include_merged '.format(right_shape))
    cubit.silent_cmd('create surface rectangle width 0.1 height {} zplane '.format(bb[5]))
    r1 = cubit.get_last_id("surface")
    cubit.silent_cmd('move Surface {} x {} y {} include_merged '.format(r1, bb[2]/2+0.05, bb[5]/2))
    cubit.silent_cmd('create surface rectangle width {} height 0.1 zplane '.format(bb[2]+0.1))
    r2 = cubit.get_last_id("surface")
    cubit.silent_cmd('move Surface {} x {} y {} include_merged '.format(r2, bb[2]/2+0.05, bb[5]+0.05))
    cubit.silent_cmd('unite body all ')
    shape = cubit.get_entities("surface")[0]

  mc.name_regions_and_nodes(shape, None)

  # mesh the geometry
  #mc.mesh_shape(shape, meshsize, meshfilename)
  # *6 yields approximately same number of nodes and elements as meshing with cubit
  mc.mesh_shape_with_gmsh(shape, meshsize*2, meshfilename)
