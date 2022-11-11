#!/usr/bin/env python
# sometimes h5py must be import first
from hdf5_tools import *
from cfs_utils import *
import matviz_3d_ortho # two scale 3d ortho stuff
from matviz_2d  import *
import matviz_io
from matviz_streamline import *
from matviz_unstructured_mesh import *
from matviz_vtk import *
from mesh_tool import *
from optimization_tools import *

import argparse
import sys
import os.path
import xml.etree.ElementTree
import xml.dom.minidom

try:
  from mpi4py import MPI
except ImportError:
  print("WARNING: Could not load mpi4py")

TWO_SCALE = {
  "hom_cross_bar" : 3,
  "hom_frame" : 2,
  "hom_framed_cross" : 4,
  "hom_ortho_3d": 3,
  "hom_rect" : 3,
  "hom_rect_mod" : 3,
  "hom_rot_cross" : 2,
  "hom_sheared_rot_cross" : 3,
  "hom_triangle": 1
  }

# we extract the main processing such that we can run it within a loop to find the optimal scaling
# @param force_scale overwrites args.scale
# @param min_bb/max_bb: min/max coordinates of bounding box
# @return volume if calculated (e.g. via --save a pixel image) otherwise None
def perform(args, h5_read, dim_2D, tensor, centers, force_scale = None, nondes = None, min_bb = None, max_bb = None, elems_in_regions = None):
  volume = None  # might never be set

  scale = force_scale if force_scale else args.scale

  coords = (centers, min_bb, max_bb, elem_dim)

  # viz is either Image or polydata
  viz = None
  if args.show in ("stream", "simp") or args.show in TWO_SCALE:
    if args.show == "simp":
      args.parametrization = "simp"

    if h5_read:
      design = matviz_io.read_design(f, args, dim_2D, centers, TWO_SCALE)
    else:
      design, coords = read_stiff_angle_matlab(args.input)

    # add angle bias, e.g. by 90 deg to correct thomas
    design['angle'] += args.angle_bias * numpy.pi / 180
    # scale angle, e.g  by -1 to correct for current standard 2D rotation direction (this is not the mathematical direction! FIXME if needed)
    design['angle'] *= -1.0
    if args.angle_factor != 1.0:
      print('scale angle by ' + str(args.angle_factor))
      design['angle'] *= args.angle_factor

#    for key, value in design.items():
#      print('unscaled {:s} in [{:s}:{:s}]'.format(key, str(numpy.min(value)), str(numpy.max(value))))

    if dim_2D:
      func = perform_2d
    else:
      func = perform_3d

    # do we have to do 1D optimization?
    if args.target_volume is None:
      viz = func(args, coords, design, scale, nondes)
    elif args.target_volume == 0.0:
      viz = func(args, coords, design, 0.0, nondes)
    else:
      if args.tv_alg == 'bisec':
        viz = match_volume_bisec(args.target_volume, func, [args, coords, design, scale, nondes], 3)
      else:
        viz = match_volume_raster(args.target_volume, func, [args, coords, design, scale, nondes], 3)

    if get_MPI_rank()==0:
      if viz is None:
        print('Error: no visualization calculated!')
      else:
        volume = matviz_io.show_or_write(viz, args)
  else:
    show_boebbale(tensor, args, coords, scale)

  return volume

def perform_2d(args, coords, design, scale, nondes = None):
  if args.show == 'simp':
    pass
  elif args.show == "hom_rect":
    if args.hom_grad == 'none':
      viz = show_frame(coords, design, args.hom_dir, args.res, args.scale)
    else:
      viz = show_frame_grad(coords, design, args.hom_grad, args.hom_dir, args.res)
  elif args.show == "hom_rect_mod":
    design['angle'] = design['angle'][:,0]
    viz = show_modified_frame(coords, design, args.res)
  elif args.show == "hom_rot_cross":
    design['angle'] = design['angle'][:,0]
    if args.hom_grad == 'none':
      viz = show_rot_cross(coords, design, args.hom_dir, args.res, scale, args.color, args.save)
    else:
      viz = show_rot_cross_grad(coords, design, args.hom_grad, args.hom_dir, args.res, scale, args.save)
  elif args.show == "hom_sheared_rot_cross":
    design['angle'] = design['angle'][:,0]
    viz = show_sheared_rot_cross(coords, design, args.hom_dir, args.res, args.scale, args.color, args.save)
  elif args.show == "hom_frame":
    viz = show_frame2(coords, design, args.res, args.scale, args.color, args.save)
  elif args.show == "hom_framed_cross":
    viz = show_framed_cross(coords, design, args.res, args.scale, args.color, args.save)
  elif args.show == "hom_cross_bar":
    design['angle'] = numpy.zeros((len(s1),1))
    viz = show_cross_bar(coords, design, args.hom_dir, args.res, args.scale, args.color, args.save)
  elif args.show == "hom_triangle":
    assert(args.hom_grad != 'none')
    if args.target_volume is None:
      viz = show_triangle_grad(coords, design, args.hom_grad, args.hom_samples, args.res, args.thres, args.save, args.access, True, args.bc_bend, args.parameter, args.repetitions)
    else:
      if args.parameter is not None:
        viz = show_triangle_grad(coords, design, args.hom_grad, args.hom_samples, args.res, args.thres, args.save, args.access, True, args.bc_bend, scale*args.parameter, args.repetitions)
      else:
        design_backup = design['s1']
        design['s1'] = scale * design['s1']
        try:
          viz = show_triangle_grad(coords, design, args.hom_grad, args.hom_samples, args.res, args.thres, args.save, args.access, True, args.bc_bend, args.parameter, args.repetitions)
        finally:
          design['s1'] = design_backup
  elif args.show == "stream":
      samples = args.hom_samples.split(',')
      design['angle'] = design['angle'][:,0]
      viz = show_streamline(coords, design, args.hom_dir, scale, args.minimal, args.stream_style, args.stream_step, float(samples[0]), float(samples[1]), args.stream_max_traces_per_cell, args.res, args.save != None, matviz_io.infoxml, args.stream_force)
  else:
    assert(False)
  return viz

def perform_3d(args, coords, design, scale, nondes = None):
  if args.hom_samples or args.cell_size:
    if args.type == "apod6":
      valid_position = valid_position_apod6
      print('Apod6 is calculated!')
    elif args.type == "robot":
      valid_position = valid_position_robot
      print('Robot is calculated!')
    elif args.type == "lufo":
      valid_position = valid_position_lufo
      print('Lufo bracket is calculated!')
    else:
      valid_position = None
      print('Warning: No type for valid_position was selected!')

    if args.show == 'hom_rect' and args.hom_grad == 'none':
      raise Exception('For hom_rect in 3D with hom_samples you need to specify hom_grad.')

    if args.cell_size:
      tmp = args.cell_size.split(',')
      csize = [float(tmp[0]),float(tmp[1]),float(tmp[2])]
      if args.mesh:
        # number of fine elements in each direction
        tmp = args.nf.split(',')
        n_f = [int(tmp[0]),int(tmp[1]),int(tmp[2])]
        if args.type == "apod6":
          valid_position = valid_position_apod6
          valid_ring_position = valid_ring_position_apod6
          print('Apod6 is validated!')
        elif args.type == "robot":
          valid_position = valid_position_robot
          #TODO: not implemented yet for robot arm
          valid_ring_position = None
          print('Robot is validated!')
        elif args.type == "lufo":
          valid_position = valid_position_lufo
          valid_ring_position = None
          print('Lufo bracket is calculated!')

        assert(nondes is not None)
        assert(nondes[0] is not None)
        assert(len(nondes[0]) == 3)

        me = create_validation_mesh(coords, nondes[0], design, None, args.hom_grad, args.hom_dir, scale, n_f, valid_position, valid_ring_position, args.type, args.thres, csize, args.show == 'simp')
        write_ansys_mesh(me, args.mesh+".mesh")
      else:
        if args.show == "hom_rot_cross":
          viz = create_3d_cross_ip(coords, design, args.hom_samples, args.hom_grad, scale, valid_position, args.thres, csize)
        elif args.show == "hom_rect":
          viz = create_3d_frame_ip(coords, design, args.hom_samples, args.hom_grad, scale, valid_position, args.thres)
        elif args.show == "hom_ortho_3d":
          raise Exception("hom_ortho_3d not possible with cell_size")
    else:
      tmp = args.hom_samples.split(',')

      samples = [int(tmp[0]),int(tmp[0]),int(tmp[0])] if len(tmp) == 1 else [int(tmp[0]),int(tmp[1]),int(tmp[2])]

      if args.show == "hom_ortho_3d" or args.mesh:
#       name = "interpretation_ortho_3d_box_varel_" + str(samples[0]) + "_" + str(samples[1]) + "_" + str(samples[2]) + "_bc_res_" + str(args.bc_res) + ".stl"
        name = args.save + ".stl" if args.save is not None else "interpretation_mc"
        # region_map: maps local node id (list idx) in given design region (e.g. "mech") to global node id in all regions
        reg_info = {"nodes":reg_nodes, "elements":elems_in_regions, "connectivity":connectivity, "region_map":reg_nodes_map}
        if nondes:
          viz = matviz_3d_ortho.create_3d_interpretation_ortho(args, reg_info, coords, min_bb, max_bb, design, scale, samples, args.hom_grad,nondes=nondes)
        else:
          viz = matviz_3d_ortho.create_3d_interpretation_ortho(args, reg_info, coords, min_bb, max_bb, design, scale, samples, args.hom_grad,nondes=None)
        if args.save:
          if args.save.endswith(".vtp"):
            name = args.save[:-4]+".stl"
          matviz_vtk.write_stl(viz, name)
        if args.type == "box_varel" or args.type == "ppbox" or args.mesh:
          if not args.save: # write surface mesh in case we haven't done it before
            matviz_vtk.write_stl(viz, name)
#          write_ansys_mesh(create_volume_mesh_from_stl(name),name[:-4]+".mesh")
          create_volume_mesh_with_gmsh(name)
          viz = None # avoid showing or writing vtp file
      else:
        viz = create_3d_cross_ip(coords, design, samples, args.hom_grad, scale, valid_position, args.thres)
  else:  # no sample
    if args.show == 'simp':
      viz = create_block(coords, design, scale, args.thres, elems_in_regions)
    elif args.show == 'hom_frame':
      if args.hom_grad == 'None':
        viz = create_3d_frame(coords, design, args.hom_dir, scale)
    elif args.unstructured:
      # dimensions of the new regular grid
      n = args.unstructured.split(',')
      nx = int(n[0])
      ny = int(n[1])
      nz = int(n[2])
      if args.nodefile:
        # load mesh from hdf5 file
        mesh = create_mesh_from_hdf5(f, ['design', 'nondesign', 'void'], 'load', 'support')
        # write design values for each element node to file
        write_node_design(args.nodefile + '_' + repr(nx) + 'x' + repr(ny) + 'x' + repr(nz), mesh, coords, design)
      else:
        assert(nondes is not None)
        assert(nondes[0] is not None)
        assert(len(nondes[0]) == 3)
        # create structured 3D mesh from given hdf5_file which can be non-structured
        mesh = create_3d_mesh_unstructured(coords, nondes[0], nondes_force, nondes_support, design, nx, ny, nz, args.hom_grad, args.scale)
        # postprocess the structured mesh for the robot arm, removal of some irregularities
        mesh = post_process_mesh(mesh, coords, nondes[0], nx, ny, nz)
    else:
      raise NotImplementedError
  return viz

def show_boebbale(tensor, args, coords, scale):
  if dim_2D:
    assert(args.show is None or args.show in ['ortho_norm', 'mono_norm', 'ortho_err'])
    if h5_read:
      tensors = get_element(f, args.tensor, args.h5_region, args.h5_step)
      tensor = numpy.zeros((len(tensors),3,3)) if tensors.shape[1] == 6 else numpy.zeros((len(tensors),6,6))
      for i, vec in enumerate(tensors):
        tensor[i] = to_mech_tensor(vec)
    angle, data = perform_rotations(tensor, int(args.sampling), args.show)

    if args.gnuplot != None:
      matviz_io.write_angle_data(args.gnuplot, angle, data)
      sys.exit() # all done
    else:
      viz = show_orientational_stiffness(coords, angle, data, args.res, scale, args.axes)
      matviz_io.show_or_write(viz, args)
  else:
    assert False, 'use tensorviz.py'
    
def match_volume_bisec(target_volume, func, params, param_idx):
# performs a bisection to match volume
# @param vol  desired volume
# @param func  function to calculate the image (show_...)
# @param params  array of parameters for func. the parameter to perform bisection with
#                can be chosen arbitrarily as we will change it anyway
# @param param_idx  index of the parameter in params to perform bisection with (zero-based)
# @return viz  Image returned by func for the optimal parameter
#
# @example  params = [coords, design, args.hom_grad, args.hom_samples, args.res, args.thres, args.save, True, 0.0, args.bc_bend]
#           viz = match_volume_bisec(args.parameter, show_triangle_grad, params, 8)
    print("match volume {:.3f} by bisection".format(target_volume))

    tol = 1e-15
    # initial values
    param1 = 0.0
    param2 = pow(2,30)
    # func should raise a ValueError, if param does not fit to other function parameters
    # e.g. show_triangle_grad raises a ValueError, if param is too large for the choosen radius
    while True:
      try:
        params[param_idx] = param1
        viz1 = func(*params)
        vol1 = 1 - numpy.average(numpy.array(viz1)) / 255 if type(viz1) is Image.Image else viz1
        if(vol1 > target_volume):
          param1 = param1 / 2
          break
        break
      except ValueError:
        if param1 == 0.0:
          param1 = param1 + .00001
        else:
          param1 = param1 * 2
    while True:
      try:
        params[param_idx] = param2
        viz2 = func(*params)
        vol2 = 1 - numpy.average(numpy.array(viz2)) / 255 if type(viz2) is Image.Image else viz2
        if(vol2 < target_volume):
          param2 = param2 * 2
        break
      except ValueError:
        param2 = param2 / 2
    assert(param2 > param1)

    vol1 = 1 - numpy.average(numpy.array(viz1)) / 255 if type(viz1) is Image.Image else viz1
    vol2 = 1 - numpy.average(numpy.array(viz2)) / 255 if type(viz2) is Image.Image else viz2
    print("bisection initial values: {}, {}".format(vol1, vol2))
    assert((vol1-target_volume)*(vol2-target_volume) < 0),'no root in initial interval [{},{}]'.format(param1,param2)

    # actual bisection
    iter = 0
    while param2-param1 > tol:
      try:
        param3 = (param1 + param2) / 2
        params[param_idx] = param3
        viz3 = func(*params)
        vol3 = 1 - numpy.average(numpy.array(viz3)) / 255 if type(viz3) is Image.Image else viz3
        # stop if there is no design change
        if np.abs(vol3-vol1) < tol or np.abs(vol2-vol1) < tol:
          break
        #print((param1,param3,param2,vol1,vol3,vol2))
        if (vol1 - target_volume)*(vol3 - target_volume) < 0:
          param2 = param3
          vol2 = vol3
        else:
          param1 = param3
          vol1 = vol3
        iter = iter + 1
        if iter % 10 == 0:
          print("bisection iter {:d}: {}, {}".format(iter, vol1, vol2))
      except ValueError:
        raise
    # check which volume is closer to the desired volume. we have to do this because we get rounding errors due to the discrete resolution of viz
    params[param_idx] = param1 if target_volume-vol1 < vol2-target_volume else param2
    viz = func(*params)
    print('best parameter: {}'.format(params[param_idx]))

    return viz

def match_volume_raster(target_volume, func, params, param_idx):
  print("match volume {:.3f} by rastering".format(target_volume))

  if args.scale > 0:
    raise Exception("Don't give --scale and --target_volume concurrently!")

  # unfortunately we are not necessary convex! therefore we may not do bisection
  lower = 0.0
  upper = 10.0
  level = 0
  best_err = 999
  err = 999
  best_s = -1
  num_intervals = 11
  if args.info is not None:
    xml_vm = xml.etree.ElementTree.SubElement(matviz_io.infoxml, "volume_match")
  while upper - lower > 1e-20 and abs(best_err) > 1e-5:
    xml_level = xml.etree.ElementTree.SubElement(xml_vm, "level")
    xml_level.set("idx", str(level))
    xml_level.set("lower", str(lower))
    xml_level.set("upper", str(upper))
    xml_level.set("delta", str(upper - lower))
    for s in numpy.linspace(lower, upper, num_intervals):
      params[param_idx] = s
      viz = func(*params)
      vol = 1 - numpy.average(numpy.array(viz)) / 255 if type(viz) is Image.Image else viz
      err = abs(vol - target_volume)
      if vol == None:
        raise Exception("volume is None")
      if(err < best_err):
        best_err = err
        best_s = s

      if args.info is not None:
        xml_search = xml.etree.ElementTree.SubElement(xml_level, "search")
        xml_search.set("target", str(target_volume))
        xml_search.set("current", str(vol))
        xml_search.set("scale_argument", str(s))
        xml_search.set("err", str(err))
      else:
        msg =  ">>>> target_volume: level=" + str(level) + " lower= " + str(lower) + " upper= " + str(upper) + " delta= " + str(upper - lower) + " best=" + str(best_err) + " test=" + str(s)
        msg += " -> " + str(vol) + " err=" + str(err) + (" *" if err == best_err else "")
        print(msg)
    xml_best = xml.etree.ElementTree.SubElement(xml_search, "best")
    xml_best.set("err", str(best_err))
    xml_best.set("scale", str(best_s))

    new_lower = numpy.max((0.0, best_s - (upper - lower) / (num_intervals-1)))
    upper = numpy.min((10.0, best_s + (upper - lower) / (num_intervals-1)))
    lower = new_lower
    level += 1

  # the last result does not mean to be the best result
  if err != best_err:
    params[param_idx] = s
    viz = func(*params)
    vol = 1 - numpy.average(numpy.array(viz)) / 255 if type(viz) is Image.Image else viz
    if args.info is not None:
      tv = xml.etree.ElementTree.SubElement(xml_vm, "best_target_volume")
      tv.set("target", str(target_volume))
      tv.set("current", str(vol))
      tv.set("scale_argument", str(best_s))
      tv.set("err", str(abs(vol - target_volume)))
    else:
      print("!!!!! best_target_volume: best_scale=" + str(best_s) + " -> " + str(vol) + " err=" + str(abs(vol - target_volume)))

  return viz

def get_MPI_rank():
  try:
    return MPI.COMM_WORLD.Get_rank()
  except NameError: # shall be name 'MPI' is not defined when from mpi4py import MPI failed
    return 0


parser = argparse.ArgumentParser()
parser.add_argument("input", help="a openCFS h5 file")
parser.add_argument("--h5_step", help="step number, too high is last (default '9999')", default=9999, type=int)
parser.add_argument("--h5_region", help="design region name (default 'mech')", default="mech")
parser.add_argument("--h5_nondes", help="non-design (solid) region name (default 'non-design')", default="None")
parser.add_argument("--h5_nondes_void", help="non-design (void) region name (default 'non-design')", default="None")
# parser.add_argument('--h5_nonreg', action='store_true', help='assume the h5 file to be nonregular')
parser.add_argument("--h5_info", action="store_true", help="dump some meta data information about the h5 file")
parser.add_argument("--hist", help="plot histograms of the results in the h5 file. Has to be used with --show", type=int, default=0)
parser.add_argument("--tensor", help="tensor name: 'mechTensor', 'piezoTensor, 'elecTensor'", default="mechTensor")
parser.add_argument("--scale", help="manual scaling factor", default=-1.0, type=float)
parser.add_argument("--target_volume", help="scale design to given volume", type=float)
parser.add_argument("--tv_alg", help="algorithm to match target volume", choices=['bisec','raster'], default='bisec')
parser.add_argument("--res", help="x-resolution (default 1000)", default=1000, type=int)
parser.add_argument("--sampling", help="sampling rate for tensor rotation (default 180)", default=180, type=float)
parser.add_argument("--show", help="mode within boebbale, hom_rect or streamline", choices=['ortho_norm', 'mono_norm', 'ortho_err', 'hom_rect', 'hom_rot_cross', 'hom_sheared_rot_cross', 'hom_frame', 'hom_framed_cross', 'hom_cross_bar', 'stream', 'hom_rect_mod', 'simp', 'hom_ortho_3d', 'hom_triangle'])
parser.add_argument("--axes", help="show axes", action='store_true', default=False)
parser.add_argument("--notation", help="mandel | voigt (default 'voigt')", choices=['mandel','voigt'], default="voigt")
parser.add_argument("--access", help="'plain' or 'smart' (default 'smart')", default="smart")
parser.add_argument("--hom_grad", help="interpolation of design: 'none', 'nearest', linear', 'cubic' (default 'linear')", default="linear", choices=['none', 'nearest', 'linear', 'cubic'])
parser.add_argument("--hom_dir", help="visualization of stiffness directions (default 'all')", default="all", choices=['all', 'horizontal', 'vertical', 'sagittal'])
parser.add_argument("--hom_samples", help="activates interpolation and the value gives samples in x,y,z direction")
parser.add_argument("--angle_factor", help="factor for angle. -1.0 turns, 0.0 disables angles", default=1.0, type=float)
parser.add_argument("--angle_bias", help="bias for the angle in deg. 90 switches s1 and s2", default=0.0, type=float)
parser.add_argument("--cell_size", help="cell size in [mm] in x,y,z direction", type=float)
parser.add_argument("--stream_style", help="select visualization", choices=['line', 'thick'], default='thick')
parser.add_argument("--stream_step", help="step length for ODE integration per macro cell", type=float, default=0.2)
parser.add_argument("--stream_max_traces_per_cell", help="maximum number of traces such that we may start a trace (>= 1)", type=int, default=1)
parser.add_argument("--stream_ode", help="method to solve the ODE", default="euler", choices=['euler', 'midpoint'])
parser.add_argument("--stream_force", help="force streamlines for special cases", choices=['right_lower', 'rhombus'])
parser.add_argument("--minimal", help="minimal stiffness to be drawn, will be scaled", type=float, default=0.0)
parser.add_argument("--parametrization", help="parametrization of the stiffness tensor", default="hom_rect", choices=['simp','hom_rect', 'trans-iso', 'ortho'])
parser.add_argument("--save", help="save 'image.png' (pixel), 'image.pdf' (vector) or VTK Poly Data file 'file.vtp'")
parser.add_argument("--cam", help="set camera (7 space-separated floats): position, focal point, roll", nargs=7, type=float, default=None)
parser.add_argument("--gnuplot", help="for single tensors: creates gnuplot file instead of image")
parser.add_argument("--penalty", help="penalty parameter for SIMP (default 5)", default=5.0)
parser.add_argument("--color", help="only for hom_rot_cross: 'black' or colormap from http://matplotlib.org/examples/color/colormaps_reference.html, default='gray'", default="grayscale")
parser.add_argument("--info", help="creates a xml file of given name with additional information")
parser.add_argument("--unstructured", help="number of structured elements per coordinate as list nx,ny,nz", default="")
parser.add_argument("--nodefile", help="name of the design to node file", default="")
parser.add_argument("--thres", help="threshold value for plot", type=float, default=0.0)
parser.add_argument("--mesh", help="create 3D mesh from optimized 2-scale result for validation", default="")
parser.add_argument("--nf", help="requires --mesh, number of fine elements in x,y,z direction")
parser.add_argument("--type", help="type of 3D object for 2-scale visualization", choices=['apod6', 'robot','ppbox','box_varel'])
# 3d ortho basecell stuff
parser.add_argument("--bc_res", help="resolution of voxelized ortho basecell", type=int)
parser.add_argument("--bc_interpolation", help="interpolation type for ortho basecell (linear or heaviside)", choices=['linear', 'heaviside'], default="heaviside")
parser.add_argument("--bc_beta", help="for heaviside interpolation (default 7.0)", type=float, default=7)
parser.add_argument("--bc_eta", help="for heaviside interpolation (default 0.5)", type=float, default=0.5)
parser.add_argument("--bc_bend", help="bending of spline (default 0.5)", type=float, default=0.5)
parser.add_argument("--bc_smooth", help="number auf Taubin smoothing steps", type=int, default=0)
parser.add_argument("--bc_thresh", help="lower and upper threshold (diameter) for ortho basecell, e.g. 1e-9,0.94")
parser.add_argument("--parameter", help="parameter for different usage", type=float, default=None)
parser.add_argument("--repetitions", help="number of repetitions of basecell", type=int, default=1)

args = parser.parse_args()

# check 3d ortho basecell parameters
if args.show == "hom_ortho_3d":
  if not args.bc_res:
    raise Exception("bc_res required")
  if not args.bc_bend:
    raise Exception("bc_bend parameter required")
  if not args.bc_interpolation:
    raise Exception("interpolation type (linear or heaviside) required")
  elif not args.bc_beta or not args.bc_eta:
    raise Exception("beta and eta values required for heaviside interpolation")
if args.type == "box_varel" or args.type == "ppbox":
  # in this case everything belongs to design domain
  args.h5_nondes = "None"

if args.info is not None:
  matviz_io.create_info_xml()
  cmd = sys.argv[0].split('/')[-1]
  for i in range(1, len(sys.argv)):
    cmd += ' ' + sys.argv[i]
  doc = xml.etree.ElementTree.SubElement(matviz_io.infoxml, "commandLine")
  doc.text = cmd

# check for tensor input
tensor = []  # becomes a single tensor or a tensor array
centers = []
dim_2D = None
h5_read = None
infoXml_read = None
elem_dim = None
min_bb = None
max_bb = None
elems_in_regions = [[None]]
# need this for interpolating cell to point data
reg_nodes = None
connectivity = None
nondes_void_elements = None
nondes_void_min = None
nondes_void_max = None
design_elems = None
design_elems_min = None
design_elems_max = None

h5_read = True
infoXml_read = False
if not os.path.exists(args.input):
  raise Exception('Error: file does not exist: ' + args.input)
f = h5py.File(args.input, 'r')
if args.h5_info:
  dump_h5_meta(f)
  sys.exit()
validate_region(f, args.h5_region)
if len(args.unstructured) != 0:
  nondes_centers, nondes_min, nondes_max, nondes_elem_dim, nondes_force, nondes_support, _, _, _, _ = centered_elements(f, args.h5_nondes, False, 'load', 'support')
  print('Reading elements from H5-file done ')
  dim_2D = nondes_min[2] == nondes_max[2]
  print('detected dimension ' + ('2D ' if dim_2D else '3D ') + "in non-design region")

if args.h5_region == 'all':
  centers = [[None, None, None]]
  min_bb = [numpy.Inf, numpy.Inf, numpy.Inf]
  max_bb = [-numpy.Inf, -numpy.Inf, -numpy.Inf]
  for region in f['/Mesh/Regions']:
    reg_centers, reg_min_bb, reg_max_bb, elem_dim, _, _, reg_elements, connectivity, reg_nodes, reg_nodes_map  = centered_elements(f, region)
    elems_in_regions.append(reg_elements)
    centers = numpy.concatenate((centers, reg_centers))
    min_bb = numpy.min([min_bb,reg_min_bb],0);
    max_bb = numpy.max([max_bb,reg_max_bb],0);
  elems_in_regions = elems_in_regions[1:]
  centers = centers[1:,:]
else:
  # similar to centers, but not centered
  centers, min_bb, max_bb, elem_dim, _, _, elems_in_regions, connectivity, reg_nodes, reg_nodes_map = centered_elements(f, args.h5_region)

design_elems = None
if get_MPI_rank() == 0:
  if args.h5_nondes != "None" or args.h5_nondes_void != "None":
    nondes_regs = args.h5_nondes
    # in case we have more than 1 non-design solid region
    if "," in args.h5_nondes:
      nondes_regs = args.h5_nondes.split(",")
    elif type(nondes_regs) == str:
      nondes_regs = [nondes_regs]

    if args.h5_nondes != "None":
      nondes_centers = []
      nondes_elements = []
      nondes_min = 999999
      nondes_max = -999999
      for nr in list(nondes_regs):
        tmp_nondes_centers, tmp_nondes_min, tmp_nondes_max, nondes_elem_dim, nondes_force, nondes_support, tmp_nondes_elements, _, _, _ = centered_elements(f, nr,centered=False)
        nondes_elements.extend(tmp_nondes_elements)
        nondes_min = numpy.minimum(tmp_nondes_min,nondes_min)
        nondes_max = numpy.maximum(tmp_nondes_max,nondes_max)

  # take centered values and interpolate to edges
  design_elems, design_elems_min, design_elems_max, _, _, _, design_elems, connectivity, reg_nodes, reg_nodes_map = centered_elements(f, args.h5_region,centered=True)

  print("args.h5_nondes_void:",args.h5_nondes_void)
  if args.h5_nondes_void != "None":
    if get_MPI_rank() == 0:
      nondes_void_centers, nondes_void_min, nondes_void_max, _, _, _, nondes_void_elements, _, _, _ = centered_elements(f, args.h5_nondes_void,centered=False)

dim_2D = min_bb[2] == max_bb[2]
print('detected dimension ' + ('2D' if dim_2D else '3D'))

if args.hist:
  if args.show == None:
    raise Exception('The option --hist can only be used in combination with --show.')
  design = matviz_io.read_design(f, dim_2D, args, TWO_SCALE)
  s1 = design['s1']
  s2 = design['s2']
  if dim_2D:
    vol = s1 + s2 - s1*s2
  else:
    s3 = design['s3']
    vol = pow(s1,2) + pow(s2,2) + pow(s3,2) - (pow(s1,2)*s2+pow(s1,2)*s3+pow(s2,2)*s1+pow(s2,2)*s3+pow(s3,2)*s1+pow(s3,2)*s2)/6.0 *2.0 # rough approximation
  matplotlib.pyplot.hist(vol,args.hist)
  matplotlib.pyplot.title('element volume')
  matplotlib.pyplot.xlim((min(vol),max(vol)))
  matplotlib.pyplot.show()
  for des in design:
    val = design[des]
    if val.ndim != 1:
      continue
  sys.exit()

if max_bb is None:
  max_bb = [1.0, 1.0] if dim_2D else [1.0, 1.0, 1.0]
if min_bb is None:
  min_bb = [0.0, 0.0] if dim_2D else [0.0, 0.0, 0.0]

# do we have to do 1D optimization?
if not args.target_volume:
  design = None
  if get_MPI_rank() == 0:
    design = (design_elems, design_elems_min, design_elems_max)
    if not infoXml_read:
      assert(len(design_elems) > 0)
  if args.h5_nondes != "None" or args.h5_nondes_void != "None":
    nondes_void = None
    nondes_solid = None
    if get_MPI_rank() == 0:
      if args.h5_nondes != "None":
        nondes_solid = (nondes_elements, nondes_min, nondes_max)
      if args.h5_nondes_void != "None":
        nondes_void = (nondes_void_elements, nondes_void_min, nondes_void_max)
    perform(args, h5_read, dim_2D, tensor, centers, None, nondes=(nondes_solid,nondes_void,design), min_bb=min_bb, max_bb=max_bb, elems_in_regions=elems_in_regions)
  else:
    perform(args, h5_read, dim_2D, tensor, centers, None, nondes=(None,None,design), min_bb=min_bb, max_bb=max_bb, elems_in_regions=elems_in_regions)
else:
  perform(args, h5_read, dim_2D, tensor, centers, min_bb=min_bb, max_bb=max_bb, elems_in_regions=elems_in_regions)

if args.info:
  matviz_io.write_info_xml(args.info)
