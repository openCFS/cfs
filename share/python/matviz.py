#!/usr/bin/env python
# sometimes h5py must be import first
from hdf5_tools import *
from matviz_vtk import *
from matviz_unstructured_mesh import *
from matviz_2d  import *
from matviz_streamline import *
from hdf5_tools import *
from mesh_tool import *
from optimization_tools import *
import argparse
import sys
import os
import io
import xml.etree.ElementTree
import xml.dom.minidom
from cfs_utils import *
import matviz_3d_ortho # two scale 3d ortho stuff
try:
  from mpi4py import MPI
except:
  print("WARNING:Could not load mpi4py!")


TWO_SCALE = (
  "hom_rot_cross",
  "hom_sheared_rot_cross",
  "hom_frame",
  "hom_framed_cross",
  "hom_rect",
  "hom_ortho_3d"
  )

## reads design_stiff*, design_shear* and design_rotAngle* for 2D and 3D. Fills other stuff by defaults 

# considers density read
# @param angle array of anglex, angley, anglez
# @return dict with design variables (s1, s2, s3, sh1, angle)  
def read_design(hdf_file, dim_2D, args):
  res = dict()
  # rot means, that we only show rotation according to rotAngle, e.g. for piezoelectric polarization
  sh1 = numpy.ones((len(centers),1)) * .5 # fix for no shearing 
  if args.parametrization == 'hom_rect':
    number_of_parameters = {
      "hom_rot_cross" : 2,
      "hom_sheared_rot_cross" : 3,
      "hom_frame" : 2,
      "hom_framed_cross" : 4,
      "hom_rect" : 3,
      "hom_ortho_3d": 3}
    if args.show in number_of_parameters:
      try:
        res['microparams'] = [get_element(f, "design_microparam{}_{}".format(i + 1, args.hom_access), args.h5_region, args.h5_step) for i in range(number_of_parameters[args.show])]
      except:
        res['s1'] = get_element(f, "design_stiff1_" + args.hom_access, args.h5_region, args.h5_step) if args.show != "rot" else numpy.ones((len(centers),1)) * .1 
        res['s2'] = get_element(f, "design_stiff2_" + args.hom_access, args.h5_region, args.h5_step) if args.show != "rot" else numpy.ones((len(centers),1)) * .1
        res['s3'] = numpy.ones((len(centers),1)) * .1 if dim_2D or args.show == "rot" else get_element(f, "design_stiff3_" + args.hom_access, args.h5_region, args.h5_step)
        res['sh1'] = get_element(f, "design_shear1_" + args.hom_access, args.h5_region, args.h5_step) if (args.show == "hom_sheared_rot_cross" or args.show == "hom_cross_bar") else sh1
  elif args.parametrization == 'trans-iso':
    s1 = get_element(f, "design_emodul-iso_" + args.hom_access, args.h5_region, args.h5_step)
    s2 = get_element(f, "design_emodul_" + args.hom_access, args.h5_region, args.h5_step)
    try:
      theta = get_element(f, "design_poisson_" + args.hom_access, args.h5_region, args.h5_step)
    except:
      theta = 0.0
      print('could not read theta (design_poisson_' + args.hom_access + '), setting to ' + str(theta))
    m = 2.0 * numpy.max([numpy.max(s1), numpy.max(s2)])
    s1 *= 1 / (m * (1 - theta))
    s2 *= 1 / (m * (1 - theta))
    res['s1'] = s1
    res['s2'] = s2
    res['s3'] = numpy.ones((len(centers), 1)) * .1  # fix for 3D
  elif args.parametrization == 'ortho':
    t11 = get_element(f, "design_tensor11_" + args.hom_access, args.h5_region, args.h5_step)
    t12 = get_element(f, "design_tensor12_" + args.hom_access, args.h5_region, args.h5_step)
    t22 = get_element(f, "design_tensor22_" + args.hom_access, args.h5_region, args.h5_step)
    t33 = get_element(f, "design_tensor33_" + args.hom_access, args.h5_region, args.h5_step)
    s1 = t11*t11+t12*t12
    s2 = t12*t12+t22*t22
    m = 2.0*numpy.max([numpy.max(s1), numpy.max(s2)])
    s1 *= 1/m
    s2 *= 1/m
    res['s1'] = s1
    res['s2'] = s2
    res['s3'] = numpy.ones((len(centers),1)) * .1 # fix for 3D
  elif args.parametrization == "hom_iso":
    # isotropic homogenized basecell e.g. lufo fuller or V7 base cell
    res['s1'] = get_element(f, "design_stiff1_" + args.hom_access, args.h5_region, args.h5_step) if args.show != "rot" else numpy.ones((len(centers),1)) * .1 
    res['s2'] = res['s1']
    res['s3'] = res['s1']
  elif args.parametrization == 'simp':
    if args.h5_region == 'all':
      s1 = [[None]]
      for region in f['/Mesh/Regions']:
        s = get_element(f, "physicalPseudoDensity", region, args.h5_step)
        s1 = numpy.concatenate((s1,s))
      res['s1'] = s1[1:]
      shape(res['s1'])
    else:
      res['s1'] = get_element(f, "physicalPseudoDensity", args.h5_region, args.h5_step)
    res['s2'] = res['s1']
    res['s3'] = res['s1']
    res['angle'] = numpy.zeros(((len(res['s1']), 3)))
    return res
  if has_element(hdf_file, "design_density_" + args.hom_access):
    print("args.h5_step:" + str(args.h5_step))
    rho = get_element(f, "design_density_" + args.hom_access, args.h5_region, args.h5_step)
    rho = pow(rho, float(args.penalty))
    res['s1'] *= rho
    res['s2'] *= rho
    res['s3'] *= rho
    print("scale stiffness values by design_density_" + args.hom_access + " with average value " + str(numpy.mean(rho)) + " and penalty " + str(args.penalty))
  
  angle = numpy.zeros(((len(list(res.values())[0]), 3)))
  
  if args.show == "hom_rot_cross" or args.show == "hom_sheared_rot_cross" or args.show == "rot" or args.show == "stream" or args.show == 'hom_rect_mod':
    try:
      if dim_2D:
      	try:
      	  angle[:,0] = get_element(f, "design_rotAngle_" + args.hom_access, args.h5_region, args.h5_step)[:,0]
      	except:
      	  print('could not read design_rotAngle_' + args.hom_access + ', trying design_rotAngle_plain')
      	  angle[:,0] = get_element(f, "design_rotAngle_plain", args.h5_region, args.h5_step)[:,0]
      else:
        angle[:, 0] = get_element(f, "design_rotAngleX_" + args.hom_access, args.h5_region, args.h5_step)[:, 0]
        angle[:, 1] = get_element(f, "design_rotAngleY_" + args.hom_access, args.h5_region, args.h5_step)[:, 0]
        angle[:, 2] = get_element(f, "design_rotAngleZ_" + args.hom_access, args.h5_region, args.h5_step)[:, 0]
    except Exception as e:
      print('could not read angle, ignore it: ', e)
  res['angle'] = angle
  
  return res

# # show or write either Image or polydata
# @param viz either Image or polydata
# @save filename for output
# @return the volume fraction if determined or None
def show_or_write(viz, args):
  assert(viz is not None)
  volume = None
  global info
  if isinstance(viz, Image.Image):
    # print 'array' + str(numpy.array(viz))
    print('average ' + str(numpy.average(numpy.array(viz))))
    volume = 1 - numpy.average(numpy.array(viz)) / 255
    print('volume fraction from image : ' + str(volume))
    if info != None:
      vol = xml.etree.ElementTree.SubElement(info, "volume")
      vol.set("imageMaterial", str(volume))  
    
    if args.save:
      viz.save(args.save)
    else:
      viz.show()
  elif isinstance(viz, tuple):
    fig = viz[0]
    sub = viz[1]
    if args.save:
      extent = sub.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
      print('write file: ' + args.save)
      fig.savefig(args.save, bbox_inches=extent)
      if args.save.split('.')[-1] != 'pdf':      
        # I was not able to render a memory image first, make an array out of the data and determine the grayness
        # So read again from file :( 
        tmp = Image.open(args.save).convert('L')  # make gray, otherwise data has the dimension x*y*4 (rgm + alpha)
        dat = numpy.array(tmp)
        volume = len(numpy.where(dat.reshape(dat.size, 1) < 128)[0]) / float(dat.size)  # cont fields below 128 which is become black with a thrshold of 0.5
        print('volume fraction from image : ' + str(volume))
        if info != None:
          vol = xml.etree.ElementTree.SubElement(info, "volume")
          vol.set("imageMaterial", str(volume))  
    else:
      matplotlib.pyplot.show()
      ax = fig.add_axes([0, 0, 1, 1])
      fig.show()  # Jannis: this is a temporary workaround as matplotlib.pyplot.show() does nothing for me

    matplotlib.pyplot.close(fig)

  else:
    show_write_vtk(viz, args.res, args.save, camera_settings=args.cam)
    
  return volume  

def plot_angle_data(file, angle, data):
  plot = open(file, "w")
  plot.write('#angle\tdata\n')
  for i in range(len(angle[0])):  # assume a single tensor, take the first element
    plot.write(str(angle[0][i]) + '\t' + str(data[0][i]) + '\n')
  plot.close()  


# we extract the main processing such that we can run it within a loop to find the optimal scaling
# @param force_scale overwrites args.scale
# @param min_bb/max_bb: min/max coordinates of bounding box
# @return volume if calculated (e.g. via --save a pixel image) otherwise None
def perform(args, h5_read, dim_2D, tensor, centers, aux_code, force_scale=None, nondes = None, min_bb = None, max_bb = None, elems_in_regions = None):
  volume = None  # might ne set
  
  scale = force_scale if force_scale else args.scale
  
  coords = (centers, min_bb, max_bb, elem_dim)
  
  # perform 2D and 3D from file
  if h5_read or dim_2D:
    # either Image or polydata  
    viz = None
    if args.show in ("hom_rect", "hom_rot_cross", "hom_sheared_rot_cross", "hom_cross_bar", "hom_frame", "hom_framed_cross", "rot", "stream", "hom_rect_mod", "simp","hom_ortho_3d"):
      if args.show == "simp":
        args.parametrization = "simp"
      if h5_read:
        design = read_design(f, dim_2D, args)
      else:
        angle, s1, s2, coords = read_stiff_angle_matlab(args.input)

      if dim_2D and args.show in TWO_SCALE:
        print("Volume for regular grid: " + str(calc_volume(design['s1'], design['s2'])))
      # add angle bias, e.g. by 90 deg to correct thomas
      design['angle'] += args.angle_bias * numpy.pi / 180
      # scale angle, e.g  by -1 to correct for current standard 2D rotation direction (this is not the mathematical direction! FIXME if needed)
      design['angle'] *= -1.0
      if args.angle_factor != 1.0:
        print('scale angle by ' + str(args.angle_factor))
        design['angle'] *= args.angle_factor
      
      for key, value in design.items():
        print('unscaled {:s} in [{:f}:{:f}]'.format(key, numpy.min(value), numpy.max(value)))
  
      # viz is either Image or polydata
      if dim_2D and not args.show == 'simp':
        if args.show == "hom_rect": 
          if args.hom_grad == 'none':
            viz = show_frame(coords, design, args.hom_dir, args.res, args.scale)
          else:
            viz = show_frame_grad(coords, design, args.hom_grad, args.hom_dir, args.res)
        elif args.show == "hom_rect_mod":
          design['angle'] = design['angle'][:,0]
          viz = show_modified_frame(coords, design, "all", args.res, scale, args.color, args.save)
        elif args.show == "hom_rot_cross" or args.show == "rot":
          # add optional angle bias
          # print 'change angle'
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
          viz = show_framed_cross(coords, microparams['microparams'], args.res, args.scale, args.color, args.save)
        elif args.show == "hom_cross_bar":
          design['angle'] = numpy.zeros((len(s1),1))
          viz = show_cross_bar(coords, design, args.hom_dir, args.res, args.scale, args.color, args.save)
        elif args.show == "stream":
            samples = args.hom_samples.split(',')
            design['angle'] = design['angle'][:,0]
            viz = show_streamline(coords, design, args.hom_dir, scale, args.minimal, args.stream_style, args.stream_step, float(samples[0]), float(samples[1]), args.stream_max_traces_per_cell, args.res, args.save != None, info, args.stream_force)            
        else:
          assert(False)
      # the 3D VTK stuff      
      else:
        if args.show == "rot":
          s1 = ones((len(s1), 1)) * 0.25
          s2 = ones((len(s1), 1)) * 0.25
          s3 = ones((len(s1), 1)) * 0.25
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
          if args.hom_grad == 'none':
            print('for hom_rect in 3D with hom_samples you need to specify hom_grad')
            exit()
          else:
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
                if args.show == 'simp':
                  me = create_validation_mesh(coords, nondes[0], design, None, args.hom_grad, args.hom_dir, scale,n_f,valid_position,valid_ring_position, args.type,args.thres,csize,args.show)
                else:
                  if args.type == "apod6" or args.type == "robot": 
                    me = create_validation_mesh(coords, nondes[0], design, None, args.hom_grad, args.hom_dir, scale,n_f,valid_position,valid_ring_position, args.type,args.thres,csize)
                write_gid_mesh(me, args.mesh+".mesh")
                exit()  
              else:
                if args.show == "hom_rot_cross":
                  viz = create_3d_cross_ip(coords, design, args.hom_samples, args.hom_grad, scale, valid_position, args.thres,csize)   
                elif args.show == "hom_rect":
                  viz = create_3d_frame_ip(coords, design, args.hom_samples, args.hom_grad, scale, valid_position, args.thres)
                elif args.show == "hom_ortho_3d":
                  print("Ohoh, I shouldn't be here...")
                  sys.exit()
            else:
              tmp = args.hom_samples.split(',')
              if len(tmp) == 1:
                samples = [int(tmp[0]),int(tmp[0]),int(tmp[0])]
              else:
                samples = [int(tmp[0]),int(tmp[1]),int(tmp[2])]
              if args.show == "hom_ortho_3d" or args.mesh:
#                 name = "interpretation_ortho_3d_box_varel_" + str(samples[0]) + "_" + str(samples[1]) + "_" + str(samples[2]) + "_bc_res_" + str(args.bc_res) + ".stl"
                name = args.save + ".stl" if args.save is not None else "interpretation_mc"
                # region_map: maps local node id (list idx) in given design region (e.g. "mech") to global node id in all regions
                reg_info = {"nodes":reg_nodes, "elements":elems_in_regions, "connectivity":connectivity, "region_map":reg_nodes_map}
                if nondes:
                  viz = matviz_3d_ortho.create_3d_interpretation_ortho(args, reg_info, coords, min_bb, max_bb, design, scale, samples, args.hom_grad,nondes=nondes)
                else:
                  viz = matviz_3d_ortho.create_3d_interpretation_ortho(args, reg_info, coords, min_bb, max_bb, design, scale, samples, args.hom_grad,nondes=None)
                me = None                
                if args.save:
                  if args.save.endswith(".vtp"):
                    name = args.save[:-4]+".stl"
                  
                  matviz_vtk.write_stl(viz, name)
#                 if args.type == "box_varel" or args.type == "ppbox" or args.mesh:    
                if args.mesh:
                  if not args.save: # write surface mesh in case we haven't done it before
                    matviz_vtk.write_stl(viz, name)
                  
                  create_volume_mesh_with_gmsh(name)
                  if not args.save:
                    viz = None # avoid showing or writing vtp file
                viz = None    
              else:
                viz = create_3d_cross_ip(coords, design, samples, args.hom_grad, scale, valid_position, args.thres)
        else:  # no sample
          if args.show == 'simp':
            viz = create_block(coords, design, scale, args.thres, elems_in_regions)
          elif args.hom_grad == 'none':
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
            print("hom_rect in 3D only for '--hom_grad none' implemented")
            exit()
              
    # no hom_rect stuff but orientational stiffness
    else:
      if args.tensor == 'mechTensor':
        print("Input data is read as " + args.notation)
      if h5_read:
        tensor = get_element(f, args.tensor, args.h5_region, args.h5_step)
      angle, data = perform_rotations(tensor, args.notation, int(args.sampling), args.tensor, args.show)
      
      if args.plot != None:
        plot_angle_data(args.plot, angle, data)
        # quit!! such we can check for viz!!
        exit() 
      else:
        viz = orientational_stiffness(coords, angle, data, args.res, scale)
    if (MPI.COMM_WORLD.Get_rank()==0):  
	    if viz is None:
	      print('Error: no visualization calculated!')
	    else:
	      volume = show_or_write(viz, args)
  
  # not from file and not 2D -> this is the single tensor with optional planes 
  else:
    angle, data, aux = perform_cfs_rotation(tensor, int(args.sampling), aux_code)
    angle_max = angle[numpy.argmax(numpy.abs(data))]
    angle_min = angle[numpy.argmin(numpy.abs(data))]

    print(" largest e-modulus: {:>13.6e}".format(numpy.max(numpy.abs(data))) + "  in direction " + str(to_vector(angle_max)))
    print("smallest e-modulus: {:>13.6e}".format(numpy.min(numpy.abs(data))) + "  in direction " + str(to_vector(angle_min)))
    if len(aux) > 0:
      print("largest " + args.show + ": " + str(numpy.max(aux)) + " smallest " + args.show + ": " + str(numpy.min(aux)))
    
    poly = create_vtk_poly_data(angle, data if args.show == None else aux)
    
    actors = []
    if not args.symmetries == "default":
      tmp = []
      if args.symmetries_mode == "minima":   
        tmp = find_minima(angle, aux)
      else:
        for i in range(len(angle)):
          tmp.append((angle[i], aux[i]))
      # sort by value  
      tmp2 = sorted(tmp, key=lambda entry: entry[1])
      # max
      tmp3 = tmp2[0:int(args.symmetries_max)]
      # handle threshold
      mins = [] 
      for i in range(len(tmp3)):
        if tmp3[i][1] <= float(args.symmetries_threshold):
          mins.append(tmp3[i]) 
      actors = create_symmetry_planes(mins, 1.2 * numpy.max(data if args.show == None else aux), not args.symmetries_planes == "false")
  
    show_write_vtk(poly, args.res, args.save, actors, args.axes, camera_settings=args.cam)  
  
  return volume

parser = argparse.ArgumentParser()
parser.add_argument("input", help="a cfs++ h5 file or a tensor \"[e11, ...]\" with 11/22/33/32/31/21 for 2D and 11/12/22/13/23/... for 3D or a '.info.xml' file or a .mat file including a matrix from matlab (2sc)")
parser.add_argument("--h5_step", help="step number, too high is last (default '9999')", default=9999, type=int)
parser.add_argument("--h5_region", help="design region name (default 'mech')", default="mech")
parser.add_argument("--h5_nondes", help="non-design (solid) region name (default 'non-design')", default="None")
parser.add_argument("--h5_nondes_void", help="non-design (void) region name (default 'non-design')", default="None")
# parser.add_argument('--h5_nonreg', action='store_true', help='assume the h5 file to be nonregular')
parser.add_argument('--h5_info', action='store_true', help='dump some meta data information about the h5 file')
parser.add_argument("--tensor", help="tensor name: 'mechTensor', 'piezoTensor, 'elecTensor'", default="mechTensor")
parser.add_argument("--scale", help="manual scaling factor", default=-1.0, type=float)
parser.add_argument("--target_volume", help="find optimal scaling. Makes only sense for streamline", type=float)
parser.add_argument("--res", help="x-resolution (default 1000)", default=800, type=int)
parser.add_argument("--sampling", help="sampling rate (default 180", default=180, type=float)
parser.add_argument("--show", help="mode within boebbale, hom_rect or streamline", choices=['ortho_norm', 'mono_norm', 'ortho_err', 'hom_rect', 'hom_rot_cross', 'hom_sheared_rot_cross', 'hom_frame', 'hom_framed_cross', 'hom_cross_bar', 'rot', 'stream', 'hom_rect_mod', 'simp', 'hom_ortho_3d'])
parser.add_argument("--axes", help="show axes", action='store_true', default=False)
parser.add_argument("--notation", help="mandel | voigt (default 'voigt')", default="voigt")
parser.add_argument("--symmetries", help="same options as for shows", default="default")
parser.add_argument("--symmetries_max", help="maximum number of symmetries (default 999)", default=999)
parser.add_argument("--symmetries_threshold", help="threshold value for symmetries (default 9999)", default=9999)
parser.add_argument("--symmetries_mode", help="'minima' or 'all' subject to max and threshold (default 'minima'", default="minima")
parser.add_argument("--symmetries_planes", help="'true' or 'false' for 3D also show planes to normals (default 'false')", default="false")
parser.add_argument("--hom_access", help="the 'plain ' or 'smart' hom values (default 'smart')", default="smart")
parser.add_argument("--hom_grad", help="interpolation of design: 'none', 'nearest', linear', 'cubic' (default 'linear')", default="linear", choices=['none', 'nearest', 'linear', 'cubic'])
parser.add_argument("--hom_dir", help="visualization of stiffness directions (default 'all')", default="all", choices=['all', 'horizontal', 'vertical', 'sagittal'])
parser.add_argument("--angle_factor", help="factor for angle. -1.0 turns, 0.0 disables angles", default=1.0, type=float)
parser.add_argument("--angle_bias", help="bias for the angle in deg. 90 switches s1 and s2", default=0.0, type=float)
parser.add_argument("--hom_samples", help="activates interpolation and the value gives samples in x,y,z direction")
parser.add_argument("--cell_size", help="cell size in [mm] in x,y,z direction")
parser.add_argument("--stream_style", help="select visualization", choices=['line', 'thick'], default='thick')
parser.add_argument("--stream_step", help="step length for ODE integration per macro cell", type=float, default=0.2)
parser.add_argument("--stream_max_traces_per_cell", help="maximum number of traces such that we may start a trace (>= 1)", type=int, default=1)
parser.add_argument("--stream_ode", help="method to solve the ODE", default="euler", choices=['euler', 'midpoint'])
parser.add_argument("--stream_force", help="force streamlines for special cases", choices=['right_lower', 'rhombus'])
parser.add_argument("--minimal", help="minimal stiffness to be drawn, will be scaled", type=float, default=0.0)
parser.add_argument("--parametrization", help="parametrization of the stiffness tensor", default="hom_rect", choices=['simp','hom_rect', 'trans-iso', 'ortho'])
parser.add_argument("--save", help="save 'image.png' (pixel), 'image.pdf' (vector) or VTK Poly Data file 'file.vtp'")
parser.add_argument("--cam", help="set camera (7 space-separated floats): position, focal point, roll", nargs=7, type=float, default=None)
parser.add_argument("--plot", help="for single tensors: creates gnuplot file instead of image")
parser.add_argument("--penalty", help="penalty parameter for SIMP (default 5)", default=5.0)
parser.add_argument("--color", help="only for hom_rot_cross: 'black' or colormap from http://matplotlib.org/examples/color/colormaps_reference.html, default='gray'", default="grayscale")
parser.add_argument("--info", help="creates a xml file of given name with additional information")
parser.add_argument("--unstructured", help="number of structured elements per coordinate as list nx,ny,nz", default="")
parser.add_argument("--nodefile", help="name of the design to node file", default="")
parser.add_argument("--thres", help="threshold value for 3D VTK plot", type=float, default=0.0)
parser.add_argument("--mesh", help="create 3D mesh from optimized 2-scale result for validation", default="")
parser.add_argument("--nf", help="requires --mesh, number of fine elements in x,y,z direction")
parser.add_argument("--type", help="type of 3D object for 2-scale visualization",choices=['apod6', 'robot','ppbox','box_varel'])
# 3d ortho basecell stuff
parser.add_argument("--bc_res", help="resolution of voxelized ortho basecell", type=int)
parser.add_argument("--bc_interpolation", help="interpolation type for ortho basecell (linear or heaviside)",choices=['linear', 'heaviside'])
parser.add_argument("--bc_beta", help="for heaviside interpolation (default 7.0)", type=float,default=7)
parser.add_argument("--bc_eta", help="for heaviside interpolation (default 0.5)", type=float,default=0.5)
parser.add_argument("--bc_bend", help="bending of spline (default 0.5)", type=float,default=0.5)
parser.add_argument("--bc_smooth", help="number auf Taubin smoothing steps", type=int,default=40)
parser.add_argument("--bc_thresh", help="lower and upper threshold (diameter) for ortho basecell, e.g. 1e-9,0.94")
# print sys.argv

args = parser.parse_args()

# check ans postproc arguments
if not args.symmetries == "default" and not args.show == None and not args.symmetries == args.show:
  print("'show' and 'symmetries' do not match")
  sys.exit(1)
aux_code = args.show if not args.show == None else args.symmetries  # might still be default

# check 3d ortho basecell parameters
if args.show == "hom_ortho_3d":
  if not args.bc_res:
    print("bc_res required")
    sys.exit(1)
  if not args.bc_bend:
    print("bc_bend parameter required")
    sys.exit(1) 
  if not args.bc_interpolation:
    print("interpolation type (linear or heaviside) required")
    sys.exit(1)
  elif not args.bc_beta or not args.bc_eta:
    print("beta and eta values required for heaviside interpolation")  
    sys.exit(1)
if args.type == "box_varel" or args.type == "ppbox":
  # in this case everything belongs to design domain
  args.h5_nondes = "None"

# in this global variable we can store meta-information to be exported as xml file 
info = None
if args.info is not None:
  info = xml.etree.ElementTree.Element("matviz")
  cmd = sys.argv[0].split('/')[-1]
  for i in range(1, len(sys.argv)):
    cmd += ' ' + sys.argv[i]
  doc = xml.etree.ElementTree.SubElement(info, "commandLine")
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

# check if we read data from command line instead from an h5 file or a info.xml was given
if args.input.startswith('[') or args.input.endswith(".info.xml") or args.input.endswith(".mat"):
    
  h5_read = False
  dim_2D = None
  input = None
  
  if args.input.startswith('['):
    input = eval(args.input)
    print(input)
    if len(input) != 21 and len(input) != 6:
      print("the input has " + str(len(input)) + " coefficients but requires 6 (2D) or 21 (3D)")
      sys.exit(1)
  
    dim_2D = len(input) != 21
  elif args.input.endswith(".info.xml"):
    xml = open_xml(args.input)
    dim = xpath(xml, "//domain/@dimensions")
    matrix = xpath(xml, "//iteration[last()]/homogenizedTensor/tensor/real/text()") # "text()" must be added due to lxml, otherwise matrix is just a string <Element real ...>
    res = list(map(float, matrix.split())) # convert list with string elements to list with float elements
    res = np.asarray(res)            # convert list to array
    if dim == '2':
      res = res.reshape(3,3)         # reshaping array
      input = [res[0][0],res[0][1],res[1][1],res[0][2],res[1][2],res[2][2]]
    else:
      assert(dim == '3')
      res = res.reshape(6,6)         # reshaping array
      input = [res[0][0],res[0][1],res[1][1],res[0][2],res[1][2],res[2][2],res[0][3],res[1][3],res[2][3],res[3][3],res[0][4],res[1][4],res[2][4],res[3][4],res[4][4],res[0][5],res[1][5],res[2][5],res[3][5],res[4][5],res[5][5]]
  else:
    #data from matlab file
    assert(args.input.endswith(".mat"))
    input = args.input
    args.tensor = 'matlab'
  if not args.tensor == 'matlab' and args.tensor == 'mechTensor':
    print(input)
    tensor = to_mech_tensor(input)
    tensor = HillMandel2Voigt(tensor) if args.notation == "mandel" else tensor
    print("Voigt notation of input tensor:")
  if not args.tensor == 'matlab' and args.tensor == 'piezoTensor':
    tensor = to_piezo_tensor(input)
  if not args.tensor == 'matlab':
    dump_tensor(tensor)

  if not args.tensor == 'matlab' and (len(tensor) == 3 or len(tensor) == 2):
    assert((len(tensor) == 3 and args.tensor == 'mechTensor') or (len(tensor) != 3 and args.tensor != 'mechTensor'))
    vec = None
    
    if len(tensor) == 3:
      # convert back to hill mandel as we simulate a h5 file
      vec = to_mech_vector(Voigt2HillMandel(tensor), as_array=True)
    else:
      vec = to_piezo_vector(tensor, as_array=True)
    tensor = []
    tensor.append(vec)
    centers.append([0.5, 0.5, 0.0])

# read 2D or 3D data from h5 file     
else:   
  h5_read = True
  infoXml_read = False
  if not os.path.exists(args.input):
    print('Error: file does not exist: ' + args.input)
    sys.exit(1)
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
  if args.h5_nondes != "None" or args.h5_nondes_void != "None":
    if (MPI.COMM_WORLD.Get_rank()==0):
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
      if (MPI.COMM_WORLD.Get_rank()==0): 
        nondes_void_centers, nondes_void_min, nondes_void_max, _, _, _, nondes_void_elements, _, _, _ = centered_elements(f, args.h5_nondes_void,centered=False)
         
  dim_2D = min_bb[2] == max_bb[2]
  print('detected dimension ' + ('2D' if dim_2D else '3D'))
  
# do we have to do 1D optimization? 
if not args.target_volume:
  if args.mesh and (args.h5_nondes != "None" or args.h5_nondes_void != "None"):
    nondes_void = None
    nondes_solid = None
    design = None
    if (MPI.COMM_WORLD.Get_rank()==0):
      if args.h5_nondes != "None":
        nondes_solid = (nondes_elements, nondes_min, nondes_max)
      if args.h5_nondes_void != "None":   
        nondes_void = (nondes_void_elements, nondes_void_min, nondes_void_max)
      design = (design_elems, design_elems_min, design_elems_max)
      
    perform(args, h5_read, dim_2D, tensor, centers, aux_code,None,nondes=(nondes_solid,nondes_void,design),min_bb=min_bb,max_bb=max_bb,elems_in_regions=elems_in_regions)
  else:
    perform(args, h5_read, dim_2D, tensor, centers, aux_code,min_bb=min_bb,max_bb=max_bb,elems_in_regions=elems_in_regions)
else:
  if args.scale > 0:
    print("Error: don't give --scale and --target_volume concurrently!")
    sys.exit(1) 
    # coords for non-design region
    if args.unstructured:
      nondes_coords = (nondes_centers, nondes_min, nondes_max, nondes_elem_dim)
  # unfortunately we are not necessary convex! therefore we may not do bisection
  lower = 1e-4
  upper = 2.0
  level = 0  
  best_err = 999
  err = 999
  best_s = -1
  while upper - lower > 0.000001 and abs(best_err) > 0.0001:
    for s in numpy.arange(lower, upper, (upper - lower) / 5.):
      vol = perform(args, h5_read, dim_2D, tensor, centers, aux_code, s)
      err = abs(vol - args.target_volume)
      if vol == None:
        sys.exit(1)
  
      tv = xml.etree.ElementTree.SubElement(info, "target_volume")
      tv.set("target", str(args.target_volume))
      tv.set("current", str(vol))
      tv.set("scale_argument", str(s))
      tv.set("err", str(err))
       
      tvs = xml.etree.ElementTree.SubElement(tv, "search")
      tvs.set("level", str(level))
      tvs.set("lower", str(lower))
      tvs.set("upper", str(upper))
      tvs.set("delta", str(upper - lower))
      tvs.set("best_err", str(best_err))
      if(err < best_err):
        best_err = err
        best_s = s
      msg =  ">>>> target_volume: level=" + str(level) + " lower= " + str(lower) + " upper= " + str(upper) + " delta= " + str(upper - lower) + " best=" + str(best_err) + " test=" + str(s)
      msg += " -> " + str(vol) + " err=" + str(err) + (" *" if err == best_err else "")
      print(msg)
            
    lower = numpy.max((1e-4, best_s - 0.5 * (upper - lower) / 5.)) 
    upper = numpy.min((2.0, best_s + 0.5 * (upper - lower) / 5.))
    level += 1

  # the last result does not mean to be the best result
  if err != best_err:
    vol = perform(args, h5_read, dim_2D, tensor, centers, aux_code, best_s)
    print("!!!!! best_target_volume: best_scale=" + str(best_s) + " -> " + str(vol) + " err=" + str(abs(vol - args.target_volume)))
    tv = xml.etree.ElementTree.SubElement(info, "best_target_volume")
    tv.set("target", str(args.target_volume))
    tv.set("current", str(vol))
    tv.set("scale_argument", str(best_s))
    tv.set("err", str(abs(vol - args.target_volume)))



if args.info:
  print('write info xml file: ' + args.info)
  out = open(args.info, "w")
  out.write(xml.dom.minidom.parseString(xml.etree.ElementTree.tostring(info)).toprettyxml())
  out.close()
  
