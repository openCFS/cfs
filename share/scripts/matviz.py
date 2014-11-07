#!/usr/bin/env python
from matviz_vtk import *
from matviz_unstructured_mesh import *
from matviz_2d  import *
from matviz_streamline import *
from hdf5_tools import *
from mesh_tool import *
import argparse
import sys
import os
import StringIO
import xml.etree.ElementTree
import xml.dom.minidom

# # reads design_stiff* and design_rotAngle* for 2D and 3D. Fills other stuff by defaults
# considers density read
# @param angle array of anglex, angley, anglez
# @return s1, s2, s3, angle 
def read_stiff_angle(hdf_file, dim_2D, args):
  # rot means, that we only show rotation according to rotAngle, e.g. for piezoelectric polarization
  if args.parametrization == 'hom_rect':
    s1 = get_element(f, "design_stiff1_" + args.hom_access, args.h5_region, args.h5_step) if args.show <> "rot" else numpy.ones((len(centers), 1)) * .1 
    s2 = get_element(f, "design_stiff2_" + args.hom_access, args.h5_region, args.h5_step) if args.show <> "rot" else numpy.ones((len(centers), 1)) * .1
    s3 = numpy.ones((len(centers), 1)) * .1 if dim_2D or args.show == "rot" else get_element(f, "design_stiff3_" + args.hom_access, args.h5_region, args.h5_step)
  elif args.parametrization == 'trans-iso':
    s1 = get_element(f, "design_emodul-iso_" + args.hom_access, args.h5_region, args.h5_step)
    s2 = get_element(f, "design_emodul_" + args.hom_access, args.h5_region, args.h5_step)
    try:
      theta = get_element(f, "design_poisson_" + args.hom_access, args.h5_region, args.h5_step)
    except:
      theta = 0.0
      print 'could not read theta (design_poisson_' + args.hom_access + '), setting to ' + str(theta)
    m = 2.0 * numpy.max([numpy.max(s1), numpy.max(s2)])
    s1 *= 1 / (m * (1 - theta))
    s2 *= 1 / (m * (1 - theta))
    s3 = numpy.ones((len(centers), 1)) * .1  # fix for 3D
  elif args.parametrization == 'ortho':
    t11 = get_element(f, "design_tensor11_" + args.hom_access, args.h5_region, args.h5_step)
    t12 = get_element(f, "design_tensor12_" + args.hom_access, args.h5_region, args.h5_step)
    t22 = get_element(f, "design_tensor22_" + args.hom_access, args.h5_region, args.h5_step)
    t33 = get_element(f, "design_tensor33_" + args.hom_access, args.h5_region, args.h5_step)
    s1 = t11 * t11 + t12 * t12
    s2 = t12 * t12 + t22 * t22
    m = 2.0 * numpy.max([numpy.max(s1), numpy.max(s2)])
    s1 *= 1 / m
    s2 *= 1 / m
    s3 = numpy.ones((len(centers), 1)) * .1  # fix for 3D
    
  if has_element(hdf_file, "design_density_" + args.hom_access):
    rho = get_element(f, "design_density_" + args.hom_access, args.h5_region, args.h5_step)
    rho = pow(rho, args.penalty)
    s1 *= rho
    s2 *= rho
    s3 *= rho
    print "scale stiffness values by design_density_" + args.hom_access + " with average value " + str(numpy.mean(rho))  
  
  angle = numpy.zeros(((len(s1), 3)))
  
  if args.show == "hom_rot_cross" or args.show == "rot" or args.show == "stream":
    try:
      if dim_2D:
      	try:
      	  angle[:, 0] = get_element(f, "design_rotAngle_" + args.hom_access, args.h5_region, args.h5_step)[:, 0]
      	except:
      	  print 'could not read design_rotAngle_' + args.hom_access + ', trying design_rotAngle_plain'
      	  angle[:, 0] = get_element(f, "design_rotAngle_plain", args.h5_region, args.h5_step)[:, 0]
      else:
        angle[:, 0] = get_element(f, "design_rotAngleX_" + args.hom_access, args.h5_region, args.h5_step)[:, 0]
        angle[:, 1] = get_element(f, "design_rotAngleY_" + args.hom_access, args.h5_region, args.h5_step)[:, 0]
        angle[:, 2] = get_element(f, "design_rotAngleZ_" + args.hom_access, args.h5_region, args.h5_step)[:, 0]
    except Exception, e:
      print 'could not read angle, ignore it: ', e   
  return s1, s2, s3, angle

# # show or write either Image or polydata
# @param viz eithe Image or polydata
# @save filename for output
def show_or_write(viz, args):
  assert(viz <> None)
  global info
  if isinstance(viz, Image.Image):
    frac = 1 - numpy.average(numpy.array(viz)) / 255
    print 'volume fraction from image : ' + str(frac)
    if info <> None:
      vol = xml.etree.ElementTree.SubElement(info, "volume")
      vol.set("imageMaterial", str(frac))  
    
    if args.save:
      viz.save(args.save)
    else:
      viz.show()
  elif isinstance(viz, tuple):
    fig = viz[0]
    sub = viz[1]
    
    if args.save:
      extent = sub.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
      print 'write file: ' + args.save
      fig.savefig(args.save, bbox_inches=extent)

      if args.save.split('.')[-1] <> 'pdf':      
        # I war not able to  render a memory image first, make an array out of the data and determine the grayness
        # So read again from file :( 
        tmp = Image.open(args.save)
        frac = 1 - numpy.average(numpy.array(tmp)) / 255
        print 'volume fraction from image : ' + str(frac)
        if info <> None:
          vol = xml.etree.ElementTree.SubElement(info, "volume")
          vol.set("imageMaterial", str(frac))  
    else:
      matplotlib.pyplot.show()
  else:
    show_write_vtk(viz, args.res, args.save)

def plot_angle_data(file, angle, data):
  plot = open(file, "w")
  plot.write('#angle\tdata\n')
  for i in range(len(angle[0])):  # assume a single tensor, take the first element
    plot.write(str(angle[0][i]) + '\t' + str(data[0][i]) + '\n')
  plot.close()  

parser = argparse.ArgumentParser()
parser.add_argument("input", help="a cfs++ h5 file or a tensor \"[e11, ...]\" with 11/22/33/32/31/21 for 2D and 11/12/22/13/23/... for 3D")
parser.add_argument("--h5_step", help="step number, too high is last (default '9999')", default=9999, type=int)
parser.add_argument("--h5_region", help="region name (default 'mech')", default="mech")
# parser.add_argument('--h5_nonreg', action='store_true', help='assume the h5 file to be nonregular')
parser.add_argument('--h5_info', action='store_true', help='dump some meta data information about the h5 file')
parser.add_argument("--tensor", help="tensor name: 'mechTensor', 'piezoTensor, 'elecTensor'", default="mechTensor")
parser.add_argument("--scale", help="manual scaling factor", default=-1.0, type=float)
parser.add_argument("--res", help="x-resolution (default 1000)", default=1000, type=int)
parser.add_argument("--sampling", help="sampling rate (default 180", default=180, type=float)
parser.add_argument("--show", help="mode within boebbale, hom_rect or streamline", choices=['ortho_norm', 'mono_norm', 'ortho_err', 'hom_rect', 'hom_rot_cross', 'rot', 'stream'])
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
parser.add_argument("--hom_samples", help="activates interpolation and, the value gives samples in x-direction", type=int)
parser.add_argument("--stream_style", help="select visualization", choices=['line', 'thick'], default='thick')
parser.add_argument("--stream_step", help="step length for ODE integration per macro cell", type=float, default=0.2)
parser.add_argument("--stream_s2_samples", help="sampling of s2 if not given hom_samples applies", type=int)
parser.add_argument("--minimal", help="minimal stiffness to be drawn, will be scaled", type=float, default=0.0)
parser.add_argument("--parametrization", help="parametrization of the stiffness tensor", default="hom_rect", choices=['hom_rect', 'trans-iso', 'ortho'])
parser.add_argument("--save", help="save 'image.png' (pixel), 'image.pdf' (vector) or VTK Poly Data file 'file.vtp'")
parser.add_argument("--plot", help="for single tensors: creates gnuplot file instead of image")
parser.add_argument("--penalty", help="penalty parameter for SIMP (default 5)", default=5.0)
parser.add_argument("--color", help="only for hom_rot_cross: black or grayscale", default="grayscale")
parser.add_argument("--info", help="creates a xml file of given name with additional information")
parser.add_argument("--unstructured", help="number of structured elements per coordinate as list nx,ny,nz", default="")
parser.add_argument("--nodefile", help="name of the design to node file", default="")
parser.add_argument("--thres", help="threshold value for 3D VTK plot", type=float, default=0.0)

args = parser.parse_args()

# check ans postproc arguments
if not args.symmetries == "default" and not args.show == "default" and not args.symmetries == args.show:
  print "'show' and 'symmetries' do not match"
  sys.exit()
aux_code = args.show if not args.show == "default" else args.symmetries  # might still be default

# in this global varibale we can store meta-information to be exported as xml file 
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
# check if we read data from command line instead from an h5 file
if args.input.startswith('['):
  h5_read = False
  input = eval(args.input)
  if len(input) <> 21 and len(input) <> 6:
    print "the input has " + str(len(input)) + " coefficients but requires 6 (2D) or 21 (3D)"
    sys.exit()

  dim_2D = len(input) <> 21
  
  if args.tensor == 'mechTensor':  
    tensor = to_mech_tensor(input)
    tensor = HillMandel2Voigt(tensor) if args.notation == "mandel" else tensor
    print "Voigt notation of input tensor:"
  if args.tensor == 'piezoTensor':
    tensor = to_piezo_tensor(input)
  
  dump_tensor(tensor)

  if len(tensor) == 3 or len(tensor) == 2:
    assert((len(tensor) == 3 and args.tensor == 'mechTensor') or (len(tensor) <> 3 and args.tensor <> 'mechTensor'))
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
  if not os.path.exists(args.input):
    print 'Error: file does not exist: ' + args.input
    sys.exit()
  f = h5py.File(args.input, 'r')
  if args.h5_info:
    dump_h5_meta(f)
    sys.exit()   
  validate_region(f, args.h5_region)
  if len(args.unstructured) <> 0:
    nondes_centers, nondes_min, nondes_max, nondes_elem_dim, nondes_force, nondes_support = centered_elements(f, 'nondesign', False, 'load', 'support')
    print 'Reading elements from H5-file done '
    dim_2D = nondes_min[2] == nondes_max[2]
    print 'detected dimension ' + ('2D ' if dim_2D else '3D ') + "in non-design region" 
  centers, min, max, elem_dim, _, _ = centered_elements(f, args.h5_region)
  dim_2D = min[2] == max[2]
  print 'detected dimension ' + ('2D' if dim_2D else '3D')

# perform 2D and 3D from file
if h5_read or dim_2D:
  # either Image or polydata  
  viz = None
  if args.show == "hom_rect" or args.show == "hom_rot_cross" or args.show == "rot" or args.show == 'stream' or args.unstructured:
    s1, s2, s3, angle = read_stiff_angle(f, dim_2D, args)
    # add angle bias
    if args.angle_factor <> 1.0:
      print 'scale angle by ' + str(args.angle_factor)
      angle *= args.angle_factor  
    
    print 'unscaled s1 in [' + str(numpy.min(s1)) + ':' + str(numpy.max(s1)) + '] s2 in [' + str(numpy.min(s2)) + ':' + str(numpy.max(s2)) + ']'
    # coords for non-design region
    if args.unstructured:
      nondes_coords = (nondes_centers, nondes_min, nondes_max, nondes_elem_dim)
    coords = (centers, min, max, elem_dim)

    # viz is either Image or polydata
    if dim_2D:
      if args.show == "hom_rect": 
        if args.hom_grad == 'none':
          viz = show_frame(coords, s1, s2, args.hom_dir, args.res)
        else:
          viz = show_frame_grad(coords, s1, s2, args.hom_grad, args.hom_dir, args.res)
      elif args.show == "hom_rot_cross" or args.show == "rot":
        # add optional angle bias
        print 'change angle'
        if args.hom_grad == 'none':
          viz = show_rot_cross(coords, s1, s2, angle[:, 0], args.hom_dir, args.res, args.scale, args.color, args.save)
        else:
          viz = show_rot_cross_grad(coords, s1, s2, angle[:, 0], args.hom_grad, args.hom_dir, args.res, args.scale)
      elif args.show == "stream":
          viz = show_streamline(coords, s1, s2, angle[:, 0], args.hom_dir, args.scale, args.minimal, args.stream_style, args.stream_step, args.hom_samples, args.stream_s2_samples, args.res, args.save <> None, info)            
      else:
        assert(False)
    # the 3D VTK stuff or 3D mesh for unstructured design and non-design case      
    else:       
      if args.show == "rot":
        s1 = ones((len(s1), 1)) * 0.25
        s2 = ones((len(s1), 1)) * 0.25
        s3 = ones((len(s1), 1)) * 0.25
      if args.hom_samples:
        if args.hom_grad == 'none':
          print 'for hom_rect in 3D with hom_samples you need to specify hom_grad'
          exit()
        else:
          viz = create_3d_frame_ip(coords, s1, s2, s3, angle, args.hom_samples, args.hom_grad, args.hom_dir, args.scale, args.thres)
      else:  # no sample        
        if args.hom_grad == 'none':
            viz = create_3d_frame(coords, s1, s2, s3, angle, args.hom_dir, args.scale)
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
            write_node_design(args.nodefile + '_' + repr(nx) + 'x' + repr(ny) + 'x' + repr(nz), mesh, coords, s1, s2, s3, angle)
          else:
            # create structured 3D mesh from given hdf5_file which can be non-structured
            mesh = create_3d_mesh_unstructured(coords, nondes_coords, nondes_force, nondes_support, s1, s2, s3, angle, nx, ny, nz, args.hom_grad, args.scale)
            # postprocess the structured mesh for the robot arm, removal of some irregularities 
            mesh = post_process_mesh(mesh, coords, nondes_coords, nx, ny, nz)


        else:
          print 'hom_rect in 3D only for --hom_grad none implemented'
          exit()
            
  # no hom_rect stuff but orientational stiffness
  else:
    if args.tensor == 'mechTensor':
      print "Input data is read as as " + args.notation
    tensor = get_element(f, args.tensor, args.h5_region, args.h5_step)
    angle, data = perform_rotations(tensor, args.notation, int(args.sampling), args.tensor, args.show)
    
    if args.plot <> None:
      plot_angle_data(args.plot, angle, data)
      # quit!! such we can check for viz!!
      exit() 
    else:
      viz = orientational_stiffness(centers, angle, data, args.res, args.scale)

  if viz == None:
    if args.unstructured:
      print 'Structured mesh is created and saved in output_' + repr(nx) + 'x' + repr(ny) + 'x' + repr(nz) + '.mesh'
      if not(args.nodefile):
        write_gid_mesh(mesh, 'output_' + repr(nx) + 'x' + repr(ny) + 'x' + repr(nz) + '.mesh')
    else:  
      print 'Error: no visualization calculated!'
  else:
    show_or_write(viz, args)      


# not from file and not 2D -> this is the single tensor with optional planes 
else:
  angle, data, aux = perform_cfs_rotation(tensor, int(args.sampling), aux_code)

  print "largest stiffness: " + str(numpy.max(data)) + " smallest stiffness: " + str(numpy.min(data))
  if len(aux) > 0:
    print "largest " + args.show + ": " + str(numpy.max(aux)) + " smallest " + args.show + ": " + str(numpy.min(aux))
  
  poly = create_vtk_poly_data(angle, data if args.show == "default" else aux)
  
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
    actors = create_symmety_planes(mins, 1.2 * numpy.max(data if args.show == "default" else aux), not args.symmetries_planes == "false")

  show_write_vtk(poly, args.res, args.save, actors)  

if args.info:
  print 'write info xml file: ' + args.info
  out = open(args.info, "w")
  out.write(xml.dom.minidom.parseString(xml.etree.ElementTree.tostring(info)).toprettyxml())
  
