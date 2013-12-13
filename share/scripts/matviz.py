#!/usr/bin/env python
from matviz_vtk import *
from matviz_2d  import *
from hdf5_tools import *
import argparse
import sys
import os

## reads design_stiff* and design_rotAngle* for 2D and 3D. Fills other stuff by defaults
# considers density read
# @return s1, s2, s3, angle1, angle2 
def read_stiff_angle(hdf_file, dim_2D, args):
  # rot means, that we only show rotAngle, e.g. for piezoelectric polarization
  s1 = get_element(f, "design_stiff2_" + args.hom_access, args.h5_region, args.h5_step) if args.show <> "rot" else numpy.ones((len(centers),1)) * .1 
  s2 = get_element(f, "design_stiff1_" + args.hom_access, args.h5_region, args.h5_step) if args.show <> "rot" else numpy.ones((len(centers),1)) * .1
  s3 = numpy.ones((len(centers),1)) * .1 if dim_2D or args.show == "rot" else get_element(f, "design_stiff3_" + args.hom_access, args.h5_region, args.h5_step) 
  if args.density_scale:
    rho = get_element(f, "design_density_" + args.hom_access, args.h5_region, args.h5_step)
    s1 *= rho
    s2 *= rho
    s3 *= rho
  
  
  a1 = numpy.zeros((len(s1),1))
  a2 = numpy.zeros((len(s1),1))
  
  if args.show == "hom_rot_cross" or args.show == "rot":
    try:
      if dim_2D:
        a1 = get_element(f, "design_rotAngle_" + args.hom_access, args.h5_region, args.h5_step)
      else:
        a1 = get_element(f, "design_rotAngleX_" + args.hom_access, args.h5_region, args.h5_step)
        a2 = get_element(f, "design_rotAngleY_" + args.hom_access, args.h5_region, args.h5_step)
    except:
      print 'could not read angle, ignore it'

  return s1, s2, s3, a1, a2

## show or write either Image or polydata
# @param viz eithe Image or polydata
# @save filename for output
def show_or_write(viz, args):
  assert(viz <> None)
  
  if isinstance(viz, Image.Image):
    if args.save:
      viz.save(args.save)
    else:
      viz.show()
  else:
    show_write_vtk(viz, args.res, args.save)

def plot_angle_data(file, angle, data):
  plot = open(file, "w")
  plot.write('#angle\tdata\n')
  for i in range(len(angle[0])): # assume a single tensor, take the first element
    plot.write(str(angle[0][i]) + '\t' + str(data[0][i]) + '\n')
  plot.close()  

parser = argparse.ArgumentParser()
parser.add_argument("input", help="a cfs++ h5 file or a tensor \"[e11, ...]\" with 11/22/33/32/31/21 for 2D and 11/12/22/13/23/... for 3D")
parser.add_argument("--h5_step", help="step number, too high is last (default '9999')", default=9999, type=int)
parser.add_argument("--h5_region", help="region name (default 'mech')", default="mech")
parser.add_argument('--h5_info', action='store_true', help='dump some meta data information about the h5 file')
parser.add_argument("--tensor", help="tensor name: 'mechTensor', 'piezoTensor, 'elecTensor'", default="mechTensor")
parser.add_argument("--scale", help="manual scaling factor", default=-1.0, type=float)
parser.add_argument("--res", help="x-resolution (default 1000)", default=1000, type=int)
parser.add_argument("--sampling", help="sampling rate (default 180", default=180, type=float)
parser.add_argument("--show", help="default | ortho_norm | mono_norm (3D) | ortho_err | e21_normed (2D) | hom_rect | hom_rot_cross | rot | shear", default="default", choices=['ortho_norm', 'mono_norm', 'ortho_err', 'hom_rect', 'hom_rot_cross', 'rot', 'shear'])
parser.add_argument("--notation", help="mandel | voigt (default 'voigt')", default="voigt")
parser.add_argument("--symmetries", help="same options as for shows", default="default")
parser.add_argument("--symmetries_max", help="maximum number of symmetries (default 999)", default=999)
parser.add_argument("--symmetries_threshold", help="threshold value for symmetries (default 9999)", default=9999)
parser.add_argument("--symmetries_mode", help="'minima' or 'all' subject to max and threshold (default 'minima'", default="minima")
parser.add_argument("--symmetries_planes", help="'true' or 'false' for 3D also show planes to normals (default 'false')", default="false")
parser.add_argument("--hom_access", help="the 'plain ' or 'smart' hom values (default 'smart')", default = "smart")
parser.add_argument("--hom_grad", help="interpolation of design: 'none', 'nearest', linear', 'cubic' (default 'linear')", default="linear", choices=['none', 'nearest', 'linear', 'cubic'] )
parser.add_argument("--hom_dir", help="visualization of stiffness directions (default 'all')", default="all", choices=['all', 'horizontal', 'vertical', 'sagittal'] )
parser.add_argument("--hom_angle", help="bias added to the angle in grad!", default=0.0, type=float )
parser.add_argument('--density_scale', action='store_true', help='scale by <hom_access> density')
parser.add_argument("--save", help="save 'image.png' or VTK Poly Data file 'file.vtp'")
parser.add_argument("--plot", help="for single tensors: creates gnuplot file instead of image")
args = parser.parse_args()

# check ans postproc arguments
if not args.symmetries == "default" and not args.show == "default" and not args.symmetries == args.show:
  print "'show' and 'symmetries' do not match"
  sys.exit()
aux_code = args.show if not args.show == "default" else args.symmetries # might still be default

# check for tensor input
tensor  = [] # becomes a single tensor or a tensor array
centers = []
dim_2D  = None
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
    centers.append([0.5,0.5,0.0])

# read 2D or 3D data from h5 file     
else:   
  h5_read = True
  if not os.path.exists(args.input):
    print 'Error: file does not exist: ' + args.input
    sys.exit()
  f = h5py.File(args.input)
  
  if args.h5_info:
    dump_h5_meta(f)
    sys.exit()
     
  validate_region(f, args.h5_region)
  centers, min, max, elem_dim  = centered_elements(f, args.h5_region)
  tensor = get_element(f, args.tensor, args.h5_region, args.h5_step)

  dim_2D = min[2] == max[2]
  print 'detected dimension ' + ('2D' if dim_2D else '3D')  
  
  
#perform 2D and 3D from file
if h5_read or dim_2D:
  # either Image or polydata  
  viz = None
  if args.show == "hom_rect" or args.show == "hom_rot_cross" or args.show == "rot":

    s1, s2, s3, a1, a2 = read_stiff_angle(f, dim_2D, args)
    coords = (centers, min, max, elem_dim)

    # viz is either Image or polydata
    if args.show == "hom_rot_cross" or args.show == "rot":
      # add optional angle bias
      print 'change angle'
      if args.hom_grad == 'none':
        viz = show_rot_cross(coords, s1, s2, a1, args.hom_dir, args.res, args.scale)
      else:
        viz = show_rot_cross_grad(coords, s1, s2, a1, args.hom_grad, args.hom_dir, args.res, args.scale)          
    else:
      if args.hom_grad == 'none':
        if dim_2D:
          viz = show_frame(coords, s1, s2, args.hom_dir, args.res)
        else:
          viz = create_3d_frame(coords, s1, s2, s3, args.hom_dir, args.scale)
      else:
        if dim_2D:
          viz = show_frame_grad(coords, s1, s2, args.hom_grad, args.hom_dir, args.res)
        else:
          print 'graded hom_rect in 3D not yet implemented'
          exit()
  # no hom_rect stuff but orientational stiffness
  else:
    if args.tensor == 'mechTensor':
      print "Input data is read as as " + args.notation
    angle, data = perform_rotations(tensor, args.notation, int(args.sampling), args.tensor, args.show)
    
    if args.plot <> None:
      plot_angle_data(args.plot, angle, data)
      # quit!! such we can check for viz!!
      exit() 
    else:
      viz = orientational_stiffness(centers, angle, data, args.res, args.scale)

  if viz == None:
    print 'Error: no visualization calculated!'
  else:
    show_or_write(viz, args)      


# not from file and not 2D -> this is the single tensor with optional planes 
else:
  angle, data, aux = perform_cfs_rotation(tensor, int(args.sampling), aux_code)

  print "largest stiffness: " + str(numpy.max(data)) + " smallest stiffness: " + str(numpy.min(data))
  if len(aux) > 0:
    print "largest " + args.show +": " + str(numpy.max(aux)) + " smallest " + args.show + ": " + str(numpy.min(aux))
  
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
