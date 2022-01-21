#!/usr/bin/env python
from mesh_tool import *
from optimization_tools import *
import argparse
import types
import numpy

def load_matrix_from_file(f):
  """
  This function is to load an ascii format matrix (float numbers separated by
  whitespace characters and newlines) into a numpy matrix object.
  f is a file object or a file path.
  """

  if type(f) == bytes:
    fo = open(f, 'r')
    matrix = load_matrix_from_file(fo)
    fo.close()
    return matrix
  elif type(f) == types.FileType:
    file_content = f.read().strip()
    file_content = file_content.replace('\r\n', ';')
    file_content = file_content.replace('\n', ';')
    file_content = file_content.replace('\r', ';')
    return numpy.matrix(file_content)
    raise TypeError('f must be a file object or a file name.') 


parser = argparse.ArgumentParser()
parser.add_argument("input", help="a (grayscale) image (any format, use gif when png makes problems) or .xml file or .txt file or h5 file")
parser.add_argument("--densemesh", help="writes a dense two region mesh (void/mech) with given name")
parser.add_argument("--sparsemesh", help="writes a sparse mesh with region mech with given name")
parser.add_argument('--density', help="write dense .density.xml (for all regions)")
parser.add_argument('--sparsedensity', help="write .density.xml only for 'mech' region elements")
parser.add_argument("--threshold", help="threshold for void material with 0 and 1 (default 0.5)", default=0.5, type=float)
parser.add_argument('--scale', help="scales by width w.r.t. one meter (default 1.0)", type=float, default=1.0)
parser.add_argument('--rhomin', help="maps pure white in tne image (default 0.001)", default=1e-3, type=float)
parser.add_argument('--showbinary', action='store_true', help='shows only a binary pop-up image')
parser.add_argument('--showsize', help="pixels in x direction for pop-up (default 800)", default=800, type=int)
parser.add_argument('--noshow', action='store_true', help='do not pop up the image window')
parser.add_argument('--colorregion', help="interpret colors as other regions", action='store_true')
parser.add_argument('--attribute', help='attribute to be read from .xml file (design / physical)', default='design')

mesh = None

args = parser.parse_args()
if not os.path.exists(args.input):
  print('input file not found: ' + args.input)
  sys.exit() 
# do it to generate statistical output of what would happen
if '.xml' in args.input:
  d = read_density(args.input, args.attribute)
  mesh = create_2d_mesh_from_array(d, args.threshold)
  mesh.scale(args.scale)
elif '.txt' in args.input:
  d = load_matrix_from_file(args.input)
  mesh = create_2d_mesh_from_array(d, args.threshold)
  mesh.scale(args.scale)
elif '.h5' in args.input:
  f = h5py.File(args.input, 'r')
  mesh = create_mesh_from_hdf5(f, ['mech'],['bottom','top','left','right'], threshold = args.threshold)
else:
  # read the png into a list
  img = Image.open(args.input)
  print("original image mode: " + img.mode)
  if img.mode == 'I':
    print("Warning: mode is stupid, may give unusable results!")
  if not args.colorregion:
    img = img.convert("L")
  else:
    img = img.convert("RGB")
  
  mesh = create_2d_mesh_from_image(img, args.threshold)
  mesh.scale(args.scale)

if not '.h5' in args.input and not args.noshow:
  show_dense_mesh_image(mesh, args.showbinary, args.showsize)
  
if args.densemesh:
  write_ansys_mesh(mesh, args.densemesh)
  print("save dense mesh: " + args.densemesh)

if args.sparsemesh:
  sparse = convert_to_sparse_mesh(mesh)
  print("save sparse mesh: " + args.sparsemesh)
  write_ansys_mesh(sparse, args.sparsemesh)
  mesh = sparse

if args.density != None:
  assert(mesh.nx * mesh.ny == len(mesh.elements))
  data = numpy.zeros((mesh.ny, mesh.nx))
  
  for x in range(mesh.nx):
    for y in range(mesh.ny):
      # the image was transposed    
      data[y,x] = mesh.elements[x * mesh.ny + y].density
  print("save dense density file '" + args.density + "'")
  write_density_file(args.density, data)
   

if args.sparsedensity != None:
  densities = []
  enr = []
  for i in range(len(mesh.elements)):
    if mesh.elements[i].region == 'mech': 
      densities.append(mesh.elements[i].density)
      enr.append(i + 1)
  data = numpy.zeros((len(densities), 1))
  data[:, 0] = densities    
  write_multi_design_file(args.sparsedensity, data, ['density'], enr)
