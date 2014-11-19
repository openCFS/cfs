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

  if type(f) == types.StringType:
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
parser.add_argument("input", help="a greyscale image (any format, use gif when png makes problems)")
parser.add_argument("--densemesh", help="writes a dense two region mesh (void/mech) with given name")
parser.add_argument("--sparsemesh", help="writes a sparse mesh with region mech with given name")
parser.add_argument("--threshold", help="threshold for void material with 0 and 1 (default 0.5)", default=0.5, type=float)
parser.add_argument('--scale', help="scales by width w.r.t. one meter (default 1.0)", type=float, default=1.0)
parser.add_argument('--rhomin', help="maps pure white in tne image (default 0.001)", default=1e-7)
parser.add_argument('--showbinary', action='store_true', help='shows only a binary pop-up image')
parser.add_argument('--showsize', help="pixels in x direction for pop-up (default 800)", default=800, type=int)
parser.add_argument('--noshow', action='store_true', help='do not pop up the image window')
parser.add_argument('--densfile', help="output .density.xml with only 'mech' densities")
parser.add_argument('--multi_d', help="number of design variables in density.xml file", type=int, default=1)
parser.add_argument('--shearangle', help="shearing angle of mesh in degree", type=float, default=0.0)
parser.add_argument('--colorregion', help="interpret colors as other regions", action='store_true')
parser.add_argument('--pressure', help='sets region for pressure in the meshfile')

mesh = Mesh()

args = parser.parse_args()
if not os.path.exists(args.input):
  print 'input file not found: ' + args.input
  sys.exit() 
print args.multi_d
# do it to generate statistical output of what would happen
if args.multi_d == 3:
  multi_d = read_multi_design(args.input, 'stiff1', 'stiff2', 'rotAngle', matrix=True)
  create_dense_mesh_density(multi_d, mesh, args.threshold, args.scale, args.rhomin, 3)
elif args.multi_d == 2:
  multi_d = read_multi_design(args.input, 'stiff1', 'stiff2', matrix=True)
  create_dense_mesh_density(multi_d, mesh, args.threshold, args.scale, args.rhomin, 2)
elif '.xml' in args.input:
  d = read_density(args.input, "design")
  create_dense_mesh_density(d, mesh, args.threshold, args.scale, args.rhomin)
elif '.txt' in args.input:
  d = load_matrix_from_file(args.input)
  create_dense_mesh_density(d, mesh, args.threshold, args.scale, args.rhomin)
else:
    # read the png into a list
  img = Image.open(args.input).transpose(Image.FLIP_TOP_BOTTOM)
  print "original image mode: " + img.mode
  if img.mode == 'I':
    print "Warning: mode is stupid, may give unusable results!"
  if not args.colorregion:
    img = img.convert("L") 
  if args.pressure:
    create_dense_mesh_img(img, mesh, float(args.threshold), float(args.scale), float(args.rhomin), float(args.shearangle))
  else:
    create_dense_mesh_img(img, mesh, float(args.threshold), float(args.scale), float(args.rhomin), float(args.shearangle), True)

if not args.noshow:
  dimension = None
  if args.multi_d > 1:
    dimensions = (multi_d.shape[0], multi_d.shape[1])
  elif '.xml' in args.input or '.txt' in args.input:
    dimensions = d.shape
  else:
    dimensions = img.size

  show_dense_mesh_image(mesh, dimensions, args.showbinary, args.showsize)
  
if args.densemesh:
  write_gid_mesh(mesh, args.densemesh)  
  print "save dense mesh: " + args.densemesh

if args.sparsemesh:
  sparse = convert_to_sparse_mesh(mesh)
  print "save sparse mesh: " + args.sparsemesh
  write_gid_mesh(sparse, args.sparsemesh)
  mesh = sparse

if args.densfile <> None:
  if args.multi_d == 1:
    densities = []
    enr = []
    for i in range(len(mesh.elements)):
      if mesh.elements[i].region == 'mech': 
        densities.append(mesh.elements[i].density)
        enr.append(i + 1)
    data = numpy.zeros((len(densities), 1))
    data[:, 0] = densities    
    write_multi_design_file(args.densfile, data, ['density'], enr)
  else:
    stiff1 = []
    stiff2 = []
    if args.multi_d == 3:
      rotAngle = []
    enr = []
    for i in range(len(mesh.elements)):
      if mesh.elements[i].region == 'mech': 
        stiff1.append(mesh.elements[i].stiff1)
        stiff2.append(mesh.elements[i].stiff2)
        if args.multi_d == 3:
          rotAngle.append(mesh.elements[i].rotAngle)
        enr.append(i + 1)
    if args.multi_d == 2:
      data = numpy.zeros((len(stiff1), 2))
      data[:, 0] = stiff1
      data[:, 1] = stiff2
      write_multi_design_file(args.densfile, data, ['stiff1', 'stiff2'], enr) 
    else:
      data = numpy.zeros((len(stiff1), 3))
      data[:, 0] = stiff1
      data[:, 1] = stiff2
      data[:, 2] = rotAngle  
      write_multi_design_file(args.densfile, data, ['stiff1', 'stiff2', 'rotAngle'], enr) 

  
