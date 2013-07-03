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
parser.add_argument("--threshold", help="threshold for void material with 0 and 1 (default 0.5)", default=0.5)
parser.add_argument('--scale', help="scales by width w.r.t. one meter (default 1.0)", type = float, default=1.0)
parser.add_argument('--rhomin', help="maps pure white in tne image (default 0.001)", default=0.0001)
parser.add_argument('--showbinary', action='store_true', help='shows only a binary pop-up image')
parser.add_argument('--showsize', help="pixels in x direction for pop-up (default 800)", default=800)
parser.add_argument('--noshow', action='store_true', help='do not pop up the image window')
parser.add_argument('--densfile', help=".density.xml with only 'mech' densities")


args = parser.parse_args()
if not os.path.exists(args.input):
  print 'input file not found: ' + args.input
  sys.exit() 

# do it to generate statistical output of what would happen
if '.xml' in args.input:
  d  = read_density(args.input, "physical")
  create_dense_mesh_density(d, mesh, args.threshold, args.scale, args.rhomin)
elif '.txt' in args.input:
  d = load_matrix_from_file(args.input)
  create_dense_mesh_density(d, mesh, args.threshold, args.scale, args.rhomin)
else:
    # read the png into a list
  input_img = Image.open(args.input)
  print "original image mode: " + input_img.mode
  if input_img.mode == 'I':
    print "Warning: mode is stupid, may give unusable results!"
  
  input_img = input_img.convert("L") #.transpose(Image.FLIP_TOP_BOTTOM)
  create_dense_mesh_img(input_img, mesh, float(args.threshold), float(args.scale), float(args.rhomin))

if not args.noshow:
  if '.xml' in args.input or '.txt' in args.input:
    show_dense_mesh_image(mesh, d.shape, args.showbinary, int(args.showsize))
  else:
    show_dense_mesh_image(mesh, input_img.size, args.showbinary, int(args.showsize))

if args.densemesh:
  write_gid_mesh(mesh, args.densemesh)  
  print "save dense mesh: " + args.densemesh

if args.sparsemesh:
  sparse = convert_to_sparse_mesh(mesh)
  print "save sparse mesh: " + args.sparsemesh
  write_gid_mesh(sparse, args.sparsemesh)
  mesh = sparse

if args.densfile <> None:
  densities = []
  enr = []
  for i in range(len(mesh.elements)):
    if mesh.elements[i].region == 'mech': 
      densities.append(mesh.elements[i].density)
      enr.append(i+1)
  data = numpy.zeros((len(densities),1))
  data[:,0] = densities    
  print data.shape
  write_multi_design_file(args.densfile, data, ['design'], enr)
