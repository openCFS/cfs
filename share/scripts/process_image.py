#!/usr/bin/env python
from mesh_tool import *
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("input", help="a greyscale image (any format, use gif when png makes problems)")
parser.add_argument("--densemesh", help="writes a dense two region mesh (void/mech) with given name")
parser.add_argument("--sparsemesh", help="writes a sparse mesh with region mech with given name")
parser.add_argument("--threshold", help="threshold for void material with 0 and 1 (default 0.5)", default=0.5)
parser.add_argument('--scale', help="scales by width w.r.t. one meter (default 1.0)", default=1.0)
parser.add_argument('--rhomin', help="maps pure white in tne image (default 0.001)", default=0.001)
parser.add_argument('--showbinary', action='store_true', help='shows only a binary pop-up image')
parser.add_argument('--showsize', help="pixels in x direction for pop-up (default 800)", default=800)
parser.add_argument('--noshow', action='store_true', help='do not pop up the image window')

args = parser.parse_args()

if not os.path.exists(args.input):
  print 'input file not found: ' + args.input
  sys.exit() 

# read the png into a list
input_img = Image.open(args.input)
print "original image mode: " + input_img.mode
if input_img.mode == 'I':
  print "Warning: mode is stupid, may give unusable results!"

input_img = input_img.convert("L") #.transpose(Image.FLIP_TOP_BOTTOM)

# do it to generate statistical output of what would happen
create_dense_mesh(input_img, mesh, float(args.threshold), float(args.scale), float(args.rhomin))

if not args.noshow:
  show_dense_mesh_image(mesh, input_img.size, args.showbinary, int(args.showsize))

if args.densemesh:
  write_gid_mesh(mesh, args.densemesh)  
  print "save dense mesh: " + args.densemesh

if args.sparsemesh:
  sparse = convert_to_sparse_mesh(mesh)
  print "save sparse mesh: " + args.sparsemesh
  write_gid_mesh(sparse, args.sparsemesh)

