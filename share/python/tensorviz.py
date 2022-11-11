#!/usr/bin/env python

## This is the frontent for single tensor visualization and rotation features. 
# It uses matviz_rot, ... as backed. Before, the features where part of matviz.py but were
# extracted to clean matviz.py

import matviz_2d
import matviz_rot
import matviz_io 
import matviz_vtk

import numpy as np
import argparse
import sys
import os.path

msg  = "This tools vizalizes a single 2d or 3d tensor. Either give the coefficients or .info.xml with a tensor\n" 
msg += "tensor coefficients are filled with zero. Give the style as hint\n"
msg += "t2d = 11 22 33 23 13 12\n"
msg += "t3d = 11 22 33 44 55 66 23 13 12 34 24 14 45 35 25 15 56 46 36 26 16\n"
msg += "c3d = 11 12 22 13 23 33 14 24 34 44 15 25 35 45 55 16 26 36 46 56 66"

parser = argparse.ArgumentParser(description = msg, formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("input", nargs='+', help="tensor coefficients 11 22 33 ... following --style or .info.xml")
# tensorviz.py has trace for default but matviz_rot.to_mech_tensor has column as default for 3d  
parser.add_argument("--style", help="default=trace style 11 22 33 ... see help. Filled with zeros", choices=['t2d', 't3d', 'c3d'], required=False)
parser.add_argument("--notation", help="mandel | voigt (default 'voigt')", choices=['mandel','voigt'], default="voigt")
parser.add_argument("--rotate", nargs='+', help="rotate tensor before sampling for max 3 angles, see axes", type=float)
parser.add_argument("--unit", help="unit of angle (default 'deg')", default='deg', choices=['deg','rad'])
parser.add_argument("--axes", help="rotation axis in 3d read from right: xyz is first around z", default='xyz',choices=['zxz','zyz','yzy','yxy','xyx','xzx','xyz','yxz','xzy','zxy','zyx','yzx'])
parser.add_argument("--res", help="x-resolution (default 1000)", default=1000, type=int)
parser.add_argument("--sampling", help="sampling rate for tensor rotation (default 180)", default=90, type=float)
parser.add_argument("--system", help="show polar plot in 2d and coord system in 3d", action='store_true', default=False)
parser.add_argument("--save", help="save 'image.png' (pixel), 'image.pdf' (vector) or VTK Poly Data file 'file.vtp'")
parser.add_argument("--cam", help="set camera (7 space-separated floats): position, focal point, roll", nargs=7, type=float, default=None)
parser.add_argument("--plot", help="write tabular data to given file name instead of showing image")

args = parser.parse_args()

dim_2D = None
input = None

# .mat is removed, add it again on need
if len(args.input) == 1 and args.input[0].endswith(".info.xml"):
  input, dim_2D = matviz_io.read_hom_tensor_from_info_xml(args.input[0])
else:
  tmp = [float(v) for v in args.input] # convert coefficient list to float list
  if len(tmp) <= 6 and args.style in [None, 't2d']:
    # fill vector of 6 elements
    input = [0] * 6 
  else:
    input = [0] * 21
  input[0:len(tmp)] = tmp

  dim_2D = len(input) != 21

style = None
if not dim_2D and args.style != 'c3d':
  style = 'trace' 
# tensorviz.py has trace for default but matviz_rot.to_mech_tensor has column as default for 3d  
tensor = matviz_rot.to_mech_tensor(input, style = style)
tensor = HillMandel2Voigt(tensor) if args.notation == "mandel" else tensor
print("Input data is read as " + args.notation + ".")
print("Voigt notation of input tensor:")
matviz_2d.dump_tensor(tensor)

if args.rotate:
  Q = None
  if args.unit == 'deg':
    args.rotate = np.radians(args.rotate)

  if dim_2D:
    assert len(args.rotate) == 1
    print('rotate in rad around z axis by',args.rotate[0])
    Q = matviz_rot.get_rot_3x3(args.rotate[0])
  else:
    assert len(args.rotate) >= 1 and len(args.rotate) <= 3
    angles = [0.0] * 3 # fill up angles to three values
    angles[0:len(args.rotate)] = args.rotate
    assert len(args.axes) == 3
    print('rotate in rad around',args.axes[2],'by',angles[2],'then around',args.axes[1],'by',angles[1],'finally around',args.axes[0],'by',angles[0])
    Q = matviz_rot.get_rot_6x6(angles[0], angles[1], angles[2], args.axes)
    
  tensor = Q @ (tensor @ Q.transpose())
  print("Voigt notation of rotated tensor:")
  matviz_rot.dump_tensor(tensor)

angle, data, _ = matviz_rot.perform_voigt_tensor_samping(tensor, int(args.sampling))

angle_max = angle[np.argmax(np.abs(data))]
angle_min = angle[np.argmin(np.abs(data))]

print(" largest entry: {:>13.6e}".format(np.max(np.abs(data)).round(15)) + "  in direction " + str(matviz_2d.to_vector(angle_max).round(15)))
print("smallest entry: {:>13.6e}".format(np.min(np.abs(data)).round(15)) + "  in direction " + str(matviz_2d.to_vector(angle_min).round(15)))

if args.plot != None:
  matviz_io.write_angle_data(args.plot, angle, data)
  sys.exit() # all done
else:
  if dim_2D:
    scale = -1.0
    coords = ([[0.5, 0.5, 0.0]], [0.0, 0.0], [1.0, 1.0], None)
    viz = matviz_2d.show_orientational_stiffness(coords, angle, data, args.res, scale, args.system)
    matviz_io.show_or_write(viz, args)
  else:
    poly = matviz_vtk.create_vtk_poly_data(angle, data)
    actors = []    
    matviz_vtk.show_write_vtk(poly, args.res, args.save, actors, args.system, camera_settings=args.cam)



