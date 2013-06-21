#!/usr/bin/env python
from mesh_tool import *


parser = argparse.ArgumentParser()
parser.add_argument("--res", help="long side resolution of mesh if action is 'mesh", type=int, required = True )
parser.add_argument('--type', help="mesh type: cantilever2d, ", choices=['cantilever2d', 'cantilever2d_reinforced'], required = True)
parser.add_argument('--file', help="optional give output file name")

args = parser.parse_args()

mesh = create_cantilever2d_mesh(args.type, args.res)

file = args.type + '_' + str(args.res) + '.mesh' if args.file == None else args.file 

write_gid_mesh(mesh, file)
print "created file '" + file + "' with " + str(len(mesh.elements)) + " elements"