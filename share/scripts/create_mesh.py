#!/usr/bin/env python
from mesh_tool import *

## creates a mesh of predefined geometry
def create_cantilever2d_mesh(type, resolution):
  width = 3.0
  height = 2.0
  
  nx = resolution
  ny = int(nx * (2./3.))
  
  dx = width / nx
  dy = height / ny 

  for y in range(ny + 1):
    for x in range(nx + 1):
      mesh.nodes.append((x * dx, y * dy))
 
  # print mesh.nodes 
  for y in range(ny):
    for x in range(nx):
      e = Element()
      e.density = 1.0
      if type == 'cantilever2d_reinforced' and float(x) >= (28./30. * nx):
        e.region = 'reinforce'
      else:
        e.region = 'mech'

      # assign nodes
      ll = (nx+1) * y + x  # lowerleft
      e.nodes = ((ll, ll+1, ll+1+nx+1, ll+nx+1))
            
      mesh.elements.append(e)
  
  mesh.bc.append(("south", range(0, nx+1)))
  mesh.bc.append(("north", range((nx+1)*ny, (nx+1)*(ny+1))))
  mesh.bc.append(("west", range(0, (nx+1)*ny+1, nx+1)))
  mesh.bc.append(("east", range(nx, (nx+1)*(ny+1), nx+1)))

  mesh.bc.append(("left_lower", [0]))
  mesh.bc.append(("right_lower", [nx]))
  mesh.bc.append(("left_upper", [(nx+1)*ny]))
  mesh.bc.append(("right_upper", [(nx+1)*(ny+1)-1]))


parser = argparse.ArgumentParser()
parser.add_argument("--res", help="long side resolution of mesh if action is 'mesh", type=int, required = True )
parser.add_argument('--type', help="mesh type: cantilever2d, ", choices=['cantilever2d', 'cantilever2d_reinforced'], required = True)
parser.add_argument('--file', help="optional give output file name")

args = parser.parse_args()

create_cantilever2d_mesh(args.type, args.res)

file = args.type + '_' + str(args.res) + '.mesh' if args.file == None else args.file 

write_gid_mesh(mesh, file)
print "created file '" + file + "' with " + str(len(mesh.elements)) + " elements"