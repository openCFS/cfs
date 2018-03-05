#!/usr/bin/env python
from basecell import getConnectivity, taubin_smoothing
import matviz_vtk 
import vtk
import pymesh
import argparse

if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument('input', help="name of input stl file")
  parser.add_argument('--vtp', help="write output as vtp",required=False,action='store_true')
  parser.add_argument('--lamb', help="factor between 0 and 1 changing speed of smoothing (default: 0.4)", default=0.4, required=False)
  parser.add_argument('--maxiter', help="max number of iterations(default: 20)", default=20, required=False)
  
  args = parser.parse_args()
  
  if not (0 <= float(args.lamb) <=1):
    print("lamba must be between 0 and 1!")
    sys.exit()
  
  filebase = args.input
  if args.input.endswith(".stl"):
    filebase = args.input[:-4]
    
  mesh = pymesh.load_mesh(filebase+".stl")
  connectivity = getConnectivity(mesh.vertices,mesh.faces)
  verts = taubin_smoothing(mesh.vertices,connectivity,bounds=None,niter=int(args.maxiter),lamb=float(args.lamb))

  pd = matviz_vtk.fill_vtk_polydata(verts,mesh.faces)
  matviz_vtk.write_stl(pd,filebase+"_smoothed.stl")
  
  if args.vtp:
    matviz_vtk.show_write_vtk(pd,0,filebase+"_smoothed.vtp")
    


