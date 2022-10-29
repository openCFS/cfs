#!/usr/bin/env python
import vtk
import numpy as np
import os
import sys
import mesh_tool
import argparse

parser = argparse.ArgumentParser(description='Convert .vtp file to volume mesh for openCFS')
parser.add_argument("name", help="specify base name of in-/ and output files without format ending .vtp")
parser.add_argument("--type", help="set b.c. for which problem type (default box_varel)",default="box_varel",choices=["box_varel","ppbox","wwb"])

args = parser.parse_args()

if not args.type in ["box_varel","ppbox"]:
  print("Works currently only for box_varel!")
  sys.exit(-1)

name = args.name

if not os.path.exists(args.name+".stl"):
  if args.name.endswith(".vtp"):
    name = args.name[:-4]

  reader = vtk.vtkXMLPolyDataReader()
  reader.SetFileName(name+".vtp")
  reader.Update()

  # Write the stl file to disk
  stlWriter = vtk.vtkSTLWriter()
  stlWriter.SetFileName(name+".stl")
  stlWriter.SetInputData(reader.GetOutput())
  stlWriter.Write()

  print("Converted .vtp to .stl file")

if not os.path.exists(args.name+".1.ele"):
  os.system("tetgen -pk " + name + ".stl")
  print("Created volume mesh with tetgen")

mesh = mesh_tool.create_mesh_from_tetgen(name,"mech")
bounds = mesh_tool.calc_min_max_coords(mesh)

if args.type == "box_varel":
  mesh = mesh_tool.add_bc_for_box_varel(mesh,bounds)
elif args.type == "ppbox":
  mesh = mesh_tool.add_bc_for_ppbox(mesh,bounds)

#mesh = mesh_tool.name_bb_faces(mesh,xmin,xmax,ymin,ymax,zmin,zmax)
mesh_tool.write_ansys_mesh(mesh,name+".mesh")

print("Converted tetgen mesh to cfs mesh and set b.c. for " + args.type)
