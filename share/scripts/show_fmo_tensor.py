#!/usr/bin/env python

from hdf5_tools import *
import argparse
import vtk
from vtk.util.colors import *
import sys

## creates 3D data to vtkPolyData
def create_vtk_poly_data(angle, data):
  
  assert(len(angle) == len(data))
  samples = int(sqrt(len(angle)))
  assert(samples**2 == len(angle))
  
  pts = vtk.vtkPoints()
  
  for i in range(len(angle)):
    radius = data[i]
    
    # we use not the usual conversion for the spherical coordinate system but the first row of the rotation matrix!
    p = radius * to_vector(angle[i])
    
    # print "phi="  + str(phi) + " theta=" + str(theta) + " radius=" + str(radius) + " -> x=" + str(x) + " y=" + str(y) + " z=" + str(z)
    pts.InsertNextPoint(p[0], p[1], p[2])

  pdo = vtk.vtkPolyData()
  pdo.SetPoints(pts)
  #Allocate memory for triangles
  pdo.Allocate((samples-1) * (samples-1) * 2)
  
  #Generate two triangles per quad
  for i in range(0, samples-1):
    for j in range(0, samples-1):
      tri = vtk.vtkTriangle()
      tri.GetPointIds().SetNumberOfIds(3)
      tri.GetPointIds().SetId(0, i * samples + j)
      tri.GetPointIds().SetId(1, i * samples + j+1)
      tri.GetPointIds().SetId(2, (i+1) * samples + j)
      # print "i=" + str(i) + " j=" + str(j) + " ids=" + str(tri.GetPointIds())
      pdo.InsertNextCell(tri.GetCellType(), tri.GetPointIds())
      
      tri = vtk.vtkTriangle()
      tri.GetPointIds().SetNumberOfIds(3)
      tri.GetPointIds().SetId(0, (i+1) * samples + j+1)
      tri.GetPointIds().SetId(1, (i+1) * samples + j)
      tri.GetPointIds().SetId(2, i * samples + j+1)
      pdo.InsertNextCell(tri.GetCellType(), tri.GetPointIds())
  
  return pdo

## create a list of vtk actors displaying symmetry planesfrom vtk.util.colors import
# @param minim
def create_symmety_planes(minima, scale, add_planes):
  # code source: http://code.google.com/p/pythonxy/source/browse/src/python/vtk/DOC/Examples/Rendering/Cylinder.py
 
  minima = []
  minima.append(((0, 0), 1))
  minima.append(((0, numpy.pi/2), 1))
  minima.append(((numpy.pi/2, numpy.pi/2), 1))
 
  actors = [] 
   
  for i in range(len(minima)):
     
    # Create cylinder
    cylinder = vtk.vtkCylinderSource()
    cylinder.SetRadius(scale / 200)
    cylinder.SetCenter( 0, 0, 0 )
    cylinder.SetHeight(2 * scale)
    cylinder.SetResolution(10)
    cylinder.Update()
     
    # The mapper is responsible for pushing the geometry into the graphics
    # library. It may also do color mapping, if scalars or other attributes are defined.
    mapper_c = vtk.vtkPolyDataMapper()
    mapper_c.SetInputConnection(cylinder.GetOutputPort())

    # The actor is a grouping mechanism: besides the geometry (mapper), it
    # also has a property, transformation matrix, and/or texture map.
    actor_c = vtk.vtkActor()
    actor_c.SetMapper(mapper_c)
        
    actor_c.GetProperty().SetColor(black)
    angle = minima[i][0]
    phi   = angle[0]
    theta = angle[1]

    
    print "angle: " + str(angle) + " -> " + str(to_vector(angle)) + " = " + str(minima[i][1])
    
    #actor_c.RotateX(90)
    #actor_c.RotateY(angle[0] * 180/numpy.pi)
    #actor_c.RotateZ(angle[1] * 180/numpy.pi)


    # cylinders!!!
    actor_c.RotateZ(90)
    actor_c.RotateY(angle[0] * 180/numpy.pi)
    actor_c.RotateX(angle[1] * 180/numpy.pi)


    actors.append(actor_c)

    if add_planes:
      disk = vtk.vtkDiskSource()
      disk.SetInnerRadius(0.1)
      disk.SetOuterRadius(0.9 * scale)
      disk.SetRadialResolution(30)
      disk.SetCircumferentialResolution(60)
      disk.Update()
  
      mapper_d = vtk.vtkPolyDataMapper()
      mapper_d.SetInputConnection(disk.GetOutputPort())
      
      actor_d = vtk.vtkActor()
      actor_d.SetMapper(mapper_d)

      actor_d.GetProperty().SetColor(banana)
      actor_d.GetProperty().SetOpacity(0.5)

      #actor_d.RotateX(angle[0] * 180/numpy.pi)
      #actor_d.RotateY(angle[1] * 180/numpy.pi)
      #actor_d.RotateZ(0)

      # plane!
      actor_d.RotateY(90)
      actor_d.RotateZ(90)
      actor_d.RotateY(angle[0] * 180/numpy.pi)
      actor_d.RotateX(angle[1] * 180/numpy.pi)
      

      
      actors.append(actor_d)
    
  return actors


# show the data on the screen
# @planes list of vtk actors containing symmetry planes 
def show_vtk(polydata,  planes,res):
  #Create a mapper and actor
  mapper = vtk.vtkPolyDataMapper()
  mapper.SetInput(polydata)
  
  actor = vtk.vtkActor()
  actor.SetMapper(mapper)
   
  # Setup a renderer, render window, and interactor
  renderer = vtk.vtkRenderer()
  renderWindow = vtk.vtkRenderWindow()
  #renderWindow.SetWindowName("Test")
   
  renderWindow.AddRenderer(renderer);
  renderWindowInteractor = vtk.vtkRenderWindowInteractor()
  renderWindowInteractor.SetRenderWindow(renderWindow)
   
  #Add the actor to the scene
  for i in range(len(planes)):
    renderer.AddActor(planes[i])

  renderer.AddActor(actor)
  renderer.SetBackground(1,1,1) # Background color white
   
  #Render and interact
  renderWindow.Render()
   
  #*** SetWindowName after renderWindow.Render() is called***
  #renderWindow.SetWindowName("Test")
   
  renderWindowInteractor.Start()



parser = argparse.ArgumentParser()
parser.add_argument("input", help="a cfs++ h5 file or a tensor \"[e11, ...]\" with 11/22/33/32/31/21 for 2D and 11/12/22/13/23/... for 3D")
parser.add_argument("--h5_step", help="step number, too high is last (default '9999')", default=9999)
parser.add_argument("--h5_region", help="region name (default 'mech')", default="mech")
parser.add_argument("--tensor", help="tensor name (default 'mechTensor')", default="mechTensor")
parser.add_argument("--scale", help="manual scaling factor", default=-1)
parser.add_argument("--res", help="x-resolution (default 1500)", default=1500)
parser.add_argument("--sampling", help="sampling rate (default 180", default=180)
parser.add_argument("--show", help="default | ortho_norm | mono_norm (3D) | ortho_err | e21_normed (2D)  hom_rect", default="default")
parser.add_argument("--notation", help="mandel | voigt (default 'mandel')", default="mandel")
parser.add_argument("--symmetries", help="same options as for shows", default="default")
parser.add_argument("--symmetries_max", help="maximum number of symmetries (default 999)", default=999)
parser.add_argument("--symmetries_threshold", help="threshold value for symmetries (default 9999)", default=9999)
parser.add_argument("--symmetries_mode", help="'minima' or 'all' subject to max and threshold (default 'minima'", default="minima")
parser.add_argument("--symmetries_planes", help="'true' or 'false' for 3D also show planes to normals (default 'false')", default="false")
parser.add_argument("--hom_access", help="the 'plain ' or 'smart' hom values (default 'smart')", default = "smart")
parser.add_argument("--hom_grad", help="interpolation of design: 'none', 'linear' (default 'linear')", default="linear")
parser.add_argument("--save", help="save 'image.png' or VTK Poly Data file 'file.vtp'")
args = parser.parse_args()

# check ans postproc arguments
if not args.symmetries == "default" and not args.show == "default" and not args.symmetries == args.show:
  print "'show' and 'symmetries' do not match"
  sys.exit()
aux_code = args.show if not args.show == "default" else args.symmetries # might still be default

# check for tensor input
tensor  = [] # becomes a single tensor or a tensor array
centers = []
dim_2D  = True
if args.input.startswith('['):
  
  input = eval(args.input)
  if len(input) <> 21 and len(input) <> 6:
    print "the input has " + str(len(input)) + " coefficients but requires 6 (2D) or 21 (3D)"
    sys.exit()
  if len(input) == 21:
    dim_2D = False  
  tensor = to_mech_tensor(input)
  
  tensor = HillMandel2Voigt(tensor) if args.notation == "mandel" else tensor
  
  print "Voigt notation of input tensor:"
  dump_tensor(tensor)
  
  if len(tensor) == 3:
    # convert back to hill mandel as we simulate a h5 file
    vec = to_mech_vector(Voigt2HillMandel(tensor), as_array=True)
    tensor = []
    tensor.append(vec)
    centers.append([0.5,0.5,0.0])
    
else:   
  # read 2D CFS optimization result 
  f = h5py.File(args.input)
  centers, min, max, elem_dim  = centered_elements(f)
  tensor = get_element(f, args.tensor, args.h5_region, int(args.h5_step))
  
#perform 2D and 3D
if dim_2D:  
  im = 0
  if args.show == "hom_rect":
    s1    = get_element(f, "design_stiff1_" + args.hom_access, args.h5_region)
    s2    = get_element(f, "design_stiff2_" + args.hom_access, args.h5_region)
    angle = get_element(f, "design_rotAngle_plain", args.h5_region)
    im = show_rot_frame((centers, min, max, elem_dim), s1, s2, angle, args.hom_grad, int(args.res), float(args.scale))	
  else:
    angle, data = perform_rotations(tensor, int(args.sampling), args.tensor, args.show)
    im = orientational_stiffness(centers, angle, data, int(args.res), float(args.scale))
  if args.save:
    im.save(args.save)
  else:
    im.show()
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
  if args.save:
    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetInput(poly)
    writer.SetFileName(args.save)
    writer.Write()    
  else:
    show_vtk(poly, actors, int(args.res))  

