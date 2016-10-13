#!/usr/bin/env python
import numpy as np
#from PIL import Image
# from mesh_tool import *
import sys
import math
import vtk
from matviz_vtk import *
from matplotlib import pyplot as plt
from scipy import interpolate
from mpl_toolkits.mplot3d import Axes3D
from scipy.spatial import Delaunay
from sympy.physics.quantum.circuitplot import pyplot
# from base_structures_images import height

def dirToString(dir):
  assert(dir > 0 and dir < 4)
  res = 'x'
  if dir == 2:
    res = 'y'
  if dir == 3:
    res = 'z'
   
  return res  

def profileLin(x1,x2,n):
  x = np.linspace(0, 1.0,n)
  return ((x2-x1) * x + x1) / 2.0

def profileCirc(x1,res):
  x = np.linspace(0, 1.0, res)
  r = 0.5-x1/2.0 # radius for circular holes not stiffness
  r2 = r**2
  cornerx = np.sqrt(r2/2.0)
  cornery = 1 - cornerx
  
  vec = np.zeros(res)
  
  for i in range(res):
    if x[i] <= cornerx:
      vec[i] = 0.5 - np.sqrt(r2 - x[i]**2)
    elif x[i] < 1 - cornerx:
      vec[i] = cornery - 0.5
    else:
      vec[i] = 0.5 - np.sqrt(r2 - (1-x[i])**2)
  
  return vec

def f03(t):
  return (1-t)**3
def f13(t):
  return 3*t*(1-t)**2
def f23(t):
  return 3*t**2*(1-t)
def f33(t):
  return t**3

def findGradOne(vec):
  idx = 0 # idx where grad is about 1
  stop = False
  while not stop and idx < len(vec):
    if vec[idx] <= 1:
      idx += 1
    else:
      stop = True
      
  if not stop:
    idx = len(vec) - 1
  
  # the idx before 1 might be closer to 1
  if abs(1-vec[idx-1]) < abs(1-vec[idx]):
    idx = idx-1 
   
  return idx

def contains_point(id_x,y,z,map):
  
  valx = (y-0.5)**2 + (z-0.5)**2
  phi = np.arccos(0.5*(y-0.5)/(0.5*np.sqrt(valx))) # angle between (0.5,0.5) and (y,z)
  
  if y < 0.5:
    phi = 2.0 * np.pi - phi
          
  r = map[int(phi/np.pi*180),id_x]
  
  if (valx <= r*r):
    return True
  
  return False

def spline_curve(x1, y1, res, bend, infoXml=None):
  rx = 0.5 - x1/2.0 # radius for center (0,1)
  ry = 0.5 - y1/2.0 # radius for center (0,1)
  
  P = np.transpose(np.array([[0,1-rx],[ry*bend,1-rx],[ry,1-rx*bend],[ry,1]]))
  t = np.linspace(0, 1.0, 2*res) # double over-sampling
  C = np.multiply.outer(P[:,0],f03(t)) + np.multiply.outer(P[:,1],f13(t)) + np.multiply.outer(P[:,2],f23(t)) + np.multiply.outer(P[:,3],f33(t))
  G = np.diff(C[1]) / np.diff(C[0]) # numerical finite differences
  G.resize(len(C[0]))
  
#   plt.plot(np.transpose(C[0,:]),np.transpose(C[1,:]), linewidth=5.0)
#   plt.plot(np.transpose(P[0,:]),np.transpose(P[1,:]), 'r', marker='o', markersize=15, linewidth=5.0)
#   plt.rcParams.update({'font.size': 18})
#   plt.savefig("spline_bend_" + str(bend) + ".png")
#   plt.show()
  
  # interpolate for regular spacing  
  fc = interpolate.interp1d(C[0],C[1])
  fg = interpolate.interp1d(C[0],G)
  
  v = np.linspace(0,ry,ry*res)
  c = fc(v)
  g = fg(v)
  
  idx = findGradOne(g)
  
  return c, idx, (v[idx], c[idx]), g[idx], P


def profileSpline(x1,y1,res,bend,dir,infoXml=None):
  assert(bend <= 1 and bend >= 0)

  vec = np.zeros(res)
  c, idx, point, g, cp = spline_curve(x1,y1,res,bend,infoXml)
#   print 'profileSpline',point
  vec[0:idx] = c[0:idx]
  vec[idx:res-idx] = c[idx] # constant value at position where grad is one
  vec[res-idx:res] = c[0:idx][::-1]
  
  # write info on x1,y1, control polygon to file
  if infoXml <> None:
    strDir = dirToString(dir)
      
    infoXml.write('  <profile type="b-spline" dir="' + strDir + '">\n')
    infoXml.write('    <bSpline type="spline" rad1 = "' + str(x1) + '" rad2="' + str(y1) + '" bend="' + str(bend) + '">\n')
    infoXml.write('      <controlPolygon>\n')
    infoXml.write('        <P1 x="' + str(cp[0][0]) + '" y="' + str(cp[1][0]) + '"/>\n')
    infoXml.write('        <P2 x="' + str(cp[0][1]) + '" y="' + str(cp[1][1]) + '"/>\n')
    infoXml.write('        <P3 x="' + str(cp[0][2]) + '" y="' + str(cp[1][2]) + '"/>\n')
    infoXml.write('        <P4 x="' + str(cp[0][3]) + '" y="' + str(cp[1][3]) + '"/>\n')
    infoXml.write('      </controlPolygon>\n')
    infoXml.write('    </bSpline>\n')
    infoXml.write('  </profile>\n\n')
    
  return vec-0.5

# @param height is second return value from profileSpline
# @param verbose 'bisec' to plot all cases
def profileSplineBisec(x1,y1,z1,res,bend,verbose,dir,infoXml):
  assert(bend <= 1 and bend >= 0)
  type = None
  ###### case 1 : bisec + quadratic --> biqua ###################
  # we have a left part from a=(0,x1) which is the average of the splines x1,y1 and x1,z1
  # where the curve has grad=1 we have point b
  # the right part is from b to x=0.5 where the heigt comes from the spline y1,z1 with grad=1 at point p
  # from the point p we determine the angle phi of this bisec profile function 
  
  # left part is same spline as orhtogonal spline up to point
  left, idx, b, gb, cp = spline_curve(x1, 0.5*(y1+z1), res, bend)
  assert(b[0] <= 0.5 and b[1] >= 0.5)

  # search for point p
  t, t, p, t, cp = spline_curve(y1,z1,res,bend)
  assert(p[0] <= 0.5 and p[1] >= 0.5)
  #height = p[1]
  height = np.sqrt((p[0]-0.5)**2+(p[1]-0.5)**2) + 0.5
  
  gamma = np.arccos((0.5-p[0])/(np.sqrt((0.5-p[0])**2 + (p[1]-0.5)**2)))
  phi = np.pi - gamma
  
#   print "x1: ",x1," y1: ",y1, " z1:",z1
#   print "p: ", p, " b: ",b
  
  # polynomial interpolation for right part from b to p
  lx = b[0]
  ly = b[1] 
  rx = 0.5
  A = np.array([
      [1, lx, lx**2, lx**3],
      [0, 1,  2*lx, 3*lx**2],
      [1, rx, rx**2, rx**3],
      [0, 1,  2*rx, 3*rx**2]
      ]) 
  rhs = np.array([ly, gb, height, 0])
  
  
  sol = np.linalg.solve(A, rhs)
  
  poly = np.poly1d(sol[::-1])
  v = np.linspace(lx,0.5,res/2-idx)
  right = poly(v)
  
  v = np.linspace(lx,0.5,res/2-idx)
  
  biqua = np.zeros(res)
  
  biqua[0:idx] = left[0:idx]
  biqua[idx:res/2] = right
  biqua[res/2:res] = biqua[0:res/2][::-1]
  
#   biqua -= 0.5
  
  #### case 2: b-spline --> bsp #############
  # if b with grad gb (approx 1) is too high for p such that the curve has a maximum within b and p we need to fallback to a b-spline from a to p
  
#     print 'need to fallback',np.amax(right)
    
  height = np.sqrt((p[0]-0.5)**2+(p[1]-0.5)**2) + 0.5
    
  #P = np.transpose(np.array([[0,0.5+x1/2.0],[0.5*bend,0.5+x1/2.0],[0.5-0.5*bend,height],[0.5,height]]))
  P = np.transpose(np.array([[0,0.5+x1/2.0],[0.5*bend,0.5+x1/2.0],[0.5-0.5*bend,height],[0.5,height]])) 
  t = np.linspace(0, 1.0, res) # double over-sampling
  C = np.multiply.outer(P[:,0],f03(t)) + np.multiply.outer(P[:,1],f13(t)) + np.multiply.outer(P[:,2],f23(t)) + np.multiply.outer(P[:,3],f33(t))
  # interpolate for regular spacing
  fc = interpolate.interp1d(C[0],C[1])
  v = np.linspace(0,0.5,res/2)
  c = fc(v)
  
  bsp = np.zeros(res)
  bsp[0:res/2] = c
  bsp[res/2:res] = c[::-1]

  #### case 3: linear --> lin ###########
  # in case undershooting for x1=0.9, y1=0.1, z1=0.1
  lin = profileLin(x1, x1, res)
  
  result = None
  
  if verbose == 'bisec':
    print "amax right: ",np.amax(right), " p: ",p[1], " amin vec: ",np.amin(biqua)
    plt.gcf().clear() 
#     plt.figure(figsize=(10,8))
    plt.plot(t,biqua-0.5,label='bspline-quadratic',linewidth=5.0)
    plt.plot(t,bsp-0.5,label='bspline',linewidth=5.0)
    plt.plot(t,lin,label='linear',linewidth=5.0)
#     plt.rcParams.update({'font.size': 18})
    plt.legend(loc='upper left', shadow=True,prop={'size':20})
#     plt.savefig("bisec_0.2.png")
    plt.show()
  
  if np.abs(np.amax(right) - height) < 1e-3:
    type = "biquadratic"
    if verbose == 'bisec':
      print "bisec: ",np.amax(right),height
      
    result = biqua - 0.5
  
  # in case we have undershooting for biqua and not for spline
  elif np.abs(np.amin(biqua) - x1/2.0) > 1e-3 and np.abs(np.amin(bsp - 0.5 - x1/2.0)) < 1e-3:
    type = "bSpline"
    if verbose == 'bisec':
      print "bspline: ",np.amin(biqua),x1/2.0
    result = bsp - 0.5
   
  # in case we have undershooting for biqua AND for spline
  else:
    type = "linear"
    if verbose == 'bisec':  
      print 'return lin'
    result = lin
  
  assert(type <> None)
  if infoXml <> None:
    strDir = dirToString(dir)
      
    infoXml.write('  <profile type="bisection" dir="' + strDir + '">\n')
    infoXml.write('    <bisectionSpline type="' + type + '" angle="' + str(phi*180/np.pi) + '">\n')
    infoXml.write('      <biquadratic coeff0="' + str(sol[0]) + '" coeff1="' + str(sol[1]) + '" coeff2="' + str(sol[2]) + '" coeff3="' + str(sol[3]) +'"/>\n')
    infoXml.write('      <bSpline type="spline" rad1 = "' + str(x1) + '" rad2="' + str(y1) + '" bend="' + str(bend) + '">\n')
    infoXml.write('        <controlPolygon>\n')
    infoXml.write('          <P1 x="' + str(P[0][0]) + '" y="' + str(P[1][0]) + '"/>\n')
    infoXml.write('          <P2 x="' + str(P[0][1]) + '" y="' + str(P[1][1]) + '"/>\n')
    infoXml.write('          <P3 x="' + str(P[0][2]) + '" y="' + str(P[1][2]) + '"/>\n')
    infoXml.write('          <P4 x="' + str(P[0][3]) + '" y="' + str(P[1][3]) + '"/>\n')
    infoXml.write('        </controlPolygon>\n')
    infoXml.write('      </bSpline>\n')
    infoXml.write('      <linear xStart="' + str(x1) + '" xEnd="' + str(x1) + '"/>\n')
    infoXml.write('    </bisectionSpline>\n')
    infoXml.write('  </profile>\n\n')
    
  return result, phi
# @return vector with profile or list of vectors
def profile(args,dir,infoXml=None):
  if (infoXml <> None):
    assert(args.profile == 'spline')
    
  if args.profile == 'linear':
    if dir == 1:
      return profileLin(args.x1, args.x2, args.res)
    if dir == 2:   
      return profileLin(args.y1, args.y2, args.res)
    if dir == 3:
      return profileLin(args.z1, args.z2, args.res)
    
  if args.profile == 'circular':
    if dir == 1:
      return profileCirc(args.x1, args.res)
    if dir == 2:   
      return profileCirc(args.y1, args.res)
    if dir == 3:
      return profileCirc(args.z1, args.res)  
    
  if args.profile == 'spline':
    if dir == 1:
      vec1 = profileSpline(args.x1, args.y1, args.res, args.bend, dir, infoXml)
      vec3 = profileSpline(args.x1, args.z1, args.res, args.bend, dir, infoXml)
      vec2,phi = profileSplineBisec(args.x1, args.y1, args.z1, args.res, args.bend, args.verbose, dir, infoXml)
      
#       t = np.linspace(0, 1.0, args.res)
#       plt.plot(t,vec1,label='x1y1', linewidth=5.0)
#       plt.plot(t,vec2,label='bisec', linewidth=5.0)
#       plt.plot(t,vec3,label='x1z1', linewidth=5.0)
#       plt.legend(loc='upper left', shadow=True)
#       plt.rcParams.update({'font.size': 18})
#       plt.savefig("3splines.png")
#       plt.show()

      return ((vec1,0), (vec2,phi), (vec3,np.pi/2.0))
    if dir == 2:
      vec1 = profileSpline(args.y1, args.x1, args.res, args.bend, dir, infoXml)
      vec2,phi = profileSplineBisec(args.y1, args.z1,args.x1, args.res, args.bend, args.verbose, dir, infoXml)
      vec3 = profileSpline(args.y1, args.z1, args.res, args.bend, dir, infoXml)

#       t = np.linspace(0, 1.0, args.res)
#       plt.plot(t,vec1,label='x1y1')
#       plt.plot(t,vec3,label='dir2 y1z1')
#       plt.ylim([0,1])
# #       plt.plot(t,vec3,label='x1z1')
# #       plt.plot(t,vec4,label='y1z1')
#       plt.legend(loc='upper left', shadow=True)
#       plt.show()
      
      return ((vec1,0), (vec2,phi), (vec3,np.pi/2.0))
    if dir == 3:
      vec1 = profileSpline(args.z1, args.y1, args.res, args.bend, dir, infoXml)
      vec2, phi = profileSplineBisec(args.z1, args.x1, args.y1, args.res, args.bend, args.verbose, dir, infoXml)
      vec3 = profileSpline(args.z1, args.x1, args.res, args.bend, dir, infoXml)
      return ((vec1,0), (vec2,phi), (vec3,np.pi/2.0))

# return information on profiles 
def create_profiles(args,infoXml=None):
  profiles = []
  if not args.skip_x:
    vec = profile(args,1,infoXml)
    profiles.append(vec)
  else:
    profiles.append(None)
  if not args.skip_y:
    vec = profile(args,2,infoXml)
    profiles.append(vec)
  else:
    profiles.append(None)
  if not args.skip_z:
    vec = profile(args,3,infoXml)
    profiles.append(vec)
  else:
    profiles.append(None)
    
  if args.verbose == "all_profiles":
    plt.gcf().clear()
    plt.gcf().subplots_adjust(bottom=0.15)
    t = np.linspace(0, 1.0, args.res)
    for vec in profiles[0:1]:
      if vec <> None:
        for v in vec:
          if type(v) == tuple:
            plt.plot(t,v[0],label=str(v[1]),linewidth=5.0)

    plt.rcParams.update({'font.size': 18})
    plt.ylim([0,0.5])
    plt.xlabel("x",labelpad=5)
    plt.ylabel("radius",labelpad=20)
    plt.savefig("profile_functions_" + str(args.stiffness) + ".png")

    plt.show()
    
  return profiles

# @param map: profile map
# @param numLines: number of lines on surface that we want to plot
# @param dir: direction (x,y or z) as (1,2,3)
def get_surface_lines(map,numLines,dir):
  assert(dir == 1  or dir == 2 or dir == 3)
  lenx = np.size(map,1)
  nodes = np.zeros((numLines, lenx, 3))
  for i,grad in enumerate(np.arange(0,360,360.0/numLines)):
    radii = map[grad,:]
    rad = grad/180.0*np.pi
     
    #transformation to cartesian coordinates
    X = radii*np.cos(rad) + 0.5*np.ones(lenx)
    Y = radii*np.sin(rad) + 0.5*np.ones(lenx)
    if dir == 1:
      nodes[i,:,0] = np.linspace(0.0,1.0,lenx)
      nodes[i,:,1] = X
      nodes[i,:,2] = Y
    if dir == 2:
      nodes[i,:,0] = Y
      nodes[i,:,1] = np.linspace(0.0,1.0,lenx)
      nodes[i,:,2] = X
       
    if dir == 3:
      nodes[i,:,0] = X
      nodes[i,:,1] = Y
      nodes[i,:,2] = np.linspace(0.0,1.0,lenx)
    
  return nodes

# for direction 'dir' with nodes 'nodes', find surface nodes on the left and right
# each surface point is also given a unique id
# @param pointId indicates the point id that already exists
def find_points_on_surface(nodes,dir,otherMap1,otherMap2, pointId):
  nodesLeft = []
  nodesRight = []
  pointId += 1
  
  if dir == 1:
    id0 = 0
    id1 = 2
    id2 = 0
    id3 = 1
  elif dir == 2:
    id0 = 1
    id1 = 2
    id2 = 0
    id3 = 1
  else: #dir == 3
    id0 = 1
    id1 = 2
    id2 = 0
    id3 = 2
        
  assert(dir == 1 or dir == 2 or dir == 3)
  res = nodes.shape[1]
  for numLine,prof in enumerate(nodes):
    lineLeft = [] # store points on one surface line temporally
    lineRight = []
    for i,p in enumerate(prof):
      if not contains_point(i,p[id0],p[id1],otherMap1) and not contains_point(i, p[id2], p[id3], otherMap2):
        if p[dir-1] < 0.5: # using plane through origin (0.5,0.5) to decide if surface point is on left or right side of structure
          lineLeft.append((p,pointId))
        else:
          lineRight.append((p,pointId))
        pointId += 1
        
    nodesLeft.append(lineLeft)
    nodesRight.append(lineRight[::-1])   
    
  return nodesLeft, nodesRight, pointId
# list is list of lists contains surface lines and all respective points
# base used for setting right ids
def define_triangles(list,cells):
  prevLine = None
  for i,line in enumerate(list):
    if i < len(list) - 1:
      prevLine = list[i+1]
    else:
      prevLine = list[0] # last line in list
      
    for j in range(0,len(line),1): # get coordinates of each point
      if (j+1 < np.size(line, 0)) and (j < np.size(prevLine, 0)):
        tri = vtk.vtkTriangle()
        tri.GetPointIds().SetId(0, line[j][1]) 
        tri.GetPointIds().SetId(1, prevLine[j][1])
        tri.GetPointIds().SetId(2, line[j+1][1])
        cells.InsertNextCell(tri)
         
#         print "triangle:",line[j][1],prevLine[j][1],line[j+1][1]
#         print "point",line[j][1], ": ",line[j][0][0],line[j][0][1],line[j][0][2]
#         print "point",prevLine[j][1], ": ",prevLine[j][0][0],prevLine[j][0][1],prevLine[j][0][2]
#         print "point",line[j+1][1], ": ",line[j+1][0][0],line[j+1][0][1],line[j+1][0][2]
      if (j+1 < np.size(line, 0)) and (j+1 < np.size(prevLine, 0)):
        tri = vtk.vtkTriangle()
        tri.GetPointIds().SetId(0, line[j+1][1])
        tri.GetPointIds().SetId(1, prevLine[j][1])
        tri.GetPointIds().SetId(2, prevLine[j+1][1])
        cells.InsertNextCell(tri)
    
#         print "triangle:",line[j+1][1],prevLine[j][1],prevLine[j+1][1]
       
#         print "point",line[j+1][1], ": ",line[j+1][0][0],line[j+1][0][1],line[j+1][0][2]
#         print "point",prevLine[j][1], ": ",prevLine[j][0][0],prevLine[j][0][1],prevLine[j][0][2]
#         print "point",prevLine[j+1][1], ": ",prevLine[j+1][0][0],prevLine[j+1][0][1],prevLine[j+1][0][2]

#     if i == 9:
#       break  
  return cells
  
def create_profiles_array(args,infoXml):
  res = args.res
  array = np.ones((res,res,res)) * (-1)
  vec = None
  t = np.linspace(0, 1.0, args.res)
  
  hf = plt.figure()
  ha = hf.add_subplot(111, projection='3d')
  ha.set_xlabel('X')
  ha.set_ylabel('Y')
  ha.set_zlabel('Z')
  
  profiles = create_profiles(args,infoXml)
  assert(len(profiles) == 3)
  
  surfNodesXLeft = []
  surfNodesYLeft = []
  surfNodesZLeft = []
  surfNodesXRight = []
  surfNodesYRight = []
  surfNodesZRight = []
  
  id = 0 # assign an id to each surface point 
  
  if args.target == "surface_mesh":
      assert(not args.skip_x and not args.skip_y and not args.skip_z)
      map_x = create_profile_map(profiles[0], args.res) if profiles[0] <> None else None
      map_y = create_profile_map(profiles[1], args.res) if profiles[1] <> None else None
      map_z = create_profile_map(profiles[2], args.res) if profiles[2] <> None else None
      
      idx = -1 # enumerating all surface nodes in order to create triangles
      
      # first dimension of nodes : surface lines
      # second dimension of nodes: resolution of unit cube
      # third dimension: tuple with x,y,z coordinate
      nodes = get_surface_lines(map_x, args.res_surf_lines, 1)
      
      #ha.scatter(nodes[:,:,0],nodes[:,:,1],nodes[:,:,2])
      
      #ha.scatter(nodes[0:4,:,0],nodes[0:4,:,1],nodes[0:4,:,2],color='red')
      
      surfNodesXLeft, surfNodesXRight, idx = find_points_on_surface(nodes, 1, map_y, map_z, idx)
      
#       for i,line in enumerate(surfNodesXRight):
#         for tuple in line:
#           ha.scatter(tuple[0][0],tuple[0][1],tuple[0][2])
#           ha.text(tuple[0][0],tuple[0][1],tuple[0][2],tuple[1])
#         if i == 10:
#           break
#       plt.show()
      
      out = open("surfNodesXLeft.coords","w")
      for i,line in enumerate(surfNodesXLeft):
        for tuple in line:
          out.write("id=" + str(tuple[1]) + ": x=" + str(tuple[0][0]) + " y=" + str(tuple[0][1]) + " z=" + str(tuple[0][2]) + "\n")
      out.close()
      nodes = get_surface_lines(map_y, args.res_surf_lines, 2)
      surfNodesYLeft, surfNodesYRight, idx = find_points_on_surface(nodes, 2, map_x, map_z, idx)  
            
      nodes = get_surface_lines(map_z, args.res_surf_lines, 3)
      surfNodesZLeft, surfNodesZRight, idx = find_points_on_surface(nodes, 3, map_x, map_y, idx)
      
      surfNodes = [surfNodesXLeft, surfNodesXRight, surfNodesYLeft, surfNodesYRight, surfNodesZLeft, surfNodesZRight]
      
      if args.show:
        for list in surfNodes:
          for line in list:
            for tuple in line: # tuple consists of on array with 3 coordinate components and 1 point id
              ha.scatter(tuple[0][0],tuple[0][1],tuple[0][2])
            break

        plt.show()
      
      # create vtk cells and points
      cells = vtk.vtkCellArray()
      points = vtk.vtkPoints()
      points.SetNumberOfPoints(idx)
      
      cells = define_triangles(surfNodesXLeft,cells)
      cells = define_triangles(surfNodesXRight,cells)
      cells = define_triangles(surfNodesYLeft,cells)
      cells = define_triangles(surfNodesYRight,cells)
      cells = define_triangles(surfNodesZLeft,cells)
      cells = define_triangles(surfNodesZRight,cells)
      
      for part in surfNodes:
        for line in part:
          for tuple in line:
            p = tuple[0]
            # set point with respectived ids as set in find_points_on_surface()
            points.SetPoint(tuple[1], (tuple[0][0], tuple[0][1], tuple[0][2]))
          
#       out = open("vtk_point_ids.txt","w")
#       for i in range(points.GetNumberOfPoints()):
#         p = points.GetPoint(i)
#         out.write("id=" + str(i) + ": x=" + str(p[0]) + " y=" + str(p[1]) + " z=" + str(p[2]) + "\n")
#       
#       out.close()  
            
      polydata = vtk.vtkPolyData()
      polydata.SetPoints(points)
      polydata.SetPolys(cells)
      
      if args.show:
        show_vtk(polydata, 1000, [], True)
        
      show_write_vtk(polydata,1000,"surface.vtp")
  else:
    for i in range(0,3):
      if profiles[i] == None:
        continue
      if args.target == "volume_mesh" or args.target == "volume_vtk":
        write_profile_to_array(array, profiles[i], i+1)
      if args.target == "3dlines":
        plot_3dlines(profiles[i], res, args.res_surf_lines, i+1, ha)
      if args.verbose == 'profile_map':
        create_profile_map(profiles[i], res, args.verbose, ha)
  
  if args.target == '3dlines':
    plt.show()
  
  return array

# creates map with info on profile depending on radius
# profile contains list of tuples with vector,angle and idx where constant part begins (bisec: res/2, orthogonal: grad is 1)
def create_profile_map(profile,res,verbose=None,ha=None):
  map = np.zeros((360, res))
  for i in range(0,res):
    f = give_interpolate_radius(profile,i)
    for alpha in range(0,360):
      rad = np.pi/180. * alpha
      map[alpha,i] = f(rad)
  
  if verbose == 'profile_map':
    X,Y = np.meshgrid(range(res),range(360))
    ha.plot_surface(X, Y, map)
    plt.show()
  
  return map

def plot_3dlines(profile,res,numLines,dir,ha):
  Z = np.linspace(0,2*np.pi,360)
  
  nodes = []
  
  map = create_profile_map(profile, res)
        
#   for ii in range(0,res,10):
#     radii = map[:,ii]
#     X = radii*np.cos(Z)+.5*np.ones(np.size(radii))
#     Y = radii*np.sin(Z)+.5*np.ones(np.size(radii))
#     if dir == 1:
#       ha.plot(ii*np.ones(np.size(X))/res,X,Y,'r')
#     if dir == 2:
#       ha.plot(Y,ii*np.ones(np.size(X))/res,X,'g')
#     if dir == 3:
#       ha.plot(X,Y,ii*np.ones(np.size(X))/res,'b')
      
  nodes = get_surface_lines(map, numLines, dir)
  for i in range(numLines):
    ha.plot(nodes[i,:,0],nodes[i,:,1],nodes[i,:,2],'r')
 
## give an interpolation for the radius within 0..2pi angle for the radius vectors at index idx
def give_interpolate_radius(vec, idx):
  if type(vec) == tuple:
    assert(len(vec) == 3)
    val = []
    rad = []
    
    val.append(vec[0][0][idx]) # 0 deg
    assert(vec[0][1] == 0.0)
    rad.append(0.0) 
    
    val.append(vec[1][0][idx]) # 45 de
    rad.append(vec[1][1])
    
    val.append(vec[2][0][idx]) # 90
    assert(vec[2][1] == np.pi/2)
    rad.append(np.pi/2)
    
    val.append(vec[1][0][idx]) # 90+45
    rad.append(np.pi-vec[1][1]) #np.pi/2+phi)
    
    val.append(vec[0][0][idx]) # 90+90=180
    rad.append(np.pi)
    
    val.append(vec[1][0][idx]) # 180+45
    rad.append(np.pi+vec[1][1])
    
    val.append(vec[2][0][idx]) # 180+90=270
    rad.append(6.0/4.0*np.pi)
    
    val.append(vec[1][0][idx]) # 270+45
    rad.append(2*np.pi-vec[1][1])
    
    val.append(vec[0][0][idx]) # 270+90=360=0 (full circle) 
    rad.append(2*np.pi)
      
    f = interpolate.interp1d(rad, val)
    return f
  else:
    assert(type(vec) <> tuple) # but a single vector
    f = interpolate.interp1d((0, 7),(vec[idx], vec[idx])) # 7 is close to 2*pi
    return f                            
                                       

# @param vec can be one vector or list of vector
def write_profile_to_array(array,vec,dir):
  res = array.shape[0]
  h = 1.0/res
  assert(dir >=1 and dir <=3)
  
  map = create_profile_map(vec, res)
#   plt.savefig("spline_bend_" + str(bend) + ".png")
      
  for i in range(0,res):
    for j in range(0,res):
      for k in range(0,res):
        y = j * h + h / 2.0
        z = k * h + h / 2.0
        valx = (y-0.5)**2 + (z-0.5)**2
        
        phi = np.arccos(0.5*(y-0.5)/(0.5*np.sqrt(valx))) # angle between (0.5,0.5) and (y,z)
        if y < 0.5:
          phi = 2*np.pi - phi
          
#         p = f(phi)
        p = map[int(phi/np.pi*180),i]
        if (valx <= p*p):
          if dir == 1:
            array[i,j,k] = dir
          if dir == 2:
            array[j,i,k] = dir
          if dir == 3:
            array[k,j,i] = dir
