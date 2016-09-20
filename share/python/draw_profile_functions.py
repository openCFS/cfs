#!/usr/bin/env python
import numpy as np
#from PIL import Image
# from mesh_tool import *
import sys
import math
from matplotlib import pyplot as plt
from scipy import interpolate
from mpl_toolkits.mplot3d import Axes3D
# from base_structures_images import height

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

def spline_curve(x1, y1, res, bend):
  rx = 0.5 - x1/2.0 # radius for center (0,1)
  ry = 0.5 - y1/2.0 # radius for center (0,1)
  
  P = np.transpose(np.array([[0,1-rx],[ry*bend,1-rx],[ry,1-rx*bend],[ry,1]]))
  t = np.linspace(0, 1.0, 2*ry*res) # double over-sampling
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
  return c, idx, (v[idx], c[idx]), g[idx]


def profileSpline(x1,y1,res,bend):
  assert(bend <= 1 and bend >= 0)

  vec = np.zeros(res)
  c, idx, point, g = spline_curve(x1,y1,res,bend)
#   print 'profileSpline',point
  vec[0:idx] = c[0:idx]
  vec[idx:res-idx] = c[idx] # constant value at position where grad is one
  vec[res-idx:res] = c[0:idx][::-1]
  
  return vec-0.5

# @param height is second return value from profileSpline
# @param verbose 'bisec' to plot all cases
def profileSplineBisec(x1,y1,z1,res,bend,verbose):
  assert(bend <= 1 and bend >= 0)
  ###### case 1 : bisec + quadratic --> biqua ###################
  # we have a left part from a=(0,x1) which is the average of the splines x1,y1 and x1,z1
  # where the curve has grad=1 we have point b
  # the right part is from b to x=0.5 where the heigt comes from the spline y1,z1 with grad=1 at point p
  # from the point p we determine the angle phi of this bisec profile function 
  
  # left part is same spline as orhtogonal spline up to point
  left, idx, b, gb = spline_curve(x1, 0.5*(y1+z1), res, bend)
  assert(b[0] <= 0.5 and b[1] >= 0.5)

  # serch for point p
  t, t, p, t = spline_curve(y1,z1,res,bend)
  assert(p[0] <= 0.5 and p[1] >= 0.5)
  height = p[1]
  
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
  
  if verbose == 'bisec':
    print "amax right: ",np.amax(right), " p: ",p[1], " amin vec: ",np.amin(biqua) 
    plt.plot(t,biqua-0.5,label='bspline-quadratic')
    plt.plot(t,bsp-0.5,label='bspline')
    plt.plot(t,lin,label='linear')
    plt.legend(loc='upper left', shadow=True)
    plt.show()
    

  if np.amax(right) <= p[1]:
    if verbose == 'bisec':
      print "bisec: ",np.amax(right),p[1]
    return biqua - 0.5,phi
  
  if np.amin(biqua) >= x1/2.0:
    if verbose == 'bisec':
      print "bspline: ",np.amin(biqua)
    return bsp - 0.5,phi
  
  if verbose == 'bisec':  
    print 'return lin'
  return lin, phi
# @return vector with profile or list of vectors
def profile(args,dir):
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
      vec1 = profileSpline(args.x1, args.y1, args.res, args.bend)
      vec3 = profileSpline(args.x1, args.z1, args.res, args.bend)
      vec2,phi = profileSplineBisec(args.x1, args.y1, args.z1, args.res, args.bend, args.verbose)
      
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
      vec1 = profileSpline(args.y1, args.x1, args.res, args.bend)
      vec2,phi = profileSplineBisec(args.y1, args.z1,args.x1, args.res, args.bend, args.verbose)
      vec3 = profileSpline(args.y1, args.z1, args.res, args.bend)

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
      vec1 = profileSpline(args.z1, args.y1, args.res, args.bend)
      vec2, phi = profileSplineBisec(args.z1, args.x1, args.y1, args.res, args.bend, args.verbose)
      vec3 = profileSpline(args.z1, args.x1, args.res, args.bend)
      return ((vec1,0), (vec2,phi), (vec3,np.pi/2.0))

# return information on profiles 
def create_profiles(args):
  profiles = []
  if not args.skip_x:
    vec = profile(args,1)
    profiles.append(vec)
  else:
    profiles.append(None)
  if not args.skip_y:
    vec = profile(args,2)
    profiles.append(vec)
  else:
    profiles.append(None)
  if not args.skip_z:
    vec = profile(args,3)
    profiles.append(vec)
  else:
    profiles.append(None)
  
  if args.verbose == "all_profiles":
    t = np.linspace(0, 1.0, args.res)
    for vec in profiles:
      if vec <> None:
        for v in vec:
          if type(v) == tuple:
            plt.plot(t,v[0],label=str(v[1]))
            
    plt.show()
    
  return profiles

def create_profiles_array(args):
  res = args.res
  array = np.ones((res,res,res)) * (-1)
  vec = None
  t = np.linspace(0, 1.0, args.res)
  
  hf = plt.figure()
  ha = hf.add_subplot(111, projection='3d')
  
  profiles = create_profiles(args)
  assert(len(profiles) == 3)
  for i in range(0,3):
    if profiles[i] == None:
      continue
    if args.target == "volume_mesh" or args.target == "volume_vtk":
      write_profile_to_array(array, profiles[i], i+1)
    if args.target == "3dlines":
      plot_3dlines(profiles[i], res, i+1, ha)
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

def plot_3dlines(profile,res,dir,ha):
  Z = np.linspace(0,2*np.pi,360)
  
  #ax = plt.axes(polar=True)
  #plt.plot(Z,o(Z))
  #plt.show()
  
  map = create_profile_map(profile, res)
        
  for ii in range(0,res,2):
    radii = map[:,ii]
    X,Y = radii*np.cos(Z),radii*np.sin(Z)
    if dir == 1:
      ha.plot(X,Y,ii*np.ones(np.size(X))/res-.5*np.ones(np.size(X)),'b')
    if dir == 2:
      ha.plot(ii*np.ones(np.size(X))/res-.5*np.ones(np.size(X)),X,Y,'r')
    if dir == 3:
      ha.plot(Y,ii*np.ones(np.size(X))/res-.5*np.ones(np.size(X)),X,'g')
  
  for grad in range(0,360,5):
    radii = map[grad,:]
    rad = grad/180.0*np.pi
    #transformation to cartesian coordinates
    X,Y = radii*np.cos(rad),radii*np.sin(rad)
    if dir == 1:
      ha.plot(X,Y,np.linspace(-.5,.5,np.size(X)),'b')
    if dir == 2:
      ha.plot(np.linspace(-.5,.5,np.size(X)),X,Y,'r')
    if dir == 3:
      ha.plot(Y,np.linspace(-.5,.5,np.size(X)),X,'g')
 
    #  plt.figure(figsize=(10,10))
    #  ax = plt.axes(polar=True)
    #  theta = np.linspace(0, 2*np.pi,360)
    #  plt.plot(theta,map[:,0],linewidth=5.0)
    #  plt.rcParams.update({'font.size': 18})
    #  plt.savefig("xpart_x_0.png")
    #  plt.plot(theta,map[:,res/4],linewidth=5.0)
    #  plt.plot(theta,map[:,res/2],linewidth=5.0)
    #  plt.savefig("xpart_all.png")

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
        
        phi = np.arccos(0.5*(y-0.5)/(0.5*np.sqrt(valx))) # angle between (0.5,0) and (y,z)
        if y < 0:
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
  