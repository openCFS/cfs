import numpy as np
import matviz_vtk, vtk
import pymesh
import sys, math

edgeTable = None
triTable = None
  # the coordinate system is a left-hand sided coorddinate system with y pointing towards the viewer  
  # 
  #  z ^  
  #    |
  #    | 
  #    | 
  #  y(^)---> x 
  #  
  # This are the node numbers if we have only one element. The .mesh file will be transformed to 1-based               
  # x is the fastet variable, z is the slowest variable 
  # 
  #             e4 
  #       4 ----------- 5          
  #      /|            / | 
  #   e7/ |   e6    e5/  | 
  #   6/-------------/7  |e9 
  #   |   |e8        |   | 
  # e11|   |      e10|   | 
  #   |   0 ---------|---/1
  #   |  /    e0     |  /
  #   | /e3          | / e1  
  #   |/             |/
  #  3 -------------2/
  #         e2 
def marching_cubes(voxels,spacing,points,triangles,normals,thresh = 0.5,cube_size=2,offset=None):
  from draw_profile_functions import grid_to_cartesian_coords
  hx, hy, hz = spacing
  nx, ny, nz = voxels.shape
  print("\ncube_size:",cube_size)
  print("hx,hy,hz:",hx,hy,hz)
  
  if offset is None:
    offset = np.array([0,0,0])
  
  global edgeTable, triTable
  #create lookup tables once
  if edgeTable is None and triTable is None:
    edgeTable, triTable = create_lookup()
  
  cube_size = int(cube_size)
#   print("cube_size:",cube_size)
#   
  assert (0 < hx and 0 < hy and 0 < hz)
  for i in range(0,nx-cube_size,cube_size):
    for j in range(0,ny-cube_size,cube_size):
      for k in range(0,nz-cube_size,cube_size):
        
        # bottom/top face:
        # NW  ----------- NE  
        #  |              |    ------>i
        #  |              |   |
        #  |              |   |
        #  |              |   | j
        #  |              |   
        # SW  ----------- SE   
        bottom_NW = (i,j,k)
        bottom_NE = (i+cube_size,j,k)
        bottom_SE = (i+cube_size,j+cube_size,k)
        bottom_SW = (i,j+cube_size,k)
        top_NW = (i,j,k+cube_size)
        top_NE = (i+cube_size,j,k+cube_size)
        top_SE = (i+cube_size,j+cube_size,k+cube_size)
        top_SW = (i,j+cube_size,k+cube_size)
        
#         print(bottom_NW,bottom_NE,bottom_SE,bottom_SW)
        
        local_verts = []
        local_verts.append(Vertex(0,bottom_NW,voxels[bottom_NW]))
        local_verts.append(Vertex(1,bottom_NE,voxels[bottom_NE]))
        local_verts.append(Vertex(2,bottom_SE,voxels[bottom_SE]))
        local_verts.append(Vertex(3,bottom_SW,voxels[bottom_SW]))
        local_verts.append(Vertex(4,top_NW,voxels[top_NW]))
        local_verts.append(Vertex(5,top_NE,voxels[top_NE]))
        local_verts.append(Vertex(6,top_SE,voxels[top_SE]))
        local_verts.append(Vertex(7,top_SW,voxels[top_SW]))
        
        # example: vertex 3 was below the isosurface, cubeindex would equal 0000 1000 or 8
        cube_idx = 0
        cube_idx |= 1 if local_verts[0].value > thresh else 0
        cube_idx |= 2 if local_verts[1].value > thresh else 0
        cube_idx |= 4 if local_verts[2].value > thresh else 0
        cube_idx |= 8 if local_verts[3].value > thresh else 0
        cube_idx |= 16 if local_verts[4].value > thresh else 0
        cube_idx |= 32 if local_verts[5].value > thresh else 0
        cube_idx |= 64 if local_verts[6].value > thresh else 0
        cube_idx |= 128 if local_verts[7].value > thresh else 0
        
        edge_idx = edgeTable[cube_idx]
        
        # we have 12 edges, find vertices where surface intersects marching cube
        if cube_idx == 0 or cube_idx == 255:
          continue
        
#         create_triangles(local_verts,thresh,points,triangles,normals,voxels,cube_size)
        verts, faces, norms = create_triangles2(voxels[i:i+cube_size+1,j:j+cube_size+1,k:k+cube_size+1],(float(hx),float(hy),float(hz)),cube_size)
        
        angles = []
        for f in faces:
          n0 = norms[f[0]]
          n1 = norms[f[1]]
          n2 = norms[f[2]]
          
          angles.append(calc_angle(n0, n1))
          angles.append(calc_angle(n0, n2))
          angles.append(calc_angle(n1, n2))
        
        if cube_size == 1 or np.any(np.array(angles) < 50):
          offset_p = grid_to_cartesian_coords((i+offset[0],j+offset[1],k+offset[2]),None,(hx,hy,hz))
          if cube_size == 1:
            offset_p -= np.array([0,hy/2,0])  
          offset_f = len(points)
          verts = [v + np.array(offset_p) for v in verts]
          points.extend(verts)
          faces = [v + np.array(offset_f) for v in faces]
          triangles.extend(faces)
          normals.extend(norms)
        else:
#           print("refining ", cube_size, " angles:",angles)
          marching_cubes(voxels[i:i+cube_size+1,j:j+cube_size+1,k:k+cube_size+1], (hx,hy,hz), points, triangles, normals, thresh, int(cube_size/2),offset=(i,j,k))

  if cube_size == 2:
#     for i in range(len(points)):
#       points[i] = grid_to_cartesian_coords(points[i],None,(hx,hy,hz))
  #     #normals[i] *= np.array([hx,hy,hz])  
    pd = matviz_vtk.fill_vtk_polydata(points,triangles)
  #   
    pointNormalsArray = vtk.vtkDoubleArray()
    pointNormalsArray.SetNumberOfComponents(3)
    pointNormalsArray.SetName("normals")
    #pointNormalsArray.SetNumberOfTuples(len(points))
     
    for n in normals:
      pointNormalsArray.InsertNextTuple3(-n[0], -n[1], -n[2])
     
    pd.GetPointData().AddArray(pointNormalsArray)  
    pd.GetPointData().SetActiveScalars("normals")
    matviz_vtk.show_write_vtk(pd, 10, "marching_cubes.vtp")   
  
#   print(triangles)
#   mesh = pymesh.form_mesh(np.asarray(points),np.asarray(triangles))
#   mesh, _ = pymesh.remove_duplicated_vertices(mesh)
#   print("#points:",len(mesh.vertices)," #faces:",len(mesh.faces))
#   print("normals:",normals)

#   pd = matviz_vtk.fill_vtk_polydata(mesh.vertices,mesh.faces)
#   matviz_vtk.show_write_vtk(pd, 10, "marching_cubes.vtp")
  
  return points, triangles
#   return mesh.vertices, mesh.faces, compute_vertex_normals(mesh.vertices,mesh.faces)
  #return 
  # assume each voxel center is a lattice node
  # -> marching cubes grid is shifted by hx/2,hy/2,hz/2 and 1 elem smaller in each dim
  # for each of the 8 lattice node, find out which ones belong to the structure
  # -> 8 bit code

def VertexInterp2(thresh,p1,p2):
  valp1 = p1.value
  valp2 = p2.value
  if np.isclose(np.abs(thresh-valp1),0):
    return p1
  if np.isclose(np.abs(thresh-valp2),0):
    return p2
  if np.isclose(np.abs(valp1-valp2),0):
    return p1
  
  mu = (thresh - valp1) / (valp2 - valp1)
  ncoords = p1.coords + mu * (p2.coords-p1.coords)
  
  return ncoords, mu

# for given (local) edge number, give (local) vertex ids defining this edge 
def edgeToVertices(edge):
  assert(type(edge) is int)
  if edge == 0:
    return 0,1
  if edge == 1:
    return 1,2 
  if edge == 2:
    return  2,3
  if edge == 3:
    return 3,0
  if edge == 4:
    return 4,5
  if edge == 5:
    return 5,6
  if edge == 6:
    return 6,7
  if edge == 7:
    return 7,4
  if edge == 8:
    return 0,4
  if edge == 9:
    return 1,5
  if edge == 10:
    return 2,6
  if edge == 11:
    return 3,7

# linearly interpolate intersection points with isovalue 'thresh'
# p1 and p2 are vertices of an edge and v1 and v2 scalar values at respective vertex
# intersection point P = p1 + (thresh - v1) * (p2-p1) / (v2-v1)
def VertexInterp(thresh,p1,p2):
  
  if smaller(p2,p1):
    # swap p1 and p2
    p1, p2 = p2, p1
  
  valp1 = p1.value
  valp2 = p2.value
  p1 = np.asarray(p1.coords)
  p2 = np.asarray(p2.coords)  
    
  if np.isclose(valp2-valp1,0):
    return p1
  else:
    mu = (thresh - valp1) / (valp2 - valp1)  
    return p1 + mu * (p2 - p1), mu   
    
# kind of comparing point with descending priority for x to z 
def smaller(p1,p2):
  p1 = p1.coords
  p2 = p2.coords
  if p1[0] < p2[0]:
    return True
  elif p1[0] > p2[0]:
    return False
  
  if p1[1] < p2[1]:
    return True
  elif p1[1] > p2[1]:
    return False
  
  if p1[2] < p2[2]:
    return True
  elif p1[2] > p2[2]:
    return False
  
  return False  

def create_triangles(vertices,thresh,points,triangles,normals,voxels,cube_size):
  hx = 1/voxels.shape[0]
  hy = 1/voxels.shape[1]
  hz = 1/voxels.shape[2]
  # 8 vertices of the Marching cube
  assert(len(vertices) == 8)
  assert(edgeTable is not None and triTable is not None)
  assert(cube_size > 0)
  tmp_triangles = []
  
  # example: vertex 3 was below the isosurface, cubeindex would equal 0000 1000 or 8
  cube_idx = 0
  cube_idx |= 1 if vertices[0].value > thresh else 0
  cube_idx |= 2 if vertices[1].value > thresh else 0
  cube_idx |= 4 if vertices[2].value > thresh else 0
  cube_idx |= 8 if vertices[3].value > thresh else 0
  cube_idx |= 16 if vertices[4].value > thresh else 0
  cube_idx |= 32 if vertices[5].value > thresh else 0
  cube_idx |= 64 if vertices[6].value > thresh else 0
  cube_idx |= 128 if vertices[7].value > thresh else 0
  
  # list of interpolated vertices (vertices of new triangles)
  vertlist = [None] * 12
  # code binary number with 12 digits (for 12 edges)
  # digit is 1 if triangle vertex lies on respective edge, else 0
  # e.g. edges 2,3,11 : edge_idx = 1000 0000 1100 and  
  edge_idx = edgeTable[cube_idx]
  
  # we have 12 edges, find vertices where surface intersects marching cube
  if edge_idx == 0:
    return
  
  Gx = Gy = Gz = 0

  for c in range(12):
    if edge_idx & 2**c:
      edge0, edge1 = edgeToVertices(c)
      interp, mu = VertexInterp(thresh,vertices[edge0],vertices[edge1])
      
      Gx0 = Gx1 = 0
      e = np.array([1,0,0])
      if vertices[edge0].coords[0] < voxels.shape[0]-1 and vertices[edge0].coords[0] > 0:
        Gx0 = voxels[tuple(vertices[edge0].coords + e)] - voxels[tuple(vertices[edge0].coords - e)]
      if vertices[edge1].coords[0] < voxels.shape[0]-1 and vertices[edge1].coords[0] > 0:
        Gx1 = voxels[tuple(vertices[edge1].coords + e)] - voxels[tuple(vertices[edge1].coords - e)]
      Gx = mu * Gx0 + (1-mu) * Gx1
      
      Gy0 = Gy1 = 0
      e = np.array([0,1,0])
      if vertices[edge0].coords[1] < voxels.shape[1]-1 and vertices[edge0].coords[1] > 0:
        Gy0 = voxels[tuple(vertices[edge0].coords + e)] - voxels[tuple(vertices[edge0].coords - e)]
      if vertices[edge1].coords[1] < voxels.shape[1]-1 and vertices[edge1].coords[1] > 0:
        Gy1 = voxels[tuple(vertices[edge1].coords + e)] - voxels[tuple(vertices[edge1].coords - e)]
      Gy = mu * Gy0 + (1-mu) * Gy1
      
      Gz0 = Gz1 = 0
      e = np.array([0,0,1])
      if vertices[edge0].coords[2] < voxels.shape[2]-1 and vertices[edge0].coords[2] > 0:
        Gz0 = voxels[tuple(vertices[edge0].coords + e)] - voxels[tuple(vertices[edge0].coords - e)]
      if vertices[edge1].coords[2] < voxels.shape[2]-1 and vertices[edge1].coords[2] > 0:
        Gz1 = voxels[tuple(vertices[edge1].coords + e)] - voxels[tuple(vertices[edge1].coords - e)]
      Gz = mu * Gz0 + (1-mu) * Gz1
      
      vertlist[c] = (interp,(Gx,Gy,Gz))
      
  t = 0
  tmp_triangles = []
  tmp_points = []
  tmp_normals = []
  angles = []
  not_refined = True
  while (triTable[cube_idx][t] != -1):
    n0 = vertlist[triTable[cube_idx][t]][1] * np.array([hx,hy,hz])
    n1 = vertlist[triTable[cube_idx][t+1]][1] * np.array([hx,hy,hz])
    n2 = vertlist[triTable[cube_idx][t+2]][1] * np.array([hx,hy,hz])
    
    angles.append(calc_angle(n0, n1))
    angles.append(calc_angle(n0, n2))
    angles.append(calc_angle(n1, n2))
    id0 = len(points) + len(tmp_points)
    tmp_points.append(vertlist[triTable[cube_idx][t]][0])
    tmp_normals.append(n0)
     
    id1 = id0 + 1
    tmp_points.append(vertlist[triTable[cube_idx][t+1]][0])
    tmp_normals.append(n1)
     
    id2 = id0 + 2
    tmp_points.append(vertlist[triTable[cube_idx][t+2]][0])
    tmp_normals.append(n2)
    
    tmp_triangles.append((id0,id1,id2))
    
    t += 3
    
  if cube_size == 1 or np.any(np.array(angles) < 50):  
    points.extend(tmp_points)
    triangles.extend(tmp_triangles)
    normals.extend(tmp_normals)
  else:  
#       print("refining:",vertices[0]," size:",cube_size/2," ",vertices[6])
    marching_cubes(voxels,(hx,hy,hz),points,triangles,normals,cube_size=cube_size/2,bottom_NW=vertices[0].coords,top_SE=vertices[6].coords+np.array([1,1,1]))    

def calc_angle(n0,n1):
  if not sameVector(n0, n1):
    dot01 = np.dot(n0,n1)
    normed = np.linalg.norm(n0)*np.linalg.norm(n1)
    if not np.isclose(normed, 0):
      return np.degrees(math.acos(dot01/normed))
    
  return 90

def compute_vertex_normals(verts,faces,linear_weights=False):
  normals = np.zeros_like(verts)
  weight = 1 # for weighted sum
  for t in faces:
    # compute surface normal for each face and add it to vertex normals
    # vertex normal is weighted sum over normals of all adjacent faces
    # use face area as weight
    
    # first edge of a triangle
    v0 = verts[t[1]] - verts[t[0]]
    v1 = verts[t[2]] - verts[t[0]]
    n = np.cross(v0,v1)
    
    if not linear_weights:
      weight = np.linalg.norm(n)
      
    normals[t[0]] += n * weight
    normals[t[1]] += n * weight
    normals[t[2]] += n * weight
  
  # normalize sum  
  normals = [v/np.linalg.norm(v) for v in normals]
  
  return normals    

def sameVector(v1,v2):
  assert(len(v1) == len(v2))
  for i in range(len(v1)):
    if v1[i] != v2[i]:
      return False
  
  return True   
    
class Vertex():
  def __init__(self,id,coords,val):
    self.id = id
    self.coords = np.array(coords)
    self.value = val
    
  def __str__(self):
    return "vertex " + str(self.id) + ": coords=" + str(self.coords) + " val=" + str(self.value)  
    
class Edge():
  def __init__(self,id,v1,v2):
    self.id = id
    self.vertices = [v1,v2]
    
# for each 8bit config, store edge configuration
# e.g. 00000001 -> 1 means 1 lattice node is solid -> triangle has vertices on edges 0,3,8 
def create_lookup():
  edgeTable = [
    0x0  , 0x109, 0x203, 0x30a, 0x406, 0x50f, 0x605, 0x70c,
    0x80c, 0x905, 0xa0f, 0xb06, 0xc0a, 0xd03, 0xe09, 0xf00,
    0x190, 0x99 , 0x393, 0x29a, 0x596, 0x49f, 0x795, 0x69c,
    0x99c, 0x895, 0xb9f, 0xa96, 0xd9a, 0xc93, 0xf99, 0xe90,
    0x230, 0x339, 0x33 , 0x13a, 0x636, 0x73f, 0x435, 0x53c,
    0xa3c, 0xb35, 0x83f, 0x936, 0xe3a, 0xf33, 0xc39, 0xd30,
    0x3a0, 0x2a9, 0x1a3, 0xaa , 0x7a6, 0x6af, 0x5a5, 0x4ac,
    0xbac, 0xaa5, 0x9af, 0x8a6, 0xfaa, 0xea3, 0xda9, 0xca0,
    0x460, 0x569, 0x663, 0x76a, 0x66 , 0x16f, 0x265, 0x36c,
    0xc6c, 0xd65, 0xe6f, 0xf66, 0x86a, 0x963, 0xa69, 0xb60,
    0x5f0, 0x4f9, 0x7f3, 0x6fa, 0x1f6, 0xff , 0x3f5, 0x2fc,
    0xdfc, 0xcf5, 0xfff, 0xef6, 0x9fa, 0x8f3, 0xbf9, 0xaf0,
    0x650, 0x759, 0x453, 0x55a, 0x256, 0x35f, 0x55 , 0x15c,
    0xe5c, 0xf55, 0xc5f, 0xd56, 0xa5a, 0xb53, 0x859, 0x950,
    0x7c0, 0x6c9, 0x5c3, 0x4ca, 0x3c6, 0x2cf, 0x1c5, 0xcc ,
    0xfcc, 0xec5, 0xdcf, 0xcc6, 0xbca, 0xac3, 0x9c9, 0x8c0,
    0x8c0, 0x9c9, 0xac3, 0xbca, 0xcc6, 0xdcf, 0xec5, 0xfcc,
    0xcc , 0x1c5, 0x2cf, 0x3c6, 0x4ca, 0x5c3, 0x6c9, 0x7c0,
    0x950, 0x859, 0xb53, 0xa5a, 0xd56, 0xc5f, 0xf55, 0xe5c,
    0x15c, 0x55 , 0x35f, 0x256, 0x55a, 0x453, 0x759, 0x650,
    0xaf0, 0xbf9, 0x8f3, 0x9fa, 0xef6, 0xfff, 0xcf5, 0xdfc,
    0x2fc, 0x3f5, 0xff , 0x1f6, 0x6fa, 0x7f3, 0x4f9, 0x5f0,
    0xb60, 0xa69, 0x963, 0x86a, 0xf66, 0xe6f, 0xd65, 0xc6c,
    0x36c, 0x265, 0x16f, 0x66 , 0x76a, 0x663, 0x569, 0x460,
    0xca0, 0xda9, 0xea3, 0xfaa, 0x8a6, 0x9af, 0xaa5, 0xbac,
    0x4ac, 0x5a5, 0x6af, 0x7a6, 0xaa , 0x1a3, 0x2a9, 0x3a0,
    0xd30, 0xc39, 0xf33, 0xe3a, 0x936, 0x83f, 0xb35, 0xa3c,
    0x53c, 0x435, 0x73f, 0x636, 0x13a, 0x33 , 0x339, 0x230,
    0xe90, 0xf99, 0xc93, 0xd9a, 0xa96, 0xb9f, 0x895, 0x99c,
    0x69c, 0x795, 0x49f, 0x596, 0x29a, 0x393, 0x99 , 0x190,
    0xf00, 0xe09, 0xd03, 0xc0a, 0xb06, 0xa0f, 0x905, 0x80c,
    0x70c, 0x605, 0x50f, 0x406, 0x30a, 0x203, 0x109, 0x0
  ]
  
  triTable = [
    (-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
    (0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
    (0, 1, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
    (1, 8, 3, 9, 8, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
    (1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
    (0, 8, 3, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
    (9, 2, 10, 0, 2, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
    (2, 8, 3, 2, 10, 8, 10, 9, 8, -1, -1, -1, -1, -1, -1, -1),
    (3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
    (0, 11, 2, 8, 11, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
    (1, 9, 0, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
    (1, 11, 2, 1, 9, 11, 9, 8, 11, -1, -1, -1, -1, -1, -1, -1),
    (3, 10, 1, 11, 10, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
    (0, 10, 1, 0, 8, 10, 8, 11, 10, -1, -1, -1, -1, -1, -1, -1),
    (3, 9, 0, 3, 11, 9, 11, 10, 9, -1, -1, -1, -1, -1, -1, -1),
    (9, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
    (4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
    (4, 3, 0, 7, 3, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
    (0, 1, 9, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
    (4, 1, 9, 4, 7, 1, 7, 3, 1, -1, -1, -1, -1, -1, -1, -1),
    (1, 2, 10, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
    (3, 4, 7, 3, 0, 4, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1),
    (9, 2, 10, 9, 0, 2, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1),
    (2, 10, 9, 2, 9, 7, 2, 7, 3, 7, 9, 4, -1, -1, -1, -1),
    (8, 4, 7, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
    (11, 4, 7, 11, 2, 4, 2, 0, 4, -1, -1, -1, -1, -1, -1, -1),
    (9, 0, 1, 8, 4, 7, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1),
    (4, 7, 11, 9, 4, 11, 9, 11, 2, 9, 2, 1, -1, -1, -1, -1),
    (3, 10, 1, 3, 11, 10, 7, 8, 4, -1, -1, -1, -1, -1, -1, -1),
    (1, 11, 10, 1, 4, 11, 1, 0, 4, 7, 11, 4, -1, -1, -1, -1),
    (4, 7, 8, 9, 0, 11, 9, 11, 10, 11, 0, 3, -1, -1, -1, -1),
    (4, 7, 11, 4, 11, 9, 9, 11, 10, -1, -1, -1, -1, -1, -1, -1),
    (9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
    (9, 5, 4, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
    (0, 5, 4, 1, 5, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
    (8, 5, 4, 8, 3, 5, 3, 1, 5, -1, -1, -1, -1, -1, -1, -1),
    (1, 2, 10, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
    (3, 0, 8, 1, 2, 10, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1),
    (5, 2, 10, 5, 4, 2, 4, 0, 2, -1, -1, -1, -1, -1, -1, -1),
    (2, 10, 5, 3, 2, 5, 3, 5, 4, 3, 4, 8, -1, -1, -1, -1),
    (9, 5, 4, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
    (0, 11, 2, 0, 8, 11, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1),
    (0, 5, 4, 0, 1, 5, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1),
    (2, 1, 5, 2, 5, 8, 2, 8, 11, 4, 8, 5, -1, -1, -1, -1),
    (10, 3, 11, 10, 1, 3, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1),
    (4, 9, 5, 0, 8, 1, 8, 10, 1, 8, 11, 10, -1, -1, -1, -1),
    (5, 4, 0, 5, 0, 11, 5, 11, 10, 11, 0, 3, -1, -1, -1, -1),
    (5, 4, 8, 5, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1),
    (9, 7, 8, 5, 7, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
    (9, 3, 0, 9, 5, 3, 5, 7, 3, -1, -1, -1, -1, -1, -1, -1),
    (0, 7, 8, 0, 1, 7, 1, 5, 7, -1, -1, -1, -1, -1, -1, -1),
    (1, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
    (9, 7, 8, 9, 5, 7, 10, 1, 2, -1, -1, -1, -1, -1, -1, -1),
    (10, 1, 2, 9, 5, 0, 5, 3, 0, 5, 7, 3, -1, -1, -1, -1),
    (8, 0, 2, 8, 2, 5, 8, 5, 7, 10, 5, 2, -1, -1, -1, -1),
    (2, 10, 5, 2, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1),
    (7, 9, 5, 7, 8, 9, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1),
    (9, 5, 7, 9, 7, 2, 9, 2, 0, 2, 7, 11, -1, -1, -1, -1),
    (2, 3, 11, 0, 1, 8, 1, 7, 8, 1, 5, 7, -1, -1, -1, -1),
    (11, 2, 1, 11, 1, 7, 7, 1, 5, -1, -1, -1, -1, -1, -1, -1),
    (9, 5, 8, 8, 5, 7, 10, 1, 3, 10, 3, 11, -1, -1, -1, -1),
    (5, 7, 0, 5, 0, 9, 7, 11, 0, 1, 0, 10, 11, 10, 0, -1),
    (11, 10, 0, 11, 0, 3, 10, 5, 0, 8, 0, 7, 5, 7, 0, -1),
    (11, 10, 5, 7, 11, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
    (10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
    (0, 8, 3, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
    (9, 0, 1, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
    (1, 8, 3, 1, 9, 8, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1),
    (1, 6, 5, 2, 6, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
    (1, 6, 5, 1, 2, 6, 3, 0, 8, -1, -1, -1, -1, -1, -1, -1),
    (9, 6, 5, 9, 0, 6, 0, 2, 6, -1, -1, -1, -1, -1, -1, -1),
    (5, 9, 8, 5, 8, 2, 5, 2, 6, 3, 2, 8, -1, -1, -1, -1),
    (2, 3, 11, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
    (11, 0, 8, 11, 2, 0, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1),
    (0, 1, 9, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1),
    (5, 10, 6, 1, 9, 2, 9, 11, 2, 9, 8, 11, -1, -1, -1, -1),
    (6, 3, 11, 6, 5, 3, 5, 1, 3, -1, -1, -1, -1, -1, -1, -1),
    (0, 8, 11, 0, 11, 5, 0, 5, 1, 5, 11, 6, -1, -1, -1, -1),
    (3, 11, 6, 0, 3, 6, 0, 6, 5, 0, 5, 9, -1, -1, -1, -1),
    (6, 5, 9, 6, 9, 11, 11, 9, 8, -1, -1, -1, -1, -1, -1, -1),
    (5, 10, 6, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
    (4, 3, 0, 4, 7, 3, 6, 5, 10, -1, -1, -1, -1, -1, -1, -1),
    (1, 9, 0, 5, 10, 6, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1),
    (10, 6, 5, 1, 9, 7, 1, 7, 3, 7, 9, 4, -1, -1, -1, -1),
    (6, 1, 2, 6, 5, 1, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1),
    (1, 2, 5, 5, 2, 6, 3, 0, 4, 3, 4, 7, -1, -1, -1, -1),
    (8, 4, 7, 9, 0, 5, 0, 6, 5, 0, 2, 6, -1, -1, -1, -1),
    (7, 3, 9, 7, 9, 4, 3, 2, 9, 5, 9, 6, 2, 6, 9, -1),
    (3, 11, 2, 7, 8, 4, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1),
    (5, 10, 6, 4, 7, 2, 4, 2, 0, 2, 7, 11, -1, -1, -1, -1),
    (0, 1, 9, 4, 7, 8, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1),
    (9, 2, 1, 9, 11, 2, 9, 4, 11, 7, 11, 4, 5, 10, 6, -1),
    (8, 4, 7, 3, 11, 5, 3, 5, 1, 5, 11, 6, -1, -1, -1, -1),
    (5, 1, 11, 5, 11, 6, 1, 0, 11, 7, 11, 4, 0, 4, 11, -1),
    (0, 5, 9, 0, 6, 5, 0, 3, 6, 11, 6, 3, 8, 4, 7, -1),
    (6, 5, 9, 6, 9, 11, 4, 7, 9, 7, 11, 9, -1, -1, -1, -1),
    (10, 4, 9, 6, 4, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
    (4, 10, 6, 4, 9, 10, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1),
    (10, 0, 1, 10, 6, 0, 6, 4, 0, -1, -1, -1, -1, -1, -1, -1),
    (8, 3, 1, 8, 1, 6, 8, 6, 4, 6, 1, 10, -1, -1, -1, -1),
    (1, 4, 9, 1, 2, 4, 2, 6, 4, -1, -1, -1, -1, -1, -1, -1),
    (3, 0, 8, 1, 2, 9, 2, 4, 9, 2, 6, 4, -1, -1, -1, -1),
    (0, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
    (8, 3, 2, 8, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1),
    (10, 4, 9, 10, 6, 4, 11, 2, 3, -1, -1, -1, -1, -1, -1, -1),
    (0, 8, 2, 2, 8, 11, 4, 9, 10, 4, 10, 6, -1, -1, -1, -1),
    (3, 11, 2, 0, 1, 6, 0, 6, 4, 6, 1, 10, -1, -1, -1, -1),
    (6, 4, 1, 6, 1, 10, 4, 8, 1, 2, 1, 11, 8, 11, 1, -1),
    (9, 6, 4, 9, 3, 6, 9, 1, 3, 11, 6, 3, -1, -1, -1, -1),
    (8, 11, 1, 8, 1, 0, 11, 6, 1, 9, 1, 4, 6, 4, 1, -1),
    (3, 11, 6, 3, 6, 0, 0, 6, 4, -1, -1, -1, -1, -1, -1, -1),
    (6, 4, 8, 11, 6, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
    (7, 10, 6, 7, 8, 10, 8, 9, 10, -1, -1, -1, -1, -1, -1, -1),
    (0, 7, 3, 0, 10, 7, 0, 9, 10, 6, 7, 10, -1, -1, -1, -1),
    (10, 6, 7, 1, 10, 7, 1, 7, 8, 1, 8, 0, -1, -1, -1, -1),
    (10, 6, 7, 10, 7, 1, 1, 7, 3, -1, -1, -1, -1, -1, -1, -1),
    (1, 2, 6, 1, 6, 8, 1, 8, 9, 8, 6, 7, -1, -1, -1, -1),
    (2, 6, 9, 2, 9, 1, 6, 7, 9, 0, 9, 3, 7, 3, 9, -1),
    (7, 8, 0, 7, 0, 6, 6, 0, 2, -1, -1, -1, -1, -1, -1, -1),
    (7, 3, 2, 6, 7, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
    (2, 3, 11, 10, 6, 8, 10, 8, 9, 8, 6, 7, -1, -1, -1, -1),
    (2, 0, 7, 2, 7, 11, 0, 9, 7, 6, 7, 10, 9, 10, 7, -1),
    (1, 8, 0, 1, 7, 8, 1, 10, 7, 6, 7, 10, 2, 3, 11, -1),
    (11, 2, 1, 11, 1, 7, 10, 6, 1, 6, 7, 1, -1, -1, -1, -1),
    (8, 9, 6, 8, 6, 7, 9, 1, 6, 11, 6, 3, 1, 3, 6, -1),
    (0, 9, 1, 11, 6, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
    (7, 8, 0, 7, 0, 6, 3, 11, 0, 11, 6, 0, -1, -1, -1, -1),
    (7, 11, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
    (7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
    (3, 0, 8, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
    (0, 1, 9, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
    (8, 1, 9, 8, 3, 1, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1),
    (10, 1, 2, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
    (1, 2, 10, 3, 0, 8, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1),
    (2, 9, 0, 2, 10, 9, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1),
    (6, 11, 7, 2, 10, 3, 10, 8, 3, 10, 9, 8, -1, -1, -1, -1),
    (7, 2, 3, 6, 2, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
    (7, 0, 8, 7, 6, 0, 6, 2, 0, -1, -1, -1, -1, -1, -1, -1),
    (2, 7, 6, 2, 3, 7, 0, 1, 9, -1, -1, -1, -1, -1, -1, -1),
    (1, 6, 2, 1, 8, 6, 1, 9, 8, 8, 7, 6, -1, -1, -1, -1),
    (10, 7, 6, 10, 1, 7, 1, 3, 7, -1, -1, -1, -1, -1, -1, -1),
    (10, 7, 6, 1, 7, 10, 1, 8, 7, 1, 0, 8, -1, -1, -1, -1),
    (0, 3, 7, 0, 7, 10, 0, 10, 9, 6, 10, 7, -1, -1, -1, -1),
    (7, 6, 10, 7, 10, 8, 8, 10, 9, -1, -1, -1, -1, -1, -1, -1),
    (6, 8, 4, 11, 8, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
    (3, 6, 11, 3, 0, 6, 0, 4, 6, -1, -1, -1, -1, -1, -1, -1),
    (8, 6, 11, 8, 4, 6, 9, 0, 1, -1, -1, -1, -1, -1, -1, -1),
    (9, 4, 6, 9, 6, 3, 9, 3, 1, 11, 3, 6, -1, -1, -1, -1),
    (6, 8, 4, 6, 11, 8, 2, 10, 1, -1, -1, -1, -1, -1, -1, -1),
    (1, 2, 10, 3, 0, 11, 0, 6, 11, 0, 4, 6, -1, -1, -1, -1),
    (4, 11, 8, 4, 6, 11, 0, 2, 9, 2, 10, 9, -1, -1, -1, -1),
    (10, 9, 3, 10, 3, 2, 9, 4, 3, 11, 3, 6, 4, 6, 3, -1),
    (8, 2, 3, 8, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1),
    (0, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
    (1, 9, 0, 2, 3, 4, 2, 4, 6, 4, 3, 8, -1, -1, -1, -1),
    (1, 9, 4, 1, 4, 2, 2, 4, 6, -1, -1, -1, -1, -1, -1, -1),
    (8, 1, 3, 8, 6, 1, 8, 4, 6, 6, 10, 1, -1, -1, -1, -1),
    (10, 1, 0, 10, 0, 6, 6, 0, 4, -1, -1, -1, -1, -1, -1, -1),
    (4, 6, 3, 4, 3, 8, 6, 10, 3, 0, 3, 9, 10, 9, 3, -1),
    (10, 9, 4, 6, 10, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
    (4, 9, 5, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
    (0, 8, 3, 4, 9, 5, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1),
    (5, 0, 1, 5, 4, 0, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1),
    (11, 7, 6, 8, 3, 4, 3, 5, 4, 3, 1, 5, -1, -1, -1, -1),
    (9, 5, 4, 10, 1, 2, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1),
    (6, 11, 7, 1, 2, 10, 0, 8, 3, 4, 9, 5, -1, -1, -1, -1),
    (7, 6, 11, 5, 4, 10, 4, 2, 10, 4, 0, 2, -1, -1, -1, -1),
    (3, 4, 8, 3, 5, 4, 3, 2, 5, 10, 5, 2, 11, 7, 6, -1),
    (7, 2, 3, 7, 6, 2, 5, 4, 9, -1, -1, -1, -1, -1, -1, -1),
    (9, 5, 4, 0, 8, 6, 0, 6, 2, 6, 8, 7, -1, -1, -1, -1),
    (3, 6, 2, 3, 7, 6, 1, 5, 0, 5, 4, 0, -1, -1, -1, -1),
    (6, 2, 8, 6, 8, 7, 2, 1, 8, 4, 8, 5, 1, 5, 8, -1),
    (9, 5, 4, 10, 1, 6, 1, 7, 6, 1, 3, 7, -1, -1, -1, -1),
    (1, 6, 10, 1, 7, 6, 1, 0, 7, 8, 7, 0, 9, 5, 4, -1),
    (4, 0, 10, 4, 10, 5, 0, 3, 10, 6, 10, 7, 3, 7, 10, -1),
    (7, 6, 10, 7, 10, 8, 5, 4, 10, 4, 8, 10, -1, -1, -1, -1),
    (6, 9, 5, 6, 11, 9, 11, 8, 9, -1, -1, -1, -1, -1, -1, -1),
    (3, 6, 11, 0, 6, 3, 0, 5, 6, 0, 9, 5, -1, -1, -1, -1),
    (0, 11, 8, 0, 5, 11, 0, 1, 5, 5, 6, 11, -1, -1, -1, -1),
    (6, 11, 3, 6, 3, 5, 5, 3, 1, -1, -1, -1, -1, -1, -1, -1),
    (1, 2, 10, 9, 5, 11, 9, 11, 8, 11, 5, 6, -1, -1, -1, -1),
    (0, 11, 3, 0, 6, 11, 0, 9, 6, 5, 6, 9, 1, 2, 10, -1),
    (11, 8, 5, 11, 5, 6, 8, 0, 5, 10, 5, 2, 0, 2, 5, -1),
    (6, 11, 3, 6, 3, 5, 2, 10, 3, 10, 5, 3, -1, -1, -1, -1),
    (5, 8, 9, 5, 2, 8, 5, 6, 2, 3, 8, 2, -1, -1, -1, -1),
    (9, 5, 6, 9, 6, 0, 0, 6, 2, -1, -1, -1, -1, -1, -1, -1),
    (1, 5, 8, 1, 8, 0, 5, 6, 8, 3, 8, 2, 6, 2, 8, -1),
    (1, 5, 6, 2, 1, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
    (1, 3, 6, 1, 6, 10, 3, 8, 6, 5, 6, 9, 8, 9, 6, -1),
    (10, 1, 0, 10, 0, 6, 9, 5, 0, 5, 6, 0, -1, -1, -1, -1),
    (0, 3, 8, 5, 6, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
    (10, 5, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
    (11, 5, 10, 7, 5, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
    (11, 5, 10, 11, 7, 5, 8, 3, 0, -1, -1, -1, -1, -1, -1, -1),
    (5, 11, 7, 5, 10, 11, 1, 9, 0, -1, -1, -1, -1, -1, -1, -1),
    (10, 7, 5, 10, 11, 7, 9, 8, 1, 8, 3, 1, -1, -1, -1, -1),
    (11, 1, 2, 11, 7, 1, 7, 5, 1, -1, -1, -1, -1, -1, -1, -1),
    (0, 8, 3, 1, 2, 7, 1, 7, 5, 7, 2, 11, -1, -1, -1, -1),
    (9, 7, 5, 9, 2, 7, 9, 0, 2, 2, 11, 7, -1, -1, -1, -1),
    (7, 5, 2, 7, 2, 11, 5, 9, 2, 3, 2, 8, 9, 8, 2, -1),
    (2, 5, 10, 2, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1),
    (8, 2, 0, 8, 5, 2, 8, 7, 5, 10, 2, 5, -1, -1, -1, -1),
    (9, 0, 1, 5, 10, 3, 5, 3, 7, 3, 10, 2, -1, -1, -1, -1),
    (9, 8, 2, 9, 2, 1, 8, 7, 2, 10, 2, 5, 7, 5, 2, -1),
    (1, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
    (0, 8, 7, 0, 7, 1, 1, 7, 5, -1, -1, -1, -1, -1, -1, -1),
    (9, 0, 3, 9, 3, 5, 5, 3, 7, -1, -1, -1, -1, -1, -1, -1),
    (9, 8, 7, 5, 9, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
    (5, 8, 4, 5, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1),
    (5, 0, 4, 5, 11, 0, 5, 10, 11, 11, 3, 0, -1, -1, -1, -1),
    (0, 1, 9, 8, 4, 10, 8, 10, 11, 10, 4, 5, -1, -1, -1, -1),
    (10, 11, 4, 10, 4, 5, 11, 3, 4, 9, 4, 1, 3, 1, 4, -1),
    (2, 5, 1, 2, 8, 5, 2, 11, 8, 4, 5, 8, -1, -1, -1, -1),
    (0, 4, 11, 0, 11, 3, 4, 5, 11, 2, 11, 1, 5, 1, 11, -1),
    (0, 2, 5, 0, 5, 9, 2, 11, 5, 4, 5, 8, 11, 8, 5, -1),
    (9, 4, 5, 2, 11, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
    (2, 5, 10, 3, 5, 2, 3, 4, 5, 3, 8, 4, -1, -1, -1, -1),
    (5, 10, 2, 5, 2, 4, 4, 2, 0, -1, -1, -1, -1, -1, -1, -1),
    (3, 10, 2, 3, 5, 10, 3, 8, 5, 4, 5, 8, 0, 1, 9, -1),
    (5, 10, 2, 5, 2, 4, 1, 9, 2, 9, 4, 2, -1, -1, -1, -1),
    (8, 4, 5, 8, 5, 3, 3, 5, 1, -1, -1, -1, -1, -1, -1, -1),
    (0, 4, 5, 1, 0, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
    (8, 4, 5, 8, 5, 3, 9, 0, 5, 0, 3, 5, -1, -1, -1, -1),
    (9, 4, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
    (4, 11, 7, 4, 9, 11, 9, 10, 11, -1, -1, -1, -1, -1, -1, -1),
    (0, 8, 3, 4, 9, 7, 9, 11, 7, 9, 10, 11, -1, -1, -1, -1),
    (1, 10, 11, 1, 11, 4, 1, 4, 0, 7, 4, 11, -1, -1, -1, -1),
    (3, 1, 4, 3, 4, 8, 1, 10, 4, 7, 4, 11, 10, 11, 4, -1),
    (4, 11, 7, 9, 11, 4, 9, 2, 11, 9, 1, 2, -1, -1, -1, -1),
    (9, 7, 4, 9, 11, 7, 9, 1, 11, 2, 11, 1, 0, 8, 3, -1),
    (11, 7, 4, 11, 4, 2, 2, 4, 0, -1, -1, -1, -1, -1, -1, -1),
    (11, 7, 4, 11, 4, 2, 8, 3, 4, 3, 2, 4, -1, -1, -1, -1),
    (2, 9, 10, 2, 7, 9, 2, 3, 7, 7, 4, 9, -1, -1, -1, -1),
    (9, 10, 7, 9, 7, 4, 10, 2, 7, 8, 7, 0, 2, 0, 7, -1),
    (3, 7, 10, 3, 10, 2, 7, 4, 10, 1, 10, 0, 4, 0, 10, -1),
    (1, 10, 2, 8, 7, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
    (4, 9, 1, 4, 1, 7, 7, 1, 3, -1, -1, -1, -1, -1, -1, -1),
    (4, 9, 1, 4, 1, 7, 0, 8, 1, 8, 7, 1, -1, -1, -1, -1),
    (4, 0, 3, 7, 4, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
    (4, 8, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
    (9, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
    (3, 0, 9, 3, 9, 11, 11, 9, 10, -1, -1, -1, -1, -1, -1, -1),
    (0, 1, 10, 0, 10, 8, 8, 10, 11, -1, -1, -1, -1, -1, -1, -1),
    (3, 1, 10, 11, 3, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
    (1, 2, 11, 1, 11, 9, 9, 11, 8, -1, -1, -1, -1, -1, -1, -1),
    (3, 0, 9, 3, 9, 11, 1, 2, 9, 2, 11, 9, -1, -1, -1, -1),
    (0, 2, 11, 8, 0, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
    (3, 2, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
    (2, 3, 8, 2, 8, 10, 10, 8, 9, -1, -1, -1, -1, -1, -1, -1),
    (9, 10, 2, 0, 9, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
    (2, 3, 8, 2, 8, 10, 0, 1, 8, 1, 10, 8, -1, -1, -1, -1),
    (1, 10, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
    (1, 3, 8, 9, 1, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
    (0, 9, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
    (0, 3, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1),
    (-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1)
  ]
  
  return edgeTable, triTable  

def write_vtp(points,cells,hx,hy,hz):    
  from draw_profile_functions import grid_to_cartesian_coords
  from matviz_vtk import show_write_vtk
  tmp_points = np.copy(points)
  for i in range(len(tmp_points)):
    tmp_points[i] = grid_to_cartesian_coords(tmp_points[i],None,(hx,hy,hz))
    
  pd = matviz_vtk.fill_vtk_polydata(points,cells)
  show_write_vtk(pd, 10, "marching_cubes.vtp")
  
def create_triangles2(array,spacing,step):
  
  from skimage import measure
  verts, faces, normals, values = measure.marching_cubes(array,spacing=spacing,step_size=step,allow_degenerate=False)
  
  return verts, faces, normals