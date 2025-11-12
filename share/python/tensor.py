# this file provides tools to compute elasticity tensors
# the formulas are based on 
# Dissertation Mike Richter, Entwicklung mechanischer Modelle zur analytischen Beschreibung der Materialeigenschaften von textilbewehrtem Feinbeton, 2005

import numpy as np

# apply symmetry from upper right triangle to lower left triangle in place! 2D and 3D
def apply_symmetry(mat):
  N = len(mat)
  for i in range(N):
    for j in range(i+1,N):
      mat[j,i] = mat[i,j]

# calculate the transversal isotropic compliance tensor based the inverse of the stiffness tensor
# https://en.wikipedia.org/wiki/Transverse_isotropy
# axial direction is z-direction!
# E_3 usually the stiffest direction  
# G_3: G_xz = G_yz
def transversely_isotropic(E_iso, E_3, nu_iso, nu_3, G_3):
   
    # Compliance matrix (S)
    S = np.zeros((6, 6))
    S[0, 0] = 1 / E_iso
    S[1, 1] = 1 / E_iso
    S[2, 2] = 1 / E_3
    S[0, 1] = -nu_iso / E_iso
    S[1, 0] = S[0, 1]
    S[0, 2] = -nu_3 / E_3
    S[2, 0] = S[0, 2]
    S[1, 2] = -nu_3 / E_3
    S[2, 1] = S[1, 2]
    S[3, 3] = 1 / G_3
    S[4, 4] = 1 / G_3
    S[5, 5] = (2 * (1 + nu_iso) / E_iso)

    # Invert to get stiffness matrix
    # print('transversely_isotropic S\n',S)
    C = np.linalg.inv(S)
    # print('transversely_isotropic inv(S)=C:\n',C)
    return C


# based on (2.50), the inverse tensor (2.52) would be easier to implement
def orthotrope(E1, E2, E3, nu_12, nu_21, nu_13, nu_31, nu_23, nu_32, G12, G23, G31):
  commn = (1-nu_12*nu_21-nu_23*nu_32-nu_31*nu_13-2*nu_21*nu_32*nu_13)/(E1*E2*E3)  
  ret = np.zeros((6,6))
  ret[0,0] = (1-nu_23*nu_32)/(E2*E3*commn)
  ret[0,1] = (nu_21+nu_31*nu_23)/(E2*E3*commn)
  ret[0,2] = (nu_31+nu_21*nu_32)/(E2*E3*commn)
  ret[1,1] = (1-nu_13*nu_31)/(E1*E3*commn)
  ret[1,2] = (nu_32+nu_12*nu_31)/(E1*E3*commn)
  ret[2,2] = (1-nu_12*nu_21)/(E1*E2*commn)
  ret[3,3] = G23
  ret[4,4] = G31
  ret[5,5] = G12
  #print(ret)
  apply_symmetry(ret)
  #print('orthotrope',ret)
  return ret

# cubic orthotrope ("iso-orthotrope" in cfs)
def cubic_orthotrope(E, nu, G):
  mat = orthotrope(E,E,E,nu,nu,nu,nu,nu,nu,G,G,G)
  #print('cubic orthotrope',mat)
  return mat

# simplified orthotrope tensor with E2=E3 as form of quasi-transversal-isotropic
# this is validate for isotropic 99lines material (use cfs -d) x_oriented_orthotrope(1, 1,0.3,0.3,1.0/(2*(1+.3)))
def x_oriented_orthotrope(E1,E2,nu_12,nu_23,G12):
  mat = orthotrope(E1, E2, E2, nu_12, nu_12, nu_23, nu_23, nu_23, nu_23, G12, G12, G12)
  return mat

# simplified orthotrope tensor with E1=E3 as form of quasi-transversal-isotropic
def y_oriented_orthotrope(E2,E1,nu_21,nu_13,G2):
  # def orthotrope(E1, E2, E3, nu_12, nu_21, nu_13, nu_31, nu_23, nu_32, G12, G23, G31):
  mat = orthotrope(E1, E2, E1, nu_21, nu_21, nu_13, nu_13, nu_21, nu_21, G2, G2, G2)
  return mat


# simplified orthotrope tensor with E1=E2 as form of quasi-transversal-isotropic
def z_oriented_orthotrope(E3,E1,nu_31,nu_12,G13):
  # def orthotrope(E1, E2, E3, nu_12, nu_21, nu_13, nu_31, nu_23, nu_32, G12, G23, G31):
  mat = orthotrope(E1, E1, E3, nu_12, nu_12, nu_31, nu_31, nu_31, nu_31, G13, G13, G13)
  return mat


# Reorder stiffness matrix C so that the TI axis is along x, y, or z.
# Default is z (no change).
# code is from ChatGPT but seams right
def rotate_transversal_isotropic_axis(C, axis='z'):
    # Voigt index mapping for axes swap (x=1, y=2, z=3)
    # mapping defines permutation of indices for (1,2,3)-> new order
    if axis == 'z':
        return C.copy()
    elif axis == 'x':
        # swap x <-> z => (1,2,3) -> (3,2,1)
        perm = [2,1,0]
    elif axis == 'y':
        # swap y <-> z => (1,2,3) -> (1,3,2)
        perm = [0,2,1]
    else:
        raise ValueError("axis must be 'x', 'y', or 'z'")
    
    # Full 6×6 permutation for Voigt notation
    voigt_pairs = [(0,0),(1,1),(2,2),(1,2),(0,2),(0,1)]
    P = np.zeros((6,6))
    for i, (a,b) in enumerate(voigt_pairs):
        new_a = perm[a]
        new_b = perm[b]
        # find new Voigt index
        if new_a == new_b:
            new_i = new_a
        else:
            new_pair = sorted([new_a,new_b])
            if new_pair == [1,2]:
                new_i = 3
            elif new_pair == [0,2]:
                new_i = 4
            elif new_pair == [0,1]:
                new_i = 5
        P[i, new_i] = 1.0
    
    # Apply transformation C' = P C P^T
    return P @ C @ P.T

# create a string with a cfs mat.xml material entry
# @param density in kg/m^3 is 1e-3 times g/cm^3
def cfs_material_tensor(name, tensor, density):
    
  t = ""
  for ri, r in enumerate(tensor):
    for c in r:
      t += format(c, f"<16.13g")
    if ri < len(tensor)-1:
      t += '\n              '
  coefficients = t

  template = f"""  <material name="{name}">
    <mechanical>
      <density>
        <linear>
          <real>{density}</real>
        </linear>
      </density>
      <elasticity>
        <linear>
          <tensor dim2="6" dim1="6">
            <real>
              {coefficients}
            </real>
          </tensor>
        </linear> 
      </elasticity>
    </mechanical>
  </material>"""

  return template



# transversal_isotropic(6.5, 0.3, 0.22, 0.5, 0.4, 0.22)  

# scaled data from FIONA Abschlussbericht E1=65Gpa, E2=5GPa, G12=2.2GPa, nu_12=0.3, nu_23=0.4, rho=1.55g/cm^3
#mat = y_oriented_orthotrope(6.5, 0.5,0.4,0.03,0.22)
#print(cfs_material_tensor('c-cfrtp_scaled', mat, 0.00155))

# a very mild ansiostropic variant for 99-lines material
#print(cfs_material_tensor('mild_aniso', y_oriented_orthotrope(2, 0.5,0.3,0.3,1.0/(2*(1+.3))),0.00155))

#print(cfs_material_tensor('strong_aniso', x_oriented_orthotrope(2, 0.5,0.3,0.3,0.25),0.00155))
#                                   def orthotrope(E1, E2, E3, nu_12, nu_21, nu_13, nu_31, nu_23, nu_32, G12, G23, G31):

# print(cfs_material_tensor('ortho_test', orthotrope(1, 2, 3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.12, 0.23, 0.31),0.00155))

# x_oriented_orthotrope(E1,E2,nu_12,nu_23,G12)
#print(cfs_material_tensor('cfrtp_scaled_3d',x_oriented_orthotrope(E1 = 6.5, E2=0.5, nu_12=0.03, nu_23=0.3, G12=0.22),0.00155))

# transversely_isotropic(E_iso, E_3, nu_iso, nu_3, G_3):
print(cfs_material_tensor('trans', transversely_isotropic(E_iso=0.5, E_3=6.5, nu_iso=0.4, nu_3=0.3, G_3=0.22),0.00155))

C = transversely_isotropic(E_iso=0.5, E_3=6.5, nu_iso=0.4, nu_3=0.3, G_3=0.22)
Cy = rotate_transversal_isotropic_axis(C, axis='y')
print(cfs_material_tensor('crftp_for_2d', Cy,0.00155))
