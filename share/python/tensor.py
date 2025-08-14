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

# based on (2.55), note that this might be simplified: nu_31=nu_32, G_31=G_23=G_32=G_13=G_3, G_12=G_21=G, E_11=E_22=E
# the axial direction is z-direction!
def transversal_isotropic(E, nu, G, E_3, nu_3, G_3):
  #print(G)
  #print(E/(2*(1+nu)))
  #assert np.isclose(G, E/(2*(1+nu))) # (2.54b)
  denom = (1+nu)/E * (1-nu-2*nu_3*nu_3*E/E_3) # (2.57)
  ret = np.zeros((6,6))
  ret[0,0] = (1-nu_3*nu_3*E/E_3)/denom
  ret[0,1] = (nu-nu_3*nu_3*E/E_3)/denom
  ret[0,2] = (1+nu)*nu_3/denom
  ret[1,1] = ret[0,0]
  ret[1,2] = ret[0,2]
  ret[2,2] = ((1-nu*nu)*E/E_3)/denom
  ret[3,3] = G_3
  ret[4,4] = ret[3,3]
  ret[5,5] = G

  print(ret)
  apply_symmetry(ret)
  print('symmetric',ret)
  return ret

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

# simplified orthotrope tensor with E1=E2 as form of quasi-transversal-isotropic
def y_oriented_orthotrope(E2,E1,nu_21,nu_13,G2):
  # def orthotrope(E1, E2, E3, nu_12, nu_21, nu_13, nu_31, nu_23, nu_32, G12, G23, G31):
  mat = orthotrope(E1, E2, E1, nu_21, nu_21, nu_13, nu_13, nu_21, nu_21, G2, G2, G2)
  return mat


# simplified orthotrope tensor with E1=E2 as form of quasi-transversal-isotropic
def z_oriented_orthotrope(E3,E1,nu_31,nu_12,G13):
  # def orthotrope(E1, E2, E3, nu_12, nu_21, nu_13, nu_31, nu_23, nu_32, G12, G23, G31):
  mat = orthotrope(E1, E1, E3, nu_12, nu_12, nu_31, nu_31, nu_31, nu_31, G13, G13, G13)
  return mat


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
mat = y_oriented_orthotrope(6.5, 0.5,0.3,0.4,0.22)
#print(cfs_material_tensor('c-cfrtp_scaled', mat, 0.00155))

# a very mild ansiostropic variant for 99-lines material
print(cfs_material_tensor('mild_aniso', y_oriented_orthotrope(2, 0.5,0.3,0.3,1.0/(2*(1+.3))),0.00155))

#print(cfs_material_tensor('strong_aniso', x_oriented_orthotrope(2, 0.5,0.3,0.3,0.25),0.00155))
#                                   def orthotrope(E1, E2, E3, nu_12, nu_21, nu_13, nu_31, nu_23, nu_32, G12, G23, G31):

# print(cfs_material_tensor('ortho_test', orthotrope(1, 2, 3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.12, 0.23, 0.31),0.00155))