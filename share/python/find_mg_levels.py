import numpy as np 
from optimization_tools import *
import math 

def find_levels(mesh_size,num_cores,nlvls):
  # This method will help to figure out for the given combination of the mesh 
  # size and processor the max level that is feasible
  [nx,ny,nz] = mesh_size
  divisor = math.pow(2.0,nlvls-1.0);
  
    nlvls=nlvls-1;
    find_levels(mesh_size, num_cores,nlvls);
  
  else:
    return nlvls;
  
  
  



nlvls=[]
for mesh in range(100):
  
  nlvls.append(find_levels([mesh+10,mesh+10,mesh+10], 1,10))


print(nlvls)