from h5py import File
from numpy import sum, pi

h5 = File('results_hdf5/dbg.cfs')

from sys import path
path.append('/home/flo/cfs/code/share/python')
from hdf5_tools import get_result

E_volumes = get_result(h5,'volume',region='V_8',multistep=0)
S_8_area = get_result(h5,'area',region='S_8',multistep=0)
S_xy_area = get_result(h5,'area',region='S_xy',multistep=0)


# radius of the sphere
r = 1.0
print('====== Analytic ======')
print('Radius = %f'%(r))
V = 1/8*4*pi/3*r**3
print('Vomlume of 1/8 sphere = %f'%V)
A_S8 = 1/8*4*pi*r**2
A_Sxy = 1/4*pi*r**2
print('Surface Area of 1/8 sphere = %f'%A_S8)
print('Surface Area of 1/4 circle = %f'%A_Sxy)
print('========= FE =========')
V_fe = sum(E_volumes)
print('Volume = %f (error = %.4e)'%(V_fe,(V_fe-V)/V))
A_S8_FE = sum( S_8_area )
print('Sphere Surface = %f (error = %.4e)'%(A_S8_FE,(A_S8_FE-A_S8)/A_S8))
A_Sxy_FE = sum( S_xy_area )
print('Circle Surface = %f (error = %.4e)'%(A_Sxy_FE,(A_Sxy_FE-A_Sxy)/A_Sxy))
