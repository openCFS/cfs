31.05.2013

Readme file for show_fmo_tensor.py

Authors: 
main author: Fabian Wein, fabian.wein@am.uni-erlangen.de
based on work of Bastian Schmidt, bastian.schmidt@am.uni-erlangen.de
contributions: Thomas Guess, thomas.guess@fau.de

show_fmo_tensor can either read many tensors from an hdf5 file in CFS++ format
or it can visualize individual tensors provided on the command line

easiest case for a 2D elasticity tensor, with the coefficients
e11, e22, e33, e32, e31, e21

show_fmo_tensor.py "[1,0,0,0,0,0]"

easiest case for a 3D elasticity tensor, with the coefficients
e11, e12, e22, e13, e23, e33, e14, e24, ...

show_fmo_tensor.py "[1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]"

the most important options:

show_fmo_tensor.py --help

show_fmo_tensor.py --notation voigt

show_fmo_tensor.py --save file.png or file.jpg for 2D and file.vtp for 3D to open in Paraview

show_fmo_tensor.py --res 1600 to change image size

show_fmo_tensor.py --sampling 360 for the discretization