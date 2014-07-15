The following section explains how to install show_fmo_tensor.py in Windows,
provided by Thomas Guess

Required Packages:

-python 2.7: http://www.python.org/download/releases/2.7.5/

-python-numpy: http://sourceforge.net/projects/numpy/files/NumPy/1.7.1/

-python-h5py: http://code.google.com/p/h5py/downloads/list

-python-matplotlib: http://matplotlib.org/downloads.html

-python-imaging: http://www.pythonware.com/products/pil/

-python-vtk: http://www.lfd.uci.edu/~gohlke/pythonlibs/

-python-scipy: http://sourceforge.net/projects/scipy/files/scipy/0.12.0/


Additional Packages: maybe needed in addition

-python-matplotlib-tk

-ipython


Bugfix in order to see the Image created by show_fmo_tensor.py:

Replace the following line

return "start /wait %s && del /f %s" % (file, file)

in C:\Python27\lib\site-packages\PIL\ImageShow.py by the line

return "start /wait %s && PING 127.0.0.1 -n 5 > NUL && del /f %s" % (file, file)

Visualize 3D tensor:
python show_fmo_tensor.py "[1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]"

Visualize 2D tensor:
python show_fmo_tensor.py "[1,0,0,0,0,0]"
