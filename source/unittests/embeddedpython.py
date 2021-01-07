from numpy import dtype
a = 0

print('embeddedpython.py: Hello World')

def inc_a():
  global a
  a += 1
  return a

# the cfs module is generated on runtime via cfs during embedding Python
import cfs
print("embeddedpython.py: callback cfs.val(4)", cfs.val(4))
assert cfs.val(4) == 5

import numpy

V = numpy.ones((10), dtype=float) # makes c-double
# assert cfs.vec(V, 10) == 11
r = cfs.vec(V,V,len(V))
print('vec -> ', r)
assert r == len(V)



A = numpy.array([[1, 2, 3], [4, 5, 6]], dtype = float) # type is important here, otherwise it is int

def print_A():
  print('embeddedpython.py: print_A() -> ', A)

# val1 needs to be larger val2
def many_values(val1, val2):
  assert val1 > val2
  print('embeddedpython.py: ', val1, 'indededd greater', val2)
  return numpy.ones((3)), numpy.zeros((3))
  
def print_dict(data):
  print('embeddedpython.py: dict=', data)
  
def get_vec(vec, emptyvec, n):
  print('embeddedpython.py get_vec ->', vec,emptyvec, n)  
  assert len(vec) == n
   
  