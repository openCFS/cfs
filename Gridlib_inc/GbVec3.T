// -*- C++ -*- automatisch in C++mode wechseln (emacs)
/*----------------------------------------------------------------------
|
|
| $Log$
| Revision 1.2  2002/03/21 14:58:57  elena
| new: changes in dat-file for reading tetrahedral (bugs in element connection)
|
| Revision 1.5  2001/09/12 09:17:40  prkipfer
| added comparision operators to be able to use it with STL containers
|
| Revision 1.4  2001/08/20 12:59:52  prkipfer
| introduced homogeneous 3D vector to avoid confusion with 4D vectors
|
| Revision 1.3  2001/08/16 16:53:22  prkipfer
| improved type safety for template parameter
|
| Revision 1.2  2001/06/20 16:13:47  prkipfer
| moved statics to .T file
|
| Revision 1.1  2000/06/14 15:39:13  prkipfer
| improved base classes and added funcstruct processing for mesh
|
|
|
+---------------------------------------------------------------------*/

template <class T>
const GbVec3<T>
GbVec3<T>::ZERO(T(0.0),T(0.0),T(0.0));

template <class T>
const GbVec3<T>
GbVec3<T>::UNIT_X(T(1.0),T(0.0),T(0.0));

template <class T>
const GbVec3<T>
GbVec3<T>::UNIT_Y(T(0.0),T(1.0),T(0.0));

template <class T>
const GbVec3<T>
GbVec3<T>::UNIT_Z(T(0.0),T(0.0),T(1.0));

template <class T>
T
GbVec3<T>::EPSILON(T(0.0));

template <class T>
GbVec3<T>::GbVec3() 
{
  // no initialization, use constructor below in case
}

template <class T>
GbVec3<T>::GbVec3(T s)
{ 
  xyz[0]=xyz[1]=xyz[2]=s; 
}

template <class T>
GbVec3<T>::GbVec3(T x, T y, T z)
{ 
  xyz[0]=x; xyz[1]=y; xyz[2]=z; 
}

template <class T>
GbVec3<T>::GbVec3(const T p[3])
{ 
  xyz[0]=p[0]; xyz[1]=p[1]; xyz[2]=p[2]; 
}
  
template <class T>
GbVec3<T>::GbVec3(const GbVec3<T> &p)
{ 
  xyz[0]=p.xyz[0]; xyz[1]=p.xyz[1]; xyz[2]=p.xyz[2]; 
}

template <class T>
GbVec3<T>::~GbVec3()
{
}
