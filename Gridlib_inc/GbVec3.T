// -*- C++ -*- automatisch in C++mode wechseln (emacs)
/*----------------------------------------------------------------------
|
|
| $Log$
| Revision 1.1  2002/02/22 14:47:56  elena
| new: dir Gridlib_inc
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
GbVec3<T>::ZERO(0,0,0);

template <class T>
const GbVec3<T>
GbVec3<T>::UNIT_X(1,0,0);

template <class T>
const GbVec3<T>
GbVec3<T>::UNIT_Y(0,1,0);

template <class T>
const GbVec3<T>
GbVec3<T>::UNIT_Z(0,0,1);

template <class T>
GbVec3<T>::GbVec3() 
{
  xyz[0]=0; xyz[1]=0; xyz[2]=0;
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
