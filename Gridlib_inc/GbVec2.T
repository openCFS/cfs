// -*- C++ -*- automatisch in C++mode wechseln (emacs)
/*----------------------------------------------------------------------
|
|
| $Log$
| Revision 1.1  2002/02/22 14:47:56  elena
| new: dir Gridlib_inc
|
| Revision 1.2  2001/06/26 12:27:14  prkipfer
| added missing default constructor
|
| Revision 1.1  2001/06/20 16:14:38  prkipfer
| introduced 2D Vector class
|
|
|
+---------------------------------------------------------------------*/

template <class T>
const GbVec2<T>
GbVec2<T>::ZERO(0,0);

template <class T>
const GbVec2<T> 
GbVec2<T>::UNIT_X(1,0);

template <class T>
const GbVec2<T> 
GbVec2<T>::UNIT_Y(0,1);

template <class T>
GbVec2<T>::GbVec2()
{
  xy[0] = 0;
  xy[1] = 0;
}

template <class T>
GbVec2<T>::GbVec2(T fX, T fY)
{
  xy[0] = fX;
  xy[1] = fY;
}

template <class T>
GbVec2<T>::GbVec2(const T afCoordinate[2])
{
  xy[0] = afCoordinate[0];
  xy[1] = afCoordinate[1];
}

template <class T>
GbVec2<T>::GbVec2(const GbVec2<T>& rkVector)
{
  xy[0] = rkVector[0];
  xy[1] = rkVector[1];
}

template <class T>
GbVec2<T>::~GbVec2()
{
}
