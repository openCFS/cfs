// -*- C++ -*- automatisch in C++mode wechseln (emacs)
/*----------------------------------------------------------------------
|
|
| $Log$
| Revision 1.3  2002/04/03 11:17:07  elena
| new: changes in CFS++ on SGI and new function for gridlib
|
| Revision 1.4  2001/09/12 09:17:39  prkipfer
| added comparision operators to be able to use it with STL containers
|
| Revision 1.3  2001/08/16 16:53:21  prkipfer
| improved type safety for template parameter
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
GbVec2<T>::ZERO(T(0.0),T(0.0));

template <class T>
const GbVec2<T> 
GbVec2<T>::UNIT_X(T(1.0),T(0.0));

template <class T>
const GbVec2<T> 
GbVec2<T>::UNIT_Y(T(0.0),T(1.0));

template <class T>
T
GbVec2<T>::EPSILON(T(0.0));

template <class T>
GbVec2<T>::GbVec2()
{
  // no initialization, use the following constructor if needed
}

template <class T>
GbVec2<T>::GbVec2(T f)
{
  xy[0] = xy[1] = f;
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
