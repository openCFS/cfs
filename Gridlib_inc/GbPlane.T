// -*- C++ -*- automatisch in C++mode wechseln (emacs)
/*----------------------------------------------------------------------
|
|
| $Log$
| Revision 1.2  2002/03/21 14:58:56  elena
| new: changes in dat-file for reading tetrahedral (bugs in element connection)
|
| Revision 1.2  2001/08/16 16:53:20  prkipfer
| improved type safety for template parameter
|
| Revision 1.1  2001/01/02 14:55:04  prkipfer
| introduced new classes
|
|
|
+---------------------------------------------------------------------*/

template <class T>
GbPlane<T>::GbPlane (const GbVec3<T>& rkPoint0, 
		     const GbVec3<T>& rkPoint1,
		     const GbVec3<T>& rkPoint2)
{
  GbVec3<T> kEdge1 = rkPoint1 - rkPoint0;
  GbVec3<T> kEdge2 = rkPoint2 - rkPoint0;
  normal_ = kEdge1.unitCross(kEdge2);
  constant_ = normal_ | rkPoint0;
}

template <class T>
GbPlane<T>::Side 
GbPlane<T>::whichSide (const GbVec3<T>& rkPoint) const
{
  T fDistance = distanceTo(rkPoint);

  if ( fDistance < T(0.0) )
    return GbPlane::NEGATIVE_SIDE;

  if ( fDistance > T(0.0) )
    return GbPlane::POSITIVE_SIDE;

  return GbPlane::NO_SIDE;
}
