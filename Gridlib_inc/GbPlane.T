// -*- C++ -*- automatisch in C++mode wechseln (emacs)
/*----------------------------------------------------------------------
|
|
| $Log$
| Revision 1.1  2002/02/22 14:47:56  elena
| new: dir Gridlib_inc
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
  constant_ = normal_.dot(rkPoint0);
}

template <class T>
GbPlane<T>::Side 
GbPlane<T>::whichSide (const GbVec3<T>& rkPoint) const
{
  T fDistance = distanceTo(rkPoint);

  if ( fDistance < 0.0 )
    return GbPlane::NEGATIVE_SIDE;

  if ( fDistance > 0.0 )
    return GbPlane::POSITIVE_SIDE;

  return GbPlane::NO_SIDE;
}
