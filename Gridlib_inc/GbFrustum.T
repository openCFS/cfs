// -*- C++ -*- automatisch in C++mode wechseln (emacs)
/*----------------------------------------------------------------------
|
|
| $Log$
| Revision 1.3  2002/04/03 11:17:07  elena
| new: changes in CFS++ on SGI and new function for gridlib
|
| Revision 1.2  2001/08/16 16:53:15  prkipfer
| improved type safety for template parameter
|
| Revision 1.1  2001/06/18 10:43:05  prkipfer
| introduced Frustum object
|
|
|
+---------------------------------------------------------------------*/

template <class T>
GbFrustum<T>::GbFrustum ()
  : origin_(GbVec3<T>::ZERO)
  , lVector_(GbVec3<T>::UNIT_X)
  , uVector_(GbVec3<T>::UNIT_Y)
  , dVector_(GbVec3<T>::UNIT_Z)
{
  lBound_ = T(1.0);
  uBound_ = T(1.0);
  dMin_ = T(1.0);
  dMax_ = T(2.0);

  update();
}

template <class T>
void 
GbFrustum<T>::computeVertices (GbVec3<T> akVertex[8]) const
{
  GbVec3<T> kDScaled = dMin_*dVector_;
  GbVec3<T> kLScaled = lBound_*lVector_;
  GbVec3<T> kUScaled = uBound_*uVector_;

  akVertex[0] = kDScaled - kLScaled - kUScaled;
  akVertex[1] = kDScaled - kLScaled + kUScaled;
  akVertex[2] = kDScaled + kLScaled + kUScaled;
  akVertex[3] = kDScaled + kLScaled - kUScaled;

  for (int i = 0, ip = 4; i < 4; i++, ip++) {
    akVertex[ip] = origin_ + dRatio_*akVertex[i];
    akVertex[i] += origin_;
  }
}

template <class T>
void 
GbFrustum<T>::update ()
{
  dRatio_ = dMax_/dMin_;
  mTwoLF_ = T(-2.0)*lBound_*dMax_;
  mTwoUF_ = T(-2.0)*uBound_*dMax_;
}
