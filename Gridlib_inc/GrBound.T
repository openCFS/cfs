// -*- C++ -*- automatisch in C++mode wechseln (emacs)
/*----------------------------------------------------------------------
|
|
| $Log$
| Revision 1.3  2002/04/03 11:17:08  elena
| new: changes in CFS++ on SGI and new function for gridlib
|
| Revision 1.2  2001/06/15 09:02:29  prkipfer
| added compiler hints and removed useless signedness
|
| Revision 1.1  2001/02/13 13:29:26  prkipfer
| introduced basic rendering tools
|
|
|
+---------------------------------------------------------------------*/

template <class T>
const T
GrBound<T>::TOLERANCE = 1e-06;

template <class T>
GrBound<T>& 
GrBound<T>::operator+= (const GrBound<T>& rkBound)
{
  // difference of current sphere center and input sphere center
  GbVec3<T> kDiff = center_ - rkBound.center_;
    
  T fLengthSqr = kDiff.getSquareNorm();
  T fDeltaRad = rkBound.radius_ - radius_;
  T fDeltaRadSqr = fDeltaRad*fDeltaRad;
  T fLength, fAlpha;
    
  if ( fDeltaRad >= 0.0 ) {
    if ( fDeltaRadSqr >= fLengthSqr ) {
      // rkBound's sphere encloses 'this' sphere
      center_ = rkBound.center_;
      radius_ = rkBound.radius_;
    }
    else {
      // this' sphere does not enclose rkBound sphere
      fLength = GbMath<T>::Sqrt(fLengthSqr);
      if ( fLength > GrBound<T>::TOLERANCE ) {
	fAlpha = (fLength - fDeltaRad)/(2.0*fLength);
	center_ = rkBound.center_ + fAlpha*kDiff;
      }
      radius_ = 0.5*(rkBound.radius_+fLength+radius_);
    }
  }
  else if ( fDeltaRadSqr < fLengthSqr ) {
    // this' sphere does not enclose rkBound sphere
    fLength = GbMath<T>::Sqrt(fLengthSqr);
    if ( fLength > GrBound<T>::TOLERANCE ) {
      fAlpha = (fLength - fDeltaRad)/(2.0*fLength);
      center_ = rkBound.center_ + fAlpha*kDiff;
    }
    radius_ = 0.5*(rkBound.radius_+fLength+radius_);
  }
  // else 'this' sphere encloses rkBound sphere

  return *this;
}

template <class T>
GrBound<T> 
GrBound<T>::transformBy (const GbMatrix3<T>& rkRotate,
			 const GbVec3<T>& rkTranslate, T fScale) const
{
  GrBound<T> bound;
  bound.center_ = fScale*(rkRotate*center_) + rkTranslate;
  bound.radius_ = fScale*radius_;
  return bound;
}

template <class T>
void 
GrBound<T>::computeFromData (int uiQuantity, const GbVec3<T>* akData)
{
  if ( uiQuantity <= 0 ) {
    center_ = GbVec3<T>::ZERO;
    radius_ = 0.0;
    return;
  }

  // compute the axis-aligned box containing the data
  GbVec3<T> kMin = akData[0], kMax = kMin;
  int uiI;
  for (uiI = 1; uiI < uiQuantity; uiI++) {
    if ( kMin[0] > akData[uiI][0] ) kMin[0] = akData[uiI][0];
    if ( kMin[1] > akData[uiI][1] ) kMin[1] = akData[uiI][1];
    if ( kMin[2] > akData[uiI][2] ) kMin[2] = akData[uiI][2];
    if ( kMax[0] < akData[uiI][0] ) kMax[0] = akData[uiI][0];
    if ( kMax[1] < akData[uiI][1] ) kMax[1] = akData[uiI][1];
    if ( kMax[2] < akData[uiI][2] ) kMax[2] = akData[uiI][2];
  }

  // sphere center is the axis-aligned box center
  center_ = T(0.5)*(kMin+kMax);

  // compute the radius
  T fRadiusSqr = 0.0;
  for (uiI = 0; uiI < uiQuantity; uiI++) {
    GbVec3<T> kDiff = akData[uiI] - center_;
    T fLengthSqr = kDiff.getSquareNorm();
    if ( fLengthSqr > fRadiusSqr )
      fRadiusSqr = fLengthSqr;
  }
  radius_ = GbMath<T>::Sqrt(fRadiusSqr);
}

template <class T>
typename GbPlane<T>::Side 
GrBound<T>::whichSide (const GbPlane<T>& rkPlane) const
{
  T fPseudoDistance = rkPlane.distanceTo(center_);

  if ( fPseudoDistance <= -radius_) {
    return GbPlane<T>::NEGATIVE_SIDE;
  }
  else if ( fPseudoDistance >= radius_) {
    return GbPlane<T>::POSITIVE_SIDE;
  }
  else {
    return GbPlane<T>::NO_SIDE;
  }
}
