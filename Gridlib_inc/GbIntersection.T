// -*- C++ -*- automatisch in C++mode wechseln (emacs)
/*----------------------------------------------------------------------
|
|
| $Log$
| Revision 1.3  2002/04/03 11:17:07  elena
| new: changes in CFS++ on SGI and new function for gridlib
|
| Revision 1.7  2002/03/18 10:00:16  prkipfer
| added tet-tet intersection
|
| Revision 1.6  2001/12/18 13:26:30  prkipfer
| introduced sphere object
|
| Revision 1.5  2001/09/12 09:16:13  prkipfer
| added tri-box test
|
| Revision 1.4  2001/08/16 16:53:16  prkipfer
| improved type safety for template parameter
|
| Revision 1.3  2001/03/20 09:36:10  prkipfer
| added tri-tri test
|
| Revision 1.2  2001/02/13 10:58:38  prkipfer
| exchanged fixed numeric tolerances by machine dependent ones
|
| Revision 1.1  2001/01/03 15:30:08  prkipfer
| introduced intersection, linear system solver and minimizer
|
|
|
+---------------------------------------------------------------------*/

template <class T>
const T
GbIntersection<T>::TOLERANCE = std::numeric_limits<T>::epsilon(); // max. 1e-12

template <class T>
GbBool 
GbIntersection<T>::test (const GbSegment<T>& rkSegment, const GbBox3<T>& rkBox)
{
  T fAWdU[3], fADdU[3], fAWxDdU[3], fRhs;
  GbVec3<T> kSDir = T(0.5)*rkSegment.getDirection();
  GbVec3<T> kSCen = rkSegment.getOrigin() + kSDir;

  GbVec3<T> kDiff = kSCen - rkBox.getCenter();
  
  fAWdU[0] = GbMath<T>::Abs(kSDir.dot(rkBox.getAxis(0)));
  fADdU[0] = GbMath<T>::Abs(kDiff.dot(rkBox.getAxis(0)));
  fRhs = rkBox.getExtent(0) + fAWdU[0];
  if ( fADdU[0] > fRhs ) return false;

  fAWdU[1] = GbMath<T>::Abs(kSDir.dot(rkBox.getAxis(1)));
  fADdU[1] = GbMath<T>::Abs(kDiff.dot(rkBox.getAxis(1)));
  fRhs = rkBox.getExtent(1) + fAWdU[1];
  if ( fADdU[1] > fRhs ) return false;

  fAWdU[2] = GbMath<T>::Abs(kSDir.dot(rkBox.getAxis(2)));
  fADdU[2] = GbMath<T>::Abs(kDiff.dot(rkBox.getAxis(2)));
  fRhs = rkBox.getExtent(2) + fAWdU[2];
  if ( fADdU[2] > fRhs ) return false;

  GbVec3<T> kWxD = kSDir.cross(kDiff);

  fAWxDdU[0] = GbMath<T>::Abs(kWxD.dot(rkBox.getAxis(0)));
  fRhs = rkBox.getExtent(1)*fAWdU[2] + rkBox.getExtent(2)*fAWdU[1];
  if ( fAWxDdU[0] > fRhs ) return false;

  fAWxDdU[1] = GbMath<T>::Abs(kWxD.dot(rkBox.getAxis(1)));
  fRhs = rkBox.getExtent(0)*fAWdU[2] + rkBox.getExtent(2)*fAWdU[0];
  if ( fAWxDdU[1] > fRhs ) return false;

  fAWxDdU[2] = GbMath<T>::Abs(kWxD.dot(rkBox.getAxis(2)));
  fRhs = rkBox.getExtent(0)*fAWdU[1] + rkBox.getExtent(1)*fAWdU[0];
  if ( fAWxDdU[2] > fRhs ) return false;

  return true;
}

template <class T>
GbBool 
GbIntersection<T>::test (const GbRay<T>& rkRay, const GbBox3<T>& rkBox)
{
  T fWdU[3], fAWdU[3], fDdU[3], fADdU[3], fAWxDdU[3], fRhs;

  GbVec3<T> kDiff = rkRay.getOrigin() - rkBox.getCenter();

  fWdU[0] = rkRay.getDirection().dot(rkBox.getAxis(0));
  fAWdU[0] = GbMath<T>::Abs(fWdU[0]);
  fDdU[0] = kDiff.dot(rkBox.getAxis(0));
  fADdU[0] = GbMath<T>::Abs(fDdU[0]);
  if ( fADdU[0] > rkBox.getExtent(0) && fDdU[0]*fWdU[0] >= T(0.0) )
    return false;

  fWdU[1] = rkRay.getDirection().dot(rkBox.getAxis(1));
  fAWdU[1] = GbMath<T>::Abs(fWdU[1]);
  fDdU[1] = kDiff.dot(rkBox.getAxis(1));
  fADdU[1] = GbMath<T>::Abs(fDdU[1]);
  if ( fADdU[1] > rkBox.getExtent(1) && fDdU[1]*fWdU[1] >= T(0.0) )
    return false;

  fWdU[2] = rkRay.getDirection().dot(rkBox.getAxis(2));
  fAWdU[2] = GbMath<T>::Abs(fWdU[2]);
  fDdU[2] = kDiff.dot(rkBox.getAxis(2));
  fADdU[2] = GbMath<T>::Abs(fDdU[2]);
  if ( fADdU[2] > rkBox.getExtent(2) && fDdU[2]*fWdU[2] >= T(0.0) )
    return false;

  GbVec3<T> kWxD = rkRay.getDirection().cross(kDiff);

  fAWxDdU[0] = GbMath<T>::Abs(kWxD.dot(rkBox.getAxis(0)));
  fRhs = rkBox.getExtent(1)*fAWdU[2] + rkBox.getExtent(2)*fAWdU[1];
  if ( fAWxDdU[0] > fRhs )
    return false;

  fAWxDdU[1] = GbMath<T>::Abs(kWxD.dot(rkBox.getAxis(1)));
  fRhs = rkBox.getExtent(0)*fAWdU[2] + rkBox.getExtent(2)*fAWdU[0];
  if ( fAWxDdU[1] > fRhs )
    return false;

  fAWxDdU[2] = GbMath<T>::Abs(kWxD.dot(rkBox.getAxis(2)));
  fRhs = rkBox.getExtent(0)*fAWdU[1] + rkBox.getExtent(1)*fAWdU[0];
  if ( fAWxDdU[2] > fRhs )
    return false;

  return true;
}

template <class T>
GbBool 
GbIntersection<T>::find (const GbVec3<T>& rkOrigin, const GbVec3<T>& rkDirection, 
			 const T afExtent[3], T& rfT0, T& rfT1)
{
  T fSaveT0 = rfT0, fSaveT1 = rfT1;

  GbBool bNotEntirelyClipped =
    clip(+rkDirection[0],-rkOrigin[0]-afExtent[0],rfT0,rfT1) &&
    clip(-rkDirection[0],+rkOrigin[0]-afExtent[0],rfT0,rfT1) &&
    clip(+rkDirection[1],-rkOrigin[1]-afExtent[1],rfT0,rfT1) &&
    clip(-rkDirection[1],+rkOrigin[1]-afExtent[1],rfT0,rfT1) &&
    clip(+rkDirection[2],-rkOrigin[2]-afExtent[2],rfT0,rfT1) &&
    clip(-rkDirection[2],+rkOrigin[2]-afExtent[2],rfT0,rfT1);

  return bNotEntirelyClipped && ( rfT0 != fSaveT0 || rfT1 != fSaveT1 );
}

template <class T>
GbBool 
GbIntersection<T>::find (const GbSegment<T>& rkSegment, const GbBox3<T>& rkBox,
			 int& riQuantity, GbVec3<T> akPoint[2])
{
  // convert segment to box coordinates
  GbVec3<T> kDiff = rkSegment.getOrigin() - rkBox.getCenter();
  GbVec3<T> kOrigin(kDiff.dot(rkBox.getAxis(0)),
		    kDiff.dot(rkBox.getAxis(1)),
		    kDiff.dot(rkBox.getAxis(2)) );
  GbVec3<T> kDirection(rkSegment.getDirection().dot(rkBox.getAxis(0)),
		       rkSegment.getDirection().dot(rkBox.getAxis(1)),
		       rkSegment.getDirection().dot(rkBox.getAxis(2)) );

  T fT0 = T(0.0), fT1 = T(1.0);
  GbBool bIntersects = GbIntersection<T>::find(kOrigin,kDirection,rkBox.getExtents(),fT0,fT1);

  if ( bIntersects ) {
    if ( fT0 > T(0.0) ) {
      if ( fT1 < T(1.0) ) {
	riQuantity = 2;
	akPoint[0] = rkSegment.getOrigin() + fT0*rkSegment.getDirection();
	akPoint[1] = rkSegment.getOrigin() + fT1*rkSegment.getDirection();
      }
      else {
	riQuantity = 1;
	akPoint[0] = rkSegment.getOrigin() + fT0*rkSegment.getDirection();
      }
    }
    else {  // fT0 == 0.0
      if ( fT1 < T(1.0) ) {
	riQuantity = 1;
	akPoint[0] = rkSegment.getOrigin() + fT1*rkSegment.getDirection();
      }
      else {  // fT1 == 1.0
	// segment entirely in box
	riQuantity = 0;
      }
    }
  }
  else {
    riQuantity = 0;
  }
  
  return bIntersects;
}

template <class T>
GbBool 
GbIntersection<T>::find (const GbRay<T>& rkRay, const GbBox3<T>& rkBox,
			 int& riQuantity, GbVec3<T> akPoint[2])
{
  // convert ray to box coordinates
  GbVec3<T> kDiff = rkRay.getOrigin() - rkBox.getCenter();
  GbVec3<T> kOrigin(kDiff.dot(rkBox.getAxis(0)),
		    kDiff.dot(rkBox.getAxis(1)),
		    kDiff.dot(rkBox.getAxis(2)));
  GbVec3<T> kDirection(rkRay.getDirection().dot(rkBox.getAxis(0)),
		       rkRay.getDirection().dot(rkBox.getAxis(1)),
		       rkRay.getDirection().dot(rkBox.getAxis(2)));

  T fT0 = T(0.0), fT1 = std::numeric_limits<T>::max();
  GbBool bIntersects = GbIntersection<T>::find(kOrigin,kDirection,rkBox.getExtents(),fT0,fT1);

  if ( bIntersects ) {
    if ( fT0 > T(0.0) ) {
      riQuantity = 2;
      akPoint[0] = rkRay.getOrigin() + fT0*rkRay.getDirection();
      akPoint[1] = rkRay.getOrigin() + fT1*rkRay.getDirection();
    }
    else { // fT0 == 0.0
      riQuantity = 1;
      akPoint[0] = rkRay.getOrigin() + fT1*rkRay.getDirection();
    }
  }
  else {
    riQuantity = 0;
  }

  return bIntersects;
}

template <class T>
GbBool 
GbIntersection<T>::test (const GbSegment<T>& rkSegment, const GoTriangleElement& rkTriangle)
{
  return GbDistance<T>::squareDistance(rkSegment,rkTriangle) <= GbIntersection::TOLERANCE;
}

template <class T>
GbBool 
GbIntersection<T>::find (const GbSegment<T>& rkSegment, const GoTriangleElement& rkTriangle, 
			 GbVec3<T>& rkPoint)
{
  T fSegP, dummy1, dummy2;
  if ( GbDistance<T>::squareDistance(rkSegment,rkTriangle,&fSegP,&dummy1,&dummy2) <= GbIntersection::TOLERANCE ) {
    rkPoint = rkSegment.getOrigin() + fSegP*rkSegment.getDirection();
    return true;
  }
  return false;
}

template <class T>
GbBool 
GbIntersection<T>::test (const GbRay<T>& rkRay, const GoTriangleElement& rkTriangle)
{
  return GbDistance<T>::squareDistance(rkRay,rkTriangle) <= GbIntersection::TOLERANCE;
}

template <class T>
GbBool 
GbIntersection<T>::find (const GbRay<T>& rkRay, const GoTriangleElement& rkTriangle, 
			 GbVec3<T>& rkPoint)
{
  T fRayP,dummy1,dummy2;
  if ( GbDistance<T>::squareDistance(rkRay,rkTriangle,&fRayP,&dummy1,&dummy2) <= GbIntersection::TOLERANCE ) {
    rkPoint = rkRay.getOrigin() + fRayP*rkRay.getDirection();
    return true;
  }
  return false;
}

template <class T>
GbBool 
GbIntersection<T>::test (const GbBox3<T>& rkBox0, const GbBox3<T>& rkBox1)
{
  // convenience variables
  const GbVec3<T>* akA = rkBox0.getAxes();
  const GbVec3<T>* akB = rkBox1.getAxes();
  const T* afEA = rkBox0.getExtents();
  const T* afEB = rkBox1.getExtents();

  // compute difference of box centers, D = C1-C0
  GbVec3<T> kD = rkBox1.getCenter() - rkBox0.getCenter();

  T aafC[3][3];     // matrix C = A^T B, c_{ij} = Dot(A_i,B_j)
  T aafAbsC[3][3];  // |c_{ij}|
  T afAD[3];        // Dot(A_i,D)
  T fR0, fR1, fR;   // interval radii and distance between centers
  T fR01;           // = R0 + R1
    
  // axis C0+t*A0
  aafC[0][0] = akA[0].dot(akB[0]);
  aafC[0][1] = akA[0].dot(akB[1]);
  aafC[0][2] = akA[0].dot(akB[2]);
  afAD[0] = akA[0].dot(kD);
  aafAbsC[0][0] = GbMath<T>::Abs(aafC[0][0]);
  aafAbsC[0][1] = GbMath<T>::Abs(aafC[0][1]);
  aafAbsC[0][2] = GbMath<T>::Abs(aafC[0][2]);
  fR = GbMath<T>::Abs(afAD[0]);
  fR1 = afEB[0]*aafAbsC[0][0]+afEB[1]*aafAbsC[0][1]+afEB[2]*aafAbsC[0][2];
  fR01 = afEA[0] + fR1;
  if ( fR > fR01 ) return false;

  // axis C0+t*A1
  aafC[1][0] = akA[1].dot(akB[0]);
  aafC[1][1] = akA[1].dot(akB[1]);
  aafC[1][2] = akA[1].dot(akB[2]);
  afAD[1] = akA[1].dot(kD);
  aafAbsC[1][0] = GbMath<T>::Abs(aafC[1][0]);
  aafAbsC[1][1] = GbMath<T>::Abs(aafC[1][1]);
  aafAbsC[1][2] = GbMath<T>::Abs(aafC[1][2]);
  fR = GbMath<T>::Abs(afAD[1]);
  fR1 = afEB[0]*aafAbsC[1][0]+afEB[1]*aafAbsC[1][1]+afEB[2]*aafAbsC[1][2];
  fR01 = afEA[1] + fR1;
  if ( fR > fR01 ) return false;

  // axis C0+t*A2
  aafC[2][0] = akA[2].dot(akB[0]);
  aafC[2][1] = akA[2].dot(akB[1]);
  aafC[2][2] = akA[2].dot(akB[2]);
  afAD[2] = akA[2].dot(kD);
  aafAbsC[2][0] = GbMath<T>::Abs(aafC[2][0]);
  aafAbsC[2][1] = GbMath<T>::Abs(aafC[2][1]);
  aafAbsC[2][2] = GbMath<T>::Abs(aafC[2][2]);
  fR = GbMath<T>::Abs(afAD[2]);
  fR1 = afEB[0]*aafAbsC[2][0]+afEB[1]*aafAbsC[2][1]+afEB[2]*aafAbsC[2][2];
  fR01 = afEA[2] + fR1;
  if ( fR > fR01 ) return false;

  // axis C0+t*B0
  fR = GbMath<T>::Abs(akB[0].dot(kD));
  fR0 = afEA[0]*aafAbsC[0][0]+afEA[1]*aafAbsC[1][0]+afEA[2]*aafAbsC[2][0];
  fR01 = fR0 + afEB[0];
  if ( fR > fR01 ) return false;

  // axis C0+t*B1
  fR = GbMath<T>::Abs(akB[1].dot(kD));
  fR0 = afEA[0]*aafAbsC[0][1]+afEA[1]*aafAbsC[1][1]+afEA[2]*aafAbsC[2][1];
  fR01 = fR0 + afEB[1];
  if ( fR > fR01 ) return false;

  // axis C0+t*B2
  fR = GbMath<T>::Abs(akB[2].dot(kD));
  fR0 = afEA[0]*aafAbsC[0][2]+afEA[1]*aafAbsC[1][2]+afEA[2]*aafAbsC[2][2];
  fR01 = fR0 + afEB[2];
  if ( fR > fR01 ) return false;

  // axis C0+t*A0xB0
  fR = GbMath<T>::Abs(afAD[2]*aafC[1][0]-afAD[1]*aafC[2][0]);
  fR0 = afEA[1]*aafAbsC[2][0] + afEA[2]*aafAbsC[1][0];
  fR1 = afEB[1]*aafAbsC[0][2] + afEB[2]*aafAbsC[0][1];
  fR01 = fR0 + fR1;
  if ( fR > fR01 ) return false;

  // axis C0+t*A0xB1
  fR = GbMath<T>::Abs(afAD[2]*aafC[1][1]-afAD[1]*aafC[2][1]);
  fR0 = afEA[1]*aafAbsC[2][1] + afEA[2]*aafAbsC[1][1];
  fR1 = afEB[0]*aafAbsC[0][2] + afEB[2]*aafAbsC[0][0];
  fR01 = fR0 + fR1;
  if ( fR > fR01 ) return false;

  // axis C0+t*A0xB2
  fR = GbMath<T>::Abs(afAD[2]*aafC[1][2]-afAD[1]*aafC[2][2]);
  fR0 = afEA[1]*aafAbsC[2][2] + afEA[2]*aafAbsC[1][2];
  fR1 = afEB[0]*aafAbsC[0][1] + afEB[1]*aafAbsC[0][0];
  fR01 = fR0 + fR1;
  if ( fR > fR01 ) return false;

  // axis C0+t*A1xB0
  fR = GbMath<T>::Abs(afAD[0]*aafC[2][0]-afAD[2]*aafC[0][0]);
  fR0 = afEA[0]*aafAbsC[2][0] + afEA[2]*aafAbsC[0][0];
  fR1 = afEB[1]*aafAbsC[1][2] + afEB[2]*aafAbsC[1][1];
  fR01 = fR0 + fR1;
  if ( fR > fR01 ) return false;

  // axis C0+t*A1xB1
  fR = GbMath<T>::Abs(afAD[0]*aafC[2][1]-afAD[2]*aafC[0][1]);
  fR0 = afEA[0]*aafAbsC[2][1] + afEA[2]*aafAbsC[0][1];
  fR1 = afEB[0]*aafAbsC[1][2] + afEB[2]*aafAbsC[1][0];
  fR01 = fR0 + fR1;
  if ( fR > fR01 ) return false;

  // axis C0+t*A1xB2
  fR = GbMath<T>::Abs(afAD[0]*aafC[2][2]-afAD[2]*aafC[0][2]);
  fR0 = afEA[0]*aafAbsC[2][2] + afEA[2]*aafAbsC[0][2];
  fR1 = afEB[0]*aafAbsC[1][1] + afEB[1]*aafAbsC[1][0];
  fR01 = fR0 + fR1;
  if ( fR > fR01 ) return false;

  // axis C0+t*A2xB0
  fR = GbMath<T>::Abs(afAD[1]*aafC[0][0]-afAD[0]*aafC[1][0]);
  fR0 = afEA[0]*aafAbsC[1][0] + afEA[1]*aafAbsC[0][0];
  fR1 = afEB[1]*aafAbsC[2][2] + afEB[2]*aafAbsC[2][1];
  fR01 = fR0 + fR1;
  if ( fR > fR01 ) return false;

  // axis C0+t*A2xB1
  fR = GbMath<T>::Abs(afAD[1]*aafC[0][1]-afAD[0]*aafC[1][1]);
  fR0 = afEA[0]*aafAbsC[1][1] + afEA[1]*aafAbsC[0][1];
  fR1 = afEB[0]*aafAbsC[2][2] + afEB[2]*aafAbsC[2][0];
  fR01 = fR0 + fR1;
  if ( fR > fR01 ) return false;

  // axis C0+t*A2xB2
  fR = GbMath<T>::Abs(afAD[1]*aafC[0][2]-afAD[0]*aafC[1][2]);
  fR0 = afEA[0]*aafAbsC[1][2] + afEA[1]*aafAbsC[0][2];
  fR1 = afEB[0]*aafAbsC[2][1] + afEB[1]*aafAbsC[2][0];
  fR01 = fR0 + fR1;
  if ( fR > fR01 ) return false;

  return true;
}

template <class T>
GbBool 
GbIntersection<T>::test (const GbLine<T>& rkLine, const GbBox3<T>& rkBox)
{
  T fAWdU[3], fAWxDdU[3], fRhs;

  GbVec3<T> kDiff = rkLine.getOrigin() - rkBox.getCenter();
  GbVec3<T> kWxD = rkLine.getDirection().cross(kDiff);

  fAWdU[1] = GbMath<T>::Abs(rkLine.getDirection().dot(rkBox.getAxis(1)));
  fAWdU[2] = GbMath<T>::Abs(rkLine.getDirection().dot(rkBox.getAxis(2)));
  fAWxDdU[0] = GbMath<T>::Abs(kWxD.dot(rkBox.getAxis(0)));
  fRhs = rkBox.getExtent(1)*fAWdU[2] + rkBox.getExtent(2)*fAWdU[1];
  if ( fAWxDdU[0] > fRhs )
    return false;

  fAWdU[0] = GbMath<T>::Abs(rkLine.getDirection().dot(rkBox.getAxis(0)));
  fAWxDdU[1] = GbMath<T>::Abs(kWxD.dot(rkBox.getAxis(1)));
  fRhs = rkBox.getExtent(0)*fAWdU[2] + rkBox.getExtent(2)*fAWdU[0];
  if ( fAWxDdU[1] > fRhs )
    return false;

  fAWxDdU[2] = GbMath<T>::Abs(kWxD.dot(rkBox.getAxis(2)));
  fRhs = rkBox.getExtent(0)*fAWdU[1] + rkBox.getExtent(1)*fAWdU[0];
  if ( fAWxDdU[2] > fRhs )
    return false;

  return true;
}

template <class T>
GbBool 
GbIntersection<T>::find (const GbLine<T>& rkLine, const GbBox3<T>& rkBox,
			 int& riQuantity, GbVec3<T> akPoint[2])
{
  // convert line to box coordinates
  GbVec3<T> kDiff = rkLine.getOrigin() - rkBox.getCenter();
  GbVec3<T> kOrigin(kDiff.dot(rkBox.getAxis(0)),
		    kDiff.dot(rkBox.getAxis(1)),
		    kDiff.dot(rkBox.getAxis(2)));
  GbVec3<T> kDirection(rkLine.getDirection().dot(rkBox.getAxis(0)),
		       rkLine.getDirection().dot(rkBox.getAxis(1)),
		       rkLine.getDirection().dot(rkBox.getAxis(2)));

  T fT0 = -std::numeric_limits<T>::max(), fT1 = std::numeric_limits<T>::max();
  GbBool bIntersects = GbIntersection<T>::find(kOrigin,kDirection,rkBox.getExtents(),fT0,fT1);

  if ( bIntersects ) {
    if ( fT0 != fT1 ) {
      riQuantity = 2;
      akPoint[0] = rkLine.getOrigin() + fT0*rkLine.getDirection();
      akPoint[1] = rkLine.getOrigin() + fT1*rkLine.getDirection();
    }
    else {
      riQuantity = 1;
      akPoint[0] = rkLine.getOrigin() + fT0*rkLine.getDirection();
    }
  }
  else {
    riQuantity = 0;
  }

  return bIntersects;
}

template <class T>
GbBool 
GbIntersection<T>::test (const GbLine<T>& rkLine, const GoTriangleElement& rkTriangle)
{
  return GbDistance<T>::squareDistance(rkLine,rkTriangle) <= GbIntersection<T>::TOLERANCE;
}

template <class T>
GbBool 
GbIntersection<T>::find (const GbLine<T>& rkLine, const GoTriangleElement& rkTriangle, 
			 GbVec3<T>& rkPoint)
{
  T fLinP,dummy1,dummy2;
  if (GbDistance<T>::squareDistance(rkLine,rkTriangle,&fLinP,&dummy1,&dummy2) <= GbIntersection<T>::TOLERANCE) {
    rkPoint = rkLine.getOrigin() + fLinP*rkLine.getDirection();
    return true;
  }
  return false;
}

template <class T>
GbBool 
GbIntersection<T>::test (const GbPlane<T>& rkPlane, const GbBox3<T>& rkBox)
{
  T fTmp[3] = {
    rkBox.getExtent(0)*(rkPlane.getNormal().dot(rkBox.getAxis(0))),
    rkBox.getExtent(1)*(rkPlane.getNormal().dot(rkBox.getAxis(1))),
    rkBox.getExtent(2)*(rkPlane.getNormal().dot(rkBox.getAxis(2)))
  };

  T fRadius = GbMath<T>::Abs(fTmp[0]) + GbMath<T>::Abs(fTmp[1]) + GbMath<T>::Abs(fTmp[2]);

  T fPseudoDistance = rkPlane.distanceTo(rkBox.getCenter());
  return GbMath<T>::Abs(fPseudoDistance) <= fRadius;
}

template <class T>
GbBool 
GbIntersection<T>::test (const GbPlane<T>& rkPlane0, const GbPlane<T>& rkPlane1)
{
  GbVec3<T> kCross = rkPlane0.getNormal().cross(rkPlane1.getNormal());
  T fSqrLength = kCross.getSquareNorm();
  return fSqrLength > GbIntersection<T>::TOLERANCE;
}

template <class T>
GbBool 
GbIntersection<T>::find (const GbPlane<T>& rkPlane0, const GbPlane<T>& rkPlane1,
			 GbLine<T>& rkLine)
{
  // If Cross(N0,N1) is zero, then either planes are parallel and separated
  // or the same plane.  In both cases, 'false' is returned.  Otherwise,
  // the intersection line is
  //
  //   L(t) = t*Cross(N0,N1) + c0*N0 + c1*N1
  //
  // for some coefficients c0 and c1 and for t any real number (the line
  // parameter).  Taking dot products with the normals,
  //
  //   d0 = Dot(N0,L) = c0*Dot(N0,N0) + c1*Dot(N0,N1)
  //   d1 = Dot(N1,L) = c0*Dot(N0,N1) + c1*Dot(N1,N1)
  //
  // which are two equations in two unknowns.  The solution is
  //
  //   c0 = (Dot(N1,N1)*d0 - Dot(N0,N1)*d1)/det
  //   c1 = (Dot(N0,N0)*d1 - Dot(N0,N1)*d0)/det
  //
  // where det = Dot(N0,N0)*Dot(N1,N1)-Dot(N0,N1)^2.

  T fN00 = rkPlane0.getNormal().getSquareNorm();
  T fN01 = rkPlane0.getNormal().dot(rkPlane1.getNormal());
  T fN11 = rkPlane1.getNormal().getSquareNorm();
  T fDet = fN00*fN11 - fN01*fN01;

  if ( GbMath<T>::Abs(fDet) < GbIntersection<T>::TOLERANCE )
    return false;

  T fInvDet = T(1.0)/fDet;
  T fC0 = (fN11*rkPlane0.getConstant() - fN01*rkPlane1.getConstant())*fInvDet;
  T fC1 = (fN00*rkPlane1.getConstant() - fN01*rkPlane0.getConstant())*fInvDet;

  rkLine.setDirection(rkPlane0.getNormal().cross(rkPlane1.getNormal()));
  rkLine.setOrigin(fC0*rkPlane0.getNormal() + fC1*rkPlane1.getNormal());
  return true;
}

template <class T>
GbBool 
GbIntersection<T>::test(const GoTriangleElement& triU, const GoTriangleElement& triV)
{
  GbVec3<T> kDir;
  T fUMin, fUMax, fVMin, fVMax;
  int i0, i1;

  // cache vertex positions
  GbVec3<T> akU[3] = { triU.getVertex(0)->getPosition(),
		       triU.getVertex(1)->getPosition(),
		       triU.getVertex(2)->getPosition() };
  GbVec3<T> akV[3] = { triV.getVertex(0)->getPosition(),
		       triV.getVertex(1)->getPosition(),
		       triV.getVertex(2)->getPosition() };

  // direction N
  GbVec3<T> akE[3] = { akU[1]-akU[0], akU[2]-akU[1], akU[0]-akU[2] };
  GbVec3<T> kN = akE[0].cross(akE[1]);
  T fNdU0 = kN.dot(akU[0]);
  projection(kN,akV,fVMin,fVMax);
  if ( fNdU0 < fVMin || fNdU0 > fVMax )
    return false;

  // direction M
  GbVec3<T> akF[3] = { akV[1]-akV[0], akV[2]-akV[1], akV[0]-akV[2] };
  GbVec3<T> kM = akF[0].cross(akF[1]);

  // cutoff: sin(Angle(N,M)) < TOLERANCE
  GbVec3<T> kNxM = kN.cross(kM);
  if ( kNxM.dot(kNxM) >= GbIntersection<T>::TOLERANCE*(kN.dot(kN)*kM.dot(kM)) ) {
    // triangles are not parallel
    T fMdV0 = kM.dot(akV[0]);
    projection(kM,akU,fUMin,fUMax);
    if ( fMdV0 < fUMin || fMdV0 > fUMax )
      return false;

    // directions E[i0]xF[i1]
    for (i1 = 0; i1 < 3; i1++) {
      for (i0 = 0; i0 < 3; i0++) {
	kDir = akE[i0].cross(akF[i1]);
	projection(kDir,akU,fUMin,fUMax);
	projection(kDir,akV,fVMin,fVMax);
	if ( fUMax < fVMin || fVMax < fUMin )
	  return false;
      }
    }
  }
  else {
    // triangles are parallel (and, in fact, coplanar)
    // directions NxE[i0]
    for (i0 = 0; i0 < 3; i0++) {
      kDir = kN.cross(akE[i0]);
      projection(kDir,akU,fUMin,fUMax);
      projection(kDir,akV,fVMin,fVMax);
      if ( fUMax < fVMin || fVMax < fUMin )
	return false;
    }

    // directions NxF[i1]
    for (i1 = 0; i1 < 3; i1++) {
      kDir = kM.cross(akF[i1]);
      projection(kDir,akU,fUMin,fUMax);
      projection(kDir,akV,fVMin,fVMax);
      if ( fUMax < fVMin || fVMax < fUMin )
	return false;
    }
  }

  return true;
}

template <class T>
GbBool 
GbIntersection<T>::test(const GoTriangleElement& rkTriangle, const GbBox3<T>& rkBox)
{
  T fMin0, fMax0, fMin1, fMax1;
  GbVec3<T> kD, akE[3];

  // test direction of triangle normal
  akE[0] = rkTriangle.getEdgeVec(0);
  akE[1] = rkTriangle.getEdgeVec(1);
  kD = akE[0].cross(akE[1]);
  fMin0 = kD.dot(rkTriangle.getVertex(0)->getPosition());
  fMax0 = fMin0;
  projectBox(kD,rkBox,fMin1,fMax1);
  if ( fMax1 < fMin0 || fMax0 < fMin1 )
    return false;

  // test direction of box faces
  for (int i = 0; i < 3; ++i) {
    kD = rkBox.getAxis(i);
    projectTriangle(kD,rkTriangle,fMin0,fMax0);
    T fDdC = kD.dot(rkBox.getCenter());
    fMin1 = fDdC - rkBox.getExtent(i);
    fMax1 = fDdC + rkBox.getExtent(i);
    if ( fMax1 < fMin0 || fMax0 < fMin1 )
      return false;
  }

  // test direction of triangle-box edge cross products
  akE[2] = akE[1] - akE[0];
  for (int i0 = 0; i0 < 3; ++i0) {
    for (int i1 = 0; i1 < 3; ++i1) {
      kD = akE[i0].cross(rkBox.getAxis(i1));
      projectTriangle(kD,rkTriangle,fMin0,fMax0);
      projectBox(kD,rkBox,fMin1,fMax1);
      if ( fMax1 < fMin0 || fMax0 < fMin1 )
	return false;
    }
  }

  return true;
}

template <class T>
GbBool
GbIntersection<T>::test (const GbBox3<T>& rkBox, const GbFrustum<T>& rkFrustum)
{
  // for convenience
  const GbVec3<T>* akA = rkBox.getAxes();
  const T* afE = rkBox.getExtents();

  GbVec3<T> kDiff = rkBox.getCenter() - rkFrustum.origin();

  T afA[3], afB[3], afC[3], afD[3];
  T afNA[3], afNB[3], afNC[3], afND[3];
  T afNApLC[3], afNAmLC[3], afNBpUC[3], afNBmUC[3];
  T afLC[3], afLD[3], afUC[3], afUD[3], afLBpUA[3], afLBmUA[3];
  T fDdD, fR, fP, fMin, fMax, fMTwoLF, fMTwoUF, fLB, fUA, fTmp;
  int i, j;

  // M = D
  afD[2] = kDiff.dot(rkFrustum.dVector());
  for (i = 0; i < 3; ++i)
    afC[i] = akA[i].dot(rkFrustum.dVector());
  fR = afE[0]*GbMath<T>::Abs(afC[0]) +
       afE[1]*GbMath<T>::Abs(afC[1]) +
       afE[2]*GbMath<T>::Abs(afC[2]);
  if ( afD[2] + fR < rkFrustum.dMin() || afD[2] - fR > rkFrustum.dMax() )
    return false;

  // M = n*L - l*D
  for (i = 0; i < 3; ++i) {
    afA[i] = akA[i].dot(rkFrustum.lVector());
    afLC[i] = rkFrustum.lBound()*afC[i];
    afNA[i] = rkFrustum.dMin()*afA[i];
    afNAmLC[i] = afNA[i] - afLC[i];
  }
  afD[0] = kDiff.dot(rkFrustum.lVector());
  fR = afE[0]*GbMath<T>::Abs(afNAmLC[0]) +
       afE[1]*GbMath<T>::Abs(afNAmLC[1]) +
       afE[2]*GbMath<T>::Abs(afNAmLC[2]);
  afND[0] = rkFrustum.dMin()*afD[0];
  afLD[2] = rkFrustum.lBound()*afD[2];
  fDdD = afND[0] - afLD[2];
  fMTwoLF = rkFrustum.getMTwoLF();
  if ( fDdD + fR < fMTwoLF || fDdD > fR )
    return false;

  // M = -n*L - l*D
  for (i = 0; i < 3; ++i)
    afNApLC[i] = afNA[i] + afLC[i];
  fR = afE[0]*GbMath<T>::Abs(afNApLC[0]) +
       afE[1]*GbMath<T>::Abs(afNApLC[1]) +
       afE[2]*GbMath<T>::Abs(afNApLC[2]);
  fDdD = -(afND[0] + afLD[2]);
  if ( fDdD + fR < fMTwoLF || fDdD > fR )
    return false;

  // M = n*U - u*D
  for (i = 0; i < 3; ++i) {
    afB[i] = akA[i].dot(rkFrustum.uVector());
    afUC[i] = rkFrustum.uBound()*afC[i];
    afNB[i] = rkFrustum.dMin()*afB[i];
    afNBmUC[i] = afNB[i] - afUC[i];
  }
  afD[1] = kDiff.dot(rkFrustum.uVector());
  fR = afE[0]*GbMath<T>::Abs(afNBmUC[0]) +
       afE[1]*GbMath<T>::Abs(afNBmUC[1]) +
       afE[2]*GbMath<T>::Abs(afNBmUC[2]);
  afND[1] = rkFrustum.dMin()*afD[1];
  afUD[2] = rkFrustum.uBound()*afD[2];
  fDdD = afND[1] - afUD[2];
  fMTwoUF = rkFrustum.getMTwoUF();
  if ( fDdD + fR < fMTwoUF || fDdD > fR )
    return false;

  // M = -n*U - u*D
  for (i = 0; i < 3; ++i)
    afNBpUC[i] = afNB[i] + afUC[i];
  fR = afE[0]*GbMath<T>::Abs(afNBpUC[0]) +
       afE[1]*GbMath<T>::Abs(afNBpUC[1]) +
       afE[2]*GbMath<T>::Abs(afNBpUC[2]);
  fDdD = -(afND[1] + afUD[2]);
  if ( fDdD + fR < fMTwoUF || fDdD > fR )
    return false;

  // M = A[i]
  for (i = 0; i < 3; ++i) {
    fP = rkFrustum.lBound()*GbMath<T>::Abs(afA[i]) +
         rkFrustum.uBound()*GbMath<T>::Abs(afB[i]);
    afNC[i] = rkFrustum.dMin()*afC[i];
    fMin = afNC[i] - fP;
    if ( fMin < T(0.0) )
      fMin *= rkFrustum.getDRatio();
    fMax = afNC[i] + fP;
    if ( fMax > T(0.0) )
      fMax *= rkFrustum.getDRatio();
    fDdD = afA[i]*afD[0] + afB[i]*afD[1] + afC[i]*afD[2];
    if ( fDdD + afE[i] < fMin || fDdD - afE[i] > fMax )
      return false;
  }

  // M = Cross(L,A[i])
  for (i = 0; i < 3; ++i) {
    fP = rkFrustum.uBound()*GbMath<T>::Abs(afC[i]);
    fMin = afNB[i] - fP;
    if ( fMin < T(0.0) )
      fMin *= rkFrustum.getDRatio();
    fMax = afNB[i] + fP;
    if ( fMax > T(0.0) )
      fMax *= rkFrustum.getDRatio();
    fDdD = -afC[i]*afD[1] + afB[i]*afD[2];
    fR = afE[0]*GbMath<T>::Abs(afB[i]*afC[0]-afB[0]*afC[i]) +
         afE[1]*GbMath<T>::Abs(afB[i]*afC[1]-afB[1]*afC[i]) +
         afE[2]*GbMath<T>::Abs(afB[i]*afC[2]-afB[2]*afC[i]);
    if ( fDdD + fR < fMin || fDdD - fR > fMax )
      return false;
  }

  // M = Cross(U,A[i])
  for (i = 0; i < 3; ++i) {
    fP = rkFrustum.lBound()*GbMath<T>::Abs(afC[i]);
    fMin = -afNA[i] - fP;
    if ( fMin < T(0.0) )
      fMin *= rkFrustum.getDRatio();
    fMax = -afNA[i] + fP;
    if ( fMax > T(0.0) )
      fMax *= rkFrustum.getDRatio();
    fDdD = afC[i]*afD[0] - afA[i]*afD[2];
    fR = afE[0]*GbMath<T>::Abs(afA[i]*afC[0]-afA[0]*afC[i]) +
         afE[1]*GbMath<T>::Abs(afA[i]*afC[1]-afA[1]*afC[i]) +
         afE[2]*GbMath<T>::Abs(afA[i]*afC[2]-afA[2]*afC[i]);
    if ( fDdD + fR < fMin || fDdD - fR > fMax )
      return false;
  }

  // M = Cross(n*D+l*L+u*U,A[i])
  for (i = 0; i < 3; ++i) {
    fLB = rkFrustum.lBound()*afB[i];
    fUA = rkFrustum.uBound()*afA[i];
    afLBpUA[i] = fLB + fUA;
    afLBmUA[i] = fLB - fUA;
  }
  for (i = 0; i < 3; ++i) {
    fP = rkFrustum.lBound()*GbMath<T>::Abs(afNBmUC[i]) +
         rkFrustum.uBound()*GbMath<T>::Abs(afNAmLC[i]);
    fTmp = rkFrustum.dMin()*afLBmUA[i];
    fMin = fTmp - fP;
    if ( fMin < T(0.0) )
      fMin *= rkFrustum.getDRatio();
    fMax = fTmp + fP;
    if ( fMax > T(0.0) )
      fMax *= rkFrustum.getDRatio();
    fDdD = -afD[0]*afNBmUC[i] + afD[1]*afNAmLC[i] + afD[2]*afLBmUA[i];
    fR = T(0.0);
    for (j = 0; j < 3; ++j) {
      fR += afE[j]*GbMath<T>::Abs(-afA[j]*afNBmUC[i]
				  + afB[j]*afNAmLC[i]
				  + afC[j]*afLBmUA[i]);
    }
    if ( fDdD + fR < fMin || fDdD - fR > fMax )
      return false;
  }

  // M = Cross(n*D+l*L-u*U,A[i])
  for (i = 0; i < 3; ++i) {
    fP = rkFrustum.lBound()*GbMath<T>::Abs(afNBpUC[i]) +
         rkFrustum.uBound()*GbMath<T>::Abs(afNAmLC[i]);
    fTmp = rkFrustum.dMin()*afLBpUA[i];
    fMin = fTmp - fP;
    if ( fMin < T(0.0) )
      fMin *= rkFrustum.getDRatio();
    fMax = fTmp + fP;
    if ( fMax > T(0.0) )
      fMax *= rkFrustum.getDRatio();
    fDdD = -afD[0]*afNBpUC[i] + afD[1]*afNAmLC[i] + afD[2]*afLBpUA[i];
    fR = T(0.0);
    for (j = 0; j < 3; ++j) {
      fR += afE[j]*GbMath<T>::Abs(-afA[j]*afNBpUC[i]
				  + afB[j]*afNAmLC[i]
				  + afC[j]*afLBpUA[i]);
    }
    if ( fDdD + fR < fMin || fDdD - fR > fMax )
      return false;
  }

  // M = Cross(n*D-l*L+u*U,A[i])
  for (i = 0; i < 3; ++i) {
    fP = rkFrustum.lBound()*GbMath<T>::Abs(afNBmUC[i]) +
         rkFrustum.uBound()*GbMath<T>::Abs(afNApLC[i]);
    fTmp = -rkFrustum.dMin()*afLBpUA[i];
    fMin = fTmp - fP;
    if ( fMin < T(0.0) ) 
      fMin *= rkFrustum.getDRatio();
    fMax = fTmp + fP;
    if ( fMax > T(0.0) )
      fMax *= rkFrustum.getDRatio();
    fDdD = -afD[0]*afNBmUC[i] + afD[1]*afNApLC[i] - afD[2]*afLBpUA[i];
    fR = T(0.0);
    for (j = 0; j < 3; ++j) {
      fR += afE[j]*GbMath<T>::Abs(-afA[j]*afNBmUC[i]
				  + afB[j]*afNApLC[i]
				  - afC[j]*afLBpUA[i]);
    }
    if ( fDdD + fR < fMin || fDdD - fR > fMax )
      return false;
  }

  // M = Cross(n*D-l*L-u*U,A[i])
  for (i = 0; i < 3; ++i) {
    fP = rkFrustum.lBound()*GbMath<T>::Abs(afNBpUC[i]) +
         rkFrustum.uBound()*GbMath<T>::Abs(afNApLC[i]);
    fTmp = -rkFrustum.dMin()*afLBmUA[i];
    fMin = fTmp - fP;
    if ( fMin < T(0.0) )
      fMin *= rkFrustum.getDRatio();
    fMax = fTmp + fP;
    if ( fMax > T(0.0) )
      fMax *= rkFrustum.getDRatio();
    fDdD = -afD[0]*afNBpUC[i] + afD[1]*afNApLC[i] - afD[2]*afLBmUA[i];
    fR = T(0.0);
    for (j = 0; j < 3; ++j) {
      fR += afE[j]*GbMath<T>::Abs(-afA[j]*afNBpUC[i]
				  + afB[j]*afNApLC[i]
				  - afC[j]*afLBmUA[i]);
    }
    if ( fDdD + fR < fMin || fDdD - fR > fMax )
      return false;
  }

  return true;
}

template <class T>
GbBool
GbIntersection<T>::test (const GbSegment<T>& rkSegment, const GbSphere<T>& rkSphere)
{
  T fSqrDist = GbDistance<T>::squareDistance(rkSphere.center(),rkSegment);
  return fSqrDist <= rkSphere.radius()*rkSphere.radius();
}

template <class T>
GbBool
GbIntersection<T>::test (const GbRay<T>& rkRay, const GbSphere<T>& rkSphere)
{
  T fSqrDist = GbDistance<T>::squareDistance(rkSphere.center(),rkRay);
  return fSqrDist <= rkSphere.radius()*rkSphere.radius();
}

template <class T>
GbBool
GbIntersection<T>::test (const GbLine<T>& rkLine, const GbSphere<T>& rkSphere)
{
  T fSqrDist = GbDistance<T>::squareDistance(rkSphere.center(),rkLine);
  return fSqrDist <= rkSphere.radius()*rkSphere.radius();
}

template <class T>
GbBool
GbIntersection<T>::find (const GbSegment<T>& rkSegment, const GbSphere<T>& rkSphere,
			 int& riQuantity, GbVec3<T> akPoint[2])
{
  // set up quadratic Q(t) = a*t^2 + 2*b*t + c
  GbVec3<T> kDiff = rkSegment.getOrigin() - rkSphere.center();
  T fA = rkSegment.getDirection().getSquareNorm();
  T fB = kDiff.dot(rkSegment.getDirection());
  T fC = kDiff.getSquareNorm() - rkSphere.radius()*rkSphere.radius();

  // no intersection if Q(t) has no real roots
  T afT[2];
  T fDiscr = fB*fB - fA*fC;
  if ( fDiscr < T(0.0) ) {
    riQuantity = 0;
    return false;
  }
  else if ( fDiscr > T(0.0) ) {
    T fRoot = GbMath<T>::Sqrt(fDiscr);
    T fInvA = T(1.0)/fA;
    afT[0] = (-fB - fRoot)*fInvA;
    afT[1] = (-fB + fRoot)*fInvA;

    // assert: t0 < t1 since A > 0

    if ( afT[0] > T(1.0) || afT[1] < T(0.0) ) {
      riQuantity = 0;
      return false;
    }
    else if ( afT[0] >= T(0.0) ) {
      if ( afT[1] > T(1.0) ) {
	riQuantity = 1;
	akPoint[0] = rkSegment.getOrigin()+afT[0]*rkSegment.getDirection();
	return true;
      }
      else {
	riQuantity = 2;
	akPoint[0] = rkSegment.getOrigin()+afT[0]*rkSegment.getDirection();
	akPoint[1] = rkSegment.getOrigin()+afT[1]*rkSegment.getDirection();
	return true;
      }
    }
    else { // afT[1] >= 0
      riQuantity = 1;
      akPoint[0] = rkSegment.getOrigin()+afT[1]*rkSegment.getDirection();
      return true;
    }
  }
  else {
    afT[0] = -fB/fA;
    if ( T(0.0) <= afT[0] && afT[0] <= T(1.0) ) {
      riQuantity = 1;
      akPoint[0] = rkSegment.getOrigin()+afT[0]*rkSegment.getDirection();
      return true;
    }
    else {
      riQuantity = 0;
      return false;
    }
  }
}

template <class T>
GbBool
GbIntersection<T>::find (const GbRay<T>& rkRay, const GbSphere<T>& rkSphere,
			 int& riQuantity, GbVec3<T> akPoint[2])
{
  // set up quadratic Q(t) = a*t^2 + 2*b*t + c
  GbVec3<T> kDiff = rkRay.getOrigin() - rkSphere.center();
  T fA = rkRay.getDirection().getSquareNorm();
  T fB = kDiff.dot(rkRay.getDirection());
  T fC = kDiff.getSquareNorm() - rkSphere.radius()*rkSphere.radius();

  T afT[2];
  T fDiscr = fB*fB - fA*fC;
  if ( fDiscr < T(0.0) ) {
    riQuantity = 0;
    return false;
  }
  else if ( fDiscr > T(0.0) ) {
    T fRoot = GbMath<T>::Sqrt(fDiscr);
    T fInvA = T(1.0)/fA;
    afT[0] = (-fB - fRoot)*fInvA;
    afT[1] = (-fB + fRoot)*fInvA;

    if ( afT[0] >= T(0.0) ) {
      riQuantity = 2;
      akPoint[0] = rkRay.getOrigin() + afT[0]*rkRay.getDirection();
      akPoint[1] = rkRay.getOrigin() + afT[1]*rkRay.getDirection();
      return true;
    }
    else if ( afT[1] >= T(0.0) ) {
      riQuantity = 1;
      akPoint[0] = rkRay.getOrigin() + afT[1]*rkRay.getDirection();
      return true;
    }
    else {
      riQuantity = 0;
      return false;
    }
  }
  else {
    afT[0] = -fB/fA;
    if ( afT[0] >= T(0.0) ) {
      riQuantity = 1;
      akPoint[0] = rkRay.getOrigin() + afT[0]*rkRay.getDirection();
      return true;
    }
    else {
      riQuantity = 0;
      return false;
    }
  }
}

template <class T>
GbBool
GbIntersection<T>::find (const GbLine<T>& rkLine, const GbSphere<T>& rkSphere,
			 int& riQuantity, GbVec3<T> akPoint[2])
{
  // set up quadratic Q(t) = a*t^2 + 2*b*t + c
  GbVec3<T> kDiff = rkLine.getOrigin() - rkSphere.center();
  T fA = rkLine.getDirection().getSquareNorm();
  T fB = kDiff.dot(rkLine.getDirection());
  T fC = kDiff.getSquareNorm() - rkSphere.radius()*rkSphere.radius();

  T afT[2];
  T fDiscr = fB*fB - fA*fC;
  if ( fDiscr < T(0.0) ) {
    riQuantity = 0;
    return false;
  }
  else if ( fDiscr > T(0.0) ) {
    T fRoot = GbMath<T>::Sqrt(fDiscr);
    T fInvA = T(1.0)/fA;
    riQuantity = 2;
    afT[0] = (-fB - fRoot)*fInvA;
    afT[1] = (-fB + fRoot)*fInvA;
    akPoint[0] = rkLine.getOrigin() + afT[0]*rkLine.getDirection();
    akPoint[1] = rkLine.getOrigin() + afT[1]*rkLine.getDirection();
    return true;
  }
  else {
    riQuantity = 1;
    afT[0] = -fB/fA;
    akPoint[0] = rkLine.getOrigin() + afT[0]*rkLine.getDirection();
    return true;
  }
}

template <class T>
GbBool
GbIntersection<T>::test (const GbPlane<T>& rkPlane, const GbSphere<T>& rkSphere,
			 GbBool bUnitNormal)
{
  GbVec3<T> kNormal = rkPlane.getNormal();
  T fConstant = rkPlane.getConstant();
  if ( !bUnitNormal ) {
    T fLength = kNormal.normalize();
    fConstant /= fLength;
  }

  T fPseudoDistance = kNormal.dot(rkSphere.center()) - fConstant;
  return GbMath<T>::Abs(fPseudoDistance) <= rkSphere.radius();
}

template <class T>
GbBool
GbIntersection<T>::test (const GbSphere<T>& rkSphere, const GbFrustum<T>& rkFrustum)
{
  return GbDistance<T>::squareDistance(rkSphere.center(),rkFrustum) <= rkSphere.radius()*rkSphere.radius();
}

/*! Intersection of sphere with sphere.

  Two spheres are |X-C0|^2 = R0^2 and |X-C1|^2 = R1^2.  If the centers are
  identical (C1=C0), then the spheres intersect if and only if R0=R1.  If
  there is a circle of intersection, the plane of intersection must have
  normal N = C1-C0.  Let U and V be unit vectors such that Dot(N,U) = 0,
  Dot(N,V) = 0, Dot(U,V) = 0.  The circle of intersection is given by

  X = C+R*(cos(A)*U+sin(A)*V)

  where C=C0+t*N with 0 <= t <= 1 (t is fixed) and 0 <= A < 2*pi
  (A is variable).

  Plug into first sphere equation to get
  
  R0^2 = t^2*|N|^2 + R^2
  
  Plug into second sphere equation to get
  
  R1^2 = (t-1)^2*|N|^2 + R^2
  
  Subtract these two equations and solve for t to get
  
  t = 0.5*(1+(R0^2-R1^2)/|N|^2)
  
  It is necessary that 0 <= t <= 1 for there to be an intersection since
  the center of the circle of intersection must be on the segment connecting
  C0 and C1.

  Then plug this value of t into the first of the previous equations and
  solve for
  
  R^2 = R0^2 - t^2*|N|^2
  
  As long as the right-hand side is nonnegative, the circle of intersection
  exists and the radius is R = sqrt(R0^2 - t^2*|N|^2).
*/
template <class T>
GbBool
GbIntersection<T>::test (const GbSphere<T>& rkS0, const GbSphere<T>& rkS1)
{
  GbVec3<T> kDiff = rkS1.center() - rkS0.center();
  T fSqrLen = kDiff.getSquareNorm();
  T fRSum = rkS0.radius() + rkS1.radius();
  T fRSumSqr = fRSum*fRSum;

  return fSqrLen <= fRSumSqr;
}

template <class T>
GbBool
GbIntersection<T>::find (const GbSphere<T>& rkS0, const GbSphere<T>& rkS1,
			 GbVec3<T>& rkU, GbVec3<T>& rkV, GbVec3<T>& rkC, T& rfR)
{
  // plane of intersection must have N as its normal
  GbVec3<T> kN = rkS1.center() - rkS0.center();
  T fNSqrLen = kN.getSquareNorm();
  T fRSum = rkS0.radius() + rkS1.radius();
  if ( fNSqrLen > fRSum*fRSum ) {
    // sphere centers are too far apart for intersection
    return false;
  }

  T fR0Sqr = rkS0.radius()*rkS0.radius();
  T fR1Sqr = rkS1.radius()*rkS1.radius();
  T fInvNSqrLen = T(1.0)/fNSqrLen;
  T fT = T(0.5)*(T(1.0)+(fR0Sqr-fR1Sqr)*fInvNSqrLen);
  if ( fT < T(0.0) || fT > T(1.0) )
    return false;

  T fRSqr = fR0Sqr - fT*fT*fNSqrLen;
  if ( fRSqr < T(0.0) )
    return false;

  // center and radius of circle of intersection
  rkC = rkS0.center() + fT*kN;
  rfR = GbMath<T>::Sqrt(fRSqr);

  // compute U and V for plane of circle
  kN *= GbMath<T>::Sqrt(fInvNSqrLen);
  GbVec3<T>::generateOrthonormalBasis(rkU,rkV,kN);

  return true;
}

template <class T>
GbBool
GbIntersection<T>::test (const GoTriangleElement& rkT, const GbSphere<T>& rkS)
{
  return GbDistance<T>::squareDistance(rkS.center(),rkT) < rkS.radius()*rkS.radius();
}

template <class T>
GbBool
GbIntersection<T>::find (const GoTetrahedronElement& rkT0, const GoTetrahedronElement& rkT1,
			 std::vector<GoTetrahedronElement*>& rkIntr)
{
  // build planar faces of T0
  GbPlane<T> akPlane[4];
  rkT0.getPlanes(akPlane);

  // initial object to clip is T1
  rkIntr.clear();
  rkIntr.push_back(const_cast<GoTetrahedronElement*>(&rkT1));

  // clip T1 against planes of T0
  for (int iP = 0; iP < 4; ++iP) {
    std::vector<GoTetrahedronElement*> kInside;
    for (int iT = 0; iT < rkIntr.size(); ++iT)
      splitAndDecompose(rkIntr[iT],akPlane[iP],kInside);
    rkIntr = kInside;
  }

  return !(rkIntr.empty());
}

