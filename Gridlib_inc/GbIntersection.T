// -*- C++ -*- automatisch in C++mode wechseln (emacs)
/*----------------------------------------------------------------------
|
|
| $Log$
| Revision 1.1  2002/02/22 14:47:56  elena
| new: dir Gridlib_inc
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
  GbVec3<T> kSDir = 0.5*rkSegment.getDirection();
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
  if ( fADdU[0] > rkBox.getExtent(0) && fDdU[0]*fWdU[0] >= 0.0 )
    return false;

  fWdU[1] = rkRay.getDirection().dot(rkBox.getAxis(1));
  fAWdU[1] = GbMath<T>::Abs(fWdU[1]);
  fDdU[1] = kDiff.dot(rkBox.getAxis(1));
  fADdU[1] = GbMath<T>::Abs(fDdU[1]);
  if ( fADdU[1] > rkBox.getExtent(1) && fDdU[1]*fWdU[1] >= 0.0 )
    return false;

  fWdU[2] = rkRay.getDirection().dot(rkBox.getAxis(2));
  fAWdU[2] = GbMath<T>::Abs(fWdU[2]);
  fDdU[2] = kDiff.dot(rkBox.getAxis(2));
  fADdU[2] = GbMath<T>::Abs(fDdU[2]);
  if ( fADdU[2] > rkBox.getExtent(2) && fDdU[2]*fWdU[2] >= 0.0 )
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

  T fT0 = 0.0, fT1 = 1.0;
  GbBool bIntersects = GbIntersection<T>::find(kOrigin,kDirection,rkBox.getExtents(),fT0,fT1);

  if ( bIntersects ) {
    if ( fT0 > 0.0 ) {
      if ( fT1 < 1.0 ) {
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
      if ( fT1 < 1.0 ) {
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

  T fT0 = 0.0, fT1 = std::numeric_limits<T>::infinity();
  GbBool bIntersects = GbIntersection<T>::find(kOrigin,kDirection,rkBox.getExtents(),fT0,fT1);

  if ( bIntersects ) {
    if ( fT0 > 0.0 ) {
      if ( fT1 < std::numeric_limits<T>::infinity() ) {
	riQuantity = 2;
	akPoint[0] = rkRay.getOrigin() + fT0*rkRay.getDirection();
	akPoint[1] = rkRay.getOrigin() + fT1*rkRay.getDirection();
      }
      else {
	riQuantity = 1;
	akPoint[0] = rkRay.getOrigin() + fT0*rkRay.getDirection();
      }
    }
    else { // fT0 == 0.0
      if ( fT1 < std::numeric_limits<T>::infinity() ) {
	riQuantity = 1;
	akPoint[0] = rkRay.getOrigin() + fT1*rkRay.getDirection();
      }
      else { // fT1 == INFINITY
	// assert:  should not get here
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

  T fT0 = -std::numeric_limits<T>::infinity(), fT1 = std::numeric_limits<T>::infinity();
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

  T fInvDet = 1.0/fDet;
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
