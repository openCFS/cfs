// -*- C++ -*- automatisch in C++mode wechseln (emacs)
/*----------------------------------------------------------------------
|
|
| $Log$
| Revision 1.3  2002/04/03 11:17:07  elena
| new: changes in CFS++ on SGI and new function for gridlib
|
| Revision 1.6  2001/09/12 09:10:40  prkipfer
| adapted to triangle layout
|
| Revision 1.5  2001/08/16 16:53:14  prkipfer
| improved type safety for template parameter
|
| Revision 1.4  2001/06/19 16:32:26  prkipfer
| fixed typos and Linux compiler target
|
| Revision 1.3  2001/06/18 10:43:03  prkipfer
| introduced Frustum object
|
| Revision 1.2  2001/02/13 10:58:38  prkipfer
| exchanged fixed numeric tolerances by machine dependent ones
|
| Revision 1.1  2001/01/02 15:07:14  prkipfer
| introduced new classes
|
|
|
+---------------------------------------------------------------------*/

template <class T>
const T
GbDistance<T>::TOLERANCE = std::numeric_limits<T>::epsilon(); // max. 1e-05

template <class T>
T 
GbDistance<T>::squareDistance (const GbVec3<T>& rkPoint, const GbBox3<T>& rkBox,
			       T* pfBParam0, T* pfBParam1, T* pfBParam2)
{
#ifdef DEBUG
  // The three parameters pointers are either all non-null or all null.
  if ( pfBParam0 ) {
    assert( pfBParam1 && pfBParam2 );
  }
  else {
    assert( !pfBParam1 && !pfBParam2 );
  }
#endif

  // compute coordinates of point in box coordinate system
  GbVec3<T> kDiff = rkPoint - rkBox.getCenter();
  GbVec3<T> kClosest(kDiff.dot(rkBox.getAxis(0)),
		     kDiff.dot(rkBox.getAxis(1)),
		     kDiff.dot(rkBox.getAxis(2)));

  // project test point onto box
  T fSqrDistance = T(0.0);
  T fDelta;

  if ( kClosest[0] < -rkBox.getExtent(0) ) {
    fDelta = kClosest[0] + rkBox.getExtent(0);
    fSqrDistance += fDelta*fDelta;
    kClosest[0] = -rkBox.getExtent(0);
  }
  else if ( kClosest[0] > rkBox.getExtent(0) ) {
    fDelta = kClosest[0] - rkBox.getExtent(0);
    fSqrDistance += fDelta*fDelta;
    kClosest[0] = rkBox.getExtent(0);
  }

  if ( kClosest[1] < -rkBox.getExtent(1) ) {
    fDelta = kClosest[1] + rkBox.getExtent(1);
    fSqrDistance += fDelta*fDelta;
    kClosest[1] = -rkBox.getExtent(1);
  }
  else if ( kClosest[1] > rkBox.getExtent(1) ) {
    fDelta = kClosest[1] - rkBox.getExtent(1);
    fSqrDistance += fDelta*fDelta;
    kClosest[1] = rkBox.getExtent(1);
  }

  if ( kClosest[2] < -rkBox.getExtent(2) ) {
    fDelta = kClosest[2] + rkBox.getExtent(2);
    fSqrDistance += fDelta*fDelta;
    kClosest[2] = -rkBox.getExtent(2);
  }
  else if ( kClosest[2] > rkBox.getExtent(2) ) {
    fDelta = kClosest[2] - rkBox.getExtent(2);
    fSqrDistance += fDelta*fDelta;
    kClosest[2] = rkBox.getExtent(2);
  }

  if ( pfBParam0 ) {
    *pfBParam0 = kClosest[0];
    *pfBParam1 = kClosest[1];
    *pfBParam2 = kClosest[2];
  }

  return fSqrDistance;
}

template <class T>
T 
GbDistance<T>::distance (const GbVec3<T>& rkPoint, const GbBox3<T>& rkBox,
			 T* pfBParam0, T* pfBParam1, T* pfBParam2)
{
  return GbMath<T>::Sqrt(GbDistance<T>::squareDistance(rkPoint,rkBox,pfBParam0,pfBParam1,pfBParam2));
}

template <class T>
T 
GbDistance<T>::squareDistance (const GbVec3<T>& rkPoint, const GoTriangleElement& rkTri,
			       T* pfSParam, T* pfTParam)
{
  GbVec3<T> kDiff = rkTri.getOrigin() - rkPoint;
  T fA00 = rkTri.getEdgeVec(0).getSquareNorm();
  T fA01 = rkTri.getEdgeVec(0).dot(rkTri.getEdgeVec(1));
  T fA11 = rkTri.getEdgeVec(1).getSquareNorm();
  T fB0 = kDiff.dot(rkTri.getEdgeVec(0));
  T fB1 = kDiff.dot(rkTri.getEdgeVec(1));
  T fC = kDiff.getSquareNorm();
  T fDet = GbMath<T>::Abs(fA00*fA11-fA01*fA01);
  T fS = fA01*fB1-fA11*fB0;
  T fT = fA01*fB0-fA00*fB1;
  T fSqrDist;

  if ( fS + fT <= fDet ) {
    if ( fS < T(0.0) ) {
      if ( fT < T(0.0) ) {  // region 4
	if ( fB0 < T(0.0) ) {
	  fT = T(0.0);
	  if ( -fB0 >= fA00 ) {
	    fS = T(1.0);
	    fSqrDist = fA00+T(2.0)*fB0+fC;
	  }
	  else {
	    fS = -fB0/fA00;
	    fSqrDist = fB0*fS+fC;
	  }
	}
	else {
	  fS = T(0.0);
	  if ( fB1 >= T(0.0) ) {
	    fT = T(0.0);
	    fSqrDist = fC;
	  }
	  else if ( -fB1 >= fA11 ) {
	    fT = T(1.0);
	    fSqrDist = fA11+T(2.0)*fB1+fC;
	  }
	  else {
	    fT = -fB1/fA11;
	    fSqrDist = fB1*fT+fC;
	  }
	}
      }
      else { // region 3
	fS = T(0.0);
	if ( fB1 >= T(0.0) ) {
	  fT = T(0.0);
	  fSqrDist = fC;
	}
	else if ( -fB1 >= fA11 ) {
	  fT = 1;
	  fSqrDist = fA11+T(2.0)*fB1+fC;
	}
	else {
	  fT = -fB1/fA11;
	  fSqrDist = fB1*fT+fC;
	}
      }
    }
    else if ( fT < T(0.0) ) {  // region 5
      fT = T(0.0);
      if ( fB0 >= T(0.0) ) {
	fS = T(0.0);
	fSqrDist = fC;
      }
      else if ( -fB0 >= fA00 ) {
	fS = T(1.0);
	fSqrDist = fA00+T(2.0)*fB0+fC;
      }
      else {
	fS = -fB0/fA00;
	fSqrDist = fB0*fS+fC;
      }
    }
    else {  // region 0
      // minimum at interior point
      T fInvDet = T(1.0)/fDet;
      fS *= fInvDet;
      fT *= fInvDet;
      fSqrDist = fS*(fA00*fS+fA01*fT+T(2.0)*fB0) +
	         fT*(fA01*fS+fA11*fT+T(2.0)*fB1) +
	         fC;
    }
  }
  else {
    T fTmp0, fTmp1, fNumer, fDenom;

    if ( fS < T(0.0) ) {  // region 2
      fTmp0 = fA01 + fB0;
      fTmp1 = fA11 + fB1;
      if ( fTmp1 > fTmp0 ) {
	fNumer = fTmp1 - fTmp0;
	fDenom = fA00-T(2.0)*fA01+fA11;
	if ( fNumer >= fDenom ) {
	  fS = T(1.0);
	  fT = T(0.0);
	  fSqrDist = fA00+T(2.0)*fB0+fC;
	}
	else {
	  fS = fNumer/fDenom;
	  fT = T(1.0) - fS;
	  fSqrDist = fS*(fA00*fS+fA01*fT+T(2.0)*fB0) +
	             fT*(fA01*fS+fA11*fT+T(2.0)*fB1) +
	             fC;
	}
      }
      else {
	fS = T(0.0);
	if ( fTmp1 <= T(0.0) ) {
	  fT = T(1.0);
	  fSqrDist = fA11+T(2.0)*fB1+fC;
	}
	else if ( fB1 >= T(0.0) ) {
	  fT = T(0.0);
	  fSqrDist = fC;
	}
	else {
	  fT = -fB1/fA11;
	  fSqrDist = fB1*fT+fC;
	}
      }
    }
    else if ( fT < T(0.0) ) {  // region 6
      fTmp0 = fA01 + fB1;
      fTmp1 = fA00 + fB0;
      if ( fTmp1 > fTmp0 ) {
	fNumer = fTmp1 - fTmp0;
	fDenom = fA00-T(2.0)*fA01+fA11;
	if ( fNumer >= fDenom ) {
	  fT = T(1.0);
	  fS = T(0.0);
	  fSqrDist = fA11+T(2.0)*fB1+fC;
	}
	else {
	  fT = fNumer/fDenom;
	  fS = T(1.0) - fT;
	  fSqrDist = fS*(fA00*fS+fA01*fT+T(2.0)*fB0) +
	             fT*(fA01*fS+fA11*fT+T(2.0)*fB1) +
	             fC;
	}
      }
      else {
	fT = T(0.0);
	if ( fTmp1 <= T(0.0) ) {
	  fS = T(1.0);
	  fSqrDist = fA00+T(2.0)*fB0+fC;
	}
	else if ( fB0 >= T(0.0) ) {
	  fS = T(0.0);
	  fSqrDist = fC;
	}
	else {
	  fS = -fB0/fA00;
	  fSqrDist = fB0*fS+fC;
	}
      }
    }
    else {  // region 1
      fNumer = fA11 + fB1 - fA01 - fB0;
      if ( fNumer <= T(0.0) ) {
	fS = T(0.0);
	fT = T(1.0);
	fSqrDist = fA11+T(2.0)*fB1+fC;
      }
      else {
	fDenom = fA00-T(2.0)*fA01+fA11;
	if ( fNumer >= fDenom ) {
	  fS = T(1.0);
	  fT = T(0.0);
	  fSqrDist = fA00+T(2.0)*fB0+fC;
	}
	else {
	  fS = fNumer/fDenom;
	  fT = T(1.0) - fS;
	  fSqrDist = fS*(fA00*fS+fA01*fT+T(2.0)*fB0) +
	             fT*(fA01*fS+fA11*fT+T(2.0)*fB1) +
	             fC;
	}
      }
    }
  }

  if ( pfSParam ) *pfSParam = fS;

  if ( pfTParam ) *pfTParam = fT;
  
  return GbMath<T>::Abs(fSqrDist);
}

template <class T>
T 
GbDistance<T>::distance (const GbVec3<T>& rkPoint, const GoTriangleElement& rkTri,
			 T* pfSParam, T* pfTParam)
{
  return GbMath<T>::Sqrt(GbDistance<T>::squareDistance(rkPoint,rkTri,pfSParam,pfTParam));
}

template <class T>
T 
GbDistance<T>::squareDistance (const GoTriangleElement& rkTri0, const GoTriangleElement& rkTri1, 
			       T* pfTri0P0, T* pfTri0P1, T* pfTri1P0, T* pfTri1P1)
{
  T fS, fT, fU, fV, fS0, fT0, fU0, fV0, fSqrDist, fSqrDist0;
  GbSegment<T> kSeg;

  // compare edges of tri0 against all of tri1
  kSeg.setOrigin(rkTri0.getOrigin());
  kSeg.setDirection(rkTri0.getEdgeVec(0));
  fSqrDist = GbDistance<T>::squareDistance(kSeg,rkTri1,&fS,&fU,&fV);
  fT = T(0.0);

  kSeg.setDirection(rkTri0.getEdgeVec(1));
  fSqrDist0 = GbDistance<T>::squareDistance(kSeg,rkTri1,&fT0,&fU0,&fV0);
  fS0 = T(0.0);
  if ( fSqrDist0 < fSqrDist ) {
    fSqrDist = fSqrDist0;
    fS = fS0;
    fT = fT0;
    fU = fU0;
    fV = fV0;
  }

  kSeg.setOrigin(kSeg.getOrigin() + rkTri0.getEdgeVec(0));
  kSeg.setDirection(kSeg.getDirection() - rkTri0.getEdgeVec(0));
  fSqrDist0 = GbDistance<T>::squareDistance(kSeg,rkTri1,&fT0,&fU0,&fV0);
  fS0 = T(1.0)-fT0;
  if ( fSqrDist0 < fSqrDist ) {
    fSqrDist = fSqrDist0;
    fS = fS0;
    fT = fT0;
    fU = fU0;
    fV = fV0;
  }

  // compare edges of tri1 against all of tri0
  kSeg.setOrigin(rkTri1.getOrigin());
  kSeg.setDirection(rkTri1.getEdgeVec(0));
  fSqrDist0 = GbDistance<T>::squareDistance(kSeg,rkTri0,&fU0,&fS0,&fT0);
  fV0 = T(0.0);
  if ( fSqrDist0 < fSqrDist ) {
    fSqrDist = fSqrDist0;
    fS = fS0;
    fT = fT0;
    fU = fU0;
    fV = fV0;
  }

  kSeg.setDirection(rkTri1.getEdgeVec(1));
  fSqrDist0 = GbDistance<T>::squareDistance(kSeg,rkTri0,&fV0,&fS0,&fT0);
  fU0 = T(0.0);
  if ( fSqrDist0 < fSqrDist ) {
    fSqrDist = fSqrDist0;
    fS = fS0;
    fT = fT0;
    fU = fU0;
    fV = fV0;
  }

  kSeg.setOrigin(kSeg.getOrigin() + rkTri1.getEdgeVec(0));
  kSeg.setDirection(kSeg.getDirection() - rkTri1.getEdgeVec(0));
  fSqrDist0 = GbDistance<T>::squareDistance(kSeg,rkTri0,&fV0,&fS0,&fT0);
  fU0 = T(1.0)-fV0;
  if ( fSqrDist0 < fSqrDist ) {
    fSqrDist = fSqrDist0;
    fS = fS0;
    fT = fT0;
    fU = fU0;
    fV = fV0;
  }

  if ( pfTri0P0 ) *pfTri0P0 = fS;

  if ( pfTri0P1 ) *pfTri0P1 = fT;

  if ( pfTri1P0 ) *pfTri1P0 = fU;

  if ( pfTri1P1 ) *pfTri1P1 = fV;

  return GbMath<T>::Abs(fSqrDist);
}

template <class T>
T 
GbDistance<T>::distance (const GoTriangleElement& rkTri0, const GoTriangleElement& rkTri1,
			 T* pfTri0P0, T* pfTri0P1, T* pfTri1P0, T* pfTri1P1)
{
  return GbMath<T>::Sqrt(GbDistance<T>::squareDistance(rkTri0,rkTri1,pfTri0P0,pfTri0P1,pfTri1P0,pfTri1P1));
}

template <class T>
T 
GbDistance<T>::squareDistance (const GbSegment<T>& rkSeg, const GoTriangleElement& rkTri,
			       T* pfSegP, T* pfTriP0, T* pfTriP1)
{
  GbVec3<T> kDiff = rkTri.getOrigin() - rkSeg.getOrigin();
  T fA00 = rkSeg.getDirection().getSquareNorm();
  T fA01 = -rkSeg.getDirection().dot(rkTri.getEdgeVec(0));
  T fA02 = -rkSeg.getDirection().dot(rkTri.getEdgeVec(1));
  T fA11 = rkTri.getEdgeVec(0).getSquareNorm();
  T fA12 = rkTri.getEdgeVec(0).dot(rkTri.getEdgeVec(1));
  T fA22 = rkTri.getEdgeVec(1).dot(rkTri.getEdgeVec(1));
  T fB0  = -kDiff.dot(rkSeg.getDirection());
  T fB1  = kDiff.dot(rkTri.getEdgeVec(0));
  T fB2  = kDiff.dot(rkTri.getEdgeVec(1));
  T fCof00 = fA11*fA22-fA12*fA12;
  T fCof01 = fA02*fA12-fA01*fA22;
  T fCof02 = fA01*fA12-fA02*fA11;
  T fDet = fA00*fCof00+fA01*fCof01+fA02*fCof02;

  GbSegment<T> kTriSeg;
  GbVec3<T> kPt;
  T fSqrDist, fSqrDist0, fR, fS, fT, fR0, fS0, fT0;

  if ( GbMath<T>::Abs(fDet) >= GbDistance<T>::TOLERANCE ) {
    T fCof11 = fA00*fA22-fA02*fA02;
    T fCof12 = fA02*fA01-fA00*fA12;
    T fCof22 = fA00*fA11-fA01*fA01;
    T fInvDet = T(1.0)/fDet;
    T fRhs0 = -fB0*fInvDet;
    T fRhs1 = -fB1*fInvDet;
    T fRhs2 = -fB2*fInvDet;

    fR = fCof00*fRhs0+fCof01*fRhs1+fCof02*fRhs2;
    fS = fCof01*fRhs0+fCof11*fRhs1+fCof12*fRhs2;
    fT = fCof02*fRhs0+fCof12*fRhs1+fCof22*fRhs2;

    if ( fR < T(0.0) ) {
      if ( fS+fT <= T(1.0) ) {
	if ( fS < T(0.0) ) {
	  if ( fT < T(0.0) ) {  // region 4m
	    // min on face s=0 or t=0 or r=0
	    kTriSeg.setOrigin(rkTri.getOrigin());
	    kTriSeg.setDirection(rkTri.getEdgeVec(1));
	    fSqrDist = GbDistance<T>::squareDistance(rkSeg,kTriSeg,&fR,&fT);
	    fS = T(0.0);
	    kTriSeg.setOrigin(rkTri.getOrigin());
	    kTriSeg.setDirection(rkTri.getEdgeVec(0));
	    fSqrDist0 = GbDistance<T>::squareDistance(rkSeg,kTriSeg,&fR0,&fS0);
	    fT0 = T(0.0);
	    if ( fSqrDist0 < fSqrDist ) {
	      fSqrDist = fSqrDist0;
	      fR = fR0;
	      fS = fS0;
	      fT = fT0;
	    }
	    fSqrDist0 = GbDistance<T>::squareDistance(rkSeg.getOrigin(),rkTri,&fS0,&fT0);
	    fR0 = T(0.0);
	    if ( fSqrDist0 < fSqrDist ) {
	      fSqrDist = fSqrDist0;
	      fR = fR0;
	      fS = fS0;
	      fT = fT0;
	    }
	  }
	  else  { // region 3m
	    // min on face s=0 or r=0
	    kTriSeg.setOrigin(rkTri.getOrigin());
	    kTriSeg.setDirection(rkTri.getEdgeVec(1));
	    fSqrDist = GbDistance<T>::squareDistance(rkSeg,kTriSeg,&fR,&fT);
	    fS = T(0.0);
	    fSqrDist0 = GbDistance<T>::squareDistance(rkSeg.getOrigin(),rkTri,&fS0,&fT0);
	    fR0 = T(0.0);
	    if ( fSqrDist0 < fSqrDist ) {
	      fSqrDist = fSqrDist0;
	      fR = fR0;
	      fS = fS0;
	      fT = fT0;
	    }
	  }
	}
	else if ( fT < T(0.0) ) {  // region 5m
	  // min on face t=0 or r=0
	  kTriSeg.setOrigin(rkTri.getOrigin());
	  kTriSeg.setDirection(rkTri.getEdgeVec(0));
	  fSqrDist = GbDistance<T>::squareDistance(rkSeg,kTriSeg,&fR,&fS);
	  fT = T(0.0);
	  fSqrDist0 = GbDistance<T>::squareDistance(rkSeg.getOrigin(),rkTri,&fS0,&fT0);
	  fR0 = T(0.0);
	  if ( fSqrDist0 < fSqrDist ) {
	    fSqrDist = fSqrDist0;
	    fR = fR0;
	    fS = fS0;
	    fT = fT0;
	  }
	}
	else {  // region 0m
	  // min on face r=0
	  fSqrDist = GbDistance<T>::squareDistance(rkSeg.getOrigin(),rkTri,&fS,&fT);
	  fR = T(0.0);
	}
      }
      else {
	if ( fS < T(0.0) ) {  // region 2m
	  // min on face s=0 or s+t=1 or r=0
	  kTriSeg.setOrigin(rkTri.getOrigin());
	  kTriSeg.setDirection(rkTri.getEdgeVec(1));
	  fSqrDist = GbDistance<T>::squareDistance(rkSeg,kTriSeg,&fR,&fT);
	  fS = T(0.0);
	  kTriSeg.setOrigin(rkTri.getOrigin()+rkTri.getEdgeVec(0));
	  kTriSeg.setDirection(rkTri.getEdgeVec(2));
	  fSqrDist0 = GbDistance<T>::squareDistance(rkSeg,kTriSeg,&fR0,&fT0);
	  fS0 = T(1.0)-fT0;
	  if ( fSqrDist0 < fSqrDist ) {
	    fSqrDist = fSqrDist0;
	    fR = fR0;
	    fS = fS0;
	    fT = fT0;
	  }
	  fSqrDist0 = GbDistance<T>::squareDistance(rkSeg.getOrigin(),rkTri,&fS0,&fT0);
	  fR0 = T(0.0);
	  if ( fSqrDist0 < fSqrDist ) {
	    fSqrDist = fSqrDist0;
	    fR = fR0;
	    fS = fS0;
	    fT = fT0;
	  }
	}
	else if ( fT < T(0.0) ) {  // region 6m
	  // min on face t=0 or s+t=1 or r=0
	  kTriSeg.setOrigin(rkTri.getOrigin());
	  kTriSeg.setDirection(rkTri.getEdgeVec(0));
	  fSqrDist = GbDistance<T>::squareDistance(rkSeg,kTriSeg,&fR,&fS);
	  fT = T(0.0);
	  kTriSeg.setOrigin(rkTri.getOrigin()+rkTri.getEdgeVec(0));
	  kTriSeg.setDirection(rkTri.getEdgeVec(2));
	  fSqrDist0 = GbDistance<T>::squareDistance(rkSeg,kTriSeg,&fR0,&fT0);
	  fS0 = T(1.0)-fT0;
	  if ( fSqrDist0 < fSqrDist ) {
	    fSqrDist = fSqrDist0;
	    fR = fR0;
	    fS = fS0;
	    fT = fT0;
	  }
	  fSqrDist0 = GbDistance<T>::squareDistance(rkSeg.getOrigin(),rkTri,&fS0,&fT0);
	  fR0 = T(0.0);
	  if ( fSqrDist0 < fSqrDist ) {
	    fSqrDist = fSqrDist0;
	    fR = fR0;
	    fS = fS0;
	    fT = fT0;
	  }
	}
	else { // region 1m
	  // min on face s+t=1 or r=0
	  kTriSeg.setOrigin(rkTri.getOrigin()+rkTri.getEdgeVec(0));
	  kTriSeg.setDirection(rkTri.getEdgeVec(2));
	  fSqrDist = GbDistance<T>::squareDistance(rkSeg,kTriSeg,&fR,&fT);
	  fS = T(1.0)-fT;
	  fSqrDist0 = GbDistance<T>::squareDistance(rkSeg.getOrigin(),rkTri,&fS0,&fT0);
	  fR0 = T(0.0);
	  if ( fSqrDist0 < fSqrDist ) {
	    fSqrDist = fSqrDist0;
	    fR = fR0;
	    fS = fS0;
	    fT = fT0;
	  }
	}
      }
    }
    else if ( fR <= T(1.0) ) {
      if ( fS+fT <= T(1.0) ) {
	if ( fS < T(0.0) ) {
	  if ( fT < T(0.0) ) {  // region 4
	    // min on face s=0 or t=0
	    kTriSeg.setOrigin(rkTri.getOrigin());
	    kTriSeg.setDirection(rkTri.getEdgeVec(1));
	    fSqrDist = GbDistance<T>::squareDistance(rkSeg,kTriSeg,&fR,&fT);
	    fS = T(0.0);
	    kTriSeg.setOrigin(rkTri.getOrigin());
	    kTriSeg.setDirection(rkTri.getEdgeVec(0));
	    fSqrDist0 = GbDistance<T>::squareDistance(rkSeg,kTriSeg,&fR0,&fS0);
	    fT0 = T(0.0);
	    if ( fSqrDist0 < fSqrDist ) {
	      fSqrDist = fSqrDist0;
	      fR = fR0;
	      fS = fS0;
	      fT = fT0;
	    }
	  }
	  else {  // region 3
	    // min on face s=0
	    kTriSeg.setOrigin(rkTri.getOrigin());
	    kTriSeg.setDirection(rkTri.getEdgeVec(1));
	    fSqrDist = GbDistance<T>::squareDistance(rkSeg,kTriSeg,&fR,&fT);
	    fS = T(0.0);
	  }
	}
	else if ( fT < T(0.0) ) {  // region 5
	  // min on face t=0
	  kTriSeg.setOrigin(rkTri.getOrigin());
	  kTriSeg.setDirection(rkTri.getEdgeVec(0));
	  fSqrDist = GbDistance<T>::squareDistance(rkSeg,kTriSeg,&fR,&fS);
	  fT = T(0.0);
	}
	else {  // region 0
	  // global minimum is interior, done
	  fSqrDist = fR*(fA00*fR+fA01*fS+fA02*fT+T(2.0)*fB0) +
	             fS*(fA01*fR+fA11*fS+fA12*fT+T(2.0)*fB1) +
	             fT*(fA02*fR+fA12*fS+fA22*fT+T(2.0)*fB2) +
	             kDiff.getSquareNorm();
	}
      }
      else {
	if ( fS < T(0.0) ) {  // region 2
	  // min on face s=0 or s+t=1
	  kTriSeg.setOrigin(rkTri.getOrigin());
	  kTriSeg.setDirection(rkTri.getEdgeVec(1));
	  fSqrDist = GbDistance<T>::squareDistance(rkSeg,kTriSeg,&fR,&fT);
	  fS = T(0.0);
	  kTriSeg.setOrigin(rkTri.getOrigin()+rkTri.getEdgeVec(0));
	  kTriSeg.setDirection(rkTri.getEdgeVec(2));
	  fSqrDist0 = GbDistance<T>::squareDistance(rkSeg,kTriSeg,&fR0,&fT0);
	  fS0 = T(1.0)-fT0;
	  if ( fSqrDist0 < fSqrDist ) {
	    fSqrDist = fSqrDist0;
	    fR = fR0;
	    fS = fS0;
	    fT = fT0;
	  }
	}
	else if ( fT < T(0.0) ) {  // region 6
	  // min on face t=0 or s+t=1
	  kTriSeg.setOrigin(rkTri.getOrigin());
	  kTriSeg.setDirection(rkTri.getEdgeVec(0));
	  fSqrDist = GbDistance<T>::squareDistance(rkSeg,kTriSeg,&fR,&fS);
	  fT = T(0.0);
	  kTriSeg.setOrigin(rkTri.getOrigin()+rkTri.getEdgeVec(0));
	  kTriSeg.setDirection(rkTri.getEdgeVec(2));
	  fSqrDist0 = GbDistance<T>::squareDistance(rkSeg,kTriSeg,&fR0,&fT0);
	  fS0 = T(1.0)-fT0;
	  if ( fSqrDist0 < fSqrDist ) {
	    fSqrDist = fSqrDist0;
	    fR = fR0;
	    fS = fS0;
	    fT = fT0;
	  }
	}
	else {  // region 1
	  // min on face s+t=1
	  kTriSeg.setOrigin(rkTri.getOrigin()+rkTri.getEdgeVec(0));
	  kTriSeg.setDirection(rkTri.getEdgeVec(2));
	  fSqrDist = GbDistance<T>::squareDistance(rkSeg,kTriSeg,&fR,&fT);
	  fS = T(1.0)-fT;
	}
      }
    }
    else {  // fR > 1
      if ( fS+fT <= T(1.0) ) {
	if ( fS < T(0.0) ) {
	  if ( fT < T(0.0) ) {  // region 4p
	    // min on face s=0 or t=0 or r=1
	    kTriSeg.setOrigin(rkTri.getOrigin());
	    kTriSeg.setDirection(rkTri.getEdgeVec(1));
	    fSqrDist = GbDistance<T>::squareDistance(rkSeg,kTriSeg,&fR,&fT);
	    fS = T(0.0);
	    kTriSeg.setOrigin(rkTri.getOrigin());
	    kTriSeg.setDirection(rkTri.getEdgeVec(0));
	    fSqrDist0 = GbDistance<T>::squareDistance(rkSeg,kTriSeg,&fR0,&fS0);
	    fT0 = T(0.0);
	    if ( fSqrDist0 < fSqrDist ) {
	      fSqrDist = fSqrDist0;
	      fR = fR0;
	      fS = fS0;
	      fT = fT0;
	    }
	    kPt = rkSeg.getOrigin()+rkSeg.getDirection();
	    fSqrDist0 = GbDistance<T>::squareDistance(kPt,rkTri,&fS0,&fT0);
	    fR0 = T(1.0);
	    if ( fSqrDist0 < fSqrDist ) {
	      fSqrDist = fSqrDist0;
	      fR = fR0;
	      fS = fS0;
	      fT = fT0;
	    }
	  }
	  else {  // region 3p
	    // min on face s=0 or r=1
	    kTriSeg.setOrigin(rkTri.getOrigin());
	    kTriSeg.setDirection(rkTri.getEdgeVec(1));
	    fSqrDist = GbDistance<T>::squareDistance(rkSeg,kTriSeg,&fR,&fT);
	    fS = T(0.0);
	    kPt = rkSeg.getOrigin()+rkSeg.getDirection();
	    fSqrDist0 = GbDistance<T>::squareDistance(kPt,rkTri,&fS0,&fT0);
	    fR0 = T(1.0);
	    if ( fSqrDist0 < fSqrDist ) {
	      fSqrDist = fSqrDist0;
	      fR = fR0;
	      fS = fS0;
	      fT = fT0;
	    }
	  }
	}
	else if ( fT < T(0.0) ) {  // region 5p
	  // min on face t=0 or r=1
	  kTriSeg.setOrigin(rkTri.getOrigin());
	  kTriSeg.setDirection(rkTri.getEdgeVec(0));
	  fSqrDist = GbDistance<T>::squareDistance(rkSeg,kTriSeg,&fR,&fS);
	  fT = T(0.0);
	  kPt = rkSeg.getOrigin()+rkSeg.getDirection();
	  fSqrDist0 = GbDistance<T>::squareDistance(kPt,rkTri,&fS0,&fT0);
	  fR0 = T(1.0);
	  if ( fSqrDist0 < fSqrDist ) {
	    fSqrDist = fSqrDist0;
	    fR = fR0;
	    fS = fS0;
	    fT = fT0;
	  }
	}
	else {  // region 0p
	  // min face on r=1
	  kPt = rkSeg.getOrigin()+rkSeg.getDirection();
	  fSqrDist = GbDistance<T>::squareDistance(kPt,rkTri,&fS,&fT);
	  fR = T(1.0);
	}
      }
      else {
	if ( fS < T(0.0) ) {  // region 2p
	  // min on face s=0 or s+t=1 or r=1
	  kTriSeg.setOrigin(rkTri.getOrigin());
	  kTriSeg.setDirection(rkTri.getEdgeVec(1));
	  fSqrDist = GbDistance<T>::squareDistance(rkSeg,kTriSeg,&fR,&fT);
	  fS = T(0.0);
	  kTriSeg.setOrigin(rkTri.getOrigin()+rkTri.getEdgeVec(0));
	  kTriSeg.setDirection(rkTri.getEdgeVec(2));
	  fSqrDist0 = GbDistance<T>::squareDistance(rkSeg,kTriSeg,&fR0,&fT0);
	  fS0 = T(1.0)-fT0;
	  if ( fSqrDist0 < fSqrDist ) {
	    fSqrDist = fSqrDist0;
	    fR = fR0;
	    fS = fS0;
	    fT = fT0;
	  }
	  kPt = rkSeg.getOrigin()+rkSeg.getDirection();
	  fSqrDist0 = GbDistance<T>::squareDistance(kPt,rkTri,&fS0,&fT0);
	  fR0 = T(1.0);
	  if ( fSqrDist0 < fSqrDist ) {
	    fSqrDist = fSqrDist0;
	    fR = fR0;
	    fS = fS0;
	    fT = fT0;
	  }
	}
	else if ( fT < T(0.0) ) {  // region 6p
	  // min on face t=0 or s+t=1 or r=1
	  kTriSeg.setOrigin(rkTri.getOrigin());
	  kTriSeg.setDirection(rkTri.getEdgeVec(0));
	  fSqrDist = GbDistance<T>::squareDistance(rkSeg,kTriSeg,&fR,&fS);
	  fT = T(0.0);
	  kTriSeg.setOrigin(rkTri.getOrigin()+rkTri.getEdgeVec(0));
	  kTriSeg.setDirection(rkTri.getEdgeVec(2));
	  fSqrDist0 = GbDistance<T>::squareDistance(rkSeg,kTriSeg,&fR0,&fT0);
	  fS0 = T(1.0)-fT0;
	  if ( fSqrDist0 < fSqrDist ) {
	    fSqrDist = fSqrDist0;
	    fR = fR0;
	    fS = fS0;
	    fT = fT0;
	  }
	  kPt = rkSeg.getOrigin()+rkSeg.getDirection();
	  fSqrDist0 = GbDistance<T>::squareDistance(kPt,rkTri,&fS0,&fT0);
	  fR0 = T(1.0);
	  if ( fSqrDist0 < fSqrDist ) {
	    fSqrDist = fSqrDist0;
	    fR = fR0;
	    fS = fS0;
	    fT = fT0;
	  }
	}
	else { // region 1p
	  // min on face s+t=1 or r=1
	  kTriSeg.setOrigin(rkTri.getOrigin()+rkTri.getEdgeVec(0));
	  kTriSeg.setDirection(rkTri.getEdgeVec(2));
	  fSqrDist = GbDistance<T>::squareDistance(rkSeg,kTriSeg,&fR,&fT);
	  fS = T(1.0)-fT;
	  kPt = rkSeg.getOrigin()+rkSeg.getDirection();
	  fSqrDist0 = GbDistance<T>::squareDistance(kPt,rkTri,&fS0,&fT0);
	  fR0 = T(1.0);
	  if ( fSqrDist0 < fSqrDist ) {
	    fSqrDist = fSqrDist0;
	    fR = fR0;
	    fS = fS0;
	    fT = fT0;
	  }
	}
      }
    }
  }
  else {
    // segment and triangle are parallel
    kTriSeg.setOrigin(rkTri.getOrigin());
    kTriSeg.setDirection(rkTri.getEdgeVec(0));
    fSqrDist = GbDistance<T>::squareDistance(rkSeg,kTriSeg,&fR,&fS);
    fT = T(0.0);

    kTriSeg.setDirection(rkTri.getEdgeVec(1));
    fSqrDist0 = GbDistance<T>::squareDistance(rkSeg,kTriSeg,&fR0,&fT0);
    fS0 = T(0.0);
    if ( fSqrDist0 < fSqrDist ) {
      fSqrDist = fSqrDist0;
      fR = fR0;
      fS = fS0;
      fT = fT0;
    }
    
    kTriSeg.setOrigin(rkTri.getOrigin() + rkTri.getEdgeVec(0));
    kTriSeg.setDirection(rkTri.getEdgeVec(2));
    fSqrDist0 = GbDistance<T>::squareDistance(rkSeg,kTriSeg,&fR0,&fT0);
    fS0 = T(1.0)-fT0;
    if ( fSqrDist0 < fSqrDist ) {
      fSqrDist = fSqrDist0;
      fR = fR0;
      fS = fS0;
      fT = fT0;
    }

    fSqrDist0 = GbDistance<T>::squareDistance(rkSeg.getOrigin(),rkTri,&fS0,&fT0);
    fR0 = T(0.0);
    if ( fSqrDist0 < fSqrDist ) {
      fSqrDist = fSqrDist0;
      fR = fR0;
      fS = fS0;
      fT = fT0;
    }

    kPt = rkSeg.getOrigin()+rkSeg.getDirection();
    fSqrDist0 = GbDistance<T>::squareDistance(kPt,rkTri,&fS0,&fT0);
    fR0 = T(1.0);
    if ( fSqrDist0 < fSqrDist ) {
      fSqrDist = fSqrDist0;
      fR = fR0;
      fS = fS0;
      fT = fT0;
    }
  }

  if ( pfSegP ) *pfSegP = fR;

  if ( pfTriP0 ) *pfTriP0 = fS;

  if ( pfTriP1 ) *pfTriP1 = fT;

  return GbMath<T>::Abs(fSqrDist);
}

template <class T>
T 
GbDistance<T>::distance (const GbSegment<T>& rkSeg, const GoTriangleElement& rkTri,
			 T* pfSegP, T* pfTriP0, T* pfTriP1)
{
  return GbMath<T>::Sqrt(GbDistance<T>::squareDistance(rkSeg,rkTri,pfSegP,pfTriP0,pfTriP1));
}

template <class T>
T 
GbDistance<T>::squareDistance (const GbSegment<T>& rkSeg0, const GbSegment<T>& rkSeg1,
			       T* pfSegP0, T* pfSegP1)
{
  GbVec3<T> kDiff = rkSeg0.getOrigin() - rkSeg1.getOrigin();
  T fA00 = rkSeg0.getDirection().getSquareNorm();
  T fA01 = -rkSeg0.getDirection().dot(rkSeg1.getDirection());
  T fA11 = rkSeg1.getDirection().getSquareNorm();
  T fB0 = kDiff.dot(rkSeg0.getDirection());
  T fC = kDiff.getSquareNorm();
  T fDet = GbMath<T>::Abs(fA00*fA11-fA01*fA01);
  T fB1, fS, fT, fSqrDist, fTmp;

  if ( fDet >= GbDistance<T>::TOLERANCE ) {
    // line segments are not parallel
    fB1 = -kDiff.dot(rkSeg1.getDirection());
    fS = fA01*fB1-fA11*fB0;
    fT = fA01*fB0-fA00*fB1;
        
    if ( fS >= T(0.0) ) {
      if ( fS <= fDet ) {
	if ( fT >= T(0.0) ) {
	  if ( fT <= fDet ) { // region 0 (interior)
	    // minimum at two interior points of 3D lines
	    T fInvDet = T(1.0)/fDet;
	    fS *= fInvDet;
	    fT *= fInvDet;
	    fSqrDist = fS*(fA00*fS+fA01*fT+T(2.0)*fB0) +
	               fT*(fA01*fS+fA11*fT+T(2.0)*fB1) +
	               fC;
	  }
	  else { // region 3 (side)
	    fT = T(1.0);
	    fTmp = fA01+fB0;
	    if ( fTmp >= T(0.0) ) {
	      fS = T(0.0);
	      fSqrDist = fA11+T(2.0)*fB1+fC;
	    }
	    else if ( -fTmp >= fA00 ) {
	      fS = T(1.0);
	      fSqrDist = fA00+fA11+fC+T(2.0)*(fB1+fTmp);
	    }
	    else {
	      fS = -fTmp/fA00;
	      fSqrDist = fTmp*fS+fA11+T(2.0)*fB1+fC;
	    }
	  }
	}
	else { // region 7 (side)
	  fT = T(0.0);
	  if ( fB0 >= T(0.0) ) {
	    fS = T(0.0);
	    fSqrDist = fC;
	  }
	  else if ( -fB0 >= fA00 ) {
	    fS = T(1.0);
	    fSqrDist = fA00+T(2.0)*fB0+fC;
	  }
	  else {
	    fS = -fB0/fA00;
	    fSqrDist = fB0*fS+fC;
	  }
	}
      }
      else {
	if ( fT >= T(0.0) ) {
	  if ( fT <= fDet ) { // region 1 (side)
	    fS = T(1.0);
	    fTmp = fA01+fB1;
	    if ( fTmp >= T(0.0) ) {
	      fT = T(0.0);
	      fSqrDist = fA00+T(2.0)*fB0+fC;
	    }
	    else if ( -fTmp >= fA11 ) {
	      fT = T(1.0);
	      fSqrDist = fA00+fA11+fC+T(2.0)*(fB0+fTmp);
	    }
	    else {
	      fT = -fTmp/fA11;
	      fSqrDist = fTmp*fT+fA00+T(2.0)*fB0+fC;
	    }
	  }
	  else { // region 2 (corner)
	    fTmp = fA01+fB0;
	    if ( -fTmp <= fA00 ) {
	      fT = T(1.0);
	      if ( fTmp >= T(0.0) ) {
		fS = T(0.0);
		fSqrDist = fA11+T(2.0)*fB1+fC;
	      }
	      else {
		fS = -fTmp/fA00;
		fSqrDist = fTmp*fS+fA11+T(2.0)*fB1+fC;
	      }
	    }
	    else {
	      fS = T(1.0);
	      fTmp = fA01+fB1;
	      if ( fTmp >= T(0.0) ) {
		fT = T(0.0);
		fSqrDist = fA00+T(2.0)*fB0+fC;
	      }
	      else if ( -fTmp >= fA11 ) {
		fT = T(1.0);
		fSqrDist = fA00+fA11+fC+T(2.0)*(fB0+fTmp);
	      }
	      else {
		fT = -fTmp/fA11;
		fSqrDist = fTmp*fT+fA00+T(2.0)*fB0+fC;
	      }
	    }
	  }
	}
	else { // region 8 (corner)
	  if ( -fB0 < fA00 ) {
	    fT = T(0.0);
	    if ( fB0 >= T(0.0) ) {
	      fS = T(0.0);
	      fSqrDist = fC;
	    }
	    else {
	      fS = -fB0/fA00;
	      fSqrDist = fB0*fS+fC;
	    }
	  }
	  else {
	    fS = T(1.0);
	    fTmp = fA01+fB1;
	    if ( fTmp >= T(0.0) ) {
	      fT = T(0.0);
	      fSqrDist = fA00+T(2.0)*fB0+fC;
	    }
	    else if ( -fTmp >= fA11 ) {
	      fT = T(1.0);
	      fSqrDist = fA00+fA11+fC+T(2.0)*(fB0+fTmp);
	    }
	    else {
	      fT = -fTmp/fA11;
	      fSqrDist = fTmp*fT+fA00+T(2.0)*fB0+fC;
	    }
	  }
	}
      }
    }
    else {
      if ( fT >= T(0.0) ) {
	if ( fT <= fDet ) { // region 5 (side)
	  fS = T(0.0);
	  if ( fB1 >= T(0.0) ) {
	    fT = T(0.0);
	    fSqrDist = fC;
	  }
	  else if ( -fB1 >= fA11 ) {
	    fT = T(1.0);
	    fSqrDist = fA11+T(2.0)*fB1+fC;
	  }
	  else {
	    fT = -fB1/fA11;
	    fSqrDist = fB1*fT+fC;
	  }
	}
	else { // region 4 (corner)
	  fTmp = fA01+fB0;
	  if ( fTmp < T(0.0) ) {
	    fT = T(1.0);
	    if ( -fTmp >= fA00 ) {
	      fS = T(1.0);
	      fSqrDist = fA00+fA11+fC+T(2.0)*(fB1+fTmp);
	    }
	    else {
	      fS = -fTmp/fA00;
	      fSqrDist = fTmp*fS+fA11+T(2.0)*fB1+fC;
	    }
	  }
	  else {
	    fS = T(0.0);
	    if ( fB1 >= T(0.0) ) {
	      fT = T(0.0);
	      fSqrDist = fC;
	    }
	    else if ( -fB1 >= fA11 ) {
	      fT = T(1.0);
	      fSqrDist = fA11+T(2.0)*fB1+fC;
	    }
	    else {
	      fT = -fB1/fA11;
	      fSqrDist = fB1*fT+fC;
	    }
	  }
	}
      }
      else {  // region 6 (corner)
	if ( fB0 < T(0.0) ) {
	  fT = T(0.0);
	  if ( -fB0 >= fA00 ) {
	    fS = T(1.0);
	    fSqrDist = fA00+T(2.0)*fB0+fC;
	  }
	  else {
	    fS = -fB0/fA00;
	    fSqrDist = fB0*fS+fC;
	  }
	}
	else {
	  fS = T(0.0);
	  if ( fB1 >= T(0.0) ) {
	    fT = T(0.0);
	    fSqrDist = fC;
	  }
	  else if ( -fB1 >= fA11 ) {
	    fT = T(1.0);
	    fSqrDist = fA11+T(2.0)*fB1+fC;
	  }
	  else {
	    fT = -fB1/fA11;
	    fSqrDist = fB1*fT+fC;
	  }
	}
      }
    }
  }
  else {
    // line segments are parallel
    if ( fA01 > T(0.0) ) {
      // direction vectors form an obtuse angle
      if ( fB0 >= T(0.0) ) {
	fS = T(0.0);
	fT = T(0.0);
	fSqrDist = fC;
      }
      else if ( -fB0 <= fA00 ) {
	fS = -fB0/fA00;
	fT = T(0.0);
	fSqrDist = fB0*fS+fC;
      }
      else {
	fB1 = -kDiff.dot(rkSeg1.getDirection());
	fS = T(1.0);
	fTmp = fA00+fB0;
	if ( -fTmp >= fA01 ) {
	  fT = T(1.0);
	  fSqrDist = fA00+fA11+fC+T(2.0)*(fA01+fB0+fB1);
	}
	else {
	  fT = -fTmp/fA01;
	  fSqrDist = fA00+T(2.0)*fB0+fC+fT*(fA11*fT+T(2.0)*(fA01+fB1));
	}
      }
    }
    else {
      // direction vectors form an acute angle
      if ( -fB0 >= fA00 ) {
	fS = T(1.0);
	fT = T(0.0);
	fSqrDist = fA00+T(2.0)*fB0+fC;
      }
      else if ( fB0 <= T(0.0) ) {
	fS = -fB0/fA00;
	fT = T(0.0);
	fSqrDist = fB0*fS+fC;
      }
      else {
	fB1 = -kDiff.dot(rkSeg1.getDirection());
	fS = T(0.0);
	if ( fB0 >= -fA01 ) {
	  fT = T(1.0);
	  fSqrDist = fA11+T(2.0)*fB1+fC;
	}
	else {
	  fT = -fB0/fA01;
	  fSqrDist = fC+fT*(T(2.0)*fB1+fA11*fT);
	}
      }
    }
  }
  
  if ( pfSegP0 ) *pfSegP0 = fS;

  if ( pfSegP1 ) *pfSegP1 = fT;

  return GbMath<T>::Abs(fSqrDist);
}

template <class T>
T 
GbDistance<T>::distance (const GbSegment<T>& rkSeg0, const GbSegment<T>& rkSeg1,
			 T* pfSegP0, T* pfSegP1)
{
  return GbMath<T>::Sqrt(GbDistance<T>::squareDistance(rkSeg0,rkSeg1,pfSegP0,pfSegP1));
}

template <class T>
T 
GbDistance<T>::squareDistance (const GbRay<T>& rkRay0, const GbRay<T>& rkRay1,
			       T* pfRayP0, T* pfRayP1)
{
  GbVec3<T> kDiff = rkRay0.getOrigin() - rkRay1.getOrigin();
  T fA00 = rkRay0.getDirection().getSquareNorm();
  T fA01 = -rkRay0.getDirection().dot(rkRay1.getDirection());
  T fA11 = rkRay1.getDirection().getSquareNorm();
  T fB0 = kDiff.dot(rkRay0.getDirection());
  T fC = kDiff.getSquareNorm();
  T fDet = GbMath<T>::Abs(fA00*fA11-fA01*fA01);
  T fB1, fS, fT, fSqrDist;

  if ( fDet >= GbDistance<T>::TOLERANCE ) {
    // rays are not parallel
    fB1 = -kDiff.dot(rkRay1.getDirection());
    fS = fA01*fB1-fA11*fB0;
    fT = fA01*fB0-fA00*fB1;

    if ( fS >= T(0.0) ) {
      if ( fT >= T(0.0) ) { // region 0 (interior)
	// minimum at two interior points of rays
	T fInvDet = T(1.0)/fDet;
	fS *= fInvDet;
	fT *= fInvDet;
	fSqrDist = fS*(fA00*fS+fA01*fT+T(2.0)*fB0) +
	           fT*(fA01*fS+fA11*fT+T(2.0)*fB1) +
	           fC;
      }
      else { // region 3 (side)
	fT = T(0.0);
	if ( fB0 >= T(0.0) ) {
	  fS = T(0.0);
	  fSqrDist = fC;
	}
	else {
	  fS = -fB0/fA00;
	  fSqrDist = fB0*fS+fC;
	}
      }
    }
    else {
      if ( fT >= T(0.0) ) { // region 1 (side)
	fS = T(0.0);
	if ( fB1 >= T(0.0) ) {
	  fT = T(0.0);
	  fSqrDist = fC;
	}
	else {
	  fT = -fB1/fA11;
	  fSqrDist = fB1*fT+fC;
	}
      }
      else { // region 2 (corner)
	if ( fB0 < T(0.0) ) {
	  fS = -fB0/fA00;
	  fT = T(0.0);
	  fSqrDist = fB0*fS+fC;
	}
	else {
	  fS = T(0.0);
	  if ( fB1 >= T(0.0) ) {
	    fT = T(0.0);
	    fSqrDist = fC;
	  }
	  else {
	    fT = -fB1/fA11;
	    fSqrDist = fB1*fT+fC;
	  }
	}
      }
    }
  }
  else {
    // rays are parallel
    if ( fA01 > T(0.0) ) {
      // opposite direction vectors
      fT = T(0.0);
      if ( fB0 >= T(0.0) ) {
	fS = T(0.0);
	fSqrDist = fC;
      }
      else {
	fS = -fB0/fA00;
	fSqrDist = fB0*fS+fC;
      }
    }
    else {
      // same direction vectors
      if ( fB0 >= T(0.0) ) {
	fB1 = -kDiff.dot(rkRay1.getDirection());
	fS = T(0.0);
	fT = -fB1/fA11;
	fSqrDist = fB1*fT+fC;
      }
      else {
	fS = -fB0/fA00;
	fT = T(0.0);
	fSqrDist = fB0*fS+fC;
      }
    }
  }
  
  if ( pfRayP0 ) *pfRayP0 = fS;
  
  if ( pfRayP1 ) *pfRayP1 = fT;

  return GbMath<T>::Abs(fSqrDist);
}

template <class T>
T 
GbDistance<T>::squareDistance (const GbRay<T>& rkRay, const GbSegment<T>& rkSeg,
			       T* pfRayP, T* pfSegP)
{
  GbVec3<T> kDiff = rkRay.getOrigin() - rkSeg.getOrigin();
  T fA00 = rkRay.getDirection().getSquareNorm();
  T fA01 = -rkRay.getDirection().dot(rkSeg.getDirection());
  T fA11 = rkSeg.getDirection().getSquareNorm();
  T fB0 = kDiff.dot(rkRay.getDirection());
  T fC = kDiff.getSquareNorm();
  T fDet = GbMath<T>::Abs(fA00*fA11-fA01*fA01);
  T fB1, fS, fT, fSqrDist, fTmp;

  if ( fDet >= GbDistance<T>::TOLERANCE ) {
    // ray and segment are not parallel
    fB1 = -kDiff.dot(rkSeg.getDirection());
    fS = fA01*fB1-fA11*fB0;
    fT = fA01*fB0-fA00*fB1;
    
    if ( fS >= T(0.0) ) {
      if ( fT >= T(0.0) ) {
	if ( fT <= fDet ) { // region 0
	  // minimum at interior points of ray and segment
	  T fInvDet = T(1.0)/fDet;
	  fS *= fInvDet;
	  fT *= fInvDet;
	  fSqrDist = fS*(fA00*fS+fA01*fT+T(2.0)*fB0) +
	             fT*(fA01*fS+fA11*fT+T(2.0)*fB1) +
	             fC;
	}
	else { // region 1
	  fT = T(1.0);
	  if ( fB0 >= -fA01 ) {
	    fS = T(0.0);
	    fSqrDist = fA11+T(2.0)*fB1+fC;
	  }
	  else {
	    fTmp = fA01 + fB0;
	    fS = -fTmp/fA00;
	    fSqrDist = fTmp*fS+fA11+T(2.0)*fB1+fC;
	  }
	}
      }
      else { // region 5
	fT = T(0.0);
	if ( fB0 >= T(0.0) ) {
	  fS = T(0.0);
	  fSqrDist = fC;
	}
	else {
	  fS = -fB0/fA00;
	  fSqrDist = fB0*fS+fC;
	}
      }
    }
    else {
      if ( fT <= T(0.0) ) { // region 4
	if ( fB0 < T(0.0) ) {
	  fS = -fB0/fA00;
	  fT = T(0.0);
	  fSqrDist = fB0*fS+fC;
	}
	else {
	  fS = T(0.0);
	  if ( fB1 >= T(0.0) ) {
	    fT = T(0.0);
	    fSqrDist = fC;
	  }
	  else if ( -fB1 >= fA11 ) {
	    fT = T(1.0);
	    fSqrDist = fA11+T(2.0)*fB1+fC;
	  }
	  else {
	    fT = -fB1/fA11;
	    fSqrDist = fB1*fT+fC;
	  }
	}
      }
      else if ( fT <= fDet ) { // region 3
	fS = T(0.0);
	if ( fB1 >= T(0.0) ) {
	  fT = T(0.0);
	  fSqrDist = fC;
	}
	else if ( -fB1 >= fA11 ) {
	  fT = T(1.0);
	  fSqrDist = fA11+T(2.0)*fB1+fC;
	}
	else {
	  fT = -fB1/fA11;
	  fSqrDist = fB1*fT+fC;
	}
      }
      else { // region 2
	fTmp = fA01+fB0;
	if ( fTmp < T(0.0) ) {
	  fS = -fTmp/fA00;
	  fT = T(1.0);
	  fSqrDist = fTmp*fS+fA11+T(2.0)*fB1+fC;
	}
	else {
	  fS = T(0.0);
	  if ( fB1 >= T(0.0) ) {
	    fT = T(0.0);
	    fSqrDist = fC;
	  }
	  else if ( -fB1 >= fA11 ) {
	    fT = T(1.0);
	    fSqrDist = fA11+2*fB1+fC;
	  }
	  else {
	    fT = -fB1/fA11;
	    fSqrDist = fB1*fT+fC;
	  }
	}
      }
    }
  }
  else {
    // ray and segment are parallel
    if ( fA01 > T(0.0) ) {
      // opposite direction vectors
      fT = T(0.0);
      if ( fB0 >= T(0.0) ) {
	fS = T(0.0);
	fSqrDist = fC;
      }
      else {
	fS = -fB0/fA00;
	fSqrDist = fB0*fS+fC;
      }
    }
    else {
      // same direction vectors
      fB1 = -kDiff.dot(rkSeg.getDirection());
      fT = T(1.0);
      fTmp = fA01+fB0;
      if ( fTmp >= T(0.0) ) {
	fS = T(0.0);
	fSqrDist = fA11+T(2.0)*fB1+fC;
      }
      else {
	fS = -fTmp/fA00;
	fSqrDist = fTmp*fS+fA11+T(2.0)*fB1+fC;
      }
    }
  }
  
  if ( pfRayP ) *pfRayP = fS;

  if ( pfSegP ) *pfSegP = fT;

  return GbMath<T>::Abs(fSqrDist);
}

template <class T>
T 
GbDistance<T>::distance (const GbRay<T>& rkRay0, const GbRay<T>& rkRay1,
			 T* pfRayP0, T* pfRayP1)
{
  return GbMath<T>::Sqrt(GbDistance<T>::squareDistance(rkRay0,rkRay1,pfRayP0,pfRayP1));
}

template <class T>
T 
GbDistance<T>::distance (const GbRay<T>& rkRay, const GbSegment<T>& rkSeg,
			 T* pfRayP, T* pfSegP)
{
  return GbMath<T>::Sqrt(GbDistance<T>::squareDistance(rkRay,rkSeg,pfRayP,pfSegP));
}

template <class T>
T 
GbDistance<T>::squareDistance (const GbRay<T>& rkRay, const GoTriangleElement& rkTri,
			       T* pfRayP, T* pfTriP0, T* pfTriP1)
{
  GbVec3<T> kDiff = rkTri.getOrigin() - rkRay.getOrigin();
  T fA00 = rkRay.getDirection().getSquareNorm();
  T fA01 = -rkRay.getDirection().dot(rkTri.getEdgeVec(0));
  T fA02 = -rkRay.getDirection().dot(rkTri.getEdgeVec(1));
  T fA11 = rkTri.getEdgeVec(0).getSquareNorm();
  T fA12 = rkTri.getEdgeVec(0).dot(rkTri.getEdgeVec(1));
  T fA22 = rkTri.getEdgeVec(1).dot(rkTri.getEdgeVec(1));
  T fB0  = -kDiff.dot(rkRay.getDirection());
  T fB1  = kDiff.dot(rkTri.getEdgeVec(0));
  T fB2  = kDiff.dot(rkTri.getEdgeVec(1));
  T fCof00 = fA11*fA22-fA12*fA12;
  T fCof01 = fA02*fA12-fA01*fA22;
  T fCof02 = fA01*fA12-fA02*fA11;
  T fDet = fA00*fCof00+fA01*fCof01+fA02*fCof02;

  GbSegment<T> kTriSeg;
  GbVec3<T> kPt;
  T fSqrDist, fSqrDist0, fR, fS, fT, fR0, fS0, fT0;
  T dummyPointer;

  if ( GbMath<T>::Abs(fDet) >= GbDistance<T>::TOLERANCE ) {
    T fCof11 = fA00*fA22-fA02*fA02;
    T fCof12 = fA02*fA01-fA00*fA12;
    T fCof22 = fA00*fA11-fA01*fA01;
    T fInvDet = T(1.0)/fDet;
    T fRhs0 = -fB0*fInvDet;
    T fRhs1 = -fB1*fInvDet;
    T fRhs2 = -fB2*fInvDet;

    fR = fCof00*fRhs0+fCof01*fRhs1+fCof02*fRhs2;
    fS = fCof01*fRhs0+fCof11*fRhs1+fCof12*fRhs2;
    fT = fCof02*fRhs0+fCof12*fRhs1+fCof22*fRhs2;

    if ( fR <= T(0.0) ) {
      if ( fS+fT <= T(1.0) ) {
	if ( fS < T(0.0) ) {
	  if ( fT < T(0.0) ) { // region 4m
	    // min on face s=0 or t=0 or r=0
	    kTriSeg.setOrigin(rkTri.getOrigin());
	    kTriSeg.setDirection(rkTri.getEdgeVec(1));
	    fSqrDist = GbDistance<T>::squareDistance(rkRay,kTriSeg,&fR,&fT);
	    fS = T(0.0);
	    kTriSeg.setOrigin(rkTri.getOrigin());
	    kTriSeg.setDirection(rkTri.getEdgeVec(0));
	    fSqrDist0 = GbDistance<T>::squareDistance(rkRay,kTriSeg,&fR0,&fS0);
	    fT0 = T(0.0);
	    if ( fSqrDist0 < fSqrDist ) {
	      fSqrDist = fSqrDist0;
	      fR = fR0;
	      fS = fS0;
	      fT = fT0;
	    }
	    fSqrDist0 = GbDistance<T>::squareDistance(rkRay.getOrigin(),rkTri,&fS0,&fT0);
	    fR0 = T(0.0);
	    if ( fSqrDist0 < fSqrDist ) {
	      fSqrDist = fSqrDist0;
	      fR = fR0;
	      fS = fS0;
	      fT = fT0;
	    }
	  }
	  else { // region 3m
	    // min on face s=0 or r=0
	    kTriSeg.setOrigin(rkTri.getOrigin());
	    kTriSeg.setDirection(rkTri.getEdgeVec(1));
	    fSqrDist = GbDistance<T>::squareDistance(rkRay,kTriSeg,&fR,&fT);
	    fS = T(0.0);
	    fSqrDist0 = GbDistance<T>::squareDistance(rkRay.getOrigin(),rkTri,&fS0,&fT0);
	    fR0 = T(0.0);
	    if ( fSqrDist0 < fSqrDist ) {
	      fSqrDist = fSqrDist0;
	      fR = fR0;
	      fS = fS0;
	      fT = fT0;
	    }
	  }
	}
	else if ( fT < T(0.0) ) { // region 5m
	  // min on face t=0 or r=0
	  kTriSeg.setOrigin(rkTri.getOrigin());
	  kTriSeg.setDirection(rkTri.getEdgeVec(0));
	  fSqrDist = GbDistance<T>::squareDistance(rkRay,kTriSeg,&fR,&fS);
	  fT = T(0.0);
	  fSqrDist0 = GbDistance<T>::squareDistance(rkRay.getOrigin(),rkTri,&fS0,&fT0);
	  fR0 = T(0.0);
	  if ( fSqrDist0 < fSqrDist ) {
	    fSqrDist = fSqrDist0;
	    fR = fR0;
	    fS = fS0;
	    fT = fT0;
	  }
	}
	else { // region 0m
	  // min on face r=0
	  fSqrDist = GbDistance<T>::squareDistance(rkRay.getOrigin(),rkTri,&fS,&fT);
	  fR = T(0.0);
	}
      }
      else {
	if ( fS < T(0.0) ) { // region 2m
	  // min on face s=0 or s+t=1 or r=0
	  kTriSeg.setOrigin(rkTri.getOrigin());
	  kTriSeg.setDirection(rkTri.getEdgeVec(1));
	  fSqrDist = GbDistance<T>::squareDistance(rkRay,kTriSeg,&fR,&fT);
	  fS = T(0.0);
	  kTriSeg.setOrigin(rkTri.getOrigin()+rkTri.getEdgeVec(0));
	  kTriSeg.setDirection(rkTri.getEdgeVec(2));
	  fSqrDist0 = GbDistance<T>::squareDistance(rkRay,kTriSeg,&fR0,&fT0);
	  fS0 = T(1.0)-fT0;
	  if ( fSqrDist0 < fSqrDist ) {
	    fSqrDist = fSqrDist0;
	    fR = fR0;
	    fS = fS0;
	    fT = fT0;
	  }
	  fSqrDist0 = GbDistance<T>::squareDistance(rkRay.getOrigin(),rkTri,&fS0,&fT0);
	  fR0 = T(0.0);
	  if ( fSqrDist0 < fSqrDist ) {
	    fSqrDist = fSqrDist0;
	    fR = fR0;
	    fS = fS0;
	    fT = fT0;
	  }
	}
	else if ( fT < T(0.0) ) { // region 6m
	  // min on face t=0 or s+t=1 or r=0
	  kTriSeg.setOrigin(rkTri.getOrigin());
	  kTriSeg.setDirection(rkTri.getEdgeVec(0));
	  fSqrDist = GbDistance<T>::squareDistance(rkRay,kTriSeg,&fR,&fS);
	  fT = T(0.0);
	  kTriSeg.setOrigin(rkTri.getOrigin()+rkTri.getEdgeVec(0));
	  kTriSeg.setDirection(rkTri.getEdgeVec(2));
	  fSqrDist0 = GbDistance<T>::squareDistance(rkRay,kTriSeg,&fR0,&fT0);
	  fS0 = T(1.0)-fT0;
	  if ( fSqrDist0 < fSqrDist ) {
	    fSqrDist = fSqrDist0;
	    fR = fR0;
	    fS = fS0;
	    fT = fT0;
	  }
	  fSqrDist0 = GbDistance<T>::squareDistance(rkRay.getOrigin(),rkTri,&fS0,&fT0);
	  fR0 = T(0.0);
	  if ( fSqrDist0 < fSqrDist ) {
	    fSqrDist = fSqrDist0;
	    fR = fR0;
	    fS = fS0;
	    fT = fT0;
	  }
	}
	else { // region 1m
	  // min on face s+t=1 or r=0
	  kTriSeg.setOrigin(rkTri.getOrigin()+rkTri.getEdgeVec(0));
	  kTriSeg.setDirection(rkTri.getEdgeVec(2));
	  fSqrDist = GbDistance<T>::squareDistance(rkRay,kTriSeg,&fR,&fT);
	  fS = T(1.0)-fT;
	  fSqrDist0 = GbDistance<T>::squareDistance(rkRay.getOrigin(),rkTri,&fS0,&fT0);
	  fR0 = T(0.0);
	  if ( fSqrDist0 < fSqrDist ) {
	    fSqrDist = fSqrDist0;
	    fR = fR0;
	    fS = fS0;
	    fT = fT0;
	  }
	}
      }
    }
    else { // fR > 0
      if ( fS+fT <= T(1.0) ) {
	if ( fS < T(0.0) ) {
	  if ( fT < T(0.0) ) { // region 4p
	    // min on face s=0 or t=0
	    kTriSeg.setOrigin(rkTri.getOrigin());
	    kTriSeg.setDirection(rkTri.getEdgeVec(1));
	    fSqrDist = GbDistance<T>::squareDistance(rkRay,kTriSeg,&dummyPointer,&fT);
	    fS = T(0.0);
	    kTriSeg.setOrigin(rkTri.getOrigin());
	    kTriSeg.setDirection(rkTri.getEdgeVec(0));
	    fSqrDist0 = GbDistance<T>::squareDistance(rkRay,kTriSeg,&dummyPointer,&fS0);
	    fT0 = T(0.0);
	    if ( fSqrDist0 < fSqrDist ) {
	      fSqrDist = fSqrDist0;
	      fS = fS0;
	      fT = fT0;
	    }
	  }
	  else { // region 3p
	    // min on face s=0
	    kTriSeg.setOrigin(rkTri.getOrigin());
	    kTriSeg.setDirection(rkTri.getEdgeVec(1));
	    fSqrDist = GbDistance<T>::squareDistance(rkRay,kTriSeg,&dummyPointer,&fT);
	    fS = T(0.0);
	  }
	}
	else if ( fT < T(0.0) ) { // region 5p
	  // min on face t=0
	  kTriSeg.setOrigin(rkTri.getOrigin());
	  kTriSeg.setDirection(rkTri.getEdgeVec(0));
	  fSqrDist = GbDistance<T>::squareDistance(rkRay,kTriSeg,&dummyPointer,&fS);
	  fT = T(0.0);
	}
	else { // region 0p
	  // ray intersects triangle
	  fSqrDist = T(0.0);
	}
      }
      else {
	if ( fS < T(0.0) ) { // region 2p
	  // min on face s=0 or s+t=1
	  kTriSeg.setOrigin(rkTri.getOrigin());
	  kTriSeg.setDirection(rkTri.getEdgeVec(1));
	  fSqrDist = GbDistance<T>::squareDistance(rkRay,kTriSeg,&dummyPointer,&fT);
	  fS = T(0.0);
	  kTriSeg.setOrigin(rkTri.getOrigin()+rkTri.getEdgeVec(0));
	  kTriSeg.setDirection(rkTri.getEdgeVec(2));
	  fSqrDist0 = GbDistance<T>::squareDistance(rkRay,kTriSeg,&dummyPointer,&fT0);
	  fS0 = T(1.0)-fT0;
	  if ( fSqrDist0 < fSqrDist ) {
	    fSqrDist = fSqrDist0;
	    fS = fS0;
	    fT = fT0;
	  }
	}
	else if ( fT < T(0.0) ) { // region 6p
	  // min on face t=0 or s+t=1
	  kTriSeg.setOrigin(rkTri.getOrigin());
	  kTriSeg.setDirection(rkTri.getEdgeVec(0));
	  fSqrDist = GbDistance<T>::squareDistance(rkRay,kTriSeg,&dummyPointer,&fS);
	  fT = T(0.0);
	  kTriSeg.setOrigin(rkTri.getOrigin()+rkTri.getEdgeVec(0));
	  kTriSeg.setDirection(rkTri.getEdgeVec(2));
	  fSqrDist0 = GbDistance<T>::squareDistance(rkRay,kTriSeg,&dummyPointer,&fT0);
	  fS0 = T(1.0)-fT0;
	  if ( fSqrDist0 < fSqrDist ) {
	    fSqrDist = fSqrDist0;
	    fS = fS0;
	    fT = fT0;
	  }
	}
	else { // region 1p
	  // min on face s+t=1
	  kTriSeg.setOrigin(rkTri.getOrigin()+rkTri.getEdgeVec(0));
	  kTriSeg.setDirection(rkTri.getEdgeVec(2));
	  fSqrDist = GbDistance<T>::squareDistance(rkRay,kTriSeg,&dummyPointer,&fT);
	  fS = T(1.0)-fT;
	}
      }
    }
  }
  else {
    // ray and triangle are parallel
    kTriSeg.setOrigin(rkTri.getOrigin());
    kTriSeg.setDirection(rkTri.getEdgeVec(0));
    fSqrDist = GbDistance<T>::squareDistance(rkRay,kTriSeg,&fR,&fS);
    fT = T(0.0);

    kTriSeg.setDirection(rkTri.getEdgeVec(1));
    fSqrDist0 = GbDistance<T>::squareDistance(rkRay,kTriSeg,&fR0,&fT0);
    fS0 = T(0.0);
    if ( fSqrDist0 < fSqrDist ) {
      fSqrDist = fSqrDist0;
      fR = fR0;
      fS = fS0;
      fT = fT0;
    }

    kTriSeg.setOrigin(rkTri.getOrigin() + rkTri.getEdgeVec(0));
    kTriSeg.setDirection(rkTri.getEdgeVec(2));
    fSqrDist0 = GbDistance<T>::squareDistance(rkRay,kTriSeg,&fR0,&fT0);
    fS0 = T(1.0)-fT0;
    if ( fSqrDist0 < fSqrDist ) {
      fSqrDist = fSqrDist0;
      fR = fR0;
      fS = fS0;
      fT = fT0;
    }

    fSqrDist0 = GbDistance<T>::squareDistance(rkRay.getOrigin(),rkTri,&fS0,&fT0);
    fR0 = T(0.0);
    if ( fSqrDist0 < fSqrDist ) {
      fSqrDist = fSqrDist0;
      fR = fR0;
      fS = fS0;
      fT = fT0;
    }
  }
  
  if ( pfRayP ) *pfRayP = fR;

  if ( pfTriP0 ) *pfTriP0 = fS;

  if ( pfTriP1 ) *pfTriP1 = fT;

  return GbMath<T>::Abs(fSqrDist);
}

template <class T>
T 
GbDistance<T>::distance (const GbRay<T>& rkRay, const GoTriangleElement& rkTri,
			 T* pfRayP, T* pfTriP0, T* pfTriP1)
{
  return GbMath<T>::Sqrt(GbDistance<T>::squareDistance(rkRay,rkTri,pfRayP,pfTriP0,pfTriP1));
}

template <class T>
T 
GbDistance<T>::squareDistance (const GbVec3<T>& rkPoint, const GbRay<T>& rkRay,
			       T* pfParam)
{
  GbVec3<T> kDiff = rkPoint - rkRay.getOrigin();
  T fT = kDiff.dot(rkRay.getDirection());

  if ( fT <= T(0.0) ) {
    fT = T(0.0);
  }
  else {
    fT /= rkRay.getDirection().getSquareNorm();
    kDiff -= fT*rkRay.getDirection();
  }

  if ( pfParam ) *pfParam = fT;

  return kDiff.getSquareNorm();
}

template <class T>
T 
GbDistance<T>::squareDistance (const GbVec3<T>& rkPoint, const GbSegment<T>& rkSegment, 
			       T* pfParam)
{
  GbVec3<T> kDiff = rkPoint - rkSegment.getOrigin();
  T fT = kDiff.dot(rkSegment.getDirection());

  if ( fT <= T(0.0) ) {
    fT = T(0.0);
  }
  else {
    T fSqrLen= rkSegment.getDirection().getSquareNorm();
    if ( fT >= fSqrLen ) {
      fT = T(1.0);
      kDiff -= rkSegment.getDirection();
    }
    else {
      fT /= fSqrLen;
      kDiff -= fT*rkSegment.getDirection();
    }
  }

  if ( pfParam ) *pfParam = fT;

  return kDiff.getSquareNorm();
}

template <class T>
T 
GbDistance<T>::distance (const GbVec3<T>& rkPoint, const GbRay<T>& rkRay,
			 T* pfParam)
{
  return GbMath<T>::Sqrt(GbDistance<T>::squareDistance(rkPoint,rkRay,pfParam));
}

template <class T>
T 
GbDistance<T>::distance (const GbVec3<T>& rkPoint, const GbSegment<T>& rkSegment,
			 T* pfParam)
{
  return GbMath<T>::Sqrt(GbDistance<T>::squareDistance(rkPoint,rkSegment,pfParam));
}

template <class T>
T 
GbDistance<T>::squareDistance (const GbRay<T>& rkRay, const GbBox3<T>& rkBox,
 			       T* pfLParam, T* pfBParam0, T* pfBParam1, T* pfBParam2)
{
#ifdef DEBUG
  // The four parameters pointers are either all non-null or all null.
  if ( pfLParam ) {
    assert( pfBParam0 && pfBParam1 && pfBParam2 );
  }
  else {
    assert( !pfBParam0 && !pfBParam1 && !pfBParam2 );
  }
#endif
  
  GbLine<T> kLine;
  kLine.setOrigin(rkRay.getOrigin());
  kLine.setDirection(rkRay.getDirection());
  
  T fLP, fBP0, fBP1, fBP2;
  T fSqrDistance = GbDistance<T>::squareDistance(kLine,rkBox,&fLP,&fBP0,&fBP1,&fBP2);
  if ( fLP >= T(0.0) ) {
    if ( pfLParam ) {
      *pfLParam = fLP;
      *pfBParam0 = fBP0;
      *pfBParam1 = fBP1;
      *pfBParam2 = fBP2;
    }
    
    return fSqrDistance;
  }
  else {
    fSqrDistance = GbDistance<T>::squareDistance(rkRay.getOrigin(),rkBox,
 						 pfBParam0,pfBParam1,pfBParam2);
    
    if ( pfLParam ) *pfLParam = T(0.0);
    
    return fSqrDistance;
  }
}

template <class T>
T 
GbDistance<T>::squareDistance (const GbSegment<T>& rkSeg, const GbBox3<T>& rkBox,
 			       T* pfLParam, T* pfBParam0, T* pfBParam1, T* pfBParam2)
{
#ifdef DEBUG
  // The four parameters pointers are either all non-null or all null.
  if ( pfLParam ) {
    assert( pfBParam0 && pfBParam1 && pfBParam2 );
  }
  else {
    assert( !pfBParam0 && !pfBParam1 && !pfBParam2 );
  }
#endif
  
  GbLine<T> kLine;
  kLine.setOrigin(rkSeg.getOrigin());
  kLine.setDirection(rkSeg.getDirection());
  
  T fLP, fBP0, fBP1, fBP2;
  T fSqrDistance = GbDistance<T>::squareDistance(kLine,rkBox,&fLP,&fBP0,&fBP1,&fBP2);
  if ( fLP >= T(0.0) ) {
    if ( fLP <= T(1.0) ) {
      if ( pfLParam ) {
 	*pfLParam = fLP;
 	*pfBParam0 = fBP0;
 	*pfBParam1 = fBP1;
 	*pfBParam2 = fBP2;
      }
      
      return fSqrDistance;
    }
    else {
      fSqrDistance = GbDistance<T>::squareDistance(rkSeg.getOrigin()+rkSeg.getDirection(),rkBox,
 						   pfBParam0,pfBParam1,pfBParam2);
      
      if ( pfLParam ) *pfLParam = T(1.0);
      
      return fSqrDistance;
    }
  }
  else {
    fSqrDistance = GbDistance<T>::squareDistance(rkSeg.getOrigin(),rkBox,
 						 pfBParam0,pfBParam1,pfBParam2);
    
    if ( pfLParam ) *pfLParam = T(0.0);
    
    return fSqrDistance;
  }
}

template <class T>
T 
GbDistance<T>::distance (const GbRay<T>& rkRay, const GbBox3<T>& rkBox,
 			 T* pfLParam, T* pfBParam0, T* pfBParam1, T* pfBParam2)
{
  return GbMath<T>::Sqrt(GbDistance<T>::squareDistance(rkRay,rkBox,
 						       pfLParam,pfBParam0,pfBParam1,pfBParam2));
}

template <class T>
T 
GbDistance<T>::distance (const GbSegment<T>& rkSeg, const GbBox3<T>& rkBox,
 			 T* pfLParam, T* pfBParam0, T* pfBParam1, T* pfBParam2)
{
  return GbMath<T>::Sqrt(GbDistance<T>::squareDistance(rkSeg,rkBox,
 						       pfLParam,pfBParam0,pfBParam1,pfBParam2));
}

template <class T>
T 
GbDistance<T>::squareDistance (const GbLine<T>& rkLine0, const GbLine<T>& rkLine1,
			       T* pfLinP0, T* pfLinP1)
{
  GbVec3<T> kDiff = rkLine0.getOrigin() - rkLine1.getOrigin();
  T fA00 = rkLine0.getDirection().getSquareNorm();
  T fA01 = -rkLine0.getDirection().dot(rkLine1.getDirection());
  T fA11 = rkLine1.getDirection().getSquareNorm();
  T fB0 = kDiff.dot(rkLine0.getDirection());
  T fC = kDiff.getSquareNorm();
  T fDet = GbMath<T>::Abs(fA00*fA11-fA01*fA01);
  T fB1, fS, fT, fSqrDist;

  if ( fDet >= GbDistance<T>::TOLERANCE ) {
    // lines are not parallel
    fB1 = -kDiff.dot(rkLine1.getDirection());
    T fInvDet = T(1.0)/fDet;
    fS = (fA01*fB1-fA11*fB0)*fInvDet;
    fT = (fA01*fB0-fA00*fB1)*fInvDet;
    fSqrDist = fS*(fA00*fS+fA01*fT+T(2.0)*fB0) +
               fT*(fA01*fS+fA11*fT+T(2.0)*fB1) +
               fC;
  }
  else {
    // lines are parallel, select any closest pair of points
    fS = -fB0/fA00;
    fT = T(0.0);
    fSqrDist = fB0*fS+fC;
  }

  if ( pfLinP0 ) *pfLinP0 = fS;

  if ( pfLinP1 ) *pfLinP1 = fT;

  return GbMath<T>::Abs(fSqrDist);
}

template <class T>
T 
GbDistance<T>::squareDistance (const GbLine<T>& rkLine, const GbRay<T>& rkRay,
			       T* pfLinP, T* pfRayP)
{
  GbVec3<T> kDiff = rkLine.getOrigin() - rkRay.getOrigin();
  T fA00 = rkLine.getDirection().getSquareNorm();
  T fA01 = -rkLine.getDirection().dot(rkRay.getDirection());
  T fA11 = rkRay.getDirection().getSquareNorm();
  T fB0 = kDiff.dot(rkLine.getDirection());
  T fC = kDiff.getSquareNorm();
  T fDet = GbMath<T>::Abs(fA00*fA11-fA01*fA01);
  T fB1, fS, fT, fSqrDist;

  if ( fDet >= GbDistance<T>::TOLERANCE ) {
    fB1 = -kDiff.dot(rkRay.getDirection());
    fT = fA01*fB0-fA00*fB1;

    if ( fT >= T(0.0) ) {
      // two interior points are closest, one on line and one on ray
      T fInvDet = T(1.0)/fDet;
      fS = (fA01*fB1-fA11*fB0)*fInvDet;
      fT *= fInvDet;
      fSqrDist = fS*(fA00*fS+fA01*fT+T(2.0)*fB0) +
	         fT*(fA01*fS+fA11*fT+T(2.0)*fB1) +
	         fC;
    }
    else {
      // end point of ray and interior point of line are closest
      fS = -fB0/fA00;
      fT = T(0.0);
      fSqrDist = fB0*fS+fC;
    }
  }
  else {
    // lines are parallel, closest pair with one point at ray origin
    fS = -fB0/fA00;
    fT = T(0.0);
    fSqrDist = fB0*fS+fC;
  }

  if ( pfLinP ) *pfLinP = fS;

  if ( pfRayP ) *pfRayP = fT;

  return GbMath<T>::Abs(fSqrDist);
}

template <class T>
T 
GbDistance<T>::squareDistance (const GbLine<T>& rkLine, const GbSegment<T>& rkSeg,
			       T* pfLinP, T* pfSegP)
{
  GbVec3<T> kDiff = rkLine.getOrigin() - rkSeg.getOrigin();
  T fA00 = rkLine.getDirection().getSquareNorm();
  T fA01 = -rkLine.getDirection().dot(rkSeg.getDirection());
  T fA11 = rkSeg.getDirection().getSquareNorm();
  T fB0 = kDiff.dot(rkLine.getDirection());
  T fC = kDiff.getSquareNorm();
  T fDet = GbMath<T>::Abs(fA00*fA11-fA01*fA01);
  T fB1, fS, fT, fSqrDist;

  if ( fDet >= GbDistance<T>::TOLERANCE ) {
    fB1 = -kDiff.dot(rkSeg.getDirection());
    fT = fA01*fB0-fA00*fB1;

    if ( fT >= T(0.0) ) {
      if ( fT <= fDet ) {
	// two interior points are closest, one on line and one on segment
	T fInvDet = T(1.0)/fDet;
	fS = (fA01*fB1-fA11*fB0)*fInvDet;
	fT *= fInvDet;
	fSqrDist = fS*(fA00*fS+fA01*fT+T(2.0)*fB0) +
	           fT*(fA01*fS+fA11*fT+T(2.0)*fB1) +
	           fC;
      }
      else {
	// end point e1 of segment and interior point of line are closest
	T fTmp = fA01+fB0;
	fS = -fTmp/fA00;
	fT = T(1.0);
	fSqrDist = fTmp*fS+fA11+T(2.0)*fB1+fC;
      }
    }
    else {
      // end point e0 of segment and interior point of line are closest
      fS = -fB0/fA00;
      fT = T(0.0);
      fSqrDist = fB0*fS+fC;
    }
  }
  else {
    // lines are parallel, closest pair with one point at segment origin
    fS = -fB0/fA00;
    fT = T(0.0);
    fSqrDist = fB0*fS+fC;
  }

  if ( pfLinP ) *pfLinP = fS;

  if ( pfSegP ) *pfSegP = fT;

  return GbMath<T>::Abs(fSqrDist);
}

template <class T>
T 
GbDistance<T>::distance (const GbLine<T>& rkLine0, const GbLine<T>& rkLine1,
			 T* pfLinP0, T* pfLinP1)
{
  return GbMath<T>::Sqrt(GbDistance<T>::squareDistance(rkLine0,rkLine1,pfLinP0,pfLinP1));
}

template <class T>
T 
GbDistance<T>::distance (const GbLine<T>& rkLine, const GbRay<T>& rkRay,
			 T* pfLinP, T* pfRayP)
{
  return GbMath<T>::Sqrt(GbDistance<T>::squareDistance(rkLine,rkRay,pfLinP,pfRayP));
}

template <class T>
T 
GbDistance<T>::distance (const GbLine<T>& rkLine, const GbSegment<T>& rkSeg,
			 T* pfLinP, T* pfSegP)
{
  return GbMath<T>::Sqrt(GbDistance<T>::squareDistance(rkLine,rkSeg,pfLinP,pfSegP));
}

template <class T>
T 
GbDistance<T>::squareDistance (const GbVec3<T>& rkPoint, const GbLine<T>& rkLine,
			       T* pfParam)
{
  GbVec3<T> kDiff = rkPoint - rkLine.getOrigin();
  T fSqrLen = rkLine.getDirection().getSquareNorm();
  T fT = kDiff.dot(rkLine.getDirection())/fSqrLen;
  kDiff -= fT*rkLine.getDirection();

  if ( pfParam ) *pfParam = fT;

  return kDiff.getSquareNorm();
}

template <class T>
T 
GbDistance<T>::distance (const GbVec3<T>& rkPoint, const GbLine<T>& rkLine,
			 T* pfParam)
{
  return GbMath<T>::Sqrt(GbDistance<T>::squareDistance(rkPoint,rkLine,pfParam));
}

template <class T>
void 
GbDistance<T>::face (int i0, int i1, int i2, 
		     GbVec3<T>& rkPnt, const GbVec3<T>& rkDir, 
		     const GbBox3<T>& rkBox, const GbVec3<T>& rkPmE,
		     T* pfLParam, T& rfSqrDistance)
{
  GbVec3<T> kPpE;
  T fLSqr, fInv, fTmp, fParam, fT, fDelta;

  kPpE[i1] = rkPnt[i1] + rkBox.getExtent(i1);
  kPpE[i2] = rkPnt[i2] + rkBox.getExtent(i2);
  if ( rkDir[i0]*kPpE[i1] >= rkDir[i1]*rkPmE[i0] ) {
    if ( rkDir[i0]*kPpE[i2] >= rkDir[i2]*rkPmE[i0] ) {
      // v[i1] >= -e[i1], v[i2] >= -e[i2] (distance = 0)
      if ( pfLParam ) {
	rkPnt[i0] = rkBox.getExtent(i0);
	fInv = T(1.0)/rkDir[i0];
	rkPnt[i1] -= rkDir[i1]*rkPmE[i0]*fInv;
	rkPnt[i2] -= rkDir[i2]*rkPmE[i0]*fInv;
	*pfLParam = -rkPmE[i0]*fInv;
      }
    }
    else {
      // v[i1] >= -e[i1], v[i2] < -e[i2]
      fLSqr = rkDir[i0]*rkDir[i0] + rkDir[i2]*rkDir[i2];
      fTmp = fLSqr*kPpE[i1] - rkDir[i1]*(rkDir[i0]*rkPmE[i0] +
					 rkDir[i2]*kPpE[i2]);
      if ( fTmp <= T(2.0)*fLSqr*rkBox.getExtent(i1) ) {
	fT = fTmp/fLSqr;
	fLSqr += rkDir[i1]*rkDir[i1];
	fTmp = kPpE[i1] - fT;
	fDelta = rkDir[i0]*rkPmE[i0] + 
	         rkDir[i1]*fTmp +
	         rkDir[i2]*kPpE[i2];
	fParam = -fDelta/fLSqr;
	rfSqrDistance += rkPmE[i0]*rkPmE[i0] + 
	                 fTmp*fTmp +
	                 kPpE[i2]*kPpE[i2] + 
	                 fDelta*fParam;

	if ( pfLParam )	{
	  *pfLParam = fParam;
	  rkPnt[i0] = rkBox.getExtent(i0);
	  rkPnt[i1] = fT - rkBox.getExtent(i1);
	  rkPnt[i2] = -rkBox.getExtent(i2);
	}
      }
      else {
	fLSqr += rkDir[i1]*rkDir[i1];
	fDelta = rkDir[i0]*rkPmE[i0] + 
	         rkDir[i1]*rkPmE[i1] +
	         rkDir[i2]*kPpE[i2];
	fParam = -fDelta/fLSqr;
	rfSqrDistance += rkPmE[i0]*rkPmE[i0] + 
	                 rkPmE[i1]*rkPmE[i1] +
	                 kPpE[i2]*kPpE[i2] + 
	                 fDelta*fParam;

	if ( pfLParam )	{
	  *pfLParam = fParam;
	  rkPnt[i0] = rkBox.getExtent(i0);
	  rkPnt[i1] = rkBox.getExtent(i1);
	  rkPnt[i2] = -rkBox.getExtent(i2);
	}
      }
    }
  }
  else {
    if ( rkDir[i0]*kPpE[i2] >= rkDir[i2]*rkPmE[i0] ) {
      // v[i1] < -e[i1], v[i2] >= -e[i2]
      fLSqr = rkDir[i0]*rkDir[i0] + rkDir[i1]*rkDir[i1];
      fTmp = fLSqr*kPpE[i2] - rkDir[i2]*(rkDir[i0]*rkPmE[i0] +
					 rkDir[i1]*kPpE[i1]);
      if ( fTmp <= T(2.0)*fLSqr*rkBox.getExtent(i2) ) {
	fT = fTmp/fLSqr;
	fLSqr += rkDir[i2]*rkDir[i2];
	fTmp = kPpE[i2] - fT;
	fDelta = rkDir[i0]*rkPmE[i0] + 
	         rkDir[i1]*kPpE[i1] +
	         rkDir[i2]*fTmp;
	fParam = -fDelta/fLSqr;
	rfSqrDistance += rkPmE[i0]*rkPmE[i0] + 
	                 kPpE[i1]*kPpE[i1] +
	                 fTmp*fTmp + 
	                 fDelta*fParam;

	if ( pfLParam )	{
	  *pfLParam = fParam;
	  rkPnt[i0] = rkBox.getExtent(i0);
	  rkPnt[i1] = -rkBox.getExtent(i1);
	  rkPnt[i2] = fT - rkBox.getExtent(i2);
	}
      }
      else {
	fLSqr += rkDir[i2]*rkDir[i2];
	fDelta = rkDir[i0]*rkPmE[i0] + 
	         rkDir[i1]*kPpE[i1] +
	         rkDir[i2]*rkPmE[i2];
	fParam = -fDelta/fLSqr;
	rfSqrDistance += rkPmE[i0]*rkPmE[i0] + 
	                 kPpE[i1]*kPpE[i1] +
	                 rkPmE[i2]*rkPmE[i2] + 
	                 fDelta*fParam;

	if ( pfLParam )	{
	  *pfLParam = fParam;
	  rkPnt[i0] = rkBox.getExtent(i0);
	  rkPnt[i1] = -rkBox.getExtent(i1);
	  rkPnt[i2] = rkBox.getExtent(i2);
	}
      }
    }
    else {
      // v[i1] < -e[i1], v[i2] < -e[i2]
      fLSqr = rkDir[i0]*rkDir[i0]+rkDir[i2]*rkDir[i2];
      fTmp = fLSqr*kPpE[i1] - rkDir[i1]*(rkDir[i0]*rkPmE[i0] +
					 rkDir[i2]*kPpE[i2]);
      if ( fTmp >= T(0.0) ) {
	// v[i1]-edge is closest
	if ( fTmp <= T(2.0)*fLSqr*rkBox.getExtent(i1) ) {
	  fT = fTmp/fLSqr;
	  fLSqr += rkDir[i1]*rkDir[i1];
	  fTmp = kPpE[i1] - fT;
	  fDelta = rkDir[i0]*rkPmE[i0] + 
	           rkDir[i1]*fTmp +
	           rkDir[i2]*kPpE[i2];
	  fParam = -fDelta/fLSqr;
	  rfSqrDistance += rkPmE[i0]*rkPmE[i0] + 
	                   fTmp*fTmp +
	                   kPpE[i2]*kPpE[i2] + 
	                   fDelta*fParam;

	  if ( pfLParam ) {
	    *pfLParam = fParam;
	    rkPnt[i0] = rkBox.getExtent(i0);
	    rkPnt[i1] = fT - rkBox.getExtent(i1);
	    rkPnt[i2] = -rkBox.getExtent(i2);
	  }
	}
	else {
	  fLSqr += rkDir[i1]*rkDir[i1];
	  fDelta = rkDir[i0]*rkPmE[i0] + 
	           rkDir[i1]*rkPmE[i1] +
	           rkDir[i2]*kPpE[i2];
	  fParam = -fDelta/fLSqr;
	  rfSqrDistance += rkPmE[i0]*rkPmE[i0] + 
	                   rkPmE[i1]*rkPmE[i1] +
	                   kPpE[i2]*kPpE[i2] + 
	                   fDelta*fParam;

	  if ( pfLParam ) {
	    *pfLParam = fParam;
	    rkPnt[i0] = rkBox.getExtent(i0);
	    rkPnt[i1] = rkBox.getExtent(i1);
	    rkPnt[i2] = -rkBox.getExtent(i2);
	  }
	}
	return;
      }

      fLSqr = rkDir[i0]*rkDir[i0] + rkDir[i1]*rkDir[i1];
      fTmp = fLSqr*kPpE[i2] - rkDir[i2]*(rkDir[i0]*rkPmE[i0] +
					 rkDir[i1]*kPpE[i1]);
      if ( fTmp >= T(0.0) ) {
	// v[i2]-edge is closest
	if ( fTmp <= T(2.0)*fLSqr*rkBox.getExtent(i2) ) {
	  fT = fTmp/fLSqr;
	  fLSqr += rkDir[i2]*rkDir[i2];
	  fTmp = kPpE[i2] - fT;
	  fDelta = rkDir[i0]*rkPmE[i0] + 
	           rkDir[i1]*kPpE[i1] +
	           rkDir[i2]*fTmp;
	  fParam = -fDelta/fLSqr;
	  rfSqrDistance += rkPmE[i0]*rkPmE[i0] + 
	                   kPpE[i1]*kPpE[i1] +
	                   fTmp*fTmp + 
	                   fDelta*fParam;

	  if ( pfLParam ) {
	    *pfLParam = fParam;
	    rkPnt[i0] = rkBox.getExtent(i0);
	    rkPnt[i1] = -rkBox.getExtent(i1);
	    rkPnt[i2] = fT - rkBox.getExtent(i2);
	  }
	}
	else {
	  fLSqr += rkDir[i2]*rkDir[i2];
	  fDelta = rkDir[i0]*rkPmE[i0] + 
	           rkDir[i1]*kPpE[i1] +
	           rkDir[i2]*rkPmE[i2];
	  fParam = -fDelta/fLSqr;
	  rfSqrDistance += rkPmE[i0]*rkPmE[i0] + 
	                   kPpE[i1]*kPpE[i1] +
	                   rkPmE[i2]*rkPmE[i2] + 
	                   fDelta*fParam;

	  if ( pfLParam ) {
	    *pfLParam = fParam;
	    rkPnt[i0] = rkBox.getExtent(i0);
	    rkPnt[i1] = -rkBox.getExtent(i1);
	    rkPnt[i2] = rkBox.getExtent(i2);
	  }
	}
	return;
      }

      // (v[i1],v[i2])-corner is closest
      fLSqr += rkDir[i2]*rkDir[i2];
      fDelta = rkDir[i0]*rkPmE[i0] + 
	       rkDir[i1]*kPpE[i1] +
	       rkDir[i2]*kPpE[i2];
      fParam = -fDelta/fLSqr;
      rfSqrDistance += rkPmE[i0]*rkPmE[i0] + 
	               kPpE[i1]*kPpE[i1] +
	               kPpE[i2]*kPpE[i2] + 
	               fDelta*fParam;

      if ( pfLParam ) {
	*pfLParam = fParam;
	rkPnt[i0] = rkBox.getExtent(i0);
	rkPnt[i1] = -rkBox.getExtent(i1);
	rkPnt[i2] = -rkBox.getExtent(i2);
      }
    }
  }
}

template <class T>
void 
GbDistance<T>::caseNoZeros (GbVec3<T>& rkPnt, const GbVec3<T>& rkDir, const GbBox3<T>& rkBox, 
			    T* pfLParam, T& rfSqrDistance)
{
  GbVec3<T> kPmE(rkPnt[0] - rkBox.getExtent(0), 
		 rkPnt[1] - rkBox.getExtent(1),
		 rkPnt[2] - rkBox.getExtent(2));

  T fProdDxPy, fProdDyPx, fProdDzPx, fProdDxPz, fProdDzPy, fProdDyPz;

  fProdDxPy = rkDir[0]*kPmE[1];
  fProdDyPx = rkDir[1]*kPmE[0];
  if ( fProdDyPx >= fProdDxPy ) {
    fProdDzPx = rkDir[2]*kPmE[0];
    fProdDxPz = rkDir[0]*kPmE[2];
    if ( fProdDzPx >= fProdDxPz ) {
      // line intersects x = e0
      GbDistance<T>::face(0,1,2,rkPnt,rkDir,rkBox,kPmE,pfLParam,rfSqrDistance);
    }
    else {
      // line intersects z = e2
      GbDistance<T>::face(2,0,1,rkPnt,rkDir,rkBox,kPmE,pfLParam,rfSqrDistance);
    }
  }
  else {
    fProdDzPy = rkDir[2]*kPmE[1];
    fProdDyPz = rkDir[1]*kPmE[2];
    if ( fProdDzPy >= fProdDyPz ) {
      // line intersects y = e1
      GbDistance<T>::face(1,2,0,rkPnt,rkDir,rkBox,kPmE,pfLParam,rfSqrDistance);
    }
    else {
      // line intersects z = e2
      GbDistance<T>::face(2,0,1,rkPnt,rkDir,rkBox,kPmE,pfLParam,rfSqrDistance);
    }
  }
}

template <class T>
void 
GbDistance<T>::case0 (int i0, int i1, int i2, 
		      GbVec3<T>& rkPnt, const GbVec3<T>& rkDir, const GbBox3<T>& rkBox, 
		      T* pfLParam, T& rfSqrDistance)
{
  T fPmE0 = rkPnt[i0] - rkBox.getExtent(i0);
  T fPmE1 = rkPnt[i1] - rkBox.getExtent(i1);
  T fProd0 = rkDir[i1]*fPmE0;
  T fProd1 = rkDir[i0]*fPmE1;
  T fDelta, fInvLSqr, fInv;

  if ( fProd0 >= fProd1 ) {
    // line intersects P[i0] = e[i0]
    rkPnt[i0] = rkBox.getExtent(i0);

    T fPpE1 = rkPnt[i1] + rkBox.getExtent(i1);
    fDelta = fProd0 - rkDir[i0]*fPpE1;
    if ( fDelta >= T(0.0) ) {
      fInvLSqr = T(1.0)/(rkDir[i0]*rkDir[i0] + rkDir[i1]*rkDir[i1]);
      rfSqrDistance += fDelta*fDelta*fInvLSqr;
      if ( pfLParam ) {
	rkPnt[i1] = -rkBox.getExtent(i1);
	*pfLParam = -(rkDir[i0]*fPmE0+rkDir[i1]*fPpE1)*fInvLSqr;
      }
    }
    else {
      if ( pfLParam ) {
	fInv = T(1.0)/rkDir[i0];
	rkPnt[i1] -= fProd0*fInv;
	*pfLParam = -fPmE0*fInv;
      }
    }
  }
  else {
    // line intersects P[i1] = e[i1]
    rkPnt[i1] = rkBox.getExtent(i1);

    T fPpE0 = rkPnt[i0] + rkBox.getExtent(i0);
    fDelta = fProd1 - rkDir[i1]*fPpE0;
    if ( fDelta >= T(0.0) ) {
      fInvLSqr = T(1.0)/(rkDir[i0]*rkDir[i0] + rkDir[i1]*rkDir[i1]);
      rfSqrDistance += fDelta*fDelta*fInvLSqr;
      if ( pfLParam ) {
	rkPnt[i0] = -rkBox.getExtent(i0);
	*pfLParam = -(rkDir[i0]*fPpE0+rkDir[i1]*fPmE1)*fInvLSqr;
      }
    }
    else {
      if ( pfLParam ) {
	fInv = T(1.0)/rkDir[i1];
	rkPnt[i0] -= fProd1*fInv;
	*pfLParam = -fPmE1*fInv;
      }
    }
  }

  if ( rkPnt[i2] < -rkBox.getExtent(i2) ) {
    fDelta = rkPnt[i2] + rkBox.getExtent(i2);
    rfSqrDistance += fDelta*fDelta;
    rkPnt[i2] = -rkBox.getExtent(i2);
  }
  else if ( rkPnt[i2] > rkBox.getExtent(i2) ) {
    fDelta = rkPnt[i2] - rkBox.getExtent(i2);
    rfSqrDistance += fDelta*fDelta;
    rkPnt[i2] = rkBox.getExtent(i2);
  }
}

template <class T>
void 
GbDistance<T>::case00 (int i0, int i1, int i2, 
		       GbVec3<T>& rkPnt, const GbVec3<T>& rkDir, const GbBox3<T>& rkBox, 
		       T* pfLParam, T& rfSqrDistance)
{
  T fDelta;

  if ( pfLParam ) *pfLParam = (rkBox.getExtent(i0) - rkPnt[i0])/rkDir[i0];

  rkPnt[i0] = rkBox.getExtent(i0);

  if ( rkPnt[i1] < -rkBox.getExtent(i1) ) {
    fDelta = rkPnt[i1] + rkBox.getExtent(i1);
    rfSqrDistance += fDelta*fDelta;
    rkPnt[i1] = -rkBox.getExtent(i1);
  }
  else if ( rkPnt[i1] > rkBox.getExtent(i1) ) {
    fDelta = rkPnt[i1] - rkBox.getExtent(i1);
    rfSqrDistance += fDelta*fDelta;
    rkPnt[i1] = rkBox.getExtent(i1);
  }

  if ( rkPnt[i2] < -rkBox.getExtent(i2) ) {
    fDelta = rkPnt[i2] + rkBox.getExtent(i2);
    rfSqrDistance += fDelta*fDelta;
    rkPnt[i1] = -rkBox.getExtent(i2);
  }
  else if ( rkPnt[i2] > rkBox.getExtent(i2) ) {
    fDelta = rkPnt[i2] - rkBox.getExtent(i2);
    rfSqrDistance += fDelta*fDelta;
    rkPnt[i2] = rkBox.getExtent(i2);
  }
}

template <class T>
void 
GbDistance<T>::case000 (GbVec3<T>& rkPnt, const GbBox3<T>& rkBox, T& rfSqrDistance)
{
  T fDelta;

  if ( rkPnt[0] < -rkBox.getExtent(0) ) {
    fDelta = rkPnt[0] + rkBox.getExtent(0);
    rfSqrDistance += fDelta*fDelta;
    rkPnt[0] = -rkBox.getExtent(0);
  }
  else if ( rkPnt[0] > rkBox.getExtent(0) ) {
    fDelta = rkPnt[0] - rkBox.getExtent(0);
    rfSqrDistance += fDelta*fDelta;
    rkPnt[0] = rkBox.getExtent(0);
  }

  if ( rkPnt[1] < -rkBox.getExtent(1) ) {
    fDelta = rkPnt[1] + rkBox.getExtent(1);
    rfSqrDistance += fDelta*fDelta;
    rkPnt[1] = -rkBox.getExtent(1);
  }
  else if ( rkPnt[1] > rkBox.getExtent(1) ) {
    fDelta = rkPnt[1] - rkBox.getExtent(1);
    rfSqrDistance += fDelta*fDelta;
    rkPnt[1] = rkBox.getExtent(1);
  }

  if ( rkPnt[2] < -rkBox.getExtent(2) ) {
    fDelta = rkPnt[2] + rkBox.getExtent(2);
    rfSqrDistance += fDelta*fDelta;
    rkPnt[2] = -rkBox.getExtent(2);
  }
  else if ( rkPnt[2] > rkBox.getExtent(2) ) {
    fDelta = rkPnt[2] - rkBox.getExtent(2);
    rfSqrDistance += fDelta*fDelta;
    rkPnt[2] = rkBox.getExtent(2);
  }
}

template <class T>
T 
GbDistance<T>::squareDistance (const GbLine<T>& rkLine, const GbBox3<T>& rkBox,
			       T* pfLParam, T* pfBParam0, T* pfBParam1, T* pfBParam2)
{
#ifdef DEBUG
  // The four parameters pointers are either all non-null or all null.
  if ( pfLParam ) {
    assert( pfBParam0 && pfBParam1 && pfBParam2 );
  }
  else {
    assert( !pfBParam0 && !pfBParam1 && !pfBParam2 );
  }
#endif

  // compute coordinates of line in box coordinate system
  GbVec3<T> kDiff = rkLine.getOrigin() - rkBox.getCenter();
  GbVec3<T> kPnt(kDiff.dot(rkBox.getAxis(0)),
		 kDiff.dot(rkBox.getAxis(1)),
		 kDiff.dot(rkBox.getAxis(2)));
  GbVec3<T> kDir(rkLine.getDirection().dot(rkBox.getAxis(0)),
		 rkLine.getDirection().dot(rkBox.getAxis(1)),
		 rkLine.getDirection().dot(rkBox.getAxis(2)));

  // Apply reflections so that direction vector has nonnegative components.
  GbBool bReflect[3];
  int i;
  for (i = 0; i < 3; i++) {
    if ( kDir[i] < T(0.0) ) {
      kPnt[i] = -kPnt[i];
      kDir[i] = -kDir[i];
      bReflect[i] = true;
    }
    else {
      bReflect[i] = false;
    }
  }

  T fSqrDistance = T(0.0);

  if ( kDir[0] > T(0.0) ) {
    if ( kDir[1] > T(0.0) ) {
      if ( kDir[2] > T(0.0) ) {
	// (+,+,+)
	GbDistance<T>::caseNoZeros(kPnt,kDir,rkBox,pfLParam,fSqrDistance);
      }
      else {
	// (+,+,0)
	GbDistance<T>::case0(0,1,2,kPnt,kDir,rkBox,pfLParam,fSqrDistance);
      }
    }
    else {
      if ( kDir[2] > T(0.0) ) {
	// (+,0,+)
	GbDistance<T>::case0(0,2,1,kPnt,kDir,rkBox,pfLParam,fSqrDistance);
      }
      else {
	// (+,0,0)
	GbDistance<T>::case00(0,1,2,kPnt,kDir,rkBox,pfLParam,fSqrDistance);
      }
    }
  }
  else {
    if ( kDir[1] > T(0.0) ) {
      if ( kDir[2] > T(0.0) ) {
	// (0,+,+)
	GbDistance<T>::case0(1,2,0,kPnt,kDir,rkBox,pfLParam,fSqrDistance);
      }
      else {
	// (0,+,0)
	GbDistance<T>::case00(1,0,2,kPnt,kDir,rkBox,pfLParam,fSqrDistance);
      }
    }
    else {
      if ( kDir[2] > T(0.0) ) {
	// (0,0,+)
	GbDistance<T>::case00(2,0,1,kPnt,kDir,rkBox,pfLParam,fSqrDistance);
      }
      else {
	// (0,0,0)
	GbDistance<T>::case000(kPnt,rkBox,fSqrDistance);
	if ( pfLParam ) *pfLParam = T(0.0);
      }
    }
  }

  if ( pfLParam ) {
    // undo reflections
    for (i = 0; i < 3; i++) {
      if ( bReflect[i] )
	kPnt[i] = -kPnt[i];
    }

    *pfBParam0 = kPnt[0];
    *pfBParam1 = kPnt[1];
    *pfBParam2 = kPnt[2];
  }

  return fSqrDistance;
}

template <class T>
T 
GbDistance<T>::distance (const GbLine<T>& rkLine, const GbBox3<T>& rkBox,
			 T* pfLParam, T* pfBParam0, T* pfBParam1, T* pfBParam2)
{
  return GbMath<T>::Sqrt(GbDistance<T>::squareDistance(rkLine,rkBox,
						       pfLParam,pfBParam0,pfBParam1,pfBParam2));
}

template <class T>
T 
GbDistance<T>::squareDistance (const GbLine<T>& rkLine, const GoTriangleElement& rkTri,
			       T* pfLinP, T* pfTriP0, T* pfTriP1)
{
  GbVec3<T> kDiff = rkTri.getOrigin() - rkLine.getOrigin();
  T fA00 = rkLine.getDirection().getSquareNorm();
  T fA01 = -rkLine.getDirection().dot(rkTri.getEdgeVec(0));
  T fA02 = -rkLine.getDirection().dot(rkTri.getEdgeVec(1));
  T fA11 = rkTri.getEdgeVec(0).getSquareNorm();
  T fA12 = rkTri.getEdgeVec(0).dot(rkTri.getEdgeVec(1));
  T fA22 = rkTri.getEdgeVec(1).dot(rkTri.getEdgeVec(1));
  T fB0  = -kDiff.dot(rkLine.getDirection());
  T fB1  = kDiff.dot(rkTri.getEdgeVec(0));
  T fB2  = kDiff.dot(rkTri.getEdgeVec(1));
  T fCof00 = fA11*fA22-fA12*fA12;
  T fCof01 = fA02*fA12-fA01*fA22;
  T fCof02 = fA01*fA12-fA02*fA11;
  T fDet = fA00*fCof00+fA01*fCof01+fA02*fCof02;

  GbSegment<T> kTriSeg;
  GbVec3<T> kPt;
  T fSqrDist, fSqrDist0, fR, fS, fT, fR0, fS0, fT0;

  if ( GbMath<T>::Abs(fDet) >= GbDistance<T>::TOLERANCE ) {
    T fCof11 = fA00*fA22-fA02*fA02;
    T fCof12 = fA02*fA01-fA00*fA12;
    T fCof22 = fA00*fA11-fA01*fA01;
    T fInvDet = T(1.0)/fDet;
    T fRhs0 = -fB0*fInvDet;
    T fRhs1 = -fB1*fInvDet;
    T fRhs2 = -fB2*fInvDet;

    fR = fCof00*fRhs0+fCof01*fRhs1+fCof02*fRhs2;
    fS = fCof01*fRhs0+fCof11*fRhs1+fCof12*fRhs2;
    fT = fCof02*fRhs0+fCof12*fRhs1+fCof22*fRhs2;

    if ( fS+fT <= T(1.0) ) {
      if ( fS < T(0.0) ) {
	if ( fT < T(0.0) ) {  // region 4
	  // min on face s=0 or t=0
	  kTriSeg.setOrigin(rkTri.getOrigin());
	  kTriSeg.setDirection(rkTri.getEdgeVec(1));
	  fSqrDist = GbDistance<T>::squareDistance(rkLine,kTriSeg,0,&fT);
	  fS = T(0.0);
	  kTriSeg.setOrigin(rkTri.getOrigin());
	  kTriSeg.setDirection(rkTri.getEdgeVec(0));
	  fSqrDist0 = GbDistance<T>::squareDistance(rkLine,kTriSeg,0,&fS0);
	  fT0 = T(0.0);
	  if ( fSqrDist0 < fSqrDist ) {
	    fSqrDist = fSqrDist0;
	    fS = fS0;
	    fT = fT0;
	  }
	}
	else { // region 3
	  // min on face s=0
	  kTriSeg.setOrigin(rkTri.getOrigin());
	  kTriSeg.setDirection(rkTri.getEdgeVec(1));
	  fSqrDist = GbDistance<T>::squareDistance(rkLine,kTriSeg,0,&fT);
	  fS = T(0.0);
	}
      }
      else if ( fT < T(0.0) ) { // region 5
	// min on face t=0
	kTriSeg.setOrigin(rkTri.getOrigin());
	kTriSeg.setDirection(rkTri.getEdgeVec(0));
	fSqrDist = GbDistance<T>::squareDistance(rkLine,kTriSeg,0,&fS);
	fT = T(0.0);
      }
      else { // region 0
	// line intersects triangle
	fSqrDist = T(0.0);
      }
    }
    else {
      if ( fS < T(0.0) ) { // region 2
	// min on face s=0 or s+t=1
	kTriSeg.setOrigin(rkTri.getOrigin());
	kTriSeg.setDirection(rkTri.getEdgeVec(1));
	fSqrDist = GbDistance<T>::squareDistance(rkLine,kTriSeg,0,&fT);
	fS = T(0.0);
	kTriSeg.setOrigin(rkTri.getOrigin()+rkTri.getEdgeVec(0));
	kTriSeg.setDirection(rkTri.getEdgeVec(2));
	fSqrDist0 = GbDistance<T>::squareDistance(rkLine,kTriSeg,0,&fT0);
	fS0 = T(1.0)-fT0;
	if ( fSqrDist0 < fSqrDist ) {
	  fSqrDist = fSqrDist0;
	  fS = fS0;
	  fT = fT0;
	}
      }
      else if ( fT < T(0.0) ) {  // region 6
	// min on face t=0 or s+t=1
	kTriSeg.setOrigin(rkTri.getOrigin());
	kTriSeg.setDirection(rkTri.getEdgeVec(0));
	fSqrDist = GbDistance<T>::squareDistance(rkLine,kTriSeg,0,&fS);
	fT = T(0.0);
	kTriSeg.setOrigin(rkTri.getOrigin()+rkTri.getEdgeVec(0));
	kTriSeg.setDirection(rkTri.getEdgeVec(2));
	fSqrDist0 = GbDistance<T>::squareDistance(rkLine,kTriSeg,0,&fT0);
	fS0 = T(1.0)-fT0;
	if ( fSqrDist0 < fSqrDist ) {
	  fSqrDist = fSqrDist0;
	  fS = fS0;
	  fT = fT0;
	}
      }
      else { // region 1
	// min on face s+t=1
	kTriSeg.setOrigin(rkTri.getOrigin()+rkTri.getEdgeVec(0));
	kTriSeg.setDirection(rkTri.getEdgeVec(2));
	fSqrDist = GbDistance<T>::squareDistance(rkLine,kTriSeg,0,&fT);
	fS = T(1.0)-fT;
      }
    }
  }
  else {
    // line and triangle are parallel
    kTriSeg.setOrigin(rkTri.getOrigin());
    kTriSeg.setDirection(rkTri.getEdgeVec(0));
    fSqrDist = GbDistance<T>::squareDistance(rkLine,kTriSeg,&fR,&fS);
    fT = T(0.0);

    kTriSeg.setDirection(rkTri.getEdgeVec(1));
    fSqrDist0 = GbDistance<T>::squareDistance(rkLine,kTriSeg,&fR0,&fT0);
    fS0 = T(0.0);
    if ( fSqrDist0 < fSqrDist ) {
      fSqrDist = fSqrDist0;
      fR = fR0;
      fS = fS0;
      fT = fT0;
    }

    kTriSeg.setOrigin(rkTri.getOrigin() + rkTri.getEdgeVec(0));
    kTriSeg.setDirection(rkTri.getEdgeVec(2));
    fSqrDist0 = GbDistance<T>::squareDistance(rkLine,kTriSeg,&fR0,&fT0);
    fS0 = T(1.0)-fT0;
    if ( fSqrDist0 < fSqrDist ) {
      fSqrDist = fSqrDist0;
      fR = fR0;
      fS = fS0;
      fT = fT0;
    }
  }

  if ( pfLinP ) *pfLinP = fR;

  if ( pfTriP0 ) *pfTriP0 = fS;

  if ( pfTriP1 ) *pfTriP1 = fT;

  return GbMath<T>::Abs(fSqrDist);
}

template <class T>
T 
GbDistance<T>::distance (const GbLine<T>& rkLine, const GoTriangleElement& rkTri,
			 T* pfLinP, T* pfTriP0, T* pfTriP1)
{
  return GbMath<T>::Sqrt(GbDistance<T>::squareDistance(rkLine,rkTri,
						       pfLinP,pfTriP0,pfTriP1));
}

template <class T>
T 
GbDistance<T>::distance (const GbVec3<T>& rkPoint, const GbVec3<T>& rkNormal,
			 const GbVec3<T>& rkOrigin, GbVec3<T>* pkClosest)
{
  GbVec3<T> kDiff = rkPoint - rkOrigin;
  T fSignedDist = rkNormal.dot(kDiff);

  if ( pkClosest )
    *pkClosest = rkPoint - fSignedDist*rkNormal;

  return GbMath<T>::Abs(fSignedDist);
}

template <class T>
T 
GbDistance<T>::distance (const GbVec3<T>& rkPoint, const GbVec3<T>& rkNormal,
			 T fConstant, GbVec3<T>* pkClosest)
{
  T fSignedDist = rkNormal.dot(rkPoint) - fConstant;

  if ( pkClosest )
    *pkClosest = rkPoint - fSignedDist*rkNormal;

  return GbMath<T>::Abs(fSignedDist);
}


// The calculations assume that the frustum is a solid, so any point inside
// the frustum has distance zero.
template <class T>
T
GbDistance<T>::squareDistance (const GbVec3<T>& rkPoint, const GbFrustum<T>& rkFrustum,
			       GbVec3<T>* pkClosest)
{
  // compute coordinates of point with respect to frustum coordinate system
  GbVec3<T> kDiff = rkPoint - rkFrustum.origin();
  GbVec3<T> kTest = GbVec3<T>(kDiff.dot(rkFrustum.lVector()),
			      kDiff.dot(rkFrustum.uVector()),
			      kDiff.dot(rkFrustum.dVector()));

  // perform calculations in octant with nonnegative L and U coordinates
  GbBool bLSignChange;
  if ( kTest[0] < T(0.0) ) {
    bLSignChange = true;
    kTest[0] = -kTest[0];
  }
  else {
    bLSignChange = false;
  }

  GbBool bUSignChange;
  if ( kTest[1] < T(0.0) ) {
    bUSignChange = true;
    kTest[1] = -kTest[1];
  }
  else {
    bUSignChange = false;
  }

  // frustum derived parameters
  T fLMin = rkFrustum.lBound(), fLMax = rkFrustum.getDRatio()*fLMin;
  T fUMin = rkFrustum.uBound(), fUMax = rkFrustum.getDRatio()*fUMin;
  T fDMin = rkFrustum.dMin(), fDMax = rkFrustum.dMax();
  T fLMinSqr = fLMin*fLMin;
  T fUMinSqr = fUMin*fUMin;
  T fDMinSqr = fDMin*fDMin;
  T fMinLDDot = fLMinSqr + fDMinSqr;
  T fMinUDDot = fUMinSqr + fDMinSqr;
  T fMinLUDDot = fLMinSqr + fMinUDDot;
  T fMaxLDDot = rkFrustum.getDRatio()*fMinLDDot;
  T fMaxUDDot = rkFrustum.getDRatio()*fMinUDDot;
  T fMaxLUDDot = rkFrustum.getDRatio()*fMinLUDDot;

  // Algorithm computes closest point in all cases by determining in which
  // Voronoi region of the vertices, edges, and faces of the frustum that
  // the test point lives.
  GbVec3<T> kClosest;
  T fLDot, fUDot, fLDDot, fUDDot, fLUDDot, fLEdgeDot, fUEdgeDot, fT;

  if ( kTest[2] >= fDMax ) {
    if ( kTest[0] <= fLMax ) {
      if ( kTest[1] <= fUMax ) {
	// F-face
	kClosest[0] = kTest[0];
	kClosest[1] = kTest[1];
	kClosest[2] = fDMax;
      }
      else {
	// UF-edge
	kClosest[0] = kTest[0];
	kClosest[1] = fUMax;
	kClosest[2] = fDMax;
      }
    }
    else {
      if ( kTest[1] <= fUMax ) {
	// LF-edge
	kClosest[0] = fLMax;
	kClosest[1] = kTest[1];
	kClosest[2] = fDMax;
      }
      else {
	// LUF-vertex
	kClosest[0] = fLMax;
	kClosest[1] = fUMax;
	kClosest[2] = fDMax;
      }
    }
  }
  else if ( kTest[2] <= fDMin ) {
    if ( kTest[0] <= fLMin ) {
      if ( kTest[1] <= fUMin ) {
	// N-face
	kClosest[0] = kTest[0];
	kClosest[1] = kTest[1];
	kClosest[2] = fDMin;
      }
      else {
	fUDDot = fUMin*kTest[1] + fDMin*kTest[2];
	if ( fUDDot >= fMaxUDDot ) {
	  // UF-edge
	  kClosest[0] = kTest[0];
	  kClosest[1] = fUMax;
	  kClosest[2] = fDMax;
	}
	else if ( fUDDot >= fMinUDDot ) {
	  // U-face
	  fUDot = fDMin*kTest[1] - fUMin*kTest[2];
	  fT = fUDot/fMinUDDot;
	  kClosest[0] = kTest[0];
	  kClosest[1] = kTest[1] - fT*fDMin;
	  kClosest[2] = kTest[2] + fT*fUMin;
	}
	else {
	  // UN-edge
	  kClosest[0] = kTest[0];
	  kClosest[1] = fUMin;
	  kClosest[2] = fDMin;
	}
      }
    }
    else {
      if ( kTest[1] <= fUMin ) {
	fLDDot = fLMin*kTest[0] + fDMin*kTest[2];
	if ( fLDDot >= fMaxLDDot ) {
	  // LF-edge
	  kClosest[0] = fLMax;
	  kClosest[1] = kTest[1];
	  kClosest[2] = fDMax;
	}
	else if ( fLDDot >= fMinLDDot ) {
	  // L-face
	  fLDot = fDMin*kTest[0] - fLMin*kTest[2];
	  fT = fLDot/fMinLDDot;
	  kClosest[0] = kTest[0] - fT*fDMin;
	  kClosest[1] = kTest[1];
	  kClosest[2] = kTest[2] + fT*fLMin;
	}
	else {
	  // LN-edge
	  kClosest[0] = fLMin;
	  kClosest[1] = kTest[1];
	  kClosest[2] = fDMin;
	}
      }
      else {
	fLUDDot = fLMin*kTest[0] + fUMin*kTest[1] + fDMin*kTest[2];
	fLEdgeDot = fUMin*fLUDDot - fMinLUDDot*kTest[1];
	if ( fLEdgeDot >= T(0.0) ) {
	  fLDDot = fLMin*kTest[0] + fDMin*kTest[2];
	  if ( fLDDot >= fMaxLDDot ) {
	    // LF-edge
	    kClosest[0] = fLMax;
	    kClosest[1] = kTest[1];
	    kClosest[2] = fDMax;
	  }
	  else if ( fLDDot >= fMinLDDot ) {
	    // L-face
	    fLDot = fDMin*kTest[0] - fLMin*kTest[2];
	    fT = fLDot/fMinLDDot;
	    kClosest[0] = kTest[0] - fT*fDMin;
	    kClosest[1] = kTest[1];
	    kClosest[2] = kTest[2] + fT*fLMin;
	  }
	  else {
	    // LN-edge
	    kClosest[0] = fLMin;
	    kClosest[1] = kTest[1];
	    kClosest[2] = fDMin;
	  }
	}
	else {
	  fUEdgeDot = fLMin*fLUDDot - fMinLUDDot*kTest[0];
	  if ( fUEdgeDot >= T(0.0) ) {
	    fUDDot = fUMin*kTest[1] + fDMin*kTest[2];
	    if ( fUDDot >= fMaxUDDot ) {
	      // UF-edge
	      kClosest[0] = kTest[0];
	      kClosest[1] = fUMax;
	      kClosest[2] = fDMax;
	    }
	    else if ( fUDDot >= fMinUDDot ) {
	      // U-face
	      fUDot = fDMin*kTest[1] - fUMin*kTest[2];
	      fT = fUDot/fMinUDDot;
	      kClosest[0] = kTest[0];
	      kClosest[1] = kTest[1] - fT*fDMin;
	      kClosest[2] = kTest[2] + fT*fUMin;
	    }
	    else {
	      // UN-edge
	      kClosest[0] = kTest[0];
	      kClosest[1] = fUMin;
	      kClosest[2] = fDMin;
	    }
	  }
	  else {
	    if ( fLUDDot >= fMaxLUDDot ) {
	      // LUF-vertex
	      kClosest[0] = fLMax;
	      kClosest[1] = fUMax;
	      kClosest[2] = fDMax;
	    }
	    else if ( fLUDDot >= fMinLUDDot ) {
	      // LU-edge
	      fT = fLUDDot/fMinLUDDot;
	      kClosest[0] = fT*fLMin;
	      kClosest[1] = fT*fUMin;
	      kClosest[2] = fT*fDMin;
	    }
	    else {
	      // LUN-vertex
	      kClosest[0] = fLMin;
	      kClosest[1] = fUMin;
	      kClosest[2] = fDMin;
	    }
	  }
	}
      }
    }
  }
  else {
    fLDot = fDMin*kTest[0] - fLMin*kTest[2];
    fUDot = fDMin*kTest[1] - fUMin*kTest[2];
    if ( fLDot <= T(0.0) ) {
      if ( fUDot <= T(0.0) ) {
	// point inside frustum
	kClosest = kTest;
      }
      else {
	fUDDot = fUMin*kTest[1] + fDMin*kTest[2];
	if ( fUDDot >= fMaxUDDot ) {
	  // UF-edge
	  kClosest[0] = kTest[0];
	  kClosest[1] = fUMax;
	  kClosest[2] = fDMax;
	}
	else {
	  // U-face
	  fUDot = fDMin*kTest[1] - fUMin*kTest[2];
	  fT = fUDot/fMinUDDot;
	  kClosest[0] = kTest[0];
	  kClosest[1] = kTest[1] - fT*fDMin;
	  kClosest[2] = kTest[2] + fT*fUMin;
	}
      }
    }
    else {
      if ( fUDot <= T(0.0) ) {
	fLDDot = fLMin*kTest[0] + fDMin*kTest[2];
	if ( fLDDot >= fMaxLDDot ) {
	  // LF-edge
	  kClosest[0] = fLMax;
	  kClosest[1] = kTest[1];
	  kClosest[2] = fDMax;
	}
	else {
	  // L-face
	  fLDot = fDMin*kTest[0] - fLMin*kTest[2];
	  fT = fLDot/fMinLDDot;
	  kClosest[0] = kTest[0] - fT*fDMin;
	  kClosest[1] = kTest[1];
	  kClosest[2] = kTest[2] + fT*fLMin;
	}
      }
      else {
	fLUDDot = fLMin*kTest[0] + fUMin*kTest[1] + fDMin*kTest[2];
	fLEdgeDot = fUMin*fLUDDot - fMinLUDDot*kTest[1];
	if ( fLEdgeDot >= T(0.0) ) {
	  fLDDot = fLMin*kTest[0] + fDMin*kTest[2];
	  if ( fLDDot >= fMaxLDDot ) {
	    // LF-edge
	    kClosest[0] = fLMax;
	    kClosest[1] = kTest[1];
	    kClosest[2] = fDMax;
	  }
	  else { // assert( fLDDot >= fMinLDDot ) from geometry
	    // L-face
	    fLDot = fDMin*kTest[0] - fLMin*kTest[2];
	    fT = fLDot/fMinLDDot;
	    kClosest[0] = kTest[0] - fT*fDMin;
	    kClosest[1] = kTest[1];
	    kClosest[2] = kTest[2] + fT*fLMin;
	  }
	}
	else {
	  fUEdgeDot = fLMin*fLUDDot - fMinLUDDot*kTest[0];
	  if ( fUEdgeDot >= T(0.0) ) {
	    fUDDot = fUMin*kTest[1] + fDMin*kTest[2];
	    if ( fUDDot >= fMaxUDDot ) {
	      // UF-edge
	      kClosest[0] = kTest[0];
	      kClosest[1] = fUMax;
	      kClosest[2] = fDMax;
	    }
	    else { // assert( fUDDot >= fMinUDDot ) from geometry
	      // U-face
	      fUDot = fDMin*kTest[1] - fUMin*kTest[2];
	      fT = fUDot/fMinUDDot;
	      kClosest[0] = kTest[0];
	      kClosest[1] = kTest[1] - fT*fDMin;
	      kClosest[2] = kTest[2] + fT*fUMin;
	    }
	  }
	  else {
	    if ( fLUDDot >= fMaxLUDDot ) {
	      // LUF-vertex
	      kClosest[0] = fLMax;
	      kClosest[1] = fUMax;
	      kClosest[2] = fDMax;
	    }
	    else { // assert( fLUDDot >= fMinLUDDot ) from geometry
	      // LU-edge
	      fT = fLUDDot/fMinLUDDot;
	      kClosest[0] = fT*fLMin;
	      kClosest[1] = fT*fUMin;
	      kClosest[2] = fT*fDMin;
	    }
	  }
	}
      }
    }
  }

  kDiff = kTest - kClosest;

  // convert back to original quadrant
  if ( bLSignChange )
    kClosest[0] = -kClosest[0];

  if ( bUSignChange )
    kClosest[1] = -kClosest[1];

  if ( pkClosest ) {
    // caller wants closest point, convert back to world coordinates
    *pkClosest = rkFrustum.origin() + 
                 kClosest[0]*rkFrustum.lVector() +
                 kClosest[1]*rkFrustum.uVector() + 
                 kClosest[2]*rkFrustum.dVector();
  }

  // compute and return squared distance
  return kDiff.getSquareNorm();
}

// The calculations assume that the frustum is a solid, so any point inside
// the frustum has distance zero.
template <class T>
T
GbDistance<T>::distance (const GbVec3<T>& rkPoint, const GbFrustum<T>& rkFrustum,
			 GbVec3<T>* pkClosest)
{
  return GbMath<T>::Sqrt(GbDistance<T>::squareDistance(rkPoint,rkFrustum,pkClosest));
}
