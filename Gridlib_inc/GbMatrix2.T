// -*- C++ -*- automatisch in C++mode wechseln (emacs)
/*----------------------------------------------------------------------
|
|
| $Log$
| Revision 1.3  2002/04/03 11:17:07  elena
| new: changes in CFS++ on SGI and new function for gridlib
|
| Revision 1.2  2001/08/16 16:53:17  prkipfer
| improved type safety for template parameter
|
| Revision 1.1  2001/06/26 08:31:22  prkipfer
| added 2D Projections
|
|
|
+---------------------------------------------------------------------*/

#ifdef LINUX
#include <cassert>
#endif

template <class T>
const T 
GbMatrix2<T>::EPSILON = std::numeric_limits<T>::epsilon(); // 1e-06f

template <class T>
GbMatrix2<T>::GbMatrix2()
{
  // no init, use constructor below if needed
}

template <class T>
GbMatrix2<T>::GbMatrix2(const T aafEntry[2][2])
{
  memcpy(entry_,aafEntry,4*sizeof(T));
}

template <class T>
GbMatrix2<T>::GbMatrix2(const T aafEntry[4])
{
  memcpy(entry_,aafEntry,4*sizeof(T));
}

template <class T>
GbMatrix2<T>::GbMatrix2(const GbMatrix2<T>& m)
{
  entry_[0][0] = m[0][0];
  entry_[0][1] = m[0][1];
  entry_[1][0] = m[1][0];
  entry_[1][1] = m[1][1];
}

template <class T>
GbMatrix2<T>::GbMatrix2(T fEntry00, T fEntry01, 
			T fEntry10, T fEntry11)
{
  entry_[0][0] = fEntry00;
  entry_[0][1] = fEntry01;
  entry_[1][0] = fEntry10;
  entry_[1][1] = fEntry11;
}

template <class T>
GbMatrix2<T>::GbMatrix2(T s)
{
  entry_[0][0] = entry_[0][1] =
  entry_[1][0] = entry_[1][1] = s;
}

template <class T>
GbMatrix2<T>::~GbMatrix2()
{
}

template <class T>
GbMatrix2<T> 
GbMatrix2<T>::transpose () const
{
  T kTranspose[4];
  int k=0;

  for (int iRow = 0; iRow < 2; ++iRow) {
    for (int iCol = 0; iCol < 2; ++iCol)
      kTranspose[k++] = entry_[iCol][iRow];
  }
  return GbMatrix2<T>(kTranspose);
}

template <class T>
GbBool 
GbMatrix2<T>::inverse (GbMatrix2<T>& rkInverse, T fTolerance) const
{
  T fDet = entry_[0][0]*entry_[1][1] - entry_[0][1]*entry_[1][0];

  if ( GbMath<T>::Abs(fDet) <= fTolerance )
    return false;

  T fInvDet = T(1.0)/fDet;
  rkInverse[0][0] =  entry_[1][1]*fInvDet;
  rkInverse[1][0] = -entry_[1][0]*fInvDet;
  rkInverse[0][1] = -entry_[0][1]*fInvDet;
  rkInverse[1][1] =  entry_[0][0]*fInvDet;

  return true;
}

template <class T>
GbMatrix2<T> 
GbMatrix2<T>::inverse (T fTolerance) const
{
  GbMatrix2<T> fInverse = GbMatrix2<T>::ZERO;
  inverse(fInverse,fTolerance);
  return fInverse;
}

template <class T>
T 
GbMatrix2<T>::determinant () const
{
  return entry_[0][0]*entry_[1][1] - entry_[0][1]*entry_[1][0];
}

template <class T>
void 
GbMatrix2<T>::singularValueDecomposition (GbMatrix2<T>& rkL,
					  GbVec2<T>& rkS, GbMatrix2<T>& rkR) const
{
  // compute M*M^t and M^t*M where M^t = Transpose(M)
  T fM00Sqr = entry_[0][0]*entry_[0][0];
  T fM01Sqr = entry_[0][1]*entry_[0][1];
  T fM10Sqr = entry_[1][0]*entry_[1][0];
  T fM11Sqr = entry_[1][1]*entry_[1][1];

  T afMMT[3], afMTM[3];
  afMMT[0] = fM00Sqr + fM01Sqr;
  afMMT[1] = entry_[0][0]*entry_[1][0] + entry_[0][1]*entry_[1][1];
  afMMT[2] = fM10Sqr + fM11Sqr;
  afMTM[0] = fM00Sqr + fM10Sqr;
  afMTM[1] = entry_[0][0]*entry_[0][1] + entry_[1][0]*entry_[1][1];
  afMTM[2] = fM01Sqr + fM11Sqr;

  // factor M*M^t = L*S^2*L^t
  T fCos, fSin, fInvLen;
  if ( GbMath<T>::Abs(afMMT[1]) > T(0.0) ) {
    T fTrace = afMMT[0] + afMMT[2];
    T fDet = afMMT[0]*afMMT[2] - afMMT[1]*afMMT[1];
    T fDiscr = fTrace*fTrace-T(4.0)*fDet;
    T fRoot = GbMath<T>::Sqrt(GbMath<T>::Abs(fDiscr));
    rkS[0] = T(0.5)*(fTrace + fRoot);
    rkS[1] = T(0.5)*(fTrace - fRoot);

    fCos = afMMT[1];
    fSin = rkS[0] - afMMT[0];
    fInvLen = T(1.0)/GbMath<T>::Sqrt(fCos*fCos + fSin*fSin);
    fCos *= fInvLen;
    fSin *= fInvLen;

    rkL[0][0] = fCos;
    rkL[1][0] = fSin;
    rkL[0][1] = -fSin;
    rkL[1][1] = fCos;
  }
  else {
    rkS[0] = afMMT[0];
    rkS[1] = afMMT[2];
    rkL = GbMatrix2<T>::IDENTITY;
  }

  // factor M^t*M = R^t*S^2*R
  if ( GbMath<T>::Abs(afMTM[1]) > T(0.0) ) {
    fCos = afMTM[1];
    fSin = rkS[0] - afMTM[0];
    fInvLen = T(1.0)/GbMath<T>::Sqrt(fCos*fCos + fSin*fSin);
    fCos *= fInvLen;
    fSin *= fInvLen;

    rkR[0][0] = fCos;
    rkR[1][0] = -fSin;
    rkR[0][1] = +fSin;
    rkR[1][1] = fCos;
  }
  else {
    rkR = GbMatrix2<T>::IDENTITY;
  }

  rkS[0] = GbMath<T>::Sqrt(GbMath<T>::Abs(rkS[0]));
  rkS[1] = GbMath<T>::Sqrt(GbMath<T>::Abs(rkS[1]));

  GbVec2<T> kDiag = rkS;

  GbMatrix2<T> kTest;
  kTest.singularValueComposition(rkL,kDiag,rkR);
  kTest = *this - kTest;
  int i = 0;
  T fMinNormSqr = kTest.L2NormSquare();

  kDiag[0] = -kDiag[0];
  kTest.singularValueComposition(rkL,kDiag,rkR);
  kTest = *this - kTest;
  T fNormSqr = kTest.L2NormSquare();
  if ( fNormSqr < fMinNormSqr ) {
    i = 1;
    fMinNormSqr = fNormSqr;
  }

  kDiag[1] = -kDiag[1];
  kTest.singularValueComposition(rkL,kDiag,rkR);
  kTest = *this - kTest;
  fNormSqr = kTest.L2NormSquare();
  if ( fNormSqr < fMinNormSqr ) {
    i = 2;
    fMinNormSqr = fNormSqr;
  }

  kDiag[0] = -kDiag[0];
  kTest.singularValueComposition(rkL,kDiag,rkR);
  kTest = *this - kTest;
  fNormSqr = kTest.L2NormSquare();
  if ( fNormSqr < fMinNormSqr ) {
    i = 3;
    fMinNormSqr = fNormSqr;
  }

  switch ( i ) {
    case 1:
      rkR[0][0] = -rkR[0][0];
      rkR[0][1] = -rkR[0][1];
      break;
    case 2:
      rkR[0][0] = -rkR[0][0];
      rkR[0][1] = -rkR[0][1];
      rkR[1][0] = -rkR[1][0];
      rkR[1][1] = -rkR[1][1];
      break;
    case 3:
      rkR[1][0] = -rkR[1][0];
      rkR[1][1] = -rkR[1][1];
      break;
  }
}

template <class T>
void 
GbMatrix2<T>::singularValueComposition (const GbMatrix2<T>& rkL,
					const GbVec2<T>& rkS, const GbMatrix2<T>& rkR)
{
  int iRow, iCol;
  GbMatrix2<T> kTmp;

  // product S*R
  for (iRow = 0; iRow < 2; ++iRow) {
    for (iCol = 0; iCol < 2; ++iCol)
      kTmp[iRow][iCol] = rkS[iRow]*rkR[iRow][iCol];
  }

  // product L*S*R
  for (iRow = 0; iRow < 2; ++iRow) {
    for (iCol = 0; iCol < 2; ++iCol) {
      entry_[iRow][iCol] = T(0.0);
      for (int iMid = 0; iMid < 2; ++iMid)
	entry_[iRow][iCol] += rkL[iRow][iMid]*kTmp[iMid][iCol];
    }
  }
}

template <class T>
void 
GbMatrix2<T>::orthonormalize ()
{
  // Algorithm uses Gram-Schmidt orthogonalization.  If 'this' matrix is
  // M = [m0|m1], then orthonormal output matrix is Q = [q0|q1],
  //
  //   q0 = m0/|m0|
  //   q1 = (m1-(q0*m1)q0)/|m1-(q0*m1)q0|
  //
  // where |V| indicates length of vector V and A*B indicates dot
  // product of vectors A and B.

  // compute q0
  T fInvLength = T(1.0)/GbMath<T>::Sqrt(entry_[0][0]*entry_[0][0] + entry_[1][0]*entry_[1][0]);

  entry_[0][0] *= fInvLength;
  entry_[1][0] *= fInvLength;

  // compute q1
  T fDot0 = entry_[0][0]*entry_[0][1] + entry_[1][0]*entry_[1][1];

  entry_[0][1] -= fDot0*entry_[0][0];
  entry_[1][1] -= fDot0*entry_[1][0];

  fInvLength = T(1.0)/GbMath<T>::Sqrt(entry_[0][1]*entry_[0][1] + entry_[1][1]*entry_[1][1]);

  entry_[0][1] *= fInvLength;
  entry_[1][1] *= fInvLength;
}

template <class T>
void 
GbMatrix2<T>::QDUDecomposition (GbMatrix2<T>& rkQ, GbVec2<T>& rkD, T& rkU) const
{
  // Factor M = QR = QDU where Q is orthogonal, D is diagonal, 
  // and U is upper triangular with ones on its diagonal.  Algorithm uses
  // Gram-Schmidt orthogonalization (the QR algorithm).
  //
  // If M = [ m0 | m1 ] and Q = [ q0 | q1 ], then
  //
  //   q0 = m0/|m0| = (a,b)
  //   q1 = (-b,a)
  //
  // where |V| indicates length of vector V and A*B indicates dot
  // product of vectors A and B.  The matrix R has entries
  //
  //   r00 = q0*m0  r01 = q0*m1
  //   r10 = 0      r11 = q1*m1
  //
  // so D = diag(r00,r11) and U has entry u01 = r01/r00.

  // Q = rotation
  // D = scaling
  // U = shear

  // D stores the two diagonal entries r00, r11
  // U stores the entries u01

  // build orthogonal matrix Q
  T fInvLength = T(1.0)/GbMath<T>::Sqrt(entry_[0][0]*entry_[0][0]);
  rkQ[0][0] = entry_[0][0]*fInvLength;
  rkQ[1][0] = entry_[1][0]*fInvLength;

  T fDot = rkQ[0][0]*entry_[0][1] + rkQ[1][0]*entry_[1][1];
  rkQ[0][1] = entry_[0][1]-fDot*rkQ[0][0];
  rkQ[1][1] = entry_[1][1]-fDot*rkQ[1][0];
  fInvLength = T(1.0)/GbMath<T>::Sqrt(rkQ[0][1]*rkQ[0][1] + rkQ[1][1]*rkQ[1][1]);
  rkQ[0][1] *= fInvLength;
  rkQ[1][1] *= fInvLength;

  // guarantee that orthogonal matrix has determinant 1 (no reflections)
  T fDet = rkQ[0][0]*rkQ[1][1] - rkQ[0][1]*rkQ[1][0];
  if ( fDet < T(0.0) ) {
    rkQ[0][1] = -rkQ[0][1];
    rkQ[1][1] = -rkQ[1][1];
  }

  // build "right" matrix R
  GbMatrix2<T> kR;
  kR[0][0] = rkQ[0][0]*entry_[0][0] + rkQ[1][0]*entry_[1][0];
  kR[0][1] = rkQ[0][0]*entry_[0][1] + rkQ[1][0]*entry_[1][1];
  kR[1][0] = T(0.0);
  kR[1][1] = rkQ[0][1]*entry_[0][1] + rkQ[1][1]*entry_[1][1];

  // the scaling component
  rkD[0] = kR[0][0];
  rkD[1] = kR[1][1];

  // the shear component
  rkU = kR[0][1]/rkD[0];
}

template <class T>
T 
GbMatrix2<T>::spectralNorm () const
{
  GbMatrix2<T> kP;
  int iRow, iCol;
  T fPmax = T(0.0);
  for (iRow = 0; iRow < 2; ++iRow) {
    for (iCol = 0; iCol < 2; ++iCol) {
      kP[iRow][iCol] = T(0.0);
      for (int iMid = 0; iMid < 2; ++iMid) {
	kP[iRow][iCol] += entry_[iMid][iRow]*entry_[iMid][iCol];
      }
      if ( kP[iRow][iCol] > fPmax )
	fPmax = kP[iRow][iCol];
    }
  }

  T fInvPmax = T(1.0)/fPmax;
  for (iRow = 0; iRow < 2; ++iRow) {
    for (iCol = 0; iCol < 3; ++iCol)
      kP[iRow][iCol] *= fInvPmax;
  }

  T fDet = kP[0][0]*kP[1][1]-kP[0][1]*kP[1][0];
  T fTrace = -(kP[0][0]+kP[1][1]);
  T fRoot = T(0.5)*(-fTrace + GbMath<T>::Sqrt(GbMath<T>::Abs(fTrace*fTrace-T(4.0)*fDet)));
  T fNorm = GbMath<T>::Sqrt(fPmax*fRoot);
  return fNorm;
}

template <class T>
T 
GbMatrix2<T>::L2NormSquare () const
{
  T fNormSqr = T(0.0);
  for (int iRow = 0; iRow < 2; ++iRow) {
    for (int iCol = 0; iCol < 2; ++iCol) {
      T fSqr = entry_[iRow][iCol]*entry_[iRow][iCol];
      fNormSqr += fSqr;
    }
  }
  return fNormSqr;
}

template <class T>
void 
GbMatrix2<T>::tridiagonal (T afDiag[2], T afSubDiag[2])
{
  // matrix is already tridiagonal, repackage for QL algorithm
  afDiag[0] = entry_[0][0];
  afDiag[1] = entry_[1][1];
  afSubDiag[0] = entry_[0][1];
  afSubDiag[1] = T(0.0);
  entry_[0][0] = T(1.0);
  entry_[0][1] = T(0.0);
  entry_[1][0] = T(0.0);
  entry_[1][1] = T(1.0);
}

template <class T>
GbBool 
GbMatrix2<T>::QLAlgorithm (T afDiag[2], T afSubDiag[2])
{
  for (int i0 = 0; i0 < 2; i0++) {
    const int iMaxIter = 32;
    int iIter;
    for (iIter = 0; iIter < iMaxIter; iIter++) {
      int i1;
      for (i1 = i0; i1 <= 0; i1++) {
	T fSum = GbMath<T>::Abs(afDiag[i1]) + GbMath<T>::Abs(afDiag[i1+1]);
	if ( GbMath<T>::Abs(afSubDiag[i1]) + fSum == fSum )
	  break;
      }
      if ( i1 == i0 )
	break;

      T fTmp0 = (afDiag[i0+1]-afDiag[i0])/(T(2.0)*afSubDiag[i0]);
      T fTmp1 = GbMath<T>::Sqrt(fTmp0*fTmp0 + T(1.0));
      if ( fTmp0 < T(0.0) )
	fTmp0 = afDiag[i1]-afDiag[i0]+afSubDiag[i0]/(fTmp0-fTmp1);
      else
	fTmp0 = afDiag[i1]-afDiag[i0]+afSubDiag[i0]/(fTmp0+fTmp1);
      T fSin = T(1.0);
      T fCos = T(1.0);
      T fTmp2 = T(0.0);
      for (int i2 = i1-1; i2 >= i0; i2--) {
	T fTmp3 = fSin*afSubDiag[i2];
	T fTmp4 = fCos*afSubDiag[i2];
	if ( GbMath<T>::Abs(fTmp3) >= GbMath<T>::Abs(fTmp0) ) {
	  fCos = fTmp0/fTmp3;
	  fTmp1 = GbMath<T>::Sqrt(fCos*fCos + T(1.0));
	  afSubDiag[i2+1] = fTmp3*fTmp1;
	  fSin = T(1.0)/fTmp1;
	  fCos *= fSin;
	}
	else {
	  fSin = fTmp3/fTmp0;
	  fTmp1 = GbMath<T>::Sqrt(fSin*fSin+T(1.0));
	  afSubDiag[i2+1] = fTmp0*fTmp1;
	  fCos = T(1.0)/fTmp1;
	  fSin *= fCos;
	}
	fTmp0 = afDiag[i2+1]-fTmp2;
	fTmp1 = (afDiag[i2]-fTmp0)*fSin + T(2.0)*fTmp4*fCos;
	fTmp2 = fSin*fTmp1;
	afDiag[i2+1] = fTmp0+fTmp2;
	fTmp0 = fCos*fTmp1-fTmp4;

	for (int iRow = 0; iRow < 2; ++iRow) {
	  fTmp3 = entry_[iRow][i2+1];
	  entry_[iRow][i2+1] = fSin*entry_[iRow][i2] + fCos*fTmp3;
	  entry_[iRow][i2] =   fCos*entry_[iRow][i2] - fSin*fTmp3;
	}
      }
      afDiag[i0] -= fTmp2;
      afSubDiag[i0] = fTmp0;
      afSubDiag[i1] = T(0.0);
    }

    if ( iIter == iMaxIter ) {
      // should not get here under normal circumstances
      return false;
    }
  }

  return true;
}

template <class T>
void 
GbMatrix2<T>::eigenSolveSymmetric (T afEigenvalue[2],
				   GbVec2<T> akEigenvector[2]) const
{
  GbMatrix2<T> kMatrix = *this;
  T afSubDiag[3];
  kMatrix.tridiagonal(afEigenvalue,afSubDiag);
  kMatrix.QLAlgorithm(afEigenvalue,afSubDiag);

  for (int i = 0; i < 2; i++) {
    akEigenvector[i][0] = kMatrix[0][i];
    akEigenvector[i][1] = kMatrix[1][i];
  }

  // make eigenvectors form a right--handed system
  T fDet = akEigenvector[0].dot(akEigenvector[1]);
  if ( fDet < T(0.0) ) {
    akEigenvector[1][0] = - akEigenvector[1][0];
    akEigenvector[1][1] = - akEigenvector[1][1];
  }
}

template <class T>
void 
GbMatrix2<T>::tensorProduct (const GbVec2<T>& rkU, const GbVec2<T>& rkV,
			     GbMatrix2<T>& rkProduct)
{
  for (int iRow = 0; iRow < 2; ++iRow) {
    for (int iCol = 0; iCol < 2; ++iCol)
      rkProduct[iRow][iCol] = rkU[iRow]*rkV[iCol];
  }
}
