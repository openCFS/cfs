// -*- C++ -*- automatisch in C++mode wechseln (emacs)
/*----------------------------------------------------------------------
|
|
| $Log$
| Revision 1.1  2002/02/22 14:47:56  elena
| new: dir Gridlib_inc
|
| Revision 1.3  2001/06/15 08:35:11  prkipfer
| added fast transpose multiplication and spherical interpolation without quaternions
|
| Revision 1.2  2001/03/20 09:38:52  prkipfer
| speedup for toEulerAngles routines
|
| Revision 1.1  2001/01/02 14:44:02  prkipfer
| introduced new classes
|
|
|
+---------------------------------------------------------------------*/


template <class T>
GbMatrix3<T>::GbMatrix3()
{
  entry_[0][0] = entry_[0][1] = entry_[0][2] = 
  entry_[1][0] = entry_[1][1] = entry_[1][2] = 
  entry_[2][0] = entry_[2][1] = entry_[2][2] = 0.0;
}

template <class T>
GbMatrix3<T>::GbMatrix3(const T entry[3][3])
{
  memcpy(entry_,entry,9*sizeof(T));
}

template <class T>
GbMatrix3<T>::GbMatrix3(const GbMatrix3<T>& m)
{
  entry_[0][0] = m[0][0];
  entry_[0][1] = m[0][1];
  entry_[0][2] = m[0][2];
  entry_[1][0] = m[1][0];
  entry_[1][1] = m[1][1];
  entry_[1][2] = m[1][2];
  entry_[2][0] = m[2][0];
  entry_[2][1] = m[2][1];
  entry_[2][2] = m[2][2];
}

template <class T>
GbMatrix3<T>::GbMatrix3(T entry00, T entry01, T entry02,
			T entry10, T entry11, T entry12, 
			T entry20, T entry21, T entry22)
{
  entry_[0][0] = entry00; entry_[0][1] = entry01; entry_[0][2] = entry02;
  entry_[1][0] = entry10; entry_[1][1] = entry11; entry_[1][2] = entry12;
  entry_[2][0] = entry20; entry_[2][1] = entry21; entry_[2][2] = entry22;
}

template <class T>
GbMatrix3<T>::~GbMatrix3()
{
}

template <class T>
GbMatrix3<T> 
GbMatrix3<T>::transpose () const
{
  GbMatrix3<T> kTranspose;
  for (int iRow = 0; iRow < 3; iRow++) {
    for (int iCol = 0; iCol < 3; iCol++)
      kTranspose[iRow][iCol] = entry_[iCol][iRow];
  }
  return kTranspose;
}

template <class T>
GbBool 
GbMatrix3<T>::inverse (GbMatrix3<T>& inv, T tolerance) const
{
  // Invert a 3x3 using cofactors.  This is about 8 times faster than
  // the Numerical Recipes code which uses Gaussian elimination.

  inv[0][0] = entry_[1][1]*entry_[2][2] - entry_[1][2]*entry_[2][1];
  inv[0][1] = entry_[0][2]*entry_[2][1] - entry_[0][1]*entry_[2][2];
  inv[0][2] = entry_[0][1]*entry_[1][2] - entry_[0][2]*entry_[1][1];
  inv[1][0] = entry_[1][2]*entry_[2][0] - entry_[1][0]*entry_[2][2];
  inv[1][1] = entry_[0][0]*entry_[2][2] - entry_[0][2]*entry_[2][0];
  inv[1][2] = entry_[0][2]*entry_[1][0] - entry_[0][0]*entry_[1][2];
  inv[2][0] = entry_[1][0]*entry_[2][1] - entry_[1][1]*entry_[2][0];
  inv[2][1] = entry_[0][1]*entry_[2][0] - entry_[0][0]*entry_[2][1];
  inv[2][2] = entry_[0][0]*entry_[1][1] - entry_[0][1]*entry_[1][0];

  T fDet = entry_[0][0]*inv[0][0] +
           entry_[0][1]*inv[1][0] +
           entry_[0][2]*inv[2][0];

  if ( GbMath<T>::Abs(fDet) <= tolerance ) return false;

  T fInvDet = 1.0/fDet;
  for (int iRow = 0; iRow < 3; iRow++) {
    for (int iCol = 0; iCol < 3; iCol++)
      inv[iRow][iCol] *= fInvDet;
  }

  return true;
}

template <class T>
GbMatrix3<T> 
GbMatrix3<T>::inverse (T tolerance) const
{
  GbMatrix3<T> kInverse = GbMatrix3<T>::ZERO;
  inverse(kInverse,tolerance);
  return kInverse;
}

template <class T>
T
GbMatrix3<T>::determinant () const
{
  T fCofactor00 = entry_[1][1]*entry_[2][2] - entry_[1][2]*entry_[2][1];
  T fCofactor10 = entry_[1][2]*entry_[2][0] - entry_[1][0]*entry_[2][2];
  T fCofactor20 = entry_[1][0]*entry_[2][1] - entry_[1][1]*entry_[2][0];

  T fDet = entry_[0][0]*fCofactor00 +
           entry_[0][1]*fCofactor10 +
           entry_[0][2]*fCofactor20;
  
  return fDet;
}

template <class T>
void 
GbMatrix3<T>::bidiagonalize (GbMatrix3<T>& kA, GbMatrix3<T>& kL, GbMatrix3<T>& kR)
{
  T afV[3], afW[3];
  T fLength, fSign, fT1, fInvT1, fT2;
  GbBool bIdentity;

  // map first column to (*,0,0)
  fLength = GbMath<T>::Sqrt(kA[0][0]*kA[0][0] + kA[1][0]*kA[1][0] + kA[2][0]*kA[2][0]);

  if ( fLength > 0.0 ) {
    fSign = (kA[0][0] > 0.0 ? 1.0 : -1.0);
    fT1 = kA[0][0] + fSign*fLength;
    fInvT1 = 1.0/fT1;
    afV[1] = kA[1][0]*fInvT1;
    afV[2] = kA[2][0]*fInvT1;

    fT2 = -2.0/(1.0+afV[1]*afV[1]+afV[2]*afV[2]);
    afW[0] = fT2*(kA[0][0]+kA[1][0]*afV[1]+kA[2][0]*afV[2]);
    afW[1] = fT2*(kA[0][1]+kA[1][1]*afV[1]+kA[2][1]*afV[2]);
    afW[2] = fT2*(kA[0][2]+kA[1][2]*afV[1]+kA[2][2]*afV[2]);
    kA[0][0] += afW[0];
    kA[0][1] += afW[1];
    kA[0][2] += afW[2];
    kA[1][1] += afV[1]*afW[1];
    kA[1][2] += afV[1]*afW[2];
    kA[2][1] += afV[2]*afW[1];
    kA[2][2] += afV[2]*afW[2];

    kL[0][0] = 1.0+fT2;
    kL[0][1] = kL[1][0] = fT2*afV[1];
    kL[0][2] = kL[2][0] = fT2*afV[2];
    kL[1][1] = 1.0+fT2*afV[1]*afV[1];
    kL[1][2] = kL[2][1] = fT2*afV[1]*afV[2];
    kL[2][2] = 1.0+fT2*afV[2]*afV[2];
    bIdentity = false;
  }
  else {
    kL = GbMatrix3<T>::IDENTITY;
    bIdentity = true;
  }

  // map first row to (*,*,0)
  fLength = GbMath<T>::Sqrt(kA[0][1]*kA[0][1]+kA[0][2]*kA[0][2]);

  if ( fLength > 0.0 ) {
    fSign = (kA[0][1] > 0.0 ? 1.0 : -1.0);
    fT1 = kA[0][1] + fSign*fLength;
    afV[2] = kA[0][2]/fT1;

    fT2 = -2.0/(1.0+afV[2]*afV[2]);
    afW[0] = fT2*(kA[0][1]+kA[0][2]*afV[2]);
    afW[1] = fT2*(kA[1][1]+kA[1][2]*afV[2]);
    afW[2] = fT2*(kA[2][1]+kA[2][2]*afV[2]);
    kA[0][1] += afW[0];
    kA[1][1] += afW[1];
    kA[1][2] += afW[1]*afV[2];
    kA[2][1] += afW[2];
    kA[2][2] += afW[2]*afV[2];
        
    kR[0][0] = 1.0;
    kR[0][1] = kR[1][0] = 0.0;
    kR[0][2] = kR[2][0] = 0.0;
    kR[1][1] = 1.0+fT2;
    kR[1][2] = kR[2][1] = fT2*afV[2];
    kR[2][2] = 1.0+fT2*afV[2]*afV[2];
  }
  else {
    kR = GbMatrix3<T>::IDENTITY;
  }

  // map second column to (*,*,0)
  fLength = GbMath<T>::Sqrt(kA[1][1]*kA[1][1]+kA[2][1]*kA[2][1]);

  if ( fLength > 0.0 ) {
    fSign = (kA[1][1] > 0.0 ? 1.0 : -1.0);
    fT1 = kA[1][1] + fSign*fLength;
    afV[2] = kA[2][1]/fT1;

    fT2 = -2.0/(1.0+afV[2]*afV[2]);
    afW[1] = fT2*(kA[1][1]+kA[2][1]*afV[2]);
    afW[2] = fT2*(kA[1][2]+kA[2][2]*afV[2]);
    kA[1][1] += afW[1];
    kA[1][2] += afW[2];
    kA[2][2] += afV[2]*afW[2];

    T fA = 1.0+fT2;
    T fB = fT2*afV[2];
    T fC = 1.0+fB*afV[2];

    if ( bIdentity ) {
      kL[0][0] = 1.0;
      kL[0][1] = kL[1][0] = 0.0;
      kL[0][2] = kL[2][0] = 0.0;
      kL[1][1] = fA;
      kL[1][2] = kL[2][1] = fB;
      kL[2][2] = fC;
    }
    else {
      for (int iRow = 0; iRow < 3; iRow++) {
	T fTmp0 = kL[iRow][1];
	T fTmp1 = kL[iRow][2];
	kL[iRow][1] = fA*fTmp0+fB*fTmp1;
	kL[iRow][2] = fB*fTmp0+fC*fTmp1;
      }
    }
  }
}

template <class T>
void 
GbMatrix3<T>::GolubKahanStep (GbMatrix3<T>& kA, GbMatrix3<T>& kL, GbMatrix3<T>& kR)
{
  T fT11 = kA[0][1]*kA[0][1]+kA[1][1]*kA[1][1];
  T fT22 = kA[1][2]*kA[1][2]+kA[2][2]*kA[2][2];
  T fT12 = kA[1][1]*kA[1][2];
  T fTrace = fT11+fT22;
  T fDiff = fT11-fT22;
  T fDiscr = GbMath<T>::Sqrt(fDiff*fDiff+4.0*fT12*fT12);
  T fRoot1 = 0.5*(fTrace+fDiscr);
  T fRoot2 = 0.5*(fTrace-fDiscr);

  // adjust right
  T fY = kA[0][0] - (GbMath<T>::Abs(fRoot1-fT22) <= GbMath<T>::Abs(fRoot2-fT22) ? fRoot1 : fRoot2);
  T fZ = kA[0][1];
  T fInvLength = 1.0/GbMath<T>::Sqrt(fY*fY+fZ*fZ);
  T fSin = fZ*fInvLength;
  T fCos = -fY*fInvLength;

  T fTmp0 = kA[0][0];
  T fTmp1 = kA[0][1];
  kA[0][0] = fCos*fTmp0-fSin*fTmp1;
  kA[0][1] = fSin*fTmp0+fCos*fTmp1;
  kA[1][0] = -fSin*kA[1][1];
  kA[1][1] *= fCos;

  int iRow;
  for (iRow = 0; iRow < 3; iRow++) {
    fTmp0 = kR[0][iRow];
    fTmp1 = kR[1][iRow];
    kR[0][iRow] = fCos*fTmp0-fSin*fTmp1;
    kR[1][iRow] = fSin*fTmp0+fCos*fTmp1;
  }

  // adjust left
  fY = kA[0][0];
  fZ = kA[1][0];
  fInvLength = 1.0/GbMath<T>::Sqrt(fY*fY+fZ*fZ);
  fSin = fZ*fInvLength;
  fCos = -fY*fInvLength;

  kA[0][0] = fCos*kA[0][0]-fSin*kA[1][0];
  fTmp0 = kA[0][1];
  fTmp1 = kA[1][1];
  kA[0][1] = fCos*fTmp0-fSin*fTmp1;
  kA[1][1] = fSin*fTmp0+fCos*fTmp1;
  kA[0][2] = -fSin*kA[1][2];
  kA[1][2] *= fCos;

  int iCol;
  for (iCol = 0; iCol < 3; iCol++) {
    fTmp0 = kL[iCol][0];
    fTmp1 = kL[iCol][1];
    kL[iCol][0] = fCos*fTmp0-fSin*fTmp1;
    kL[iCol][1] = fSin*fTmp0+fCos*fTmp1;
  }

  // adjust right
  fY = kA[0][1];
  fZ = kA[0][2];
  fInvLength = 1.0/GbMath<T>::Sqrt(fY*fY+fZ*fZ);
  fSin = fZ*fInvLength;
  fCos = -fY*fInvLength;

  kA[0][1] = fCos*kA[0][1]-fSin*kA[0][2];
  fTmp0 = kA[1][1];
  fTmp1 = kA[1][2];
  kA[1][1] = fCos*fTmp0-fSin*fTmp1;
  kA[1][2] = fSin*fTmp0+fCos*fTmp1;
  kA[2][1] = -fSin*kA[2][2];
  kA[2][2] *= fCos;

  for (iRow = 0; iRow < 3; iRow++) {
    fTmp0 = kR[1][iRow];
    fTmp1 = kR[2][iRow];
    kR[1][iRow] = fCos*fTmp0-fSin*fTmp1;
    kR[2][iRow] = fSin*fTmp0+fCos*fTmp1;
  }

  // adjust left
  fY = kA[1][1];
  fZ = kA[2][1];
  fInvLength = 1.0/GbMath<T>::Sqrt(fY*fY+fZ*fZ);
  fSin = fZ*fInvLength;
  fCos = -fY*fInvLength;

  kA[1][1] = fCos*kA[1][1]-fSin*kA[2][1];
  fTmp0 = kA[1][2];
  fTmp1 = kA[2][2];
  kA[1][2] = fCos*fTmp0-fSin*fTmp1;
  kA[2][2] = fSin*fTmp0+fCos*fTmp1;

  for (iCol = 0; iCol < 3; iCol++) {
    fTmp0 = kL[iCol][1];
    fTmp1 = kL[iCol][2];
    kL[iCol][1] = fCos*fTmp0-fSin*fTmp1;
    kL[iCol][2] = fSin*fTmp0+fCos*fTmp1;
  }
}

template <class T>
void 
GbMatrix3<T>::singularValueDecomposition (GbMatrix3<T>& kL, GbVec3<T>& kS, GbMatrix3<T>& kR) const
{
  int iRow, iCol;

  GbMatrix3<T> kA = *this;
  bidiagonalize(kA,kL,kR);

  for (int i = 0; i < svdMaxIterations_; i++) {
    T fTmp, fTmp0, fTmp1;
    T fSin0, fCos0, fTan0;
    T fSin1, fCos1, fTan1;

    GbBool bTest1 = (GbMath<T>::Abs(kA[0][1]) <= svdEpsilon_*(GbMath<T>::Abs(kA[0][0])+GbMath<T>::Abs(kA[1][1])));
    GbBool bTest2 = (GbMath<T>::Abs(kA[1][2]) <= svdEpsilon_*(GbMath<T>::Abs(kA[1][1])+GbMath<T>::Abs(kA[2][2])));
    if ( bTest1 ) {
      if ( bTest2 ) {
	kS[0] = kA[0][0];
	kS[1] = kA[1][1];
	kS[2] = kA[2][2];
	break;
      }
      else {
	// 2x2 closed form factorization
	fTmp = (kA[1][1]*kA[1][1] - kA[2][2]*kA[2][2] +
		kA[1][2]*kA[1][2])/(kA[1][2]*kA[2][2]);
	fTan0 = 0.5*(fTmp+GbMath<T>::Sqrt(fTmp*fTmp + 4.0));
	fCos0 = 1.0/GbMath<T>::Sqrt(1.0+fTan0*fTan0);
	fSin0 = fTan0*fCos0;

	for (iCol = 0; iCol < 3; iCol++) {
	  fTmp0 = kL[iCol][1];
	  fTmp1 = kL[iCol][2];
	  kL[iCol][1] = fCos0*fTmp0-fSin0*fTmp1;
	  kL[iCol][2] = fSin0*fTmp0+fCos0*fTmp1;
	}
                
	fTan1 = (kA[1][2]-kA[2][2]*fTan0)/kA[1][1];
	fCos1 = 1.0/GbMath<T>::Sqrt(1.0+fTan1*fTan1);
	fSin1 = -fTan1*fCos1;

	for (iRow = 0; iRow < 3; iRow++) {
	  fTmp0 = kR[1][iRow];
	  fTmp1 = kR[2][iRow];
	  kR[1][iRow] = fCos1*fTmp0-fSin1*fTmp1;
	  kR[2][iRow] = fSin1*fTmp0+fCos1*fTmp1;
	}

	kS[0] = kA[0][0];
	kS[1] = fCos0*fCos1*kA[1][1] - fSin1*(fCos0*kA[1][2]-fSin0*kA[2][2]);
	kS[2] = fSin0*fSin1*kA[1][1] + fCos1*(fSin0*kA[1][2]+fCos0*kA[2][2]);
	break;
      }
    }
    else {
      if ( bTest2 ) {
	// 2x2 closed form factorization 
	fTmp = (kA[0][0]*kA[0][0] + kA[1][1]*kA[1][1] -
		kA[0][1]*kA[0][1])/(kA[0][1]*kA[1][1]);
	fTan0 = 0.5*(-fTmp+GbMath<T>::Sqrt(fTmp*fTmp + 4.0));
	fCos0 = 1.0/GbMath<T>::Sqrt(1.0+fTan0*fTan0);
	fSin0 = fTan0*fCos0;

	for (iCol = 0; iCol < 3; iCol++) {
	  fTmp0 = kL[iCol][0];
	  fTmp1 = kL[iCol][1];
	  kL[iCol][0] = fCos0*fTmp0-fSin0*fTmp1;
	  kL[iCol][1] = fSin0*fTmp0+fCos0*fTmp1;
	}
                
	fTan1 = (kA[0][1]-kA[1][1]*fTan0)/kA[0][0];
	fCos1 = 1.0/GbMath<T>::Sqrt(1.0+fTan1*fTan1);
	fSin1 = -fTan1*fCos1;

	for (iRow = 0; iRow < 3; iRow++) {
	  fTmp0 = kR[0][iRow];
	  fTmp1 = kR[1][iRow];
	  kR[0][iRow] = fCos1*fTmp0-fSin1*fTmp1;
	  kR[1][iRow] = fSin1*fTmp0+fCos1*fTmp1;
	}

	kS[0] = fCos0*fCos1*kA[0][0] - fSin1*(fCos0*kA[0][1]-fSin0*kA[1][1]);
	kS[1] = fSin0*fSin1*kA[0][0] + fCos1*(fSin0*kA[0][1]+fCos0*kA[1][1]);
	kS[2] = kA[2][2];
	break;
      }
      else {
	GolubKahanStep(kA,kL,kR);
      }
    }
  }

  // positize diagonal
  for (iRow = 0; iRow < 3; iRow++) {
    if ( kS[iRow] < 0.0 ) {
      kS[iRow] = -kS[iRow];
      for (iCol = 0; iCol < 3; iCol++) 
	kR[iRow][iCol] = -kR[iRow][iCol];
    }
  }
}

template <class T>
void 
GbMatrix3<T>::singularValueComposition (const GbMatrix3<T>& kL, const GbVec3<T>& kS, const GbMatrix3<T>& kR)
{
  int iRow, iCol;
  GbMatrix3<T> kTmp;

  // product S*R
  for (iRow = 0; iRow < 3; iRow++) {
    for (iCol = 0; iCol < 3; iCol++)
      kTmp[iRow][iCol] = kS[iRow]*kR[iRow][iCol];
  }

  // product L*S*R
  for (iRow = 0; iRow < 3; iRow++) {
    for (iCol = 0; iCol < 3; iCol++) {
      entry_[iRow][iCol] = 0.0;
      for (int iMid = 0; iMid < 3; iMid++)
	entry_[iRow][iCol] += kL[iRow][iMid]*kTmp[iMid][iCol];
    }
  }
}

template <class T>
void 
GbMatrix3<T>::orthonormalize ()
{
  // Algorithm uses Gram-Schmidt orthogonalization.  If 'this' matrix is
  // M = [m0|m1|m2], then orthonormal output matrix is Q = [q0|q1|q2],
  //
  //   q0 = m0/|m0|
  //   q1 = (m1-(q0*m1)q0)/|m1-(q0*m1)q0|
  //   q2 = (m2-(q0*m2)q0-(q1*m2)q1)/|m2-(q0*m2)q0-(q1*m2)q1|
  //
  // where |V| indicates length of vector V and A*B indicates dot
  // product of vectors A and B.

  // compute q0
  T fInvLength = 1.0/GbMath<T>::Sqrt(entry_[0][0]*entry_[0][0] +
				     entry_[1][0]*entry_[1][0] +
				     entry_[2][0]*entry_[2][0]);

  entry_[0][0] *= fInvLength;
  entry_[1][0] *= fInvLength;
  entry_[2][0] *= fInvLength;

  // compute q1
  T fDot0 = entry_[0][0]*entry_[0][1] +
            entry_[1][0]*entry_[1][1] +
            entry_[2][0]*entry_[2][1];

  entry_[0][1] -= fDot0*entry_[0][0];
  entry_[1][1] -= fDot0*entry_[1][0];
  entry_[2][1] -= fDot0*entry_[2][0];

  fInvLength = 1.0/GbMath<T>::Sqrt(entry_[0][1]*entry_[0][1] +
				   entry_[1][1]*entry_[1][1] +
				   entry_[2][1]*entry_[2][1]);

  entry_[0][1] *= fInvLength;
  entry_[1][1] *= fInvLength;
  entry_[2][1] *= fInvLength;

  // compute q2
  T fDot1 = entry_[0][1]*entry_[0][2] +
            entry_[1][1]*entry_[1][2] +
            entry_[2][1]*entry_[2][2];

  fDot0 = entry_[0][0]*entry_[0][2] +
          entry_[1][0]*entry_[1][2] +
          entry_[2][0]*entry_[2][2];

  entry_[0][2] -= fDot0*entry_[0][0] + fDot1*entry_[0][1];
  entry_[1][2] -= fDot0*entry_[1][0] + fDot1*entry_[1][1];
  entry_[2][2] -= fDot0*entry_[2][0] + fDot1*entry_[2][1];

  fInvLength = 1.0/GbMath<T>::Sqrt(entry_[0][2]*entry_[0][2] +
				   entry_[1][2]*entry_[1][2] +
				   entry_[2][2]*entry_[2][2]);

  entry_[0][2] *= fInvLength;
  entry_[1][2] *= fInvLength;
  entry_[2][2] *= fInvLength;
}

template <class T>
void 
GbMatrix3<T>::QDUDecomposition (GbMatrix3<T>& kQ, GbVec3<T>& kD, GbVec3<T>& kU) const
{
  // Factor M = QR = QDU where Q is orthogonal, D is diagonal, 
  // and U is upper triangular with ones on its diagonal.  Algorithm uses
  // Gram-Schmidt orthogonalization (the QR algorithm).
  //
  // If M = [ m0 | m1 | m2 ] and Q = [ q0 | q1 | q2 ], then
  //
  //   q0 = m0/|m0|
  //   q1 = (m1-(q0*m1)q0)/|m1-(q0*m1)q0|
  //   q2 = (m2-(q0*m2)q0-(q1*m2)q1)/|m2-(q0*m2)q0-(q1*m2)q1|
  //
  // where |V| indicates length of vector V and A*B indicates dot
  // product of vectors A and B.  The matrix R has entries
  //
  //   r00 = q0*m0  r01 = q0*m1  r02 = q0*m2
  //   r10 = 0      r11 = q1*m1  r12 = q1*m2
  //   r20 = 0      r21 = 0      r22 = q2*m2
  //
  // so D = diag(r00,r11,r22) and U has entries u01 = r01/r00,
  // u02 = r02/r00, and u12 = r12/r11.

  // Q = rotation
  // D = scaling
  // U = shear

  // D stores the three diagonal entries r00, r11, r22
  // U stores the entries U[0] = u01, U[1] = u02, U[2] = u12

  // build orthogonal matrix Q
  T fInvLength = 1.0/GbMath<T>::Sqrt(entry_[0][0]*entry_[0][0] +
				     entry_[1][0]*entry_[1][0] +
				     entry_[2][0]*entry_[2][0]);
  kQ[0][0] = entry_[0][0]*fInvLength;
  kQ[1][0] = entry_[1][0]*fInvLength;
  kQ[2][0] = entry_[2][0]*fInvLength;

  T fDot = kQ[0][0]*entry_[0][1] + 
           kQ[1][0]*entry_[1][1] +
           kQ[2][0]*entry_[2][1];
  kQ[0][1] = entry_[0][1]-fDot*kQ[0][0];
  kQ[1][1] = entry_[1][1]-fDot*kQ[1][0];
  kQ[2][1] = entry_[2][1]-fDot*kQ[2][0];
  fInvLength = 1.0/GbMath<T>::Sqrt(kQ[0][1]*kQ[0][1] + 
				   kQ[1][1]*kQ[1][1] +
				   kQ[2][1]*kQ[2][1]);
  kQ[0][1] *= fInvLength;
  kQ[1][1] *= fInvLength;
  kQ[2][1] *= fInvLength;

  fDot = kQ[0][0]*entry_[0][2] + 
         kQ[1][0]*entry_[1][2] +
         kQ[2][0]*entry_[2][2];
  kQ[0][2] = entry_[0][2]-fDot*kQ[0][0];
  kQ[1][2] = entry_[1][2]-fDot*kQ[1][0];
  kQ[2][2] = entry_[2][2]-fDot*kQ[2][0];
  fDot = kQ[0][1]*entry_[0][2] + 
         kQ[1][1]*entry_[1][2] +
         kQ[2][1]*entry_[2][2];
  kQ[0][2] -= fDot*kQ[0][1];
  kQ[1][2] -= fDot*kQ[1][1];
  kQ[2][2] -= fDot*kQ[2][1];
  fInvLength = 1.0/GbMath<T>::Sqrt(kQ[0][2]*kQ[0][2] + 
				   kQ[1][2]*kQ[1][2] +
				   kQ[2][2]*kQ[2][2]);
  kQ[0][2] *= fInvLength;
  kQ[1][2] *= fInvLength;
  kQ[2][2] *= fInvLength;

  // guarantee that orthogonal matrix has determinant 1 (no reflections)
  T fDet = kQ[0][0]*kQ[1][1]*kQ[2][2] + 
           kQ[0][1]*kQ[1][2]*kQ[2][0] +
           kQ[0][2]*kQ[1][0]*kQ[2][1] - 
           kQ[0][2]*kQ[1][1]*kQ[2][0] -
           kQ[0][1]*kQ[1][0]*kQ[2][2] - 
           kQ[0][0]*kQ[1][2]*kQ[2][1];

  if ( fDet < 0.0 ) {
    for (int iRow = 0; iRow < 3; iRow++)
      for (int iCol = 0; iCol < 3; iCol++)
	kQ[iRow][iCol] = -kQ[iRow][iCol];
  }

  // build "right" matrix R
  GbMatrix3<T> kR;
  kR[0][0] = kQ[0][0]*entry_[0][0] + 
             kQ[1][0]*entry_[1][0] +
             kQ[2][0]*entry_[2][0];
  kR[0][1] = kQ[0][0]*entry_[0][1] + 
             kQ[1][0]*entry_[1][1] +
             kQ[2][0]*entry_[2][1];
  kR[1][1] = kQ[0][1]*entry_[0][1] + 
             kQ[1][1]*entry_[1][1] +
             kQ[2][1]*entry_[2][1];
  kR[0][2] = kQ[0][0]*entry_[0][2] + 
             kQ[1][0]*entry_[1][2] +
             kQ[2][0]*entry_[2][2];
  kR[1][2] = kQ[0][1]*entry_[0][2] + 
             kQ[1][1]*entry_[1][2] +
             kQ[2][1]*entry_[2][2];
  kR[2][2] = kQ[0][2]*entry_[0][2] + 
             kQ[1][2]*entry_[1][2] +
             kQ[2][2]*entry_[2][2];

  // the scaling component
  kD[0] = kR[0][0];
  kD[1] = kR[1][1];
  kD[2] = kR[2][2];

  // the shear component
  T fInvD0 = 1.0/kD[0];
  kU[0] = kR[0][1]*fInvD0;
  kU[1] = kR[0][2]*fInvD0;
  kU[2] = kR[1][2]/kD[1];
}

template <class T>
T 
GbMatrix3<T>::maxCubicRoot (T afCoeff[3])
{
  // Spectral norm is for A^T*A, so characteristic polynomial
  // P(x) = c[0]+c[1]*x+c[2]*x^2+x^3 has three positive real roots.
  // This yields the assertions c[0] < 0 and c[2]*c[2] >= 3*c[1].

  // quick out for uniform scale (triple root)
  const T fOneThird = 1.0/3.0;
  const T fEpsilon = GbMatrix3<T>::EPSILON;
  T fDiscr = afCoeff[2]*afCoeff[2] - 3.0*afCoeff[1];
  if ( fDiscr <= fEpsilon )
    return -fOneThird*afCoeff[2];

  // Compute an upper bound on roots of P(x).  This assumes that A^T*A
  // has been scaled by its largest entry.
  T fX = 1.0;
  T fPoly = afCoeff[0]+fX*(afCoeff[1]+fX*(afCoeff[2]+fX));
  if ( fPoly < 0.0 ) {
    // uses a matrix norm to find an upper bound on maximum root
    fX = GbMath<T>::Abs(afCoeff[0]);
    T fTmp = 1.0+GbMath<T>::Abs(afCoeff[1]);
    if ( fTmp > fX )
      fX = fTmp;
    fTmp = 1.0+GbMath<T>::Abs(afCoeff[2]);
    if ( fTmp > fX )
      fX = fTmp;
  }

  // Newton's method to find root
  T fTwoC2 = 2.0*afCoeff[2];
  for (int i = 0; i < 16; i++) {
    fPoly = afCoeff[0]+fX*(afCoeff[1]+fX*(afCoeff[2]+fX));
    if ( GbMath<T>::Abs(fPoly) <= fEpsilon )
      return fX;

    T fDeriv = afCoeff[1]+fX*(fTwoC2+3.0*fX);
    fX -= fPoly/fDeriv;
  }

  return fX;
}

template <class T>
T 
GbMatrix3<T>::spectralNorm () const
{
  GbMatrix3<T> kP;
  int iRow, iCol;
  T fPmax = 0.0;
  for (iRow = 0; iRow < 3; iRow++) {
    for (iCol = 0; iCol < 3; iCol++) {
      kP[iRow][iCol] = 0.0;
      for (int iMid = 0; iMid < 3; iMid++) {
	kP[iRow][iCol] += entry_[iMid][iRow]*entry_[iMid][iCol];
      }
      if ( kP[iRow][iCol] > fPmax )
	fPmax = kP[iRow][iCol];
    }
  }

  T fInvPmax = 1.0/fPmax;
  for (iRow = 0; iRow < 3; iRow++) {
    for (iCol = 0; iCol < 3; iCol++)
      kP[iRow][iCol] *= fInvPmax;
  }

  T afCoeff[3];
  afCoeff[0] = -(kP[0][0]*(kP[1][1]*kP[2][2]-kP[1][2]*kP[2][1]) +
		 kP[0][1]*(kP[2][0]*kP[1][2]-kP[1][0]*kP[2][2]) +
		 kP[0][2]*(kP[1][0]*kP[2][1]-kP[2][0]*kP[1][1]));
  afCoeff[1] = kP[0][0]*kP[1][1]-kP[0][1]*kP[1][0] +
               kP[0][0]*kP[2][2]-kP[0][2]*kP[2][0] +
               kP[1][1]*kP[2][2]-kP[1][2]*kP[2][1];
  afCoeff[2] = -(kP[0][0]+kP[1][1]+kP[2][2]);

  T fRoot = maxCubicRoot(afCoeff);
  T fNorm = GbMath<T>::Sqrt(fPmax*fRoot);
  return fNorm;
}

template <class T>
void 
GbMatrix3<T>::toAxisAngle (GbVec3<T>& rkAxis, T& rfRadians) const
{
  // Let (x,y,z) be the unit-length axis and let A be an angle of rotation.
  // The rotation matrix is R = I + sin(A)*P + (1-cos(A))*P^2 where
  // I is the identity and
  //
  //       +-        -+
  //   P = |  0 -z +y |
  //       | +z  0 -x |
  //       | -y +x  0 |
  //       +-        -+
  //
  // If A > 0, R represents a counterclockwise rotation about the axis in
  // the sense of looking from the tip of the axis vector towards the
  // origin.  Some algebra will show that
  //
  //   cos(A) = (trace(R)-1)/2  and  R - R^t = 2*sin(A)*P
  //
  // In the event that A = pi, R-R^t = 0 which prevents us from extracting
  // the axis through P.  Instead note that R = I+2*P^2 when A = pi, so
  // P^2 = (R-I)/2.  The diagonal entries of P^2 are x^2-1, y^2-1, and
  // z^2-1.  We can solve these for axis (x,y,z).  Because the angle is pi,
  // it does not matter which sign you choose on the square roots.

  T fTrace = entry_[0][0] + entry_[1][1] + entry_[2][2];
  T fCos = 0.5*(fTrace-1.0);
  rfRadians = GbMath<T>::ACos(fCos);  // in [0,PI]

  if ( rfRadians > 0.0 ) {
    if ( rfRadians < GbMath<T>::PI ) {
      rkAxis[0] = entry_[2][1]-entry_[1][2];
      rkAxis[1] = entry_[0][2]-entry_[2][0];
      rkAxis[2] = entry_[1][0]-entry_[0][1];
      rkAxis.normalize();
    }
    else {
      // angle is PI
      float fHalfInverse;
      if ( entry_[0][0] >= entry_[1][1] ) {
	// r00 >= r11
	if ( entry_[0][0] >= entry_[2][2] ) {
	  // r00 is maximum diagonal term
	  rkAxis[0] = 0.5*GbMath<T>::Sqrt(entry_[0][0] - entry_[1][1] - entry_[2][2] + 1.0);
	  fHalfInverse = 0.5/rkAxis[0];
	  rkAxis[1] = fHalfInverse*entry_[0][1];
	  rkAxis[2] = fHalfInverse*entry_[0][2];
	}
	else {
	  // r22 is maximum diagonal term
	  rkAxis[2] = 0.5*GbMath<T>::Sqrt(entry_[2][2] - entry_[0][0] - entry_[1][1] + 1.0);
	  fHalfInverse = 0.5/rkAxis[2];
	  rkAxis[0] = fHalfInverse*entry_[0][2];
	  rkAxis[1] = fHalfInverse*entry_[1][2];
	}
      }
      else {
	// r11 > r00
	if ( entry_[1][1] >= entry_[2][2] ) {
	  // r11 is maximum diagonal term
	  rkAxis[1] = 0.5*GbMath<T>::Sqrt(entry_[1][1] - entry_[0][0] - entry_[2][2] + 1.0);
	  fHalfInverse  = 0.5/rkAxis[1];
	  rkAxis[0] = fHalfInverse*entry_[0][1];
	  rkAxis[2] = fHalfInverse*entry_[1][2];
	}
	else {
	  // r22 is maximum diagonal term
	  rkAxis[2] = 0.5*GbMath<T>::Sqrt(entry_[2][2] - entry_[0][0] - entry_[1][1] + 1.0);
	  fHalfInverse = 0.5/rkAxis[2];
	  rkAxis[0] = fHalfInverse*entry_[0][2];
	  rkAxis[1] = fHalfInverse*entry_[1][2];
	}
      }
    }
  }
  else {
    // The angle is 0 and the matrix is the identity.  Any axis will
    // work, so just use the x-axis.
    rkAxis[0] = 1.0;
    rkAxis[1] = 0.0;
    rkAxis[2] = 0.0;
  }
}

template <class T>
void 
GbMatrix3<T>::fromAxisAngle (const GbVec3<T>& rkAxis, T fRadians)
{
  T fCos = GbMath<T>::Cos(fRadians);
  T fSin = GbMath<T>::Sin(fRadians);
  T fOneMinusCos = 1.0-fCos;
  T fX2 = rkAxis[0]*rkAxis[0];
  T fY2 = rkAxis[1]*rkAxis[1];
  T fZ2 = rkAxis[2]*rkAxis[2];
  T fXYM = rkAxis[0]*rkAxis[1]*fOneMinusCos;
  T fXZM = rkAxis[0]*rkAxis[2]*fOneMinusCos;
  T fYZM = rkAxis[1]*rkAxis[2]*fOneMinusCos;
  T fXSin = rkAxis[0]*fSin;
  T fYSin = rkAxis[1]*fSin;
  T fZSin = rkAxis[2]*fSin;
    
  entry_[0][0] = fX2*fOneMinusCos+fCos;
  entry_[0][1] = fXYM-fZSin;
  entry_[0][2] = fXZM+fYSin;
  entry_[1][0] = fXYM+fZSin;
  entry_[1][1] = fY2*fOneMinusCos+fCos;
  entry_[1][2] = fYZM-fXSin;
  entry_[2][0] = fXZM-fYSin;
  entry_[2][1] = fYZM+fXSin;
  entry_[2][2] = fZ2*fOneMinusCos+fCos;
}

template <class T>
GbBool 
GbMatrix3<T>::toEulerAnglesXYZ (T& rfXAngle, T& rfYAngle, T& rfZAngle)
{
  // rot =  cy*cz          -cy*sz           sy
  //        cz*sx*sy+cx*sz  cx*cz-sx*sy*sz -cy*sx
  //       -cx*cz*sy+sx*sz  cz*sx+cx*sy*sz  cx*cy

  if ( entry_[0][2] < 1.0 ) {
    if ( entry_[0][2] > -1.0 ) {
      rfXAngle = GbMath<T>::ATan2(-entry_[1][2],entry_[2][2]);
      rfYAngle = GbMath<T>::ASin(entry_[0][2]);
      rfZAngle = GbMath<T>::ATan2(-entry_[0][1],entry_[0][0]);
      return true;
    }
    else {
      // WARNING.  Not unique.  XA - ZA = -atan2(r10,r11)
      rfXAngle = -GbMath<T>::ATan2(entry_[1][0],entry_[1][1]);
      rfYAngle = -GbMath<T>::PIHALF;
      rfZAngle = 0.0;
      return false;
    }
  }
  else {
    // WARNING.  Not unique.  XAngle + ZAngle = atan2(r10,r11)
    rfXAngle = GbMath<T>::ATan2(entry_[1][0],entry_[1][1]);
    rfYAngle = GbMath<T>::PIHALF;
    rfZAngle = 0.0;
    return false;
  }
}

template <class T>
GbBool 
GbMatrix3<T>::toEulerAnglesXZY (T& rfXAngle, T& rfZAngle, T& rfYAngle)
{
  // rot =  cy*cz          -sz              cz*sy
  //        sx*sy+cx*cy*sz  cx*cz          -cy*sx+cx*sy*sz
  //       -cx*sy+cy*sx*sz  cz*sx           cx*cy+sx*sy*sz

  if ( entry_[0][1] < 1.0 ) {
    if ( entry_[0][1] > -1.0 ) {
      rfXAngle = GbMath<T>::ATan2(entry_[2][1],entry_[1][1]);
      rfZAngle = GbMath<T>::ASin(-entry_[0][1]);
      rfYAngle = GbMath<T>::ATan2(entry_[0][2],entry_[0][0]);
      return true;
    }
    else {
      // WARNING.  Not unique.  XA - YA = atan2(r20,r22)
      rfXAngle = GbMath<T>::ATan2(entry_[2][0],entry_[2][2]);
      rfZAngle = GbMath<T>::PIHALF;
      rfYAngle = 0.0;
      return false;
    }
  }
  else {
    // WARNING.  Not unique.  XA + YA = atan2(-r20,r22)
    rfXAngle = GbMath<T>::ATan2(-entry_[2][0],entry_[2][2]);
    rfZAngle = -GbMath<T>::PIHALF;
    rfYAngle = 0.0;
    return false;
  }
}

template <class T>
GbBool 
GbMatrix3<T>::toEulerAnglesYXZ (T& rfYAngle, T& rfXAngle, T& rfZAngle)
{
  // rot =  cy*cz+sx*sy*sz  cz*sx*sy-cy*sz  cx*sy
  //        cx*sz           cx*cz          -sx
  //       -cz*sy+cy*sx*sz  cy*cz*sx+sy*sz  cx*cy

  if ( entry_[1][2] < 1.0 ) {
    if ( entry_[1][2] > -1.0 ) {
      rfYAngle = GbMath<T>::ATan2(entry_[0][2],entry_[2][2]);
      rfXAngle = GbMath<T>::ASin(-entry_[1][2]);
      rfZAngle = GbMath<T>::ATan2(entry_[1][0],entry_[1][1]);
      return true;
    }
    else {
      // WARNING.  Not unique.  YA - ZA = atan2(r01,r00)
      rfYAngle = GbMath<T>::ATan2(entry_[0][1],entry_[0][0]);
      rfXAngle = GbMath<T>::PIHALF;
      rfZAngle = 0.0;
      return false;
    }
  }
  else {
    // WARNING.  Not unique.  YA + ZA = atan2(-r01,r00)
    rfYAngle = GbMath<T>::ATan2(-entry_[0][1],entry_[0][0]);
    rfXAngle = -GbMath<T>::PIHALF;
    rfZAngle = 0.0;
    return false;
  }
}

template <class T>
GbBool 
GbMatrix3<T>::toEulerAnglesYZX (T& rfYAngle, T& rfZAngle, T& rfXAngle)
{
  // rot =  cy*cz           sx*sy-cx*cy*sz  cx*sy+cy*sx*sz
  //        sz              cx*cz          -cz*sx
  //       -cz*sy           cy*sx+cx*sy*sz  cx*cy-sx*sy*sz

  if ( entry_[1][0] < 1.0 ) {
    if ( entry_[1][0] > -1.0 ) {
      rfYAngle = GbMath<T>::ATan2(-entry_[2][0],entry_[0][0]);
      rfZAngle = GbMath<T>::ASin(entry_[1][0]);
      rfXAngle = GbMath<T>::ATan2(-entry_[1][2],entry_[1][1]);
      return true;
    }
    else {
      // WARNING.  Not unique.  YA - XA = -atan2(r21,r22);
      rfYAngle = -GbMath<T>::ATan2(entry_[2][1],entry_[2][2]);
      rfZAngle = -GbMath<T>::PIHALF;
      rfXAngle = 0.0;
      return false;
    }
  }
  else {
    // WARNING.  Not unique.  YA + XA = atan2(r21,r22)
    rfYAngle = GbMath<T>::ATan2(entry_[2][1],entry_[2][2]);
    rfZAngle = GbMath<T>::PIHALF;
    rfXAngle = 0.0;
    return false;
  }
}

template <class T>
GbBool 
GbMatrix3<T>::toEulerAnglesZXY (T& rfZAngle, T& rfXAngle, T& rfYAngle)
{
  // rot =  cy*cz-sx*sy*sz -cx*sz           cz*sy+cy*sx*sz
  //        cz*sx*sy+cy*sz  cx*cz          -cy*cz*sx+sy*sz
  //       -cx*sy           sx              cx*cy

  if ( entry_[2][1] < 1.0 ) {
    if ( entry_[2][1] > -1.0 ) {
      rfZAngle = GbMath<T>::ATan2(-entry_[0][1],entry_[1][1]);
      rfXAngle = GbMath<T>::ASin(entry_[2][1]);
      rfYAngle = GbMath<T>::ATan2(-entry_[2][0],entry_[2][2]);
      return true;
    }
    else {
      // WARNING.  Not unique.  ZA - YA = -atan(r02,r00)
      rfZAngle = -GbMath<T>::ATan2(entry_[0][2],entry_[0][0]);
      rfXAngle = -GbMath<T>::PIHALF;
      rfYAngle = 0.0;
      return false;
    }
  }
  else {
    // WARNING.  Not unique.  ZA + YA = atan2(r02,r00)
    rfZAngle = GbMath<T>::ATan2(entry_[0][2],entry_[0][0]);
    rfXAngle = GbMath<T>::PIHALF;
    rfYAngle = 0.0;
    return false;
  }
}

template <class T>
GbBool 
GbMatrix3<T>::toEulerAnglesZYX (T& rfZAngle, T& rfYAngle, T& rfXAngle)
{
  // rot =  cy*cz           cz*sx*sy-cx*sz  cx*cz*sy+sx*sz
  //        cy*sz           cx*cz+sx*sy*sz -cz*sx+cx*sy*sz
  //       -sy              cy*sx           cx*cy

  if ( entry_[2][0] < 1.0 ) {
    if ( entry_[2][0] > -1.0 ) {
      rfZAngle = GbMath<T>::ATan2(entry_[1][0],entry_[0][0]);
      rfYAngle = GbMath<T>::ASin(-entry_[2][0]);
      rfXAngle = GbMath<T>::ATan2(entry_[2][1],entry_[2][2]);
      return true;
    }
    else {
      // WARNING.  Not unique.  ZA - XA = -atan2(r01,r02)
      rfZAngle = -GbMath<T>::ATan2(entry_[0][1],entry_[0][2]);
      rfYAngle = GbMath<T>::PIHALF;
      rfXAngle = 0.0;
      return false;
    }
  }
  else {
    // WARNING.  Not unique.  ZA + XA = atan2(-r01,-r02)
    rfZAngle = GbMath<T>::ATan2(-entry_[0][1],-entry_[0][2]);
    rfYAngle = -GbMath<T>::PIHALF;
    rfXAngle = 0.0;
    return false;
  }
}

template <class T>
void 
GbMatrix3<T>::fromEulerAnglesXYZ (T fYAngle, T fPAngle, T fRAngle)
{
  T fCos, fSin;

  fCos = GbMath<T>::Cos(fYAngle);
  fSin = GbMath<T>::Sin(fYAngle);
  GbMatrix3<T> kXMat(1.0,0.0,0.0,
		     0.0,fCos,-fSin,
		     0.0,fSin,fCos);

  fCos = GbMath<T>::Cos(fPAngle);
  fSin = GbMath<T>::Sin(fPAngle);
  GbMatrix3<T> kYMat(fCos,0.0,fSin,
		     0.0,1.0,0.0,
		     -fSin,0.0,fCos);

  fCos = GbMath<T>::Cos(fRAngle);
  fSin = GbMath<T>::Sin(fRAngle);
  GbMatrix3<T> kZMat(fCos,-fSin,0.0,
		     fSin,fCos,0.0,
		     0.0,0.0,1.0);

  *this = kXMat*(kYMat*kZMat);
}

template <class T>
void 
GbMatrix3<T>::fromEulerAnglesXZY (T fYAngle, T fPAngle, T fRAngle)
{
  T fCos, fSin;

  fCos = GbMath<T>::Cos(fYAngle);
  fSin = GbMath<T>::Sin(fYAngle);
  GbMatrix3<T> kXMat(1.0,0.0,0.0,
		     0.0,fCos,-fSin,
		     0.0,fSin,fCos);

  fCos = GbMath<T>::Cos(fPAngle);
  fSin = GbMath<T>::Sin(fPAngle);
  GbMatrix3<T> kZMat(fCos,-fSin,0.0,
		     fSin,fCos,0.0,
		     0.0,0.0,1.0);

  fCos = GbMath<T>::Cos(fRAngle);
  fSin = GbMath<T>::Sin(fRAngle);
  GbMatrix3<T> kYMat(fCos,0.0,fSin,
		     0.0,1.0,0.0,
		     -fSin,0.0,fCos);

  *this = kXMat*(kZMat*kYMat);
}

template <class T>
void 
GbMatrix3<T>::fromEulerAnglesYXZ (T fYAngle, T fPAngle, T fRAngle)
{
  T fCos, fSin;

  fCos = GbMath<T>::Cos(fYAngle);
  fSin = GbMath<T>::Sin(fYAngle);
  GbMatrix3<T> kYMat(fCos,0.0,fSin,
		     0.0,1.0,0.0,
		     -fSin,0.0,fCos);

  fCos = GbMath<T>::Cos(fPAngle);
  fSin = GbMath<T>::Sin(fPAngle);
  GbMatrix3<T> kXMat(1.0,0.0,0.0,
		     0.0,fCos,-fSin,
		     0.0,fSin,fCos);

  fCos = GbMath<T>::Cos(fRAngle);
  fSin = GbMath<T>::Sin(fRAngle);
  GbMatrix3<T> kZMat(fCos,-fSin,0.0,
		     fSin,fCos,0.0,
		     0.0,0.0,1.0);

  *this = kYMat*(kXMat*kZMat);
}

template <class T>
void 
GbMatrix3<T>::fromEulerAnglesYZX (T fYAngle, T fPAngle, T fRAngle)
{
  T fCos, fSin;

  fCos = GbMath<T>::Cos(fYAngle);
  fSin = GbMath<T>::Sin(fYAngle);
  GbMatrix3<T> kYMat(fCos,0.0,fSin,
		     0.0,1.0,0.0,
		     -fSin,0.0,fCos);

  fCos = GbMath<T>::Cos(fPAngle);
  fSin = GbMath<T>::Sin(fPAngle);
  GbMatrix3<T> kZMat(fCos,-fSin,0.0,
		     fSin,fCos,0.0,
		     0.0,0.0,1.0);

  fCos = GbMath<T>::Cos(fRAngle);
  fSin = GbMath<T>::Sin(fRAngle);
  GbMatrix3<T> kXMat(1.0,0.0,0.0,
		     0.0,fCos,-fSin,
		     0.0,fSin,fCos);

  *this = kYMat*(kZMat*kXMat);
}

template <class T>
void 
GbMatrix3<T>::fromEulerAnglesZXY (T fYAngle, T fPAngle, T fRAngle)
{
  T fCos, fSin;

  fCos = GbMath<T>::Cos(fYAngle);
  fSin = GbMath<T>::Sin(fYAngle);
  GbMatrix3<T> kZMat(fCos,-fSin,0.0,
		     fSin,fCos,0.0,
		     0.0,0.0,1.0);

  fCos = GbMath<T>::Cos(fPAngle);
  fSin = GbMath<T>::Sin(fPAngle);
  GbMatrix3<T> kXMat(1.0,0.0,0.0,
		     0.0,fCos,-fSin,
		     0.0,fSin,fCos);

  fCos = GbMath<T>::Cos(fRAngle);
  fSin = GbMath<T>::Sin(fRAngle);
  GbMatrix3<T> kYMat(fCos,0.0,fSin,
		     0.0,1.0,0.0,
		     -fSin,0.0,fCos);

  *this = kZMat*(kXMat*kYMat);
}

template <class T>
void 
GbMatrix3<T>::fromEulerAnglesZYX (T fYAngle, T fPAngle, T fRAngle)
{
  T fCos, fSin;

  fCos = GbMath<T>::Cos(fYAngle);
  fSin = GbMath<T>::Sin(fYAngle);
  GbMatrix3<T> kZMat(fCos,-fSin,0.0,
		     fSin,fCos,0.0,
		     0.0,0.0,1.0);

  fCos = GbMath<T>::Cos(fPAngle);
  fSin = GbMath<T>::Sin(fPAngle);
  GbMatrix3<T> kYMat(fCos,0.0,fSin,
		     0.0,1.0,0.0,
		     -fSin,0.0,fCos);

  fCos = GbMath<T>::Cos(fRAngle);
  fSin = GbMath<T>::Sin(fRAngle);
  GbMatrix3<T> kXMat(1.0,0.0,0.0,
		     0.0,fCos,-fSin,
		     0.0,fSin,fCos);

  *this = kZMat*(kYMat*kXMat);
}

template <class T>
void 
GbMatrix3<T>::tridiagonal (T afDiag[3], T afSubDiag[3])
{
  // Householder reduction T = Q^t M Q
  //   Input:   
  //     mat, symmetric 3x3 matrix M
  //   Output:  
  //     mat, orthogonal matrix Q
  //     diag, diagonal entries of T
  //     subd, subdiagonal entries of T (T is symmetric)

  T fA = entry_[0][0];
  T fB = entry_[0][1];
  T fC = entry_[0][2];
  T fD = entry_[1][1];
  T fE = entry_[1][2];
  T fF = entry_[2][2];

  afDiag[0] = fA;
  afSubDiag[2] = 0.0;
  if ( GbMath<T>::Abs(fC) >= GbMatrix3<T>::EPSILON ) {
    T fLength = GbMath<T>::Sqrt(fB*fB+fC*fC);
    T fInvLength = 1.0/fLength;
    fB *= fInvLength;
    fC *= fInvLength;
    T fQ = 2.0*fB*fE+fC*(fF-fD);
    afDiag[1] = fD+fC*fQ;
    afDiag[2] = fF-fC*fQ;
    afSubDiag[0] = fLength;
    afSubDiag[1] = fE-fB*fQ;
    entry_[0][0] = 1.0;
    entry_[0][1] = 0.0;
    entry_[0][2] = 0.0;
    entry_[1][0] = 0.0;
    entry_[1][1] = fB;
    entry_[1][2] = fC;
    entry_[2][0] = 0.0;
    entry_[2][1] = fC;
    entry_[2][2] = -fB;
  }
  else {
    afDiag[1] = fD;
    afDiag[2] = fF;
    afSubDiag[0] = fB;
    afSubDiag[1] = fE;
    entry_[0][0] = 1.0;
    entry_[0][1] = 0.0;
    entry_[0][2] = 0.0;
    entry_[1][0] = 0.0;
    entry_[1][1] = 1.0;
    entry_[1][2] = 0.0;
    entry_[2][0] = 0.0;
    entry_[2][1] = 0.0;
    entry_[2][2] = 1.0;
  }
}

template <class T>
GbBool 
GbMatrix3<T>::QLAlgorithm (T afDiag[3], T afSubDiag[3])
{
  // QL iteration with implicit shifting to reduce matrix from tridiagonal
  // to diagonal

  for (int i0 = 0; i0 < 3; i0++) {
    const int iMaxIter = 32;
    int iIter;
    for (iIter = 0; iIter < iMaxIter; iIter++) {
      int i1;
      for (i1 = i0; i1 <= 1; i1++) {
	T fSum = GbMath<T>::Abs(afDiag[i1]) + GbMath<T>::Abs(afDiag[i1+1]);
	if ( GbMath<T>::Abs(afSubDiag[i1]) + fSum == fSum )
	  break;
      }
      if ( i1 == i0 )
	break;

      T fTmp0 = (afDiag[i0+1]-afDiag[i0])/(2.0*afSubDiag[i0]);
      T fTmp1 = GbMath<T>::Sqrt(fTmp0*fTmp0+1.0);
      if ( fTmp0 < 0.0 )
	fTmp0 = afDiag[i1]-afDiag[i0]+afSubDiag[i0]/(fTmp0-fTmp1);
      else
	fTmp0 = afDiag[i1]-afDiag[i0]+afSubDiag[i0]/(fTmp0+fTmp1);
      T fSin = 1.0;
      T fCos = 1.0;
      T fTmp2 = 0.0;
      for (int i2 = i1-1; i2 >= i0; i2--) {
	T fTmp3 = fSin*afSubDiag[i2];
	T fTmp4 = fCos*afSubDiag[i2];
	if ( GbMath<T>::Abs(fTmp3) >= GbMath<T>::Abs(fTmp0) ) {
	  fCos = fTmp0/fTmp3;
	  fTmp1 = GbMath<T>::Sqrt(fCos*fCos+1.0);
	  afSubDiag[i2+1] = fTmp3*fTmp1;
	  fSin = 1.0/fTmp1;
	  fCos *= fSin;
	}
	else {
	  fSin = fTmp3/fTmp0;
	  fTmp1 = GbMath<T>::Sqrt(fSin*fSin+1.0);
	  afSubDiag[i2+1] = fTmp0*fTmp1;
	  fCos = 1.0/fTmp1;
	  fSin *= fCos;
	}
	fTmp0 = afDiag[i2+1]-fTmp2;
	fTmp1 = (afDiag[i2]-fTmp0)*fSin+2.0*fTmp4*fCos;
	fTmp2 = fSin*fTmp1;
	afDiag[i2+1] = fTmp0+fTmp2;
	fTmp0 = fCos*fTmp1-fTmp4;

	for (int iRow = 0; iRow < 3; iRow++) {
	  fTmp3 = entry_[iRow][i2+1];
	  entry_[iRow][i2+1] = fSin*entry_[iRow][i2] + fCos*fTmp3;
	  entry_[iRow][i2]   = fCos*entry_[iRow][i2] - fSin*fTmp3;
	}
      }
      afDiag[i0] -= fTmp2;
      afSubDiag[i0] = fTmp0;
      afSubDiag[i1] = 0.0;
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
GbMatrix3<T>::eigenSolveSymmetric (T afEigenvalue[3], GbVec3<T> akEigenvector[3]) const
{
  GbMatrix3<T> kMatrix = *this;
  T afSubDiag[3];
  kMatrix.tridiagonal(afEigenvalue,afSubDiag);
  kMatrix.QLAlgorithm(afEigenvalue,afSubDiag);

  for (int i = 0; i < 3; i++) {
    akEigenvector[i][0] = kMatrix[0][i];
    akEigenvector[i][1] = kMatrix[1][i];
    akEigenvector[i][2] = kMatrix[2][i];
  }

  // make eigenvectors form a right--handed system
  GbVec3<T> kCross = akEigenvector[1].cross(akEigenvector[2]);
  T fDet = akEigenvector[0].dot(kCross);
  if ( fDet < 0.0 ) {
    akEigenvector[2][0] = - akEigenvector[2][0];
    akEigenvector[2][1] = - akEigenvector[2][1];
    akEigenvector[2][2] = - akEigenvector[2][2];
  }
}

template <class T>
void 
GbMatrix3<T>::tensorProduct (const GbVec3<T>& rkU, const GbVec3<T>& rkV, GbMatrix3<T>& rkProduct)
{
  for (int iRow = 0; iRow < 3; iRow++) {
    for (int iCol = 0; iCol < 3; iCol++)
      rkProduct[iRow][iCol] = rkU[iRow]*rkV[iCol];
  }
}

template <class T>
GbMatrix3<T>
GbMatrix3<T>::slerp (T fT, const GbMatrix3<T>& rkR0, const GbMatrix3<T>& rkR1)
{
  GbVec3<T> kAxis;
  T fAngle;
  GbMatrix3<T> kProd = rkR0.transposeTimes(rkR1);
  kProd.toAxisAngle(kAxis,fAngle);

  GbMatrix3<T> kR;
  kR.fromAxisAngle(kAxis,fT*fAngle);
  
  return kR;
}

