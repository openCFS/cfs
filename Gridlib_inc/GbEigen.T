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
| Revision 1.1  2001/01/02 15:02:12  prkipfer
| introduced new classes
|
|
|
+---------------------------------------------------------------------*/

template <class T>
GbEigen<T>::GbEigen(int iSize)
{
  if ( iSize < 2 ) {
    errormsg("eigensolver only works on matrices, you specified rank "<<iSize);
    size_ = 2;
  }
  else 
    size_ = iSize;

  matrix_ = new T*[size_];
  for (int i = 0; i < size_; i++)
    matrix_[i] = new T[size_];

  diagonal_ = new T[size_];
  subdiagonal_ = new T[size_];
}

template <class T>
GbEigen<T>::~GbEigen()
{
  delete[] subdiagonal_;
  delete[] diagonal_;
  for (int i = 0; i < size_; i++)
    delete[] matrix_[i];
  delete[] matrix_;
}

template <class T>
void 
GbEigen<T>::tridiagonal2 ()
{
  // matrix is already tridiagonal
  diagonal_[0] = matrix_[0][0];
  diagonal_[1] = matrix_[1][1];
  subdiagonal_[0] = matrix_[0][1];
  subdiagonal_[1] = T(0.0);
  matrix_[0][0] = T(1.0);
  matrix_[0][1] = T(0.0);
  matrix_[1][0] = T(0.0);
  matrix_[1][1] = T(1.0);
}

template <class T>
void 
GbEigen<T>::tridiagonal3 ()
{
  T fM00 = matrix_[0][0];
  T fM01 = matrix_[0][1];
  T fM02 = matrix_[0][2];
  T fM11 = matrix_[1][1];
  T fM12 = matrix_[1][2];
  T fM22 = matrix_[2][2];

  diagonal_[0] = fM00;
  subdiagonal_[2] = T(0.0);
  if ( fM02 != T(0.0) ) {
    T fLength = GbMath<T>::Sqrt(fM01*fM01+fM02*fM02);
    T fInvLength = T(1.0)/fLength;
    fM01 *= fInvLength;
    fM02 *= fInvLength;
    T fQ = T(2.0)*fM01*fM12+fM02*(fM22-fM11);
    diagonal_[1] = fM11+fM02*fQ;
    diagonal_[2] = fM22-fM02*fQ;
    subdiagonal_[0] = fLength;
    subdiagonal_[1] = fM12-fM01*fQ;
    matrix_[0][0] = T(1.0); matrix_[0][1] = T(0.0);  matrix_[0][2] = T(0.0);
    matrix_[1][0] = T(0.0); matrix_[1][1] = fM01; matrix_[1][2] = fM02;
    matrix_[2][0] = T(0.0); matrix_[2][1] = fM02; matrix_[2][2] = -fM01;
  }
  else {
    diagonal_[1] = fM11;
    diagonal_[2] = fM22;
    subdiagonal_[0] = fM01;
    subdiagonal_[1] = fM12;
    matrix_[0][0] = T(1.0); matrix_[0][1] = T(0.0); matrix_[0][2] = T(0.0);
    matrix_[1][0] = T(0.0); matrix_[1][1] = T(1.0); matrix_[1][2] = T(0.0);
    matrix_[2][0] = T(0.0); matrix_[2][1] = T(0.0); matrix_[2][2] = T(1.0);
  }
}

template <class T>
void 
GbEigen<T>::tridiagonal4 ()
{
  // save matrix M
  T fM00 = matrix_[0][0];
  T fM01 = matrix_[0][1];
  T fM02 = matrix_[0][2];
  T fM03 = matrix_[0][3];
  T fM11 = matrix_[1][1];
  T fM12 = matrix_[1][2];
  T fM13 = matrix_[1][3];
  T fM22 = matrix_[2][2];
  T fM23 = matrix_[2][3];
  T fM33 = matrix_[3][3];

  diagonal_[0] = fM00;
  subdiagonal_[3] = T(0.0);

  matrix_[0][0] = T(1.0);
  matrix_[0][1] = T(0.0);
  matrix_[0][2] = T(0.0);
  matrix_[0][3] = T(0.0);
  matrix_[1][0] = T(0.0);
  matrix_[2][0] = T(0.0);
  matrix_[3][0] = T(0.0);
  
  float fLength, fInvLength;

  if ( fM02 != T(0.0) || fM03 != T(0.0) ) {
    T fQ11, fQ12, fQ13;
    T fQ21, fQ22, fQ23;
    T fQ31, fQ32, fQ33;

    // build column Q1
    fLength = GbMath<T>::Sqrt(fM01*fM01+fM02*fM02+fM03*fM03);
    fInvLength = T(1.0)/fLength;
    fQ11 = fM01*fInvLength;
    fQ21 = fM02*fInvLength;
    fQ31 = fM03*fInvLength;

    subdiagonal_[0] = fLength;

    // compute S*Q1
    T fV0 = fM11*fQ11+fM12*fQ21+fM13*fQ31;
    T fV1 = fM12*fQ11+fM22*fQ21+fM23*fQ31;
    T fV2 = fM13*fQ11+fM23*fQ21+fM33*fQ31;

    diagonal_[1] = fQ11*fV0+fQ21*fV1+fQ31*fV2;

    // build column Q3 = Q1x(S*Q1)
    fQ13 = fQ21*fV2-fQ31*fV1;
    fQ23 = fQ31*fV0-fQ11*fV2;
    fQ33 = fQ11*fV1-fQ21*fV0;
    fLength = GbMath<T>::Sqrt(fQ13*fQ13+fQ23*fQ23+fQ33*fQ33);
    if ( fLength > T(0.0) ) {
      fInvLength = T(1.0)/fLength;
      fQ13 *= fInvLength;
      fQ23 *= fInvLength;
      fQ33 *= fInvLength;

      // build column Q2 = Q3xQ1
      fQ12 = fQ23*fQ31-fQ33*fQ21; 
      fQ22 = fQ33*fQ11-fQ13*fQ31;
      fQ32 = fQ13*fQ21-fQ23*fQ11;

      fV0 = fQ12*fM11+fQ22*fM12+fQ32*fM13;
      fV1 = fQ12*fM12+fQ22*fM22+fQ32*fM23;
      fV2 = fQ12*fM13+fQ22*fM23+fQ32*fM33;
      subdiagonal_[1] = fQ11*fV0+fQ21*fV1+fQ31*fV2;
      diagonal_[2] = fQ12*fV0+fQ22*fV1+fQ32*fV2;
      subdiagonal_[2] = fQ13*fV0+fQ23*fV1+fQ33*fV2;

      fV0 = fQ13*fM11+fQ23*fM12+fQ33*fM13;
      fV1 = fQ13*fM12+fQ23*fM22+fQ33*fM23;
      fV2 = fQ13*fM13+fQ23*fM23+fQ33*fM33;
      diagonal_[3] = fQ13*fV0+fQ23*fV1+fQ33*fV2;
    }
    else {
      // S*Q1 parallel to Q1, choose any valid Q2 and Q3
      subdiagonal_[1] = T(0.0);

      fLength = fQ21*fQ21+fQ31*fQ31;
      if ( fLength > T(0.0) ) {
	fInvLength = T(1.0)/fLength;
	T fTmp = fQ11-T(1.0);
	fQ12 = -fQ21;
	fQ22 = T(1.0)+fTmp*fQ21*fQ21*fInvLength;
	fQ32 = fTmp*fQ21*fQ31*fInvLength;
	
	fQ13 = -fQ31;
	fQ23 = fQ32;
	fQ33 = T(1.0)+fTmp*fQ31*fQ31*fInvLength;
	
	fV0 = fQ12*fM11+fQ22*fM12+fQ32*fM13;
	fV1 = fQ12*fM12+fQ22*fM22+fQ32*fM23;
	fV2 = fQ12*fM13+fQ22*fM23+fQ32*fM33;
	diagonal_[2] = fQ12*fV0+fQ22*fV1+fQ32*fV2;
	subdiagonal_[2] = fQ13*fV0+fQ23*fV1+fQ33*fV2;

	fV0 = fQ13*fM11+fQ23*fM12+fQ33*fM13;
	fV1 = fQ13*fM12+fQ23*fM22+fQ33*fM23;
	fV2 = fQ13*fM13+fQ23*fM23+fQ33*fM33;
	diagonal_[3] = fQ13*fV0+fQ23*fV1+fQ33*fV2;
      }
      else {
	// Q1 = (+-1,0,0)
	fQ12 = T(0.0); fQ22 = T(1.0); fQ32 = T(0.0);
	fQ13 = T(0.0); fQ23 = T(0.0); fQ33 = T(1.0);

	diagonal_[2] = fM22;
	diagonal_[3] = fM33;
	subdiagonal_[2] = fM23;
      }
    }
    
    matrix_[1][1] = fQ11; matrix_[1][2] = fQ12; matrix_[1][3] = fQ13;
    matrix_[2][1] = fQ21; matrix_[2][2] = fQ22; matrix_[2][3] = fQ23;
    matrix_[3][1] = fQ31; matrix_[3][2] = fQ32; matrix_[3][3] = fQ33;
  }
  else {
    diagonal_[1] = fM11;
    subdiagonal_[0] = fM01;
    matrix_[1][1] = T(1.0);
    matrix_[2][1] = T(0.0);
    matrix_[3][1] = T(0.0); 

    if ( fM13 != T(0.0) ) {
      fLength = GbMath<T>::Sqrt(fM12*fM12+fM13*fM13);
      fInvLength = T(1.0)/fLength;
      fM12 *= fInvLength;
      fM13 *= fInvLength;
      T fQ = T(2.0)*fM12*fM23+fM13*(fM33-fM22);

      diagonal_[2] = fM22+fM13*fQ;
      diagonal_[3] = fM33-fM13*fQ;
      subdiagonal_[1] = fLength;
      subdiagonal_[2] = fM23-fM12*fQ;
      matrix_[1][2] = T(0.0);
      matrix_[1][3] = T(0.0);
      matrix_[2][2] = fM12;
      matrix_[2][3] = fM13;
      matrix_[3][2] = fM13;
      matrix_[3][3] = -fM12;
    }
    else {
      diagonal_[2] = fM22;
      diagonal_[3] = fM33;
      subdiagonal_[1] = fM12;
      subdiagonal_[2] = fM23;
      matrix_[1][2] = T(0.0);
      matrix_[1][3] = T(0.0);
      matrix_[2][2] = T(1.0);
      matrix_[2][3] = T(0.0);
      matrix_[3][2] = T(0.0);
      matrix_[3][3] = T(1.0);
    }
  }
}

template <class T>
void 
GbEigen<T>::tridiagonalN ()
{
  int i0, i1, i2, i3;

  for (i0 = size_-1, i3 = size_-2; i0 >= 1; i0--, i3--) {
    T fH = T(0.0), fScale = T(0.0);

    if ( i3 > 0 ) {
      for (i2 = 0; i2 <= i3; i2++)
	fScale += GbMath<T>::Abs(matrix_[i0][i2]);
      if ( fScale == 0 ) {
	subdiagonal_[i0] = matrix_[i0][i3];
      }
      else {
	T fInvScale = T(1.0)/fScale;
	for (i2 = 0; i2 <= i3; i2++) {
	  matrix_[i0][i2] *= fInvScale;
	  fH += matrix_[i0][i2]*matrix_[i0][i2];
	}
	T fF = matrix_[i0][i3];
	T fG = GbMath<T>::Sqrt(fH);
	if ( fF > T(0.0) )
	  fG = -fG;
	subdiagonal_[i0] = fScale*fG;
	fH -= fF*fG;
	matrix_[i0][i3] = fF-fG;
	fF = T(0.0);
	T fInvH = T(1.0)/fH;
	for (i1 = 0; i1 <= i3; i1++) {
	  matrix_[i1][i0] = matrix_[i0][i1]*fInvH;
	  fG = T(0.0);
	  for (i2 = 0; i2 <= i1; i2++)
	    fG += matrix_[i1][i2]*matrix_[i0][i2];
	  for (i2 = i1+1; i2 <= i3; i2++)
	    fG += matrix_[i2][i1]*matrix_[i0][i2];
	  subdiagonal_[i1] = fG*fInvH;
	  fF += subdiagonal_[i1]*matrix_[i0][i1];
	}
	T fHalfFdivH = T(0.5)*fF*fInvH;
	for (i1 = 0; i1 <= i3; i1++) {
	  fF = matrix_[i0][i1];
	  fG = subdiagonal_[i1] - fHalfFdivH*fF;
	  subdiagonal_[i1] = fG;
	  for (i2 = 0; i2 <= i1; i2++) {
	    matrix_[i1][i2] -= fF*subdiagonal_[i2] + fG*matrix_[i0][i2];
	  }
	}
      }
    }
    else {
      subdiagonal_[i0] = matrix_[i0][i3];
    }

    diagonal_[i0] = fH;
  }

  diagonal_[0] = subdiagonal_[0] = T(0.0);
  for (i0 = 0, i3 = -1; i0 <= size_-1; i0++, i3++) {
    if ( diagonal_[i0] ) {
      for (i1 = 0; i1 <= i3; i1++) {
	T fSum = T(0.0);
	for (i2 = 0; i2 <= i3; i2++)
	  fSum += matrix_[i0][i2]*matrix_[i2][i1];
	for (i2 = 0; i2 <= i3; i2++)
	  matrix_[i2][i1] -= fSum*matrix_[i2][i0];
      }
    }
    diagonal_[i0] = matrix_[i0][i0];
    matrix_[i0][i0] = 1;
    for (i1 = 0; i1 <= i3; i1++)
      matrix_[i1][i0] = matrix_[i0][i1] = T(0.0);
  }

  // re-ordering if GbEigen<T>::QLAlgorithm is used subsequently
  for (i0 = 1, i3 = 0; i0 < size_; i0++, i3++)
    subdiagonal_[i3] = subdiagonal_[i0];
  subdiagonal_[size_-1] = T(0.0);
}

template <class T>
GbBool 
GbEigen<T>::QLAlgorithm ()
{
  const int iMaxIter = 32;

  for (int i0 = 0; i0 < size_; i0++) {
    int i1;
    for (i1 = 0; i1 < iMaxIter; i1++) {
      int i2;
      for (i2 = i0; i2 <= size_-2; i2++) {
	T fTmp = GbMath<T>::Abs(diagonal_[i2])+GbMath<T>::Abs(diagonal_[i2+1]);
	if ( GbMath<T>::Abs(subdiagonal_[i2]) + fTmp == fTmp )
	  break;
      }
      if ( i2 == i0 )
	break;

      T fG = (diagonal_[i0+1]-diagonal_[i0])/(T(2.0)*subdiagonal_[i0]);
//      T fR = (T)sqrt(fG*fG+T(1.0));
      T fR = GbMath<T>::Sqrt(fG*fG+T(1.0));
      if ( fG < T(0.0) )
	fG = diagonal_[i2]-diagonal_[i0]+subdiagonal_[i0]/(fG-fR);
      else
	fG = diagonal_[i2]-diagonal_[i0]+subdiagonal_[i0]/(fG+fR);
      T fs = T(1.0), fc = T(1.0), fP = T(0.0);
      for (int i3 = i2-1; i3 >= i0; i3--) {
	T fF = fs*subdiagonal_[i3];
	T fB = fc*subdiagonal_[i3];
	if ( GbMath<T>::Abs(fF) >= GbMath<T>::Abs(fG) ) {
	  fc = fG/fF;
  	  fR = GbMath<T>::Sqrt(fc*fc+T(1.0));
//	  fR = (T)sqrt(fc*fc+T(1.0));
	  subdiagonal_[i3+1] = fF*fR;
	  fs = T(1.0)/fR;
	  fc *= fs;
	}
	else {
	  fs = fF/fG;
  	  fR = GbMath<T>::Sqrt(fs*fs+T(1.0));
//	  fR = (T)sqrt(fs*fs+T(1.0));
	  subdiagonal_[i3+1] = fG*fR;
	  fc = T(1.0)/fR;
	  fs *= fc;
	}
	fG = diagonal_[i3+1]-fP;
	fR = (diagonal_[i3]-fG)*fs+T(2.0)*fB*fc;
	fP = fs*fR;
	diagonal_[i3+1] = fG+fP;
	fG = fc*fR-fB;

	for (int i4 = 0; i4 < size_; i4++) {
	  fF = matrix_[i4][i3+1];
	  matrix_[i4][i3+1] = fs*matrix_[i4][i3]+fc*fF;
	  matrix_[i4][i3] = fc*matrix_[i4][i3]-fs*fF;
	}
      }
      diagonal_[i0] -= fP;
      subdiagonal_[i0] = fG;
      subdiagonal_[i2] = T(0.0);
    }
    if ( i1 == iMaxIter )
      return false;
  }
  
  return true;
}

template <class T>
void 
GbEigen<T>::decreasingSort ()
{
  // sort eigenvalues in decreasing order, e[0] >= ... >= e[size_-1]
  for (int i0 = 0, i1; i0 <= size_-2; i0++) {
    // locate maximum eigenvalue
    i1 = i0;
    T fMax = diagonal_[i1];
    int i2;
    for (i2 = i0+1; i2 < size_; i2++) {
      if ( diagonal_[i2] > fMax ) {
	i1 = i2;
	fMax = diagonal_[i1];
      }
    }

    if ( i1 != i0 ) {
      // swap eigenvalues
      diagonal_[i1] = diagonal_[i0];
      diagonal_[i0] = fMax;

      // swap eigenvectors
      for (i2 = 0; i2 < size_; i2++) {
	T fTmp = matrix_[i2][i0];
	matrix_[i2][i0] = matrix_[i2][i1];
	matrix_[i2][i1] = fTmp;
      }
    }
  }
}

template <class T>
void 
GbEigen<T>::increasingSort ()
{
  // sort eigenvalues in increasing order, e[0] <= ... <= e[size_-1]
  for (int i0 = 0, i1; i0 <= size_-2; i0++) {
    // locate minimum eigenvalue
    i1 = i0;
    T fMin = diagonal_[i1];
    int i2;
    for (i2 = i0+1; i2 < size_; i2++) {
      if ( diagonal_[i2] < fMin ) {
	i1 = i2;
	fMin = diagonal_[i1];
      }
    }
    
    if ( i1 != i0 ) {
      // swap eigenvalues
      diagonal_[i1] = diagonal_[i0];
      diagonal_[i0] = fMin;

      // swap eigenvectors
      for (i2 = 0; i2 < size_; i2++) {
	T fTmp = matrix_[i2][i0];
	matrix_[i2][i0] = matrix_[i2][i1];
	matrix_[i2][i1] = fTmp;
      }
    }
  }
}

