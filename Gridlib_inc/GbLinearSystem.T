// -*- C++ -*- automatisch in C++mode wechseln (emacs)
/*----------------------------------------------------------------------
|
|
| $Log$
| Revision 1.3  2002/04/03 11:17:07  elena
| new: changes in CFS++ on SGI and new function for gridlib
|
| Revision 1.5  2001/08/16 16:53:16  prkipfer
| improved type safety for template parameter
|
| Revision 1.4  2001/06/20 16:14:03  prkipfer
| added compiler hints
|
| Revision 1.3  2001/04/23 10:08:28  prkipfer
| fixed typo
|
| Revision 1.2  2001/02/13 10:58:39  prkipfer
| exchanged fixed numeric tolerances by machine dependent ones
|
| Revision 1.1  2001/01/03 15:30:09  prkipfer
| introduced intersection, linear system solver and minimizer
|
|
|
+---------------------------------------------------------------------*/

#ifdef LINUX
#include <cassert>
#endif

template <class T>
const T 
GbLinearSystem<T>::TOLERANCE = std::numeric_limits<T>::epsilon(); // max. 1e-06

template <class T>
T** 
GbLinearSystem<T>::newMatrix (int iSize)
{
  T** aafA = new T*[iSize];
  assert( aafA );

  for (int iRow = 0; iRow < iSize; iRow++) {
    aafA[iRow] = new T[iSize];
    memset(aafA[iRow],0,iSize*sizeof(T));
  }

  return aafA;
}

template <class T>
void 
GbLinearSystem<T>::deleteMatrix (int iSize, T** aafA)
{
  for (int iRow = 0; iRow < iSize; iRow++)
    delete[] aafA[iRow];
  delete[] aafA;
}

template <class T>
T* 
GbLinearSystem<T>::newVector (int iSize)
{
  T* afB = new T[iSize];
  assert( afB );
  memset(afB,0,iSize*sizeof(T));
  return afB;
}

template <class T>
void 
GbLinearSystem<T>::deleteVector (T* afB)
{
  delete[] afB;
}

template <class T>
GbBool 
GbLinearSystem<T>::solve2 (const T aafA[2][2], const T afB[2], T afX[2])
{
  T fDet = aafA[0][0]*aafA[1][1]-aafA[0][1]*aafA[1][0];
  if ( GbMath<T>::Abs(fDet) < GbLinearSystem<T>::TOLERANCE )
    return false;

  T fInvDet = T(1.0)/fDet;
  afX[0] = (aafA[1][1]*afB[0]-aafA[0][1]*afB[1])*fInvDet;
  afX[1] = (aafA[0][0]*afB[1]-aafA[1][0]*afB[0])*fInvDet;
  return true;
}

template <class T>
GbBool 
GbLinearSystem<T>::solve3 (const T aafA[3][3], const T afB[3], T afX[3])
{
  T aafAInv[3][3];
  aafAInv[0][0] = aafA[1][1]*aafA[2][2]-aafA[1][2]*aafA[2][1];
  aafAInv[0][1] = aafA[0][2]*aafA[2][1]-aafA[0][1]*aafA[2][2];
  aafAInv[0][2] = aafA[0][1]*aafA[1][2]-aafA[0][2]*aafA[1][1];
  aafAInv[1][0] = aafA[1][2]*aafA[2][0]-aafA[1][0]*aafA[2][2];
  aafAInv[1][1] = aafA[0][0]*aafA[2][2]-aafA[0][2]*aafA[2][0];
  aafAInv[1][2] = aafA[0][2]*aafA[1][0]-aafA[0][0]*aafA[1][2];
  aafAInv[2][0] = aafA[1][0]*aafA[2][1]-aafA[1][1]*aafA[2][0];
  aafAInv[2][1] = aafA[0][1]*aafA[2][0]-aafA[0][0]*aafA[2][1];
  aafAInv[2][2] = aafA[0][0]*aafA[1][1]-aafA[0][1]*aafA[1][0];
  T fDet = aafA[0][0]*aafAInv[0][0] + aafA[0][1]*aafAInv[1][0] + aafA[0][2]*aafAInv[2][0];
  if ( GbMath<T>::Abs(fDet) < GbLinearSystem<T>::TOLERANCE )
    return false;

  T fInvDet = T(1.0)/fDet;
  for (int iRow = 0; iRow < 3; iRow++) {
    for (int iCol = 0; iCol < 3; iCol++)
      aafAInv[iRow][iCol] *= fInvDet;
  }

  afX[0] = aafAInv[0][0]*afB[0]+aafAInv[0][1]*afB[1]+aafAInv[0][2]*afB[2];
  afX[1] = aafAInv[1][0]*afB[0]+aafAInv[1][1]*afB[1]+aafAInv[1][2]*afB[2];
  afX[2] = aafAInv[2][0]*afB[0]+aafAInv[2][1]*afB[1]+aafAInv[2][2]*afB[2];
  return true;
}

template <class T>
GbBool 
GbLinearSystem<T>::inverse (int iSize, T** aafA)
{
  int* aiColIndex = new int[iSize];
  assert( aiColIndex );

  int* aiRowIndex = new int[iSize];
  assert( aiRowIndex );

  GbBool* abPivoted = new GbBool[iSize];
  assert( abPivoted );
  memset(abPivoted,0,iSize*sizeof(GbBool));

  int i1, i2, iRow, iCol;
  T fSave;

  // elimination by full pivoting
  for (int i0 = 0; i0 < iSize; i0++) {
    // search matrix (excluding pivoted rows) for maximum absolute entry
    T fMax = T(0.0);
    for (i1 = 0; i1 < iSize; i1++) {
      if ( !abPivoted[i1] ) {
	for (i2 = 0; i2 < iSize; i2++) {
	  if ( !abPivoted[i2] ) {
	    T fAbs = GbMath<T>::Abs(aafA[i1][i2]);
	    if ( fAbs > fMax ) {
	      fMax = fAbs;
	      iRow = i1;
	      iCol = i2;
	    }
	  }
	}
      }
    }

    if ( fMax == T(0.0) ) {
      // matrix is not invertible
      delete[] aiColIndex;
      delete[] aiRowIndex;
      delete[] abPivoted;
      return false;
    }

    abPivoted[iCol] = true;

    // swap rows so that A[iCol][iCol] contains the pivot entry
    if ( iRow != iCol ) {
      T* afRowPtr = aafA[iRow];
      aafA[iRow] = aafA[iCol];
      aafA[iCol] = afRowPtr;
    }

    // keep track of the permutations of the rows
    aiRowIndex[i0] = iRow;
    aiColIndex[i0] = iCol;

    // scale the row so that the pivot entry is 1
    T fInv = T(1.0)/aafA[iCol][iCol];
    aafA[iCol][iCol] = T(1.0);
    for (i2 = 0; i2 < iSize; i2++)
      aafA[iCol][i2] *= fInv;

    // zero out the pivot column locations in the other rows
    for (i1 = 0; i1 < iSize; i1++) {
      if ( i1 != iCol ) {
	fSave = aafA[i1][iCol];
	aafA[i1][iCol] = T(0.0);
	for (i2 = 0; i2 < iSize; i2++)
	  aafA[i1][i2] -= aafA[iCol][i2]*fSave;
      }
    }
  }

  // reorder rows so that A[][] stores the inverse of the original matrix
  for (i1 = iSize-1; i1 >= 0; i1--) {
    if ( aiRowIndex[i1] != aiColIndex[i1] ) {
      for (i2 = 0; i2 < iSize; i2++) {
	fSave = aafA[i2][aiRowIndex[i1]];
	aafA[i2][aiRowIndex[i1]] = aafA[i2][aiColIndex[i1]];
	aafA[i2][aiColIndex[i1]] = fSave;
      }
    }
  }

  delete[] aiColIndex;
  delete[] aiRowIndex;
  delete[] abPivoted;
  return true;
}

template <class T>
GbBool 
GbLinearSystem<T>::solve (int iSize, T** aafA, T* afB)
{
  int* aiColIndex = new int[iSize];
  assert( aiColIndex );

  int* aiRowIndex = new int[iSize];
  assert( aiRowIndex );

  GbBool* abPivoted = new GbBool[iSize];
  assert( abPivoted );
  memset(abPivoted,0,iSize*sizeof(GbBool));

  int i1, i2, iRow, iCol;
  T fSave;

  // elimination by full pivoting
  for (int i0 = 0; i0 < iSize; i0++) {
    // search matrix (excluding pivoted rows) for maximum absolute entry
    T fMax = T(0.0);
    for (i1 = 0; i1 < iSize; i1++) {
      if ( !abPivoted[i1] ) {
	for (i2 = 0; i2 < iSize; i2++) {
	  if ( !abPivoted[i2] ) {
	    T fAbs = GbMath<T>::Abs(aafA[i1][i2]);
	    if ( fAbs > fMax ) {
	      fMax = fAbs;
	      iRow = i1;
	      iCol = i2;
	    }
	  }
	}
      }
    }

    if ( fMax == T(0.0) ) {
      // matrix is not invertible
      delete[] aiColIndex;
      delete[] aiRowIndex;
      delete[] abPivoted;
      return false;
    }

    abPivoted[iCol] = true;

    // swap rows so that A[iCol][iCol] contains the pivot entry
    if ( iRow != iCol ) {
      T* afRowPtr = aafA[iRow];
      aafA[iRow] = aafA[iCol];
      aafA[iCol] = afRowPtr;

      fSave = afB[iRow];
      afB[iRow] = afB[iCol];
      afB[iCol] = fSave;
    }

    // keep track of the permutations of the rows
    aiRowIndex[i0] = iRow;
    aiColIndex[i0] = iCol;

    // scale the row so that the pivot entry is 1
    T fInv = T(1.0)/aafA[iCol][iCol];
    aafA[iCol][iCol] = T(1.0);
    for (i2 = 0; i2 < iSize; i2++)
      aafA[iCol][i2] *= fInv;
    afB[iCol] *= fInv;

    // zero out the pivot column locations in the other rows
    for (i1 = 0; i1 < iSize; i1++) {
      if ( i1 != iCol ) {
	fSave = aafA[i1][iCol];
	aafA[i1][iCol] = T(0.0);
	for (i2 = 0; i2 < iSize; i2++)
	  aafA[i1][i2] -= aafA[iCol][i2]*fSave;
	afB[i1] -= afB[iCol]*fSave;
      }
    }
  }

  // reorder rows so that A[][] stores the inverse of the original matrix
  for (i1 = iSize-1; i1 >= 0; i1--) {
    if ( aiRowIndex[i1] != aiColIndex[i1] ) {
      for (i2 = 0; i2 < iSize; i2++) {
	fSave = aafA[i2][aiRowIndex[i1]];
	aafA[i2][aiRowIndex[i1]] = aafA[i2][aiColIndex[i1]];
	aafA[i2][aiColIndex[i1]] = fSave;
      }
    }
  }

  delete[] aiColIndex;
  delete[] aiRowIndex;
  delete[] abPivoted;
  return true;
}

template <class T>
GbBool 
GbLinearSystem<T>::solveTri (int iSize, T* afA, T* afB,
			     T* afC, T* afR, T* afU)
{
  if ( afB[0] == T(0.0) )
    return false;

  T* afD = new T[iSize-1];
  assert( afD );

  T fE = afB[0];
  T fInvE = T(1.0)/fE;
  afU[0] = afR[0]*fInvE;

  int i0, i1;
  for (i0 = 0, i1 = 1; i1 < iSize; i0++, i1++) {
    afD[i0] = afC[i0]*fInvE;
    fE = afB[i1] - afA[i0]*afD[i0];
    if ( fE == T(0.0) ) {
      delete[] afD;
      return false;
    }
    fInvE = T(1.0)/fE;
    afU[i1] = (afR[i1] - afA[i0]*afU[i0])*fInvE;
  }

  for (i0 = iSize-1, i1 = iSize-2; i1 >= 0; i0--, i1--)
    afU[i1] -= afD[i1]*afU[i0];

  delete[] afD;
  return true;
}

template <class T>
GbBool 
GbLinearSystem<T>::solveConstTri (int iSize, T fA, T fB,
				  T fC, T* afR, T* afU)
{
  if ( fB == T(0.0) )
    return false;

  T* afD = new T[iSize-1];
  assert( afD );

  T fE = fB;
  T fInvE = T(1.0)/fE;
  afU[0] = afR[0]*fInvE;

  int i0, i1;
  for (i0 = 0, i1 = 1; i1 < iSize; i0++, i1++) {
    afD[i0] = fC*fInvE;
    fE = fB - fA*afD[i0];
    if ( fE == T(0.0) ) {
      delete[] afD;
      return false;
    }
    fInvE = T(1.0)/fE;
    afU[i1] = (afR[i1] - fA*afU[i0])*fInvE;
  }
  for (i0 = iSize-1, i1 = iSize-2; i1 >= 0; i0--, i1--)
    afU[i1] -= afD[i1]*afU[i0];

  delete[] afD;
  return true;
}

template <class T>
GbBool 
GbLinearSystem<T>::solveSymmetric (int iSize, T** aafA, T* afB)
{
  // A = L D L^t decomposition with diagonal terms of L equal to 1.  The
  // algorithm stores D terms in A[i][i] and off-diagonal L terms in
  // A[i][j] for i > j.  (G. Golub and C. Van Loan, Matrix Computations)

  int i0, i1;
  T* afV = new T[iSize];
  assert( afV );

  for (i1 = 0; i1 < iSize; i1++) {
    for (i0 = 0; i0 < i1; i0++)
      afV[i0] = aafA[i1][i0]*aafA[i0][i0];

    afV[i1] = aafA[i1][i1];
    for (i0 = 0; i0 < i1; i0++)
      afV[i1] -= aafA[i1][i0]*afV[i0];

    aafA[i1][i1] = afV[i1];
    if ( GbMath<T>::Abs(afV[i1]) <= GbLinearSystem<T>::TOLERANCE ) {
      delete[] afV;
      return false;
    }

    T fInv = T(1.0)/afV[i1];
    for (i0 = i1+1; i0 < iSize; i0++) {
      for (int i2 = 0; i2 < i1; i2++)
	aafA[i0][i1] -= aafA[i0][i2]*afV[i2];
      aafA[i0][i1] *= fInv;
    }
  }
  delete[] afV;

  // Solve Ax = B.

  // Forward substitution:  Let z = DL^t x, then Lz = B.  Algorithm
  // stores z terms in B vector.
  for (i0 = 0; i0 < iSize; i0++) {
    for (i1 = 0; i1 < i0; i1++)
      afB[i0] -= aafA[i0][i1]*afB[i1];

  }

  // Diagonal division:  Let y = L^t x, then Dy = z.  Algorithm stores
  // y terms in B vector.
  for (i0 = 0; i0 < iSize; i0++) {
    if ( GbMath<T>::Abs(aafA[i0][i0]) <= GbLinearSystem<T>::TOLERANCE )
      return false;
    afB[i0] /= aafA[i0][i0];
  }

  // Back substitution:  Solve L^t x = y.  Algorithm stores x terms in
  // B vector.
  for (i0 = iSize-2; i0 >= 0; i0--) {
    for (i1 = i0+1; i1 < iSize; i1++)
      afB[i0] -= aafA[i1][i0]*afB[i1];
  }

  return true;
}

template <class T>
GbBool 
GbLinearSystem<T>::symmetricInverse (int iSize, T** aafA, T** aafAInv)
{
  // Same algorithm as SolveSymmetric, but applied simultaneously to
  // columns of identity matrix.

  T* afV = new T[iSize];
  assert( afV );

  int i0, i1;
  for (i0 = 0; i0 < iSize; i0++) {
    for (i1 = 0; i1 < iSize; i1++)
      aafAInv[i0][i1] = ( i0 != i1 ? T(0.0) : T(1.0) );
  }

  for (i1 = 0; i1 < iSize; i1++) {
    for (i0 = 0; i0 < i1; i0++)
      afV[i0] = aafA[i1][i0]*aafA[i0][i0];

    afV[i1] = aafA[i1][i1];
    for (i0 = 0; i0 < i1; i0++)
      afV[i1] -= aafA[i1][i0]*afV[i0];

    aafA[i1][i1] = afV[i1];
    for (i0 = i1+1; i0 < iSize; i0++) {
      for (int i2 = 0; i2 < i1; i2++)
	aafA[i0][i1] -= aafA[i0][i2]*afV[i2];
      aafA[i0][i1] /= afV[i1];
    }
  }
  delete[] afV;

  for (int iCol = 0; iCol < iSize; iCol++) {
    // forward substitution
    for (i0 = 0; i0 < iSize; i0++) {
      for (i1 = 0; i1 < i0; i1++)
	aafAInv[i0][iCol] -= aafA[i0][i1]*aafAInv[i1][iCol];
    }

    // diagonal division
    for (i0 = 0; i0 < iSize; i0++) {
      if ( GbMath<T>::Abs(aafA[i0][i0]) <= GbLinearSystem<T>::TOLERANCE )
	return false;
      aafAInv[i0][iCol] /= aafA[i0][i0];
    }

    // back substitution
    for (i0 = iSize-2; i0 >= 0; i0--) {
      for (i1 = i0+1; i1 < iSize; i1++)
	aafAInv[i0][iCol] -= aafA[i1][i0]*aafAInv[i1][iCol];
    }
  }

  return true;
}
