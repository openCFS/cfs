/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/

#ifndef  GBMATRIX2_HH
#define  GBMATRIX2_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include <iostream>
#include <typeinfo>

#include "GbVec2.hh"
#include "GbMath.hh"

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/

// NOTE.  Rotation matrices are of the form
//   R = cos(t) -sin(t)
//       sin(t)  cos(t)
// where t > 0 indicates a counterclockwise rotation in the xy-plane.

template <class T>
class GbMatrix2
{
public:
  // construction
  GbMatrix2 ();
  GbMatrix2 (const T aafEntry[2][2]);
  GbMatrix2 (const GbMatrix2<T>& rkMatrix);
  GbMatrix2 (T fEntry00, T fEntry01,
	     T fEntry10, T fEntry11);

  // destruction
  ~GbMatrix2();

  // member access, allows use of construct mat[r][c]
  INLINE T* operator[] (int iRow) const;
//  operator T* ();
  INLINE GbVec2<T> getColumn (int iCol) const;

  //! assignment and comparison
  INLINE GbMatrix2<T>& operator=  (const GbMatrix2<T>& rkMatrix);
  INLINE GbBool        operator== (const GbMatrix2<T>& rkMatrix) const;
  INLINE GbBool        operator!= (const GbMatrix2<T>& rkMatrix) const;

  //! arithmetic operations
  INLINE GbMatrix2<T> operator+ (const GbMatrix2<T>& rkMatrix) const;
  INLINE GbMatrix2<T> operator- (const GbMatrix2<T>& rkMatrix) const;
  INLINE GbMatrix2<T> operator* (const GbMatrix2<T>& rkMatrix) const;
  INLINE GbMatrix2<T> operator- () const;

  //! matrix * vector [2x2 * 2x1 = 2x1]
  INLINE GbVec2<T> operator* (const GbVec2<T>& rkVector) const;

  //! vector * matrix [1x2 * 2x2 = 1x2]
  friend GbVec2<T> operator* (const GbVec2<T>& rkVector,
			      const GbMatrix2<T>& rkMatrix);

  //! matrix * scalar
  INLINE GbMatrix2<T> operator* (T fScalar) const;
  INLINE GbMatrix2<T> operator* (int s) const;

  //! scalar * matrix
  friend GbMatrix2<T> operator* (T fScalar,
				 const GbMatrix2<T>& rkMatrix);

  //! utilities
  GbMatrix2<T> transpose () const;
  // tolerance max. 1e-06
  GbBool inverse (GbMatrix2<T>& rkInverse, T fTolerance = std::numeric_limits<T>::epsilon()) const;
  GbMatrix2<T> inverse (T fTolerance = std::numeric_limits<T>::epsilon()) const;
  T determinant () const;

  //! singular value decomposition
  void singularValueDecomposition (GbMatrix2<T>& rkL, GbVec2<T>& rkS, GbMatrix2<T>& rkR) const;
  void singularValueComposition (const GbMatrix2<T>& rkL, const GbVec2<T>& rkS, const GbMatrix2<T>& rkR);

  //! Gram-Schmidt orthonormalization (applied to columns of rotation matrix)
  void orthonormalize ();

  //! orthogonal Q, diagonal D, upper triangular U stored as (u01)
  void QDUDecomposition (GbMatrix2<T>& rkQ, GbVec2<T>& rkD, T& rfU) const;

  T spectralNorm () const;
  T L2NormSquare () const;

  //! matrix must be orthonormal for these methods to work
  INLINE void toAngle (T& rfRadians) const;
  INLINE void fromAngle (T fRadians);

  // eigensolver, matrix must be symmetric
  void eigenSolveSymmetric (T afEigenvalue[2], GbVec2<T> akEigenvector[2]) const;

  static void tensorProduct (const GbVec2<T>& rkU, const GbVec2<T>& rkV, GbMatrix2<T>& rkProduct);

  static const T EPSILON;
  static const GbMatrix2<T> ZERO;
  static const GbMatrix2<T> IDENTITY;

  //! This operator displays the coordinate values of the matrix
  friend std::ostream& operator<<(std::ostream&, const GbMatrix2<T>&);

private:
  //! support for eigensolver
  void tridiagonal (T afDiag[2], T afSubDiag[2]);
  GbBool QLAlgorithm (T afDiag[2], T afSubDiag[2]);

  T entry_[2][2];
};

template<class T>
GbMatrix2<T> 
operator* (T fScalar, const GbMatrix2<T>& rkMatrix)
{
  GbMatrix2<T> kProd;
  for (int iRow = 0; iRow < 2; iRow++) {
    for (int iCol = 0; iCol < 2; iCol++)
      kProd[iRow][iCol] = fScalar*rkMatrix.entry_[iRow][iCol];
  }
  return kProd;
}

template<class T>
GbVec2<T> 
operator* (const GbVec2<T>& rkVector, const GbMatrix2<T>& rkMatrix)
{
  GbVec2<T> kProd;
  for (int iRow = 0; iRow < 2; iRow++) {
    kProd[iRow] = rkVector[0]*rkMatrix[0][iRow] + rkVector[1]*rkMatrix[1][iRow];
  }
  return kProd;
}

template<class T>
std::ostream&
operator<<(std::ostream& s, const GbMatrix2<T>& v)
{
  s<<typeid(v).name()<<":";
  s<<"\n/ "<<v[0][0]<<" "<<v[1][0]<<" \\";
  s<<"\n\\ "<<v[0][1]<<" "<<v[1][1]<<" /"<<std::endl;
  return s;
}

#ifdef __GNUG__

#include "GbMatrix2.in"
#include "GbMatrix2.T"

#else

#pragma instantiate GbMatrix2<float>
#pragma instantiate GbMatrix2<double>
#pragma instantiate GbVec2<float>  operator* (const GbVec2<float>& v, const GbMatrix2<float>& m)
#pragma instantiate GbVec2<double> operator* (const GbVec2<double>& v, const GbMatrix2<double>& m)
#pragma instantiate std::ostream& operator<<(std::ostream&, const GbMatrix2<float>&)
#pragma instantiate std::ostream& operator<<(std::ostream&, const GbMatrix2<double>&)

#ifndef OUTLINE
#include "GbMatrix2.in"
#endif

#endif  // g++

#endif  // GBMATRIX2_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.1  2002/02/22 14:47:56  elena
| new: dir Gridlib_inc
|
| Revision 1.1  2001/06/26 08:31:22  prkipfer
| added 2D Projections
|
|
+---------------------------------------------------------------------*/
