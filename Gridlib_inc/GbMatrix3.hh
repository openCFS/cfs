/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/

// Simple 3x3 Matrix class
// This is intended for use with floating point values, therefore please
// note the special treatment of integer scalars

// The (x,y,z) coordinate system is assumed to be right-handed.
// Coordinate axis rotation matrices are of the form
//   RX =    1       0       0
//           0     cos(t) -sin(t)
//           0     sin(t)  cos(t)
// where t > 0 indicates a counterclockwise rotation in the yz-plane
//   RY =  cos(t)    0     sin(t)
//           0       1       0
//        -sin(t)    0     cos(t)
// where t > 0 indicates a counterclockwise rotation in the zx-plane
//   RZ =  cos(t) -sin(t)    0
//         sin(t)  cos(t)    0
//           0       0       1
// where t > 0 indicates a counterclockwise rotation in the xy-plane.

#ifndef  GBMATRIX3_HH
#define  GBMATRIX3_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include <iostream>
#include <typeinfo>

#include "GbVec3.hh"
#include "GbMath.hh"

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/

template <class T>
class GbMatrix3
{
public:
  //! Various Constructors
  GbMatrix3 ();
  GbMatrix3 (const T entry[3][3]);
  GbMatrix3 (const T entry[9]);
  GbMatrix3 (const GbMatrix3<T>& m);
  GbMatrix3 (T entry00, T entry01, T entry02,
	     T entry10, T entry11, T entry12,
	     T entry20, T entry21, T entry22);
  explicit GbMatrix3 (T s);

  //! Destructor
  ~GbMatrix3();

  //! Member access, allows use of construct mat[r][c]
  INLINE T* operator[] (int row) const;
  INLINE GbVec3<T> getColumn (int col) const;

  //! Assignment and comparison
//!!  INLINE GbMatrix3<T>& operator=  (const GbMatrix3<T>& m);
  INLINE GbBool        operator== (const GbMatrix3<T>& m) const;
  INLINE GbBool        operator!= (const GbMatrix3<T>& m) const;
  INLINE GbMatrix3<T>& operator+= (const GbMatrix3<T>& m);
  INLINE GbMatrix3<T>& operator-= (const GbMatrix3<T>& m);
  INLINE GbMatrix3<T>& operator*= (const GbMatrix3<T>& m);

  //! Arithmetic operations
  INLINE GbMatrix3<T> operator+ (const GbMatrix3<T>& m) const;
  INLINE GbMatrix3<T> operator- (const GbMatrix3<T>& m) const;
  INLINE GbMatrix3<T> operator* (const GbMatrix3<T>& m) const;
  INLINE GbMatrix3<T> operator- () const;

  //! Matrix(self) * vector [3x3 * 3x1 = 3x1]
  INLINE GbVec3<T> operator* (const GbVec3<T>& v) const;

  //! Vector * matrix [1x3 * 3x3 = 1x3]
  friend GbVec3<T> operator*(const GbVec3<T>& v, const GbMatrix3<T>& m);

  //! Matrix * scalar
  INLINE GbMatrix3<T> operator* (const T& s) const;
  INLINE GbMatrix3<T> operator* (const int& s) const;

  //! Scalar * matrix
  friend GbMatrix3<T> operator*(const T& s, const GbMatrix3<T>& m);

  //! M0.transposeTimes(M1) = M0^t*M1 where M0^t is the transpose of M0
  INLINE GbMatrix3<T> transposeTimes (const GbMatrix3<T>& rkM) const;

  //! M0.timesTranspose(M1) = M0*M1^t where M1^t is the transpose of M1
  INLINE GbMatrix3<T> timesTranspose (const GbMatrix3<T>& rkM) const;

  //! Utilities
  INLINE GbMatrix3<T> transpose () const;
  // tolerance max. 1e-06
  GbBool inverse (GbMatrix3<T>& inv, T tolerance = std::numeric_limits<T>::epsilon()) const;
  GbMatrix3<T> inverse (T tolerance = std::numeric_limits<T>::epsilon()) const;
  INLINE T determinant () const;

  //! SLERP (spherical linear interpolation) without quaternions.  Computes
  //! R(t) = R0*(transpose(R0)*R1)^t.  If Q is a rotation matrix with
  //! unit-length axis U and angle A, then Q^t is a rotation matrix with
  //! unit-length axis U and rotation angle t*A.
  static GbMatrix3<T> slerp (T fT, const GbMatrix3<T>& rkR0, const GbMatrix3<T>& rkR1);

  //! Singular value decomposition
  void singularValueDecomposition (GbMatrix3<T>& l, GbVec3<T>& s, GbMatrix3<T>& r) const;
  void singularValueComposition (const GbMatrix3<T>& l, const GbVec3<T>& s, const GbMatrix3<T>& r);

  //! Gram-Schmidt orthonormalization (applied to columns of rotation matrix)
  void orthonormalize ();

  //! Orthogonal Q, diagonal D, upper triangular U stored as (u01,u02,u12)
  void QDUDecomposition (GbMatrix3<T>& q, GbVec3<T>& d, GbVec3<T>& u) const;

  T spectralNorm () const;

  //! Matrix must be orthonormal for these methods to work
  void toAxisAngle (GbVec3<T>& axis, T& radians) const;
  void fromAxisAngle (const GbVec3<T>& axis, T radians);

  //! The matrix must be orthonormal.  The decomposition is yaw*pitch*roll
  //! where yaw is rotation about the Up vector, pitch is rotation about the
  //! Right axis, and roll is rotation about the Direction axis.
  GbBool toEulerAnglesXYZ (T& angle_y, T& angle_p, T& angle_r);
  GbBool toEulerAnglesXZY (T& angle_y, T& angle_p, T& angle_r);
  GbBool toEulerAnglesYXZ (T& angle_y, T& angle_p, T& angle_r);
  GbBool toEulerAnglesYZX (T& angle_y, T& angle_p, T& angle_r);
  GbBool toEulerAnglesZXY (T& angle_y, T& angle_p, T& angle_r);
  GbBool toEulerAnglesZYX (T& angle_y, T& angle_p, T& angle_r);
  void fromEulerAnglesXYZ (T angle_y, T angle_p, T angle_r);
  void fromEulerAnglesXZY (T angle_y, T angle_p, T angle_r);
  void fromEulerAnglesYXZ (T angle_y, T angle_p, T angle_r);
  void fromEulerAnglesYZX (T angle_y, T angle_p, T angle_r);
  void fromEulerAnglesZXY (T angle_y, T angle_p, T angle_r);
  void fromEulerAnglesZYX (T angle_y, T angle_p, T angle_r);

  //! Eigensolver, matrix must be symmetric
  void eigenSolveSymmetric (T eigenvalue[3], GbVec3<T> eigenvector[3]) const;

  static void tensorProduct (const GbVec3<T>& u, const GbVec3<T>& v, GbMatrix3<T>& product);

  static const T EPSILON;
  static const GbMatrix3<T> ZERO;
  static const GbMatrix3<T> IDENTITY;

  //! This operator displays and reads the coordinate values of the matrix
  friend std::ostream& operator<<(std::ostream&, const GbMatrix3<T>&);
  friend std::istream& operator>>(std::istream&, GbMatrix3<T>&);

private:
  //! Support for eigensolver
  void tridiagonal (T diag[3], T subDiag[3]);
  GbBool QLAlgorithm (T diag[3], T subDiag[3]);

  //! Support for singular value decomposition
  static const T svdEpsilon_;
  static const int svdMaxIterations_;
  static void bidiagonalize  (GbMatrix3<T>& A, GbMatrix3<T>& L, GbMatrix3<T>& R);
  static void GolubKahanStep (GbMatrix3<T>& A, GbMatrix3<T>& L, GbMatrix3<T>& R);

  //! Support for spectral norm
  static T maxCubicRoot (T coeff[3]);

  T entry_[3][3];
};

template<class T>
GbMatrix3<T>
operator* (const T& s, const GbMatrix3<T>& m)
{
  T kProd[9];
  int k=0;

  for (int iRow = 0; iRow < 3; ++iRow) {
    for (int iCol = 0; iCol < 3; ++iCol)
      kProd[k++] = s*m.entry_[iRow][iCol];
  }
  return GbMatrix3<T>(kProd);
}

template<class T>
GbVec3<T> 
operator* (const GbVec3<T>& v, const GbMatrix3<T>& m)
{
  T kProd[3];
  for (int iRow = 0; iRow < 3; ++iRow) {
    kProd[iRow] = v[0]*m[0][iRow] + v[1]*m[1][iRow] + v[2]*m[2][iRow];
  }
  return GbVec3<T>(kProd);
}

template<class T>
std::ostream&
operator<<(std::ostream& s, const GbMatrix3<T>& v)
{
  s<<"[("<<v[0][0]<<", "<<v[1][0]<<", "<<v[2][0]<<") ";
  s<<"("<<v[0][1]<<", "<<v[1][1]<<", "<<v[2][1]<<") ";
  s<<"("<<v[0][2]<<", "<<v[1][2]<<", "<<v[2][2]<<")]";
  return s;
}

template<class T>
std::istream&
operator>>(std::istream& s, GbMatrix3<T>& v)
{
  char c;
  char dummy[3];
  T x,y,z;
  
  s>>c>>c>>x>>dummy>>y>>dummy>>z>>c>>c;
  v[0][0]=x; v[1][0]=y; v[2][0]=z;
  s>>c>>x>>dummy>>y>>dummy>>z>>c>>c;
  v[0][1]=x; v[1][1]=y; v[2][1]=z;
  s>>c>>x>>dummy>>y>>dummy>>z>>c>>c;
  v[0][2]=x; v[1][2]=y; v[2][2]=z;
  return s;
}

#ifdef __GNUG__

#include "GbMatrix3.in"
#include "GbMatrix3.T"

#else

#pragma instantiate GbMatrix3<float>
#pragma instantiate GbMatrix3<double>
#pragma instantiate GbVec3<float>  operator* (const GbVec3<float>& v, const GbMatrix3<float>& m)
#pragma instantiate GbVec3<double> operator* (const GbVec3<double>& v, const GbMatrix3<double>& m)
#pragma instantiate std::ostream& operator<<(std::ostream&, const GbMatrix3<float>&)
#pragma instantiate std::ostream& operator<<(std::ostream&, const GbMatrix3<double>&)
#pragma instantiate std::istream& operator>>(std::istream&, GbMatrix3<float>&)
#pragma instantiate std::istream& operator>>(std::istream&, GbMatrix3<double>&)

#ifndef OUTLINE
#include "GbMatrix3.in"
#endif

#endif  // g++

#endif  // GBMATRIX3_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.3  2002/04/03 11:17:07  elena
| new: changes in CFS++ on SGI and new function for gridlib
|
| Revision 1.5  2001/09/12 09:14:35  prkipfer
| made constructor explicit to avoid implicit type conversion
|
| Revision 1.4  2001/08/16 16:53:18  prkipfer
| improved type safety for template parameter
|
| Revision 1.3  2001/06/15 08:35:10  prkipfer
| added fast transpose multiplication and spherical interpolation without quaternions
|
| Revision 1.2  2001/02/13 10:58:39  prkipfer
| exchanged fixed numeric tolerances by machine dependent ones
|
| Revision 1.1  2001/01/02 14:44:02  prkipfer
| introduced new classes
|
|
+---------------------------------------------------------------------*/
