/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/

#ifndef  GBLINEARSYSTEM_HH
#define  GBLINEAESYSTEM_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include <iostream>
#include <typeinfo>

#include "GbTypes.hh"
#include "GbMath.hh"

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/

template <class T>
class GbLinearSystem
{
public:
  GbLinearSystem () { /**/ }

  // 2x2 and 3x3 systems (avoids overhead of Gaussian elimination)
  static GbBool solve2 (const T aafA[2][2], const T afB[2], T afX[2]);
  static GbBool solve3 (const T aafA[3][3], const T afB[3], T afX[3]);

  // convenience for allocation, memory is zeroed out
  static T** newMatrix (int iSize);
  static void deleteMatrix (int iSize, T** aafA);
  static T* newVector (int iSize);
  static void deleteVector (T* afB);

  /*! Input:
      A[iSize][iSize], entries are A[row][col]
      Output:
      return value is TRUE if successful, FALSE if pivoting failed
      A[iSize][iSize], inverse matrix
  */
  static GbBool inverse (int iSize, T** aafA);

  /*! Input:
      A[iSize][iSize] coefficient matrix, entries are A[row][col]
      B[iSize] vector, entries are B[row]
      Output:
      return value is TRUE if successful, FALSE if pivoting failed
      A[iSize][iSize] is inverse matrix
      B[iSize] is solution x to Ax = B
  */
  static GbBool solve (int iSize, T** aafA, T* afB);

  /*! Input:
      Matrix is tridiagonal.
      Lower diagonal A[iSize-1]
      Main  diagonal B[iSize]
      Upper diagonal C[iSize-1]
      Right-hand side R[iSize]
      Output:
      return value is TRUE if successful, FALSE if pivoting failed
      U[iSize] is solution
  */
  static GbBool solveTri (int iSize, T* afA, T* afB, T* afC, T* afR, T* afU);

  /*! Input:
      Matrix is tridiagonal.
      Lower diagonal is constant, A
      Main  diagonal is constant, B
      Upper diagonal is constant, C
      Right-hand side Rr[iSize]
      Output:
      return value is TRUE if successful, FALSE if pivoting failed
      U[iSize] is solution
  */
  static GbBool solveConstTri (int iSize, T fA, T fB, T fC, T* afR, T* afU);

  /*! Input:
      A[iSize][iSize] symmetric matrix, entries are A[row][col]
      B[iSize] vector, entries are B[row]
      Output:
      return value is TRUE if successful, FALSE if (nearly) singular
      decomposition A = L D L^t (diagonal terms of L are all 1)
      A[i][i] = entries of diagonal D
      A[i][j] for i > j = lower triangular part of L
      B[iSize] is solution to x to Ax = B
  */
  static GbBool solveSymmetric (int iSize, T** aafA, T* afB);

  /*! Input:
      A[iSize][iSize], entries are A[row][col]
      Output:
      return value is TRUE if successful, FALSE if algorithm failed
      AInv[iSize][iSize], inverse matrix
  */
  static GbBool symmetricInverse (int iSize, T** aafA, T** aafAInv);

  //! tolerance for 2x2 and 3x3 system solving
  static T TOLERANCE;
};



#ifdef __GNUG__

#include "GbLinearSystem.in"
#include "GbLinearSystem.T"

#else

#pragma instantiate GbLinearSystem<float>
#pragma instantiate GbLinearSystem<double>

#ifndef OUTLINE
#include "GbLinearSystem.in"
#endif

#endif  // g++

#endif  // GBLINEARSYSTEM_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.1  2002/02/22 14:47:56  elena
| new: dir Gridlib_inc
|
| Revision 1.2  2001/04/23 10:08:28  prkipfer
| fixed typo
|
| Revision 1.1  2001/01/03 15:30:09  prkipfer
| introduced intersection, linear system solver and minimizer
|
|
+---------------------------------------------------------------------*/
