// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FUNCTIONS_2001
#define FUNCTIONS_2001

#include "General/environment.hh"
#include "Matrix/matrix.hh"

namespace CoupledField
{
  //! Returns the value ln[gamma(xx)] for xx > 0
  /*!
    \param real value greater zero
  */
  Double gammaln(Double xx);


  //! Determines approximative the eigenvalues of a symmetric positive matrix MAT
  /*!
    \param S.P.D. Matrix MAT
    \param error tolerance (we just approximate the EVs up to the tolerance eps)
    \param eigen - returns EVs upwards sorted, i.e. ev_1<ev_2< ... < ev_n
  */
  void eigenValues(Matrix<Double> & mat, Double eps, Vector<Double> & eigen);
  //! Performs Given Rotations of Matrix MAT, used by eigenValues
  void givensRotation(Integer ndim, Matrix<Double> & mat);
  //! This criterion is used by method eigenValues, terminates approximation of eigenvalues.
  void terminationCriterion(Integer ndim, Matrix<Double> & mat, Double a2, Double eps2, Integer *l_conv);
  //! sorting a one dimensional array
  /*!
    \param l_lort decides about sorting order, i.e. l_sort==1 -> sorting downwards
    lsort == 0 -> sorting upwards 
  */
  void sortArray(Integer ndim, Integer l_sort, Vector<Double> & d);

  //! NC_SIMON: Normalizes a vector using the L2-norm
  Double Normalize(Vector<Double>& vec);

  //! NC_SIMON: Computes the cross product of two 3-dimensional vectors.
  //! \param a (in) left operand of the cross product
  //! \param b (in) right operand of the cross product
  //! \param c (out) receives the result of the cross product
  void CrossProd(const Vector<Double>& a,
                 const Vector<Double>& b,
                 Vector<Double> &result);

  //! Returns if vectors a and b are collinear.
  bool CoLinear(const Vector<Double> &a, const Vector<Double> &b);

  //! NC_SIMON: get barycentric coordinates of a point with respect to p1,
  //! p2 and p3
  void GetBarycentricCoords(const Vector<Double>& p1,
                            const Vector<Double>& p2,
                            const Vector<Double>& p3,
                            const Vector<Double>& p,
                            Vector<Double>& b);

} // end of namespace CoupledField

#endif  // FILE_FUNCTIONS
