#ifndef FUNCTIONS_2001
#define FUNCTIONS_2001

#include "General/environment.hh"
#include "Matrix/matrix.hh"

namespace CoupledField
{
  Double gammaln(Double xx);
  //! Returns the value ln[gamma(xx)] for xx > 0
  /*!
    \param real value greater zero
  */

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


} // end of namespace CoupledField

#endif  // FILE_FUNCTIONS
