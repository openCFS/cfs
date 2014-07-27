#ifndef FUNCTIONS_2001
#define FUNCTIONS_2001

#include "General/Environment.hh"
#include "MatVec/Matrix.hh"

namespace CoupledField {

//! =======================================================================
//! List of mathematical functions, primarily for the use in mathParser
//! to generate signals
//! =======================================================================
  //! Generate a sinus burst signal
  Double SinBurst( Double freq, Double numPeriods,
                Double nPerFadeIn, Double nPerFadeOut,
                Double t);

  //! Fade-in function
  Double FadeIn( Double fadeInTime, Double mode, Double t);

  //! Generate a spike signal
  Double Spike( Double duration, Double t );

  //! Generate a band-filtered spike signal
  Double SpikeBPF( Double cutOff, Double slewRate, Double t );

  //! Generate a chirp signal
  
  //! This function generates a chirp signal, i.e. a signal where the
  //! frequency rises proportional to the time. Additionally,
  //! there can be fade-in/-out, defined by the number of periods 
  //! related to the center frequency.
  Double Chirp(Double sweepTime, Double startFreq, Double stopFreq, 
               Double nPerFadeIn, Double nPerFadeOut, Double t);

  Double CosPulseComb( Double freq, Double pulseWidth, Double t );

  //! Generate a square pulse signal
  typedef enum { UNI_POLAR = 0, BI_POLAR = 1 } PolarType;

  Double SquarePulse(  Double freq, Double numPeriods, Double biPolarType,
                       Double pulseWidth, Double riseTime, Double t );

  //! Modulo function
  Double Mod( Double x, Double m );

  //! Generate a Gauss curve (pdf-curve)
  Double Gauss( Double mue, Double sigma, Double normVal, Double t );

  //! Calculate cylindric bessel function of first kind
  Double BesselCylJ( Double x, Double v );

  //! Calculate cylindric bessel function of second kind
  Double BesselCylY( Double x, Double v );

  //! Calculate spherical bessel function of first kind
  Double BesselSphJ( Double x, Double v );

  //! Calculate spherical bessel function of second kind
  Double BesselSphY( Double x, Double v );

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

  /* Commented out due to refactoring.
   * Can be removed, if implementation in Vector class works.
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
*/

} // end of namespace CoupledField

#endif  // FILE_FUNCTIONS
