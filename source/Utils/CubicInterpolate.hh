#ifndef FILE_CUBICINTERPOLATE
#define FILE_CUBICINTERPOLATE

#include "ApproxData.hh"
#include "MatVec/Matrix.hh"
#include "Utils/StdVector.hh"

namespace CoupledField {

class CubicInterpolate : public ApproxData
{
public:
  //! constructor
  // if periodic = true, values on both boundaries have to be given (and have to match)
  CubicInterpolate(const StdVector<double>& data, const StdVector<double>& a, const bool periodic = false);

  //! constructor
  CubicInterpolate(const StdVector<double>& data, const StdVector<double>& a, const Matrix<double>& coeff, const bool periodic = false);

  //! destructor
  virtual ~ CubicInterpolate() {};

  //! computes the approximation polynom
  virtual void CalcApproximation(const bool start=true) override;

  //! returns function value at x -> f(x)
  virtual double EvaluateFunc(const double x) const override;

  virtual double EvaluateFunc(const Vector<double>& p) const override { return EvaluateFunc(p[0]); };

  //! returns first derivative -> d f(x) / dx
  virtual double EvaluateDeriv(const double x) const override;

  virtual double EvaluateDeriv(const Vector<double>& p, const int dparam) const override { return EvaluateDeriv(p[0]); };

  //! returns second derivative -> d^2 f(x) / dx^2
  virtual double EvaluateSecondDeriv(const double x) const override;

  //! returns grad f(x) = d f(x) / dx
  virtual Vector<double> EvaluatePrime(const Vector<double>& p) const override;

  ///
  int GetSize() {return numMeas_;};

private:

  typedef enum { NONE, X, XX } Derivative;

  // approximate partial derivatives by finite differences
  void ApproxPartialDeriv(StdVector<double>& dFda) const;

  // calculate coefficients of interpolating polynoms
  void CalcCoeff(StdVector<double>& coeff, const StdVector<double>& F, const StdVector<double>& Fda) const;

  // evaluation of the interpolation polynomial at point x,y
  double EvaluatePolynom(const unsigned int index, const double x, const Derivative deriv = Derivative::NONE) const;

  /** Calculate local x, y, i.e. relative to corresponding interpolation patch
   *  @return index of corresponding interpolation patch
   */
  unsigned int GetLocalValues(double x, double& xloc, double& dxloc) const;

  // sample points
  Vector<double> x_;

  // sample values
  StdVector<double> data_;

  // coefficients of interpolating polynoms
  StdVector<StdVector<double>> coeff_;

  double offsetX_;
  double scaleX_;

  bool periodic_;
};

} //end of namespace


#endif
