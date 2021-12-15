#ifndef FILE_BICUBICINTERPOLATE
#define FILE_BICUBICINTERPOLATE

#include "ApproxData.hh"
#include "MatVec/Matrix.hh"
#include "Utils/StdVector.hh"

namespace CoupledField {

class BiCubicInterpolate : public ApproxData
{
public:
  //! constructor
  // if periodic = true, values on all boundaries have to be given (and have to match pairwise)
  BiCubicInterpolate(const StdVector<double>& data, const StdVector<double>& a, const StdVector<double>& b,
      const bool periodic = false);

  //! constructor
  BiCubicInterpolate(const StdVector<double>& data, const StdVector<double>& a, const StdVector<double>& b,
      const Matrix<double>& coeff, const bool periodic = false);

  //! destructor
  virtual ~ BiCubicInterpolate() {};

  //! computes the approximation polynom
  virtual void CalcApproximation(const bool start=true) override;

  //! returns f(x,y)
  virtual double EvaluateFunc(const double x, const double y) const override;

  virtual double EvaluateFunc(const Vector<double>& p) const override { return EvaluateFunc(p[0], p[1]); };

  //! returns d f(x,y) / d param
  double EvaluateDeriv(const double x, const double y, const int dparam) const;

  virtual double EvaluateDeriv(const Vector<double>& p, const int dparam) const override { return EvaluateDeriv(p[0], p[1], dparam); };

  //! returns grad f(x,y)
  Vector<double> EvaluatePrime(const double x, const double y) const;

  ///
  int GetSize() {return numMeas_;};

private:

  typedef enum { NONE, X, Y } Derivative;

  // approximate partial derivatives by finite differences
  void ApproxPartialDeriv(StdVector<double>& dFda, StdVector<double>& dFdb, StdVector<double>& dFdadb) const;

  // calculate coefficients of interpolating polynoms
  void CalcCoeff(StdVector<double>& coeff, const StdVector<double>& F,
      const StdVector<double>& Fda, const StdVector<double>& Fdb,
      const StdVector<double>& Fdadb) const;

  // evaluation of the interpolation polynomial at point x,y
  double EvaluatePolynom(const unsigned int index, const double x, const double y, const Derivative deriv = Derivative::NONE) const;

  /** Calculate local x, y, i.e. relative to corresponding interpolation patch
   *  @return index of corresponding interpolation patch
   */
  unsigned int GetLocalValues(double x, double y, double& xloc, double& yloc, double& dxloc, double& dyloc) const;

  // convert subscript to linear index
  unsigned int sub2ind(const unsigned int ii, const unsigned int jj) const;

  // sample points
  Vector<double> x_;
  Vector<double> y_;

  // sample values
  StdVector<double> data_;

  // coefficients of interpolating polynoms
  StdVector<StdVector<double>> coeff_;

  double offsetX_;
  double offsetY_;
  double scaleX_;
  double scaleY_;

  bool periodic_;
};

} //end of namespace


#endif
