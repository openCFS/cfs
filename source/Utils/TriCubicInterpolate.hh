#ifndef FILE_TRICUBICINTERPOLATE
#define FILE_TRICUBICINTERPOLATE

#include "ApproxData.hh"
#include "MatVec/Matrix.hh"
#include "Utils/StdVector.hh"

namespace CoupledField {

class TriCubicInterpolate : public ApproxData
{
public:
  //! constructor
  // if periodic = true, values on all boundaries have to be given (and have to match pairwise)
  TriCubicInterpolate(StdVector<double>& data, StdVector<double>& a, StdVector<double>& b, StdVector<double>& c,
      bool periodic = false);

  //! constructor
  TriCubicInterpolate(StdVector<double>& data, StdVector<double>& a, StdVector<double>& b, StdVector<double>& c,
      Matrix<double>& coeff, bool periodic = false);

  //! destructor
  virtual ~ TriCubicInterpolate() {};

  //! computes the approximation polynom
  virtual void CalcApproximation(const bool start=true) override;

  //! returns f(x,y,z)
  virtual double EvaluateFunc(const double x, const double y, const double z) const override;

  virtual double EvaluateFunc(const Vector<double>& p) const override { return EvaluateFunc(p[0], p[1], p[2]); };

  //! returns d f(x,y,z) / d param
  double EvaluateDeriv(const double x, const double y, const double z, const int dparam) const;

  virtual double EvaluateDeriv(const Vector<double>& p, const int dparam) const override { return EvaluateDeriv(p[0], p[1], p[2], dparam); };

  //! returns grad f(x,y,z)
  Vector<double> EvaluatePrime(const double x, const double y, const double z) const;

  ///
  int GetSize() {return numMeas_;};

private:

  typedef enum { NONE = -1, X, Y, Z } Derivative;

  // approximate partial derivatives by finite differences
  void ApproxPartialDeriv(StdVector<double>& dFda, StdVector<double>& dFdb, StdVector<double>& dFdc,
      StdVector<double>& dFdadb, StdVector<double>& dFdadc, StdVector<double>& dFdbdc, StdVector<double>& dFdadbdc) const;

  // calculate coefficients of interpolating polynoms
  void CalcCoeff(StdVector<double>& coeff, const StdVector<double>& F,
      const StdVector<double>& Fda, const StdVector<double>& Fdb, const StdVector<double>& Fdc,
      const StdVector<double>& Fdadb, const StdVector<double>& Fdadc, const StdVector<double>& Fdbdc,
      const StdVector<double>& Fdadbdc) const;

  // evaluation of the interpolation polynomial at point x,y
  double EvaluatePolynom(const unsigned int index, const double x, const double y, const double z, const Derivative deriv = Derivative::NONE) const;

  /** Calculate local x, y, i.e. relative to corresponding interpolation patch
   *  @return index of corresponding interpolation patch
   */
  unsigned int GetLocalValues(double x, double y, const double z,
      double& xloc, double& yloc, double& zloc, double& dxloc, double& dyloc, double& dzloc) const;

  // convert subscript to linear index
  unsigned int sub2ind(const unsigned int ii, const unsigned int jj, const unsigned int kk) const;

  // sample points
  Vector<double> x_;
  Vector<double> y_;
  Vector<double> z_;

  // sample values
  StdVector<double> data_;

  // coefficients of interpolating polynoms
  StdVector<StdVector<double>> coeff_;

  double offsetX_;
  double offsetY_;
  double offsetZ_;
  double scaleX_;
  double scaleY_;
  double scaleZ_;

  bool periodic_;
};

} //end of namespace


#endif
