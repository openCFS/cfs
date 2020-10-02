#ifndef FILE_TRICUBICINTERPOLATE
#define FILE_TRICUBICINTERPOLATE

#include "ApproxData.hh"
#include "Utils/StdVector.hh"

namespace CoupledField {

class TriCubicInterpolate : public ApproxData
{
public:
  TriCubicInterpolate(StdVector<double> data, StdVector<double> a, StdVector<double> b, StdVector<double> c, bool periodic = false);

  //! destructor
  virtual ~ TriCubicInterpolate() {};

  //computes the approximation polynom
  virtual void CalcApproximation(bool start=true);

  virtual Double EvaluateFunc(double x, double y, double z);

  //! returns grad f(x,y)
  Vector<Double> EvaluatePrime(double x, double y, double z);

  ///
  virtual void CalcBestParameter() {
    EXCEPTION("not implemented");
  }

  ///
  virtual Double EvaluateFunc(Double x) {
    EXCEPTION("EF: you need to provide x, y and z for tricubic interpolation");
    return -1.0;
  }

  ///
  virtual Double EvaluatePrime(Double x) {
    EXCEPTION("EP: you need to provide x, y and z for tricubic interpolation");
    return -1.0;
  };

  ///
  virtual Double EvaluateFuncInv(Double t) {
    EXCEPTION("not implemented");
    return -1.0;
  }

  ///
  virtual Double EvaluatePrimeInv(Double t) {
    EXCEPTION("not implemented");
    return -1.0;
  }

  ///
  int GetSize() {return numMeas_;};

private:

  typedef enum { NONE, X, Y, Z } Derivative;

  // approximate partial derivatives by finite differences
  void ApproxPartialDeriv(StdVector<double>& dFda, StdVector<double>& dFdb, StdVector<double>& dFdc, StdVector<double>& dFdadb, StdVector<double>& dFdadc, StdVector<double>& dFdbdc, StdVector<double>& dFdadbdc) const;

  // calculate coefficients of interpolating polynoms
  void CalcCoeff(StdVector<double>& coeff, StdVector<double> F, StdVector<double> Fda, StdVector<double> Fdb, StdVector<double> Fdc, StdVector<double> Fdadb, StdVector<double> Fdadc, StdVector<double> Fdbdc, StdVector<double> Fdadbdc) const;

  // evaluation of the interpolation polynomial at point x,y
  double EvaluatePolynom(unsigned int index, double x, double y, double z, Derivative deriv = Derivative::NONE) const;

  /** Calculate local x, y, i.e. relative to corresponding interpolation patch
   *  @return index of corresponding interpolation patch
   */
  unsigned int GetLocalValues(double x, double y, double z, double& xloc, double& yloc, double& zloc);

  unsigned int sub2ind(unsigned int ii, unsigned int jj, unsigned int kk) const;

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
