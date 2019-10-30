#ifndef FILE_BICUBICINTERPOLATE
#define FILE_BICUBICINTERPOLATE

#include "ApproxData.hh"
#include "Utils/StdVector.hh"

namespace CoupledField {

class BiCubicInterpolate : public ApproxData
{
public:
  BiCubicInterpolate(StdVector<double> data, StdVector<double> a, StdVector<double> b, bool periodic = false);

  //! destructor
  virtual ~ BiCubicInterpolate() {};

  //computes the approximation polynom
  virtual void CalcApproximation(bool start=true);

  virtual Double EvaluateFunc(double x, double y);

  //! returns grad f(x,y)
  Vector<Double> EvaluatePrime(double x, double y);

  ///
  virtual void CalcBestParameter() {
    EXCEPTION("not implemented");
  }

  ///
  virtual Double EvaluateFunc(Double x) {
    EXCEPTION("EF: you need to provide x and y for bicubic interpolation");
    return -1.0;
  }

  ///
  virtual Double EvaluatePrime(Double x) {
    EXCEPTION("EP: you need to provide x and y for bicubic interpolation");
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

  typedef enum { NONE, X, Y } Derivative;

  // approximate partial derivatives by finite differences
  void ApproxPartialDeriv(StdVector<double>& dFda, StdVector<double>& dFdb, StdVector<double>& dFdadb) const;

  // calculate coefficients of interpolating polynoms
  void CalcCoeff(StdVector<double>& coeff, StdVector<double> F, StdVector<double> Fda, StdVector<double> Fdb, StdVector<double> Fdadb) const;

  // evaluation of the interpolation polynomial at point x,y
  double EvaluatePolynom(unsigned int index, double x, double y, Derivative deriv = Derivative::NONE) const;

  unsigned int sub2ind(unsigned int ii, unsigned int jj) const;

  // sample points
  Vector<double> x_;
  Vector<double> y_;

  // sample values
  StdVector<double> data_;

  // coefficients of interpolating polynoms
  StdVector<StdVector<double>> coeff_;

  double offsetX;
  double offsetY;
  double scaleX;
  double scaleY;

  bool periodic_;
};

} //end of namespace


#endif
