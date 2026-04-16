#include "DataInOut/Logging/LogConfigurator.hh"
#include "Utils/CubicInterpolate.hh"
#include "Utils/CubicInterpolateCoeff.hh"

namespace CoupledField
{

DEFINE_LOG(cubicappx, "cubicappx")

CubicInterpolate::CubicInterpolate(const StdVector<double>& data, const StdVector<double>& a, const bool periodic) {
  unsigned int m = a.GetSize();
  numMeas_ = m;
  assert(m > 1);
  assert(data.GetSize() == numMeas_);
  data_ = data;
  factor_ = 1.;
  coeff_.Resize( (m-1) );
  periodic_ = periodic;

  if (periodic_)
    if (std::abs(data_[m-1] - data_[0]) > 1e-12)
      EXCEPTION("CubicInterpolate: data is not periodic. "<< data[0] << "  " << data[m-1]);

  StdVector<double> x = a;
  // sort x
  // be aware, that data and/or coeff are not reordered! this might lead to evaluation of wrong pieces
  std::sort(x.Begin(), x.End());
  offsetX_ = x[0];
  scaleX_ = 1./(x[m-1]-x[0]);
  // copy from StdVector to Vector (no cast available)
  // and scale to [0,1]
  x_.Resize(x.GetSize());
  for(unsigned int i = 0; i < m; ++i)
    x_[i] = (x[i] - offsetX_) * scaleX_;


  LOG_DBG2(cubicappx) << "CI: x = " << x.ToString();
  LOG_DBG(cubicappx) << "CI: offsetX= " << offsetX_ << " scaleX= " << scaleX_;
  LOG_DBG2(cubicappx) << "CI: x_= " << x_.ToString();
  LOG_DBG3(cubicappx) << "CI: data_= " << data_.ToString();
}

CubicInterpolate::CubicInterpolate(const StdVector<double>& data, const StdVector<double>& a,
    const Matrix<double>& coeff, const bool periodic)
  : CubicInterpolate(data, a, periodic) {

  assert(coeff.GetNumCols() == 4);
  // convert Matrix<double> to StdVector<StdVector<double>>
  coeff_.Resize(coeff.GetNumRows());
  for (unsigned int i=0; i < coeff.GetNumRows(); i++) {
    Vector<double> temp;
    coeff.GetRow(temp, i);
    coeff_[i].Resize(temp.GetSize());
    std::copy(temp.GetPointer(), temp.GetPointer()+temp.GetSize(), coeff_[i].Begin());
  }
}

double CubicInterpolate::EvaluateFunc(double x) const {
  double xloc, dxloc;
  unsigned int index = GetLocalValues(x, xloc, dxloc);

  double fValue = EvaluatePolynom(index, xloc) * factor_;

  LOG_DBG2(cubicappx) << "CI: Eval monocubic interpolator at point (" << x << "): " << fValue;
  return fValue;
}

Vector<double> CubicInterpolate::EvaluatePrime(const Vector<double>& p) const {
  Vector<double> fValue(1, 0.0);

  fValue[0] = EvaluateDeriv(p[0]);

  return fValue;
}

double CubicInterpolate::EvaluateDeriv(double x) const {
  double fValue = 0.0;

  double xloc, dxloc;
  unsigned int index = GetLocalValues(x, xloc, dxloc);

  fValue = EvaluatePolynom(index, xloc, Derivative::X) / dxloc * scaleX_ * factor_;

  LOG_DBG2(cubicappx) << "CI::ED df(" << x << ")= " << fValue;
  return fValue;
}

double CubicInterpolate::EvaluateSecondDeriv(double x) const {
  double fValue = 0.0;

  double xloc, dxloc;
  unsigned int index = GetLocalValues(x, xloc, dxloc);

  fValue = EvaluatePolynom(index, xloc, Derivative::XX) / dxloc * scaleX_ / dxloc * scaleX_ * factor_;

  LOG_DBG2(cubicappx) << "CI::ESD ddf(" << x << ")= " << fValue;
  return fValue;
}

double CubicInterpolate::EvaluatePolynom(unsigned int index, double x, Derivative deriv) const {
  double ret = 0.0;
  switch(deriv) {
  case Derivative::NONE:
    for(unsigned int i = 0; i < 4; ++i) {
      ret += coeff_[index][i] * pow(x,i);
    }
    break;
  case Derivative::X:
    for(unsigned int i = 1; i < 4; ++i) {
      ret += coeff_[index][i] * i * pow(x,i-1);
    }
    break;
  case Derivative::XX:
    for(unsigned int i = 2; i < 4; ++i) {
      ret += coeff_[index][i] * i * (i-1) * pow(x,i-2);
    }
    break;
  }
  return ret;
}

void CubicInterpolate::CalcApproximation(bool start) {
  unsigned int m = x_.GetSize();

  // Calculation of the coefficients of the interpolation polynomial for all possible intervals
  StdVector<double> dFdx;
  ApproxPartialDeriv(dFdx);

  double dx;
  StdVector<double> F(2);
  StdVector<double> Fdx(2);
  for(unsigned int i = 0; i < m-1; ++i) {
    F[0]= data_[i];
    F[1]= data_[i+1];

    Fdx[0]= dFdx[i];
    Fdx[1]= dFdx[i+1];

    // each patch evaluates on [0,1]^dim
    // thus we have to scale values
    dx = x_[i+1] - x_[i];
    for (int q = 0; q < 2; ++q) {
      F[q] *= 1.0;
      Fdx[q] *= dx;
    }

    StdVector<double> coeff(4);
    CalcCoeff(coeff, F, Fdx);
    coeff_[i] = coeff;
    LOG_DBG2(cubicappx) << "CA: coeff_[" << i << "]= " << coeff_[i].ToString();
  }
}

void CubicInterpolate::ApproxPartialDeriv(StdVector<double>& dFdx) const {
  // approximation of the derivatives with finite differences at the sample points
  int m = x_.GetSize();

  dFdx.Resize(numMeas_);

  unsigned int iu, il;
  double dx;

  for(int i = 0; i < m; ++i) {
    if(periodic_) {
      iu = i+1 < m ? i+1 : 1;
      il = i-1 >= 0 ? i-1 : m-2;
    } else {
      iu = i+1 < m ? i+1 : i;
      il = i-1 >= 0 ? i-1 : 0;
    }
    LOG_DBG3(cubicappx) << "APD: iu= " << iu << " il= " << il;

    dx = x_[iu] - x_[il];
    if (periodic_ && x_[iu] < x_[il])
      // we cross the boundary of (x[0],x[m-1])
      dx += 1;

    // the scaled derivative. for the real one, this has to be multiplied with scaleX_
    // (f(i+1)-f(i-1)) / (x(i+1)-x(i-1))
    dFdx[i] = (data_[iu] - data_[il]) / dx;

    LOG_DBG3(cubicappx) << "APD: idx= " << i << " dFdx= " << dFdx[i] * scaleX_;
  }
}


void CubicInterpolate::CalcCoeff(StdVector<double>& coeff, const StdVector<double>& F, const StdVector<double>& Fdx) const {
  Vector<double> x(4, 0.0);
  // Solve a linear system in order to get the interpolation coefficients
  for(unsigned int i = 0; i < 2; ++i) {
    x[i] = F[i];
    x[2+i] = Fdx[i];
  }
  LOG_DBG3(cubicappx) << "CC: x = " << x.ToString();

  // solve equation system
  for(unsigned int i = 0; i < 4; ++i) {
    coeff[i] = 0.;
    for (unsigned int j = 0; j < 4; ++j) {
      coeff[i] += CIC[i][j] * x[j];
    }
  }
}

inline unsigned int CubicInterpolate::GetLocalValues(double x, double& xloc, double& dxloc) const {
  // scale to [0,1] like in constructor
  x = (x - offsetX_) * scaleX_;
  // for periodic functions we map values outside of data range to the inside
  if (periodic_) {
    if (x > 1)
      x -= (int)x;
    if (x < 0)
      x += (int)x + 1;
  }

  assert(0 - std::numeric_limits<double>::epsilon() <= x && x <= 1 + std::numeric_limits<double>::epsilon());

  unsigned int kx;

  // if coordinate is out of bounds, return boundary value
  if( x < x_[0] ) {
    kx = 0;
  } else {
    // upper_bound gives pointer to first x_ with x < x_. we subtract 1 to get largest x_ < x.
    kx = std::upper_bound(x_.GetPointer(), x_.GetPointer() + x_.GetSize() - 1, x) - x_.GetPointer() - 1;
  }
  LOG_DBG3(cubicappx) << "kx: " << kx;

  unsigned int index = kx;

  // find relative x,y in patch
  dxloc = x_[kx + 1] - x_[kx];
  xloc = (x - x_[kx]) / dxloc;

  return index;
}

}
