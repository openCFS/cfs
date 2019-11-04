#include "DataInOut/Logging/LogConfigurator.hh"
#include "DataInOut/Logging/log.hpp"
#include "Utils/BiCubicInterpolate.hh"
#include "Utils/BiCubicInterpolateCoeff.hh"

namespace CoupledField
{

DECLARE_LOG(bicubicappx)
DEFINE_LOG(bicubicappx, "bicubicappx")

BiCubicInterpolate::BiCubicInterpolate(StdVector<double> data, StdVector<double> a, StdVector<double> b, bool periodic) {
  unsigned int m = a.GetSize();
  unsigned int n = b.GetSize();
  numMeas_ = m*n;
  assert(data.GetSize() == numMeas_);
  data_ = data;
  factor_ = 1.;
  coeff_.Resize( (m-1) * (n-1) );
  periodic_ = periodic;

  StdVector<double> x = a;
  StdVector<double> y = b;
  // sort x and y
  std::sort(x.Begin(), x.End());
  std::sort(y.Begin(), y.End());
  offsetX_ = x[0];
  offsetY_ = y[0];
  scaleX_ = 1./x[m-1];
  scaleY_ = 1./y[n-1];
  // copy from StdVector to Vector (no cast available)
  // and scale to [0,1]
  x_.Resize(x.GetSize());
  for(unsigned int i = 0; i < m; ++i)
    x_[i] = (x[i] - offsetX_) * scaleX_;
  y_.Resize(y.GetSize());
  for(unsigned int i = 0; i < n; ++i)
    y_[i] = (y[i] - offsetY_) * scaleY_;

  LOG_DBG(bicubicappx) << "BCI: offsetX= " << offsetX_ << " scaleX= " << scaleX_
      << " offsetY= " << offsetY_ << " scaleY= " << scaleY_;
  LOG_DBG2(bicubicappx) << "BCI: x_= " << x_.ToString();
  LOG_DBG2(bicubicappx) << "BCI: y_= " << y_.ToString();
  LOG_DBG3(bicubicappx) << "BCI: data_= " << data_.ToString();
}

Double BiCubicInterpolate::EvaluateFunc(double x, double y) {
  double xloc, yloc;
  unsigned int index = GetLocalValues(x, y, xloc, yloc);

  Double fValue = EvaluatePolynom(index, xloc, yloc) * factor_;
  LOG_DBG2(bicubicappx) << "BCI: Eval bicubic interpolator at point (" << x << ", " << y << "): " << fValue;

  return fValue;
}

Vector<Double> BiCubicInterpolate::EvaluatePrime(double x, double y) {
  Vector<Double> fValue(2, 0.0);

  double xloc, yloc;
  unsigned int index = GetLocalValues(x, y, xloc, yloc);

  fValue[0] = EvaluatePolynom(index, xloc, yloc, Derivative::X) * factor_;
  fValue[1] = EvaluatePolynom(index, xloc, yloc, Derivative::Y) * factor_;

  return fValue;
}

double BiCubicInterpolate::EvaluatePolynom(unsigned int index, double x, double y, Derivative deriv) const {
  double ret = 0.0;
  switch(deriv) {
  case Derivative::NONE:
    for(unsigned int i = 0; i < 4; ++i) {
      for(unsigned int j = 0; j < 4; ++j) {
        ret += coeff_[index][i + 4*j] * pow(x,i) * pow(y,j);
      }
    }
    break;
  case Derivative::X:
    for(unsigned int i = 1; i < 4; ++i) {
      for(unsigned int j = 0; j < 4; ++j) {
        ret += coeff_[index][i + 4*j] * i * pow(x,i-1) * pow(y,j);
      }
    }
    break;
  case Derivative::Y:
    for(unsigned int i = 0; i < 4; ++i) {
      for(unsigned int j = 1; j < 4; ++j) {
        ret += coeff_[index][i + 4*j] * pow(x,i) * j * pow(y,j-1);
      }
    }
  }
  return ret;
}

void BiCubicInterpolate::CalcApproximation(bool start) {
  unsigned int m = x_.GetSize();
  unsigned int n = y_.GetSize();
  double dx = x_[1] - x_[0];
  double dy = y_[1] - y_[0];

  LOG_DBG(bicubicappx) << "CA: dx= " << dx << " dy= " << dy;

  // Calculation of the coefficients of the interpolation polynomial for all possible intervals
  StdVector<double> dFdx;
  StdVector<double> dFdy;
  StdVector<double> dFdxdy;
  ApproxPartialDeriv(dFdx, dFdy, dFdxdy);

  StdVector<double> F(4);
  StdVector<double> Fdx(4);
  StdVector<double> Fdy(4);
  StdVector<double> Fdxdy(4);
  for(unsigned int i = 0; i < m-1; ++i) {
    for(unsigned int j = 0; j < n-1; ++j) {
      F[0]= data_[sub2ind(i,j)];
      F[1]= data_[sub2ind(i+1,j)];
      F[2]= data_[sub2ind(i,j+1)];
      F[3]= data_[sub2ind(i+1,j+1)];

      Fdx[0]= dFdx[sub2ind(i,j)];
      Fdx[1]= dFdx[sub2ind(i+1,j)];
      Fdx[2]= dFdx[sub2ind(i,j+1)];
      Fdx[3]= dFdx[sub2ind(i+1,j+1)];

      Fdy[0]= dFdy[sub2ind(i,j)];
      Fdy[1]= dFdy[sub2ind(i+1,j)];
      Fdy[2]= dFdy[sub2ind(i,j+1)];
      Fdy[3]= dFdy[sub2ind(i+1,j+1)];

      Fdxdy[0]= dFdxdy[sub2ind(i,j)];
      Fdxdy[1]= dFdxdy[sub2ind(i+1,j)];
      Fdxdy[2]= dFdxdy[sub2ind(i,j+1)];
      Fdxdy[3]= dFdxdy[sub2ind(i+1,j+1)];

      // each patch evaluates on [0,1]^dim
      // thus we have to scale values
      for (int q = 0; q < 4; ++q) {
        F[q] *= 1.0;
        Fdx[q] *= dx;
        Fdy[q] *= dy;
        Fdxdy[q] *= dx * dy;
      }

      StdVector<double> coeff(16);
      CalcCoeff(coeff, F, Fdx, Fdy, Fdxdy);
      coeff_[(n-1)*i + j] = coeff;
      LOG_DBG2(bicubicappx) << "CA: coeff_[" << (n-1)*i+j << "]= " << coeff_[(n-1)*i+j].ToString();
    }
  }
}

void BiCubicInterpolate::ApproxPartialDeriv(StdVector<double>& dFdx, StdVector<double>& dFdy, StdVector<double>& dFdxdy) const {
  // approximation of the derivatives with finite differences at the sample points
  int m = x_.GetSize();
  int n = y_.GetSize();

  double dx = x_[1] - x_[0];
  double dy = y_[1] - y_[0];

  dFdx.Resize(numMeas_);
  dFdy.Resize(numMeas_);
  dFdxdy.Resize(numMeas_);

  unsigned int iu, il, ju, jl;

  for(int i = 0; i < m; ++i) {
    for(int j = 0; j < n; ++j) {
      if(periodic_) {
        iu = i+1 < m ? i+1 : 0;
        il = i-1 >= 0 ? i-1 : m-1;
        ju = j+1 < n ? j+1 : 0;
        jl = j-1 >= 0 ? j-1 : n-1;
      } else {
        iu = i+1 < m ? i+1 : i;
        il = i-1 >= 0 ? i-1 : 0;
        ju = j+1 < n ? j+1 : j;
        jl = j-1 >= 0 ? j-1 : 0;
      }
      LOG_DBG3(bicubicappx) << "APD: iu= " << iu << " il= " << il << " ju= " << ju << " jl= " << jl;
      //f(i+1,j)-f(i-1,j) / 2
      dFdx[sub2ind(i,j)] = (data_[sub2ind(iu,j)] - data_[sub2ind(il,j)]) / (2.*dx);
      dFdy[sub2ind(i,j)] = (data_[sub2ind(i,ju)] - data_[sub2ind(i,jl)]) / (2.*dy);
      //f(i+1,j+1)-f(i+1,j-1)-f(i-1,j+1)+f(i-1,j-1))/4
      dFdxdy[sub2ind(i,j)] = (data_[sub2ind(iu,ju)] - data_[sub2ind(iu,jl)]
                            - data_[sub2ind(il,ju)] + data_[sub2ind(il,jl)]) / (4. * dx * dy);

      LOG_DBG3(bicubicappx) << "APD: idx= " << sub2ind(i,j) << " dFdx= " << dFdx[sub2ind(i,j)]
        << " dFdy= " << dFdy[sub2ind(i,j)] << " dFdxdy= " << dFdxdy[sub2ind(i,j)];
    }
  }
}


void BiCubicInterpolate::CalcCoeff(StdVector<double>& coeff, StdVector<double> F, StdVector<double> Fdx, StdVector<double> Fdy, StdVector<double> Fdxdy) const {
  Vector<double> x(16, 0.0);
  // Solve a linear system in order to get the interpolation coefficients
  for(unsigned int i = 0; i < 4; ++i) {
    x[i] = F[i];
    x[4+i] = Fdx[i];
    x[8+i] = Fdy[i];
    x[12+i] = Fdxdy[i];
  }
  LOG_DBG(bicubicappx) << "CC: x = " << x.ToString();

  // solve equation system
  for(unsigned int i = 0; i < 16; ++i) {
    coeff[i] = 0.;
    for (unsigned int j = 0; j < 16; ++j) {
      coeff[i] += BCIC[i][j] * x[j];
    }
  }
}

unsigned int BiCubicInterpolate::GetLocalValues(double x, double y, double& xloc, double& yloc) {
  // scale to [0,1] like in constructor
  x = (x - offsetX_) * scaleX_;
  y = (y - offsetY_) * scaleY_;

  // get index of last element
  const unsigned int kxend = x_.GetSize() - 1;
  const unsigned int kyend = y_.GetSize() - 1;

  unsigned int kxlo, kxhi, kylo, kyhi;
  Double diffx, diffy;

  // if coordinate is out of bounds, return boundary value
  if( x < x_[0] ) {
    kxlo = 0;
  } else if( x > x_[kxend] ) {
    kxlo = kxend-1;
  } else {
    findBracketIndices(x, x_, kxlo, kxhi, diffx);
  }
  if( y < y_[0] ) {
    kylo = 0;
  } else if ( y > y_[kyend] ) {
    kylo = kyend-1;
  } else {
    findBracketIndices(y, y_, kylo, kyhi, diffy);
  }
  unsigned int index = (y_.GetSize()-1)*kxlo + kylo;

  // find relative x,y in patch
  double dx = x_[1] - x_[0];
  double dy = y_[1] - y_[0];
  xloc = (x - x_[kxlo]) / dx;
  yloc = (y - y_[kylo]) / dy;

  return index;
}

unsigned int BiCubicInterpolate::sub2ind(unsigned int ii, unsigned int jj) const {
  return jj*x_.GetSize() + ii;
}

}
