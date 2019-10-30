#include "DataInOut/Logging/LogConfigurator.hh"
#include "DataInOut/Logging/log.hpp"
#include "Utils/TriCubicInterpolate.hh"
#include "Utils/TriCubicInterpolateCoeff.hh"

namespace CoupledField
{

DECLARE_LOG(tricubicappx)
DEFINE_LOG(tricubicappx, "tricubicappx")

TriCubicInterpolate::TriCubicInterpolate(StdVector<double> data, StdVector<double> a, StdVector<double> b, StdVector<double> c, bool periodic) {
  unsigned int m = a.GetSize();
  unsigned int n = b.GetSize();
  unsigned int o = c.GetSize();
  numMeas_ = m*n*o;
  assert(data.GetSize() == numMeas_);
  data_ = data;
  factor_ = 1.;
  coeff_.Resize( (m-1) * (n-1) * (o-1) );
  periodic_ = periodic;

  StdVector<double> x = a;
  StdVector<double> y = b;
  StdVector<double> z = c;
  // sort x, y and z
  std::sort(x.Begin(), x.End());
  std::sort(y.Begin(), y.End());
  std::sort(z.Begin(), z.End());
  offsetX_ = x[0];
  offsetY_ = y[0];
  offsetZ_ = z[0];
  scaleX_ = 1./x[m-1];
  scaleY_ = 1./y[n-1];
  scaleZ_ = 1./z[o-1];
  // copy from StdVector to Vector (no cast available)
  // and scale to [0,1]
  x_.Resize(x.GetSize());
  for(unsigned int i = 0; i < m; ++i)
    x_[i] = (x[i] - offsetX_) * scaleX_;
  y_.Resize(y.GetSize());
  for(unsigned int i = 0; i < n; ++i)
    y_[i] = (y[i] - offsetY_) * scaleY_;
  z_.Resize(z.GetSize());
  for(unsigned int i = 0; i < o; ++i)
    z_[i] = (z[i] - offsetZ_) * scaleZ_;

  LOG_DBG(tricubicappx) << "TCI: offsetX= " << offsetX_ << " scaleX= " << scaleX_
      << " offsetY= " << offsetY_ << " scaleY= " << scaleY_
      << " offsetZ= " << offsetZ_ << " scaleZ= " << scaleZ_;
  LOG_DBG2(tricubicappx) << "TCI: x_= " << x_.ToString();
  LOG_DBG2(tricubicappx) << "TCI: y_= " << y_.ToString();
  LOG_DBG2(tricubicappx) << "TCI: z_= " << z_.ToString();
  LOG_DBG3(tricubicappx) << "TCI: data_= " << data_.ToString();
}

Double TriCubicInterpolate::EvaluateFunc(double x, double y, double z) {
  double xloc, yloc, zloc;
  unsigned int index = GetLocalValues(x, y, z, xloc, yloc, zloc);

  Double fValue = EvaluatePolynom(index, xloc, yloc, zloc);
  LOG_DBG2(tricubicappx) << "TCI: Eval tricubic interpolator at point (" << x << ", " << y << ", " << z << "): " << fValue*factor_;

  return fValue*factor_;
}

Vector<Double> TriCubicInterpolate::EvaluatePrime(double x, double y, double z) {
  Vector<Double> fValue(3, 0.0);

  double xloc, yloc, zloc;
  unsigned int index = GetLocalValues(x, y, z, xloc, yloc, zloc);

  fValue[0] = EvaluatePolynom(index, xloc, yloc, zloc, Derivative::X);
  fValue[1] = EvaluatePolynom(index, xloc, yloc, zloc, Derivative::Y);
  fValue[2] = EvaluatePolynom(index, xloc, yloc, zloc, Derivative::Z);

  return fValue*factor_;
}

double TriCubicInterpolate::EvaluatePolynom(unsigned int index, double x, double y, double z, Derivative deriv) const {
  double ret = 0.0;
  switch(deriv) {
  case Derivative::NONE:
    for(unsigned int i = 0; i < 4; ++i) {
      for(unsigned int j = 0; j < 4; ++j) {
        for(unsigned int k = 0; k < 4; ++k) {
          ret += coeff_[index][i + 4*j + 16*k] * pow(x,i) * pow(y,j) * pow(z,k);
        }
      }
    }
    break;
  case Derivative::X:
    for(unsigned int i = 1; i < 4; ++i) {
      for(unsigned int j = 0; j < 4; ++j) {
        for(unsigned int k = 0; k < 4; ++k) {
          ret += coeff_[index][i + 4*j + 16*k] * i * pow(x,i-1) * pow(y,j) * pow(z,k);
        }
      }
    }
    break;
  case Derivative::Y:
    for(unsigned int i = 0; i < 4; ++i) {
      for(unsigned int j = 1; j < 4; ++j) {
        for(unsigned int k = 0; k < 4; ++k) {
          ret += coeff_[index][i + 4*j + 16*k] * pow(x,i) * j * pow(y,j-1) * pow(z,k);
        }
      }
    }
    break;
  case Derivative::Z:
    for(unsigned int i = 0; i < 4; ++i) {
      for(unsigned int j = 1; j < 4; ++j) {
        for(unsigned int k = 0; k < 4; ++k) {
          ret += coeff_[index][i + 4*j + 16*k] * pow(x,i) * pow(y,j) * k * pow(z,k-1);
        }
      }
    }
  }
  return ret;
}

void TriCubicInterpolate::CalcApproximation(bool start) {
  unsigned int m = x_.GetSize();
  unsigned int n = y_.GetSize();
  unsigned int o = z_.GetSize();
  double dx = x_[1] - x_[0];
  double dy = y_[1] - y_[0];
  double dz = z_[1] - z_[0];

  LOG_DBG(tricubicappx) << "CA: dx= " << dx << " dy= " << dy << " dz= " << dz;

  // Calculation of the coefficients of the interpolation polynomial for all possible intervals
  StdVector<double> dFdx;
  StdVector<double> dFdy;
  StdVector<double> dFdz;
  StdVector<double> dFdxdy;
  StdVector<double> dFdxdz;
  StdVector<double> dFdydz;
  StdVector<double> dFdxdydz;
  ApproxPartialDeriv(dFdx, dFdy, dFdz, dFdxdy, dFdxdz, dFdydz, dFdxdydz);

  StdVector<double> F(8);
  StdVector<double> Fdx(8);
  StdVector<double> Fdy(8);
  StdVector<double> Fdz(8);
  StdVector<double> Fdxdy(8);
  StdVector<double> Fdxdz(8);
  StdVector<double> Fdydz(8);
  StdVector<double> Fdxdydz(8);
  for(unsigned int i = 0; i < m-1; ++i) {
    for(unsigned int j = 0; j < n-1; ++j) {
      for(unsigned int k = 0; k < o-1; ++k) {
        F[0]= data_[sub2ind(i,j,k)];
        F[1]= data_[sub2ind(i+1,j,k)];
        F[2]= data_[sub2ind(i,j+1,k)];
        F[3]= data_[sub2ind(i+1,j+1,k)];
        F[4]= data_[sub2ind(i,j,k+1)];
        F[5]= data_[sub2ind(i+1,j,k+1)];
        F[6]= data_[sub2ind(i,j+1,k+1)];
        F[7]= data_[sub2ind(i+1,j+1,k+1)];

        Fdx[0]= dFdx[sub2ind(i,j,k)];
        Fdx[1]= dFdx[sub2ind(i+1,j,k)];
        Fdx[2]= dFdx[sub2ind(i,j+1,k)];
        Fdx[3]= dFdx[sub2ind(i+1,j+1,k)];
        Fdx[4]= dFdx[sub2ind(i,j,k+1)];
        Fdx[5]= dFdx[sub2ind(i+1,j,k+1)];
        Fdx[6]= dFdx[sub2ind(i,j+1,k+1)];
        Fdx[7]= dFdx[sub2ind(i+1,j+1,k+1)];

        Fdy[0]= dFdy[sub2ind(i,j,k)];
        Fdy[1]= dFdy[sub2ind(i+1,j,k)];
        Fdy[2]= dFdy[sub2ind(i,j+1,k)];
        Fdy[3]= dFdy[sub2ind(i+1,j+1,k)];
        Fdy[4]= dFdy[sub2ind(i,j,k+1)];
        Fdy[5]= dFdy[sub2ind(i+1,j,k+1)];
        Fdy[6]= dFdy[sub2ind(i,j+1,k+1)];
        Fdy[7]= dFdy[sub2ind(i+1,j+1,k+1)];

        Fdxdy[0]= dFdxdy[sub2ind(i,j,k)];
        Fdxdy[1]= dFdxdy[sub2ind(i+1,j,k)];
        Fdxdy[2]= dFdxdy[sub2ind(i,j+1,k)];
        Fdxdy[3]= dFdxdy[sub2ind(i+1,j+1,k)];
        Fdxdy[4]= dFdxdy[sub2ind(i,j,k+1)];
        Fdxdy[5]= dFdxdy[sub2ind(i+1,j,k+1)];
        Fdxdy[6]= dFdxdy[sub2ind(i,j+1,k+1)];
        Fdxdy[7]= dFdxdy[sub2ind(i+1,j+1,k+1)];

        // each patch evaluates on [0,1]^dim
        // thus we have to scale values
        for (int q = 0; q < 8; ++q) {
          F[q] *= 1.0;
          Fdx[q] *= dx;
          Fdy[q] *= dy;
          Fdxdy[q] *= dx * dy;
        }

        StdVector<double> coeff(64);
        CalcCoeff(coeff, F, Fdx, Fdy, Fdz, Fdxdy, Fdxdz, Fdydz, Fdxdydz);
        coeff_[(m-1)*(n-1)*i + (n-1)*j + k] = coeff;
      }
    }
  }
}

void TriCubicInterpolate::ApproxPartialDeriv(StdVector<double>& dFdx, StdVector<double>& dFdy, StdVector<double>& dFdz,
    StdVector<double>& dFdxdy, StdVector<double>& dFdxdz, StdVector<double>& dFdydz, StdVector<double>& dFdxdydz) const {
  // approximation of the derivatives with finite differences at the sample points
  int m = x_.GetSize();
  int n = y_.GetSize();
  int o = z_.GetSize();

  double dx = x_[1] - x_[0];
  double dy = y_[1] - y_[0];
  double dz = z_[1] - z_[0];

  dFdx.Resize(numMeas_);
  dFdy.Resize(numMeas_);
  dFdz.Resize(numMeas_);
  dFdxdy.Resize(numMeas_);
  dFdxdz.Resize(numMeas_);
  dFdydz.Resize(numMeas_);
  dFdxdydz.Resize(numMeas_);

  unsigned int iu, il, ju, jl, ku, kl;

  for(int i = 0; i < m; ++i) {
    for(int j = 0; j < n; ++j) {
      for(int k = 0; k < o; ++k) {
        if(periodic_) {
          iu = i+1 < m ? i+1 : 0;
          il = i-1 >= 0 ? i-1 : m-1;
          ju = j+1 < n ? j+1 : 0;
          jl = j-1 >= 0 ? j-1 : n-1;
          ku = k+1 < o ? k+1 : 0;
          kl = k-1 >= 0 ? k-1 : k-1;
        } else {
          iu = i+1 < m ? i+1 : i;
          il = i-1 >= 0 ? i-1 : 0;
          ju = j+1 < n ? j+1 : j;
          jl = j-1 >= 0 ? j-1 : 0;
          ku = k+1 < o ? k+1 : k;
          kl = k-1 >= 0 ? k-1 : 0;
        }
        LOG_DBG3(tricubicappx) << "APD: iu= " << iu << " il= " << il << " ju= " << ju << " jl= " << jl << " ku= " << ku << " kl= " << kl;
        //f(i+1,j,k)-f(i-1,j,k) / 2
        dFdx[sub2ind(i,j,k)] = (data_[sub2ind(iu,j,k)] - data_[sub2ind(il,j,k)]) / (2.*dx);
        dFdy[sub2ind(i,j,k)] = (data_[sub2ind(i,ju,k)] - data_[sub2ind(i,jl,k)]) / (2.*dy);
        dFdz[sub2ind(i,j,k)] = (data_[sub2ind(i,j,ku)] - data_[sub2ind(i,j,kl)]) / (2.*dy);
        //f(i+1,j+1,k)-f(i+1,j-1,k)-f(i-1,j+1,k)+f(i-1,j-1,k))/4
        dFdxdy[sub2ind(i,j,k)] = (data_[sub2ind(iu,ju,k)] - data_[sub2ind(iu,jl,k)]
                              - data_[sub2ind(il,ju,k)] + data_[sub2ind(il,jl,k)]) / (4. * dx * dy);
        dFdxdz[sub2ind(i,j,k)] = (data_[sub2ind(iu,j,ku)] - data_[sub2ind(iu,j,kl)]
                              - data_[sub2ind(il,j,ku)] + data_[sub2ind(il,j,kl)]) / (4. * dx * dy);
        dFdydz[sub2ind(i,j,k)] = (data_[sub2ind(i,ju,ku)] - data_[sub2ind(i,ju,kl)]
                              - data_[sub2ind(i,jl,ku)] + data_[sub2ind(i,jl,kl)]) / (4. * dx * dy);
        //(f(i+1,j+1,z+1) - f(i+1,j+1,k-1) - f(i+1,j-1,k+1) + f(i+1,j-1,k-1)
        //- f(i-1,j+1,k+1) + f(i-1,j+1,k-1) + f(i-1,j-1,k+1) - f(i-1,j-1,k-1))/8
        dFdxdydz[sub2ind(i,j,k)] = (data_[sub2ind(iu,ju,ku)] - data_[sub2ind(iu,ju,kl)] - data_[sub2ind(iu,jl,ku)] + data_[sub2ind(iu,jl,kl)]
                              - data_[sub2ind(il,ju,ku)] + data_[sub2ind(il,ju,kl)] + data_[sub2ind(il,jl,ku)] - data_[sub2ind(il,jl,kl)])
                              / (8. * dx * dy * dz);
      }
    }
  }
}


void TriCubicInterpolate::CalcCoeff(StdVector<double>& coeff, StdVector<double> F, StdVector<double> Fdx, StdVector<double> Fdy, StdVector<double> Fdz,
    StdVector<double> Fdxdy, StdVector<double> Fdxdz, StdVector<double> Fdydz, StdVector<double> Fdxdydz) const {
  Vector<double> x(64, 0.0);
  // Solve a linear system in order to get the interpolation coefficients
  for(unsigned int i = 0; i < 8; ++i) {
    x[i] = F[i];
    x[8+i] = Fdx[i];
    x[16+i] = Fdy[i];
    x[24+i] = Fdz[i];
    x[32+i] = Fdxdy[i];
    x[40+i] = Fdxdz[i];
    x[48+i] = Fdydz[i];
    x[56+i] = Fdxdydz[i];
  }

  // solve equation system
  for(unsigned int i = 0; i < 64; ++i) {
    coeff[i] = 0.;
    for (unsigned int j = 0; j < 64; ++j) {
      coeff[i] += TCIC[i][j] * x[j];
    }
  }
}

unsigned int TriCubicInterpolate::GetLocalValues(double x, double y, double z, double& xloc, double& yloc, double& zloc) {
  // scale to [0,1] like in constructor
  x = (x - offsetX_) * scaleX_;
  y = (y - offsetY_) * scaleY_;
  z = (y - offsetY_) * scaleZ_;

  // get index of last element
  const unsigned int kxend = x_.GetSize() - 1;
  const unsigned int kyend = y_.GetSize() - 1;
  const unsigned int kzend = z_.GetSize() - 1;

  unsigned int kxlo, kxhi, kylo, kyhi, kzlo, kzhi;
  Double diffx, diffy, diffz;

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
  if( z < z_[0] ) {
    kzlo = 0;
  } else if ( z > z_[kzend] ) {
    kzlo = kzend-1;
  } else {
    findBracketIndices(z, z_, kzlo, kzhi, diffz);
  }
  unsigned int index = (z_.GetSize()-1)*(y_.GetSize()-1)*kxlo + (z_.GetSize()-1)*kylo + kzlo;

  // find relative x,y in patch
  double dx = x_[1] - x_[0];
  double dy = y_[1] - y_[0];
  double dz = z_[1] - z_[0];
  xloc = (x - x_[kxlo]) / dx;
  yloc = (y - y_[kylo]) / dy;
  zloc = (z - z_[kzlo]) / dz;

  return index;
}

unsigned int TriCubicInterpolate::sub2ind(unsigned int ii, unsigned int jj, unsigned kk) const {
  return kk*y_.GetSize()*x_.GetSize() + jj*x_.GetSize() + ii;
}

}
