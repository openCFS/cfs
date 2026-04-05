#include <algorithm>
#include "DataInOut/Logging/LogConfigurator.hh"
#include "Utils/TriCubicInterpolate.hh"
#include "Utils/TriCubicInterpolateCoeff.hh"

namespace CoupledField
{

DEFINE_LOG(tricubicappx, "tricubicappx")

TriCubicInterpolate::TriCubicInterpolate(StdVector<double>& data, StdVector<double>& a, StdVector<double>& b,
    StdVector<double>& c, bool periodic) {
  unsigned int m = a.GetSize();
  unsigned int n = b.GetSize();
  unsigned int o = c.GetSize();
  numMeas_ = m*n*o;
  assert(m > 1);
  assert(n > 1);
  assert(o > 1);
  assert(data.GetSize() == numMeas_);
  data_ = data;
  factor_ = 1.;
  coeff_.Resize( (m-1) * (n-1) * (o-1) );
  periodic_ = periodic;

  if (periodic_) {
    for(unsigned int i = 0; i < m-1; ++i)
      for(unsigned int j = 0; j < n-1; ++j)
        if (std::abs(data_[sub2ind(i,j,0)] - data_[sub2ind(i,j,o-1)]) > 1e-12)
          EXCEPTION("CubicInterpolate: data is not periodic.");
    for(unsigned int i = 0; i < m-1; ++i)
      for(unsigned int k = 0; k < o-1; ++k)
        if (std::abs(data_[sub2ind(i,0,k)] - data_[sub2ind(i,n-1,k)]) > 1e-12)
          EXCEPTION("CubicInterpolate: data is not periodic.");
    for(unsigned int j = 0; j < n-1; ++j)
      for(unsigned int k = 0; k < o-1; ++k)
        if (std::abs(data_[sub2ind(0,j,k)] - data_[sub2ind(m-1,j,k)]) > 1e-12)
          EXCEPTION("CubicInterpolate: data is not periodic.");
  }

  StdVector<double> x = a;
  StdVector<double> y = b;
  StdVector<double> z = c;
  // sort x, y and z
  // be aware, that data and/or coeff are not reordered! this might lead to evaluation of wrong pieces
  std::sort(x.Begin(), x.End());
  std::sort(y.Begin(), y.End());
  std::sort(z.Begin(), z.End());
  offsetX_ = x[0];
  offsetY_ = y[0];
  offsetZ_ = z[0];
  scaleX_ = 1./(x[m-1] - x[0]);
  scaleY_ = 1./(y[n-1] - y[0]);
  scaleZ_ = 1./(z[o-1] - z[0]);
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

  LOG_DBG2(tricubicappx) << "TCI: x = " << x.ToString();
  LOG_DBG2(tricubicappx) << "TCI: y = " << y.ToString();
  LOG_DBG2(tricubicappx) << "TCI: z = " << z.ToString();
  LOG_DBG(tricubicappx) << "TCI: offsetX= " << offsetX_ << " scaleX= " << scaleX_
      << " offsetY= " << offsetY_ << " scaleY= " << scaleY_
      << " offsetZ= " << offsetZ_ << " scaleZ= " << scaleZ_;
  LOG_DBG2(tricubicappx) << "TCI: x_= " << x_.ToString();
  LOG_DBG2(tricubicappx) << "TCI: y_= " << y_.ToString();
  LOG_DBG2(tricubicappx) << "TCI: z_= " << z_.ToString();
  LOG_DBG3(tricubicappx) << "TCI: data_= " << data_.ToString();
}

TriCubicInterpolate::TriCubicInterpolate(StdVector<double>& data, StdVector<double>& a,
    StdVector<double>& b, StdVector<double>& c, Matrix<double>& coeff, bool periodic)
  : TriCubicInterpolate(data, a, b, c, periodic) {

  assert(coeff.GetNumCols() == 64);
  // convert Matrix<double> to StdVector<StdVector<double>>
  coeff_.Resize(coeff.GetNumRows());
  for (unsigned int i=0; i < coeff.GetNumRows(); i++) {
    Vector<double> temp;
    coeff.GetRow(temp, i);
    coeff_[i].Resize(temp.GetSize());
    std::copy(temp.GetPointer(), temp.GetPointer()+temp.GetSize(), coeff_[i].Begin());
  }
}

double TriCubicInterpolate::EvaluateFunc(double x, double y, double z) const {
  double xloc, yloc, zloc, dxloc, dyloc, dzloc;
  unsigned int index = GetLocalValues(x, y, z, xloc, yloc, zloc, dxloc, dyloc, dzloc);

  double fValue = EvaluatePolynom(index, xloc, yloc, zloc) * factor_;

  LOG_DBG2(tricubicappx) << "TCI: Eval tricubic interpolator at point ("
      << x << ", " << y << ", " << z << ") with index " << index << ": " << fValue;
  return fValue;
}

Vector<double> TriCubicInterpolate::EvaluatePrime(double x, double y, double z) const {
  Vector<double> fValue(3, 0.0);

  double xloc, yloc, zloc, dxloc, dyloc, dzloc;
  unsigned int index = GetLocalValues(x, y, z, xloc, yloc, zloc, dxloc, dyloc, dzloc);

  fValue[0] = EvaluatePolynom(index, xloc, yloc, zloc, Derivative::X) / dxloc * scaleX_ * factor_;
  fValue[1] = EvaluatePolynom(index, xloc, yloc, zloc, Derivative::Y) / dyloc * scaleY_ * factor_;
  fValue[2] = EvaluatePolynom(index, xloc, yloc, zloc, Derivative::Z) / dzloc * scaleZ_ * factor_;

  LOG_DBG2(tricubicappx) << "TCI::EP df(" << x << ", " << y << ", " << z << "): " << fValue.ToString();
  return fValue;
}

double TriCubicInterpolate::EvaluateDeriv(double x, double y, double z, int dparam) const {
  double fValue = 0.0;

  double xloc, yloc, zloc, dxloc, dyloc, dzloc;
  unsigned int index = GetLocalValues(x, y, z, xloc, yloc, zloc, dxloc, dyloc, dzloc);

  switch (dparam) {
  case 0:
    fValue = EvaluatePolynom(index, xloc, yloc, zloc, Derivative::X) / dxloc * scaleX_ * factor_;
    break;
  case 1:
    fValue = EvaluatePolynom(index, xloc, yloc, zloc, Derivative::Y) / dyloc * scaleY_ * factor_;
    break;
  case 2:
    fValue = EvaluatePolynom(index, xloc, yloc, zloc, Derivative::Z) / dzloc * scaleZ_ * factor_;
    break;
  }
  LOG_DBG2(tricubicappx) << "TCI::ED deriv=" << dparam << ", df(" << x << ", " << y << ", " << z << ")= " << fValue;
  return fValue;
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
      for(unsigned int j = 0; j < 4; ++j) {
        for(unsigned int k = 1; k < 4; ++k) {
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

  // Calculation of the coefficients of the interpolation polynomial for all possible intervals
  StdVector<double> dFdx;
  StdVector<double> dFdy;
  StdVector<double> dFdz;
  StdVector<double> dFdxdy;
  StdVector<double> dFdxdz;
  StdVector<double> dFdydz;
  StdVector<double> dFdxdydz;
  ApproxPartialDeriv(dFdx, dFdy, dFdz, dFdxdy, dFdxdz, dFdydz, dFdxdydz);

  double dx, dy, dz;
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

        Fdz[0]= dFdz[sub2ind(i,j,k)];
        Fdz[1]= dFdz[sub2ind(i+1,j,k)];
        Fdz[2]= dFdz[sub2ind(i,j+1,k)];
        Fdz[3]= dFdz[sub2ind(i+1,j+1,k)];
        Fdz[4]= dFdz[sub2ind(i,j,k+1)];
        Fdz[5]= dFdz[sub2ind(i+1,j,k+1)];
        Fdz[6]= dFdz[sub2ind(i,j+1,k+1)];
        Fdz[7]= dFdz[sub2ind(i+1,j+1,k+1)];

        Fdxdy[0]= dFdxdy[sub2ind(i,j,k)];
        Fdxdy[1]= dFdxdy[sub2ind(i+1,j,k)];
        Fdxdy[2]= dFdxdy[sub2ind(i,j+1,k)];
        Fdxdy[3]= dFdxdy[sub2ind(i+1,j+1,k)];
        Fdxdy[4]= dFdxdy[sub2ind(i,j,k+1)];
        Fdxdy[5]= dFdxdy[sub2ind(i+1,j,k+1)];
        Fdxdy[6]= dFdxdy[sub2ind(i,j+1,k+1)];
        Fdxdy[7]= dFdxdy[sub2ind(i+1,j+1,k+1)];

        Fdxdz[0]= dFdxdz[sub2ind(i,j,k)];
        Fdxdz[1]= dFdxdz[sub2ind(i+1,j,k)];
        Fdxdz[2]= dFdxdz[sub2ind(i,j+1,k)];
        Fdxdz[3]= dFdxdz[sub2ind(i+1,j+1,k)];
        Fdxdz[4]= dFdxdz[sub2ind(i,j,k+1)];
        Fdxdz[5]= dFdxdz[sub2ind(i+1,j,k+1)];
        Fdxdz[6]= dFdxdz[sub2ind(i,j+1,k+1)];
        Fdxdz[7]= dFdxdz[sub2ind(i+1,j+1,k+1)];

        Fdydz[0]= dFdydz[sub2ind(i,j,k)];
        Fdydz[1]= dFdydz[sub2ind(i+1,j,k)];
        Fdydz[2]= dFdydz[sub2ind(i,j+1,k)];
        Fdydz[3]= dFdydz[sub2ind(i+1,j+1,k)];
        Fdydz[4]= dFdydz[sub2ind(i,j,k+1)];
        Fdydz[5]= dFdydz[sub2ind(i+1,j,k+1)];
        Fdydz[6]= dFdydz[sub2ind(i,j+1,k+1)];
        Fdydz[7]= dFdydz[sub2ind(i+1,j+1,k+1)];

        Fdxdydz[0]= dFdxdydz[sub2ind(i,j,k)];
        Fdxdydz[1]= dFdxdydz[sub2ind(i+1,j,k)];
        Fdxdydz[2]= dFdxdydz[sub2ind(i,j+1,k)];
        Fdxdydz[3]= dFdxdydz[sub2ind(i+1,j+1,k)];
        Fdxdydz[4]= dFdxdydz[sub2ind(i,j,k+1)];
        Fdxdydz[5]= dFdxdydz[sub2ind(i+1,j,k+1)];
        Fdxdydz[6]= dFdxdydz[sub2ind(i,j+1,k+1)];
        Fdxdydz[7]= dFdxdydz[sub2ind(i+1,j+1,k+1)];


        // each patch evaluates on [0,1]^dim
        // thus we have to scale values
        dx = x_[i+1] - x_[i];
        dy = y_[j+1] - y_[j];
        dz = z_[k+1] - z_[k];
        for (int q = 0; q < 8; ++q) {
          F[q] *= 1.0;
          Fdx[q] *= dx;
          Fdy[q] *= dy;
          Fdz[q] *= dz;
          Fdxdy[q] *= dx * dy;
          Fdxdz[q] *= dx * dz;
          Fdydz[q] *= dy * dz;
          Fdxdydz[q] *= dx * dy * dz;
        }

        StdVector<double> coeff(64);
        CalcCoeff(coeff, F, Fdx, Fdy, Fdz, Fdxdy, Fdxdz, Fdydz, Fdxdydz);
        coeff_[(o-1)*(n-1)*i + (o-1)*j + k] = coeff;
//        LOG_DBG3(tricubicappx) << "CA: coeff_[" << (o-1)*(n-1)*i + (o-1)*j + k << "] = " << coeff_[(o-1)*(n-1)*i+(o-1)*j+k].ToString();
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

  dFdx.Resize(numMeas_);
  dFdy.Resize(numMeas_);
  dFdz.Resize(numMeas_);
  dFdxdy.Resize(numMeas_);
  dFdxdz.Resize(numMeas_);
  dFdydz.Resize(numMeas_);
  dFdxdydz.Resize(numMeas_);

  unsigned int iu, il, ju, jl, ku, kl;
  double dx, dy, dz;

  for(int i = 0; i < m; ++i) {
    for(int j = 0; j < n; ++j) {
      for(int k = 0; k < o; ++k) {
        if(periodic_) {
          iu = i+1 < m ? i+1 : 1;
          il = i-1 >= 0 ? i-1 : m-2;
          ju = j+1 < n ? j+1 : 1;
          jl = j-1 >= 0 ? j-1 : n-2;
          ku = k+1 < o ? k+1 : 1;
          kl = k-1 >= 0 ? k-1 : o-2;
        } else {
          iu = i+1 < m ? i+1 : i;
          il = i-1 >= 0 ? i-1 : 0;
          ju = j+1 < n ? j+1 : j;
          jl = j-1 >= 0 ? j-1 : 0;
          ku = k+1 < o ? k+1 : k;
          kl = k-1 >= 0 ? k-1 : 0;
        }
        LOG_DBG3(tricubicappx) << "APD: iu= " << iu << " il= " << il << " ju= " << ju << " jl= " << jl << " ku= " << ku << " kl= " << kl;

        dx = x_[iu] - x_[il];
        if (periodic_ && x_[iu] < x_[il])
          // we cross the boundary of (x[0],x[m-1])
          dx += 1;
        dy = y_[ju] - y_[jl];
        if (periodic_ && y_[ju] < y_[jl])
          // we cross the boundary of (y[0],y[n-1])
          dy += 1;
        dz = z_[ku] - z_[kl];
        if (periodic_ && z_[ku] < z_[kl])
          // we cross the boundary of (z[0],z[o-1])
          dz += 1;

        // the scaled derivative. for the real one, this has to be multiplied with scaleX_, scaleY_ or scaleY_
        //f(i+1,j,k)-f(i-1,j,k) / 2
        dFdx[sub2ind(i,j,k)] = (data_[sub2ind(iu,j,k)] - data_[sub2ind(il,j,k)]) / dx;
        dFdy[sub2ind(i,j,k)] = (data_[sub2ind(i,ju,k)] - data_[sub2ind(i,jl,k)]) / dy;
        dFdz[sub2ind(i,j,k)] = (data_[sub2ind(i,j,ku)] - data_[sub2ind(i,j,kl)]) / dz;
        //f(i+1,j+1,k)-f(i+1,j-1,k)-f(i-1,j+1,k)+f(i-1,j-1,k))/4
        dFdxdy[sub2ind(i,j,k)] = (data_[sub2ind(iu,ju,k)] - data_[sub2ind(iu,jl,k)]
                              - data_[sub2ind(il,ju,k)] + data_[sub2ind(il,jl,k)]) / (dx * dy);
        dFdxdz[sub2ind(i,j,k)] = (data_[sub2ind(iu,j,ku)] - data_[sub2ind(iu,j,kl)]
                              - data_[sub2ind(il,j,ku)] + data_[sub2ind(il,j,kl)]) / (dx * dz);
        dFdydz[sub2ind(i,j,k)] = (data_[sub2ind(i,ju,ku)] - data_[sub2ind(i,ju,kl)]
                              - data_[sub2ind(i,jl,ku)] + data_[sub2ind(i,jl,kl)]) / (dy * dz);
        //(f(i+1,j+1,z+1) - f(i+1,j+1,k-1) - f(i+1,j-1,k+1) + f(i+1,j-1,k-1)
        //- f(i-1,j+1,k+1) + f(i-1,j+1,k-1) + f(i-1,j-1,k+1) - f(i-1,j-1,k-1))/8
        dFdxdydz[sub2ind(i,j,k)] = (data_[sub2ind(iu,ju,ku)] - data_[sub2ind(iu,ju,kl)] - data_[sub2ind(iu,jl,ku)] + data_[sub2ind(iu,jl,kl)]
                              - data_[sub2ind(il,ju,ku)] + data_[sub2ind(il,ju,kl)] + data_[sub2ind(il,jl,ku)] - data_[sub2ind(il,jl,kl)])
                              / (dx * dy * dz);

        LOG_DBG3(tricubicappx) << "APD: i= " << i << ", j= " << j << ", k= " << k << ", idx= " << sub2ind(i,j,k)
            << " dFdx= " << dFdx[sub2ind(i,j,k)] * scaleX_
            << " dFdy= " << dFdy[sub2ind(i,j,k)] * scaleY_
            << " dFdz= " << dFdz[sub2ind(i,j,k)] * scaleZ_
            << " dFdxdy= " << dFdxdy[sub2ind(i,j,k)] * scaleX_ * scaleY_
            << " dFdxdz= " << dFdxdz[sub2ind(i,j,k)] * scaleX_ * scaleY_
            << " dFdydz= " << dFdydz[sub2ind(i,j,k)] * scaleX_ * scaleY_
            << " dFdxdydz= " << dFdxdydz[sub2ind(i,j,k)] * scaleX_ * scaleY_ * scaleZ_;
      }
    }
  }
}


void TriCubicInterpolate::CalcCoeff(StdVector<double>& coeff, const StdVector<double>& F,
    const StdVector<double>& Fdx, const StdVector<double>& Fdy, const StdVector<double>& Fdz,
    const StdVector<double>& Fdxdy, const StdVector<double>& Fdxdz, const StdVector<double>& Fdydz,
    const StdVector<double>& Fdxdydz) const {
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
  LOG_DBG3(tricubicappx) << "CC: x = " << x.ToString();

  // solve equation system
  for(unsigned int i = 0; i < 64; ++i) {
    coeff[i] = 0.;
    for (unsigned int j = 0; j < 64; ++j) {
      coeff[i] += TCIC[i][j] * x[j];
    }
  }
}

inline unsigned int TriCubicInterpolate::GetLocalValues(double x, double y, double z,
    double& xloc, double& yloc, double& zloc, double& dxloc, double& dyloc, double& dzloc) const {
  // scale to [0,1] like in constructor
  x = (x - offsetX_) * scaleX_;
  y = (y - offsetY_) * scaleY_;
  z = (z - offsetZ_) * scaleZ_;
  // for periodic functions we map values outside of data range to the inside
  if (periodic_) {
    if (x > 1)
      x -= (int)x;
    if (x < 0)
      x += (int)x + 1;
    if (y > 1)
      y -= (int)y;
    if (y < 0)
      y += (int)y + 1;
    if (z > 1)
      z -= (int)z;
    if (z < 0)
      z += (int)z + 1;
  }
  assert(0 <= x && x <= 1);
  assert(0 <= y && y <= 1);
  assert(0 <= z && z <= 1);

  unsigned int kx, ky, kz;

  // if coordinate is out of bounds, return boundary value
  if( x < x_[0] ) {
    kx = 0;
  } else {
    // upper_bound gives pointer to first x_ with x < x_. we subtract 1 to get largest x_ < x.
    kx = std::upper_bound(x_.GetPointer(), x_.GetPointer() + x_.GetSize() - 1, x) - x_.GetPointer() - 1;
  }
  if( y < y_[0] ) {
    ky = 0;
  } else {
    ky = std::upper_bound(y_.GetPointer(), y_.GetPointer() + y_.GetSize() - 1, y) - y_.GetPointer() - 1;
  }
  if( z < z_[0] ) {
    kz = 0;
  } else {
    kz = std::upper_bound(z_.GetPointer(), z_.GetPointer() + z_.GetSize() - 1, z) - z_.GetPointer() - 1;
  }
  LOG_DBG3(tricubicappx) << "kx: " << kx << ", ky: " << ky << ", kz: " << kz;

  unsigned int index = (z_.GetSize()-1)*(y_.GetSize()-1)*kx + (z_.GetSize()-1)*ky + kz;

  // find relative x,y in patch
  dxloc = x_[kx + 1] - x_[kx];
  dyloc = y_[ky + 1] - y_[ky];
  dzloc = z_[kz + 1] - z_[kz];
  xloc = (x - x_[kx]) / dxloc;
  yloc = (y - y_[ky]) / dyloc;
  zloc = (z - z_[kz]) / dzloc;

  return index;
}

unsigned int TriCubicInterpolate::sub2ind(unsigned int ii, unsigned int jj, unsigned kk) const {
  return kk*y_.GetSize()*x_.GetSize() + jj*x_.GetSize() + ii;
}

}
