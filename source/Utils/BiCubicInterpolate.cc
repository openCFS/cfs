#include <algorithm>
#include "DataInOut/Logging/LogConfigurator.hh"
#include "Utils/BiCubicInterpolate.hh"
#include "Utils/BiCubicInterpolateCoeff.hh"

#include <fstream>

namespace CoupledField
{

DEFINE_LOG(bicubicappx, "bicubicappx")

BiCubicInterpolate::BiCubicInterpolate(const StdVector<double>& data, const StdVector<double>& a,
    const StdVector<double>& b, const bool periodic) {
  unsigned int m = a.GetSize();
  unsigned int n = b.GetSize();
  numMeas_ = m*n;
  assert(m > 1);
  assert(n > 1);
  assert(data.GetSize() == numMeas_);
  data_ = data;
  factor_ = 1.;
  coeff_.Resize( (m-1) * (n-1) );
  periodic_ = periodic;

  if (periodic_) {
    for(unsigned int i = 0; i < m-1; ++i)
      if (std::abs(data_[sub2ind(i,0)] - data_[sub2ind(i,n-1)]) > 1e-12)
        EXCEPTION("CubicInterpolate: data is not periodic.");
    for(unsigned int j = 0; j < n-1; ++j)
      if (std::abs(data_[sub2ind(0,j)] - data_[sub2ind(m-1,j)]) > 1e-12)
        EXCEPTION("CubicInterpolate: data is not periodic.");
  }

  StdVector<double> x = a;
  StdVector<double> y = b;
  // sort x and y
  // be aware, that data and/or coeff are not reordered! this might lead to evaluation of wrong pieces
  std::sort(x.Begin(), x.End());
  std::sort(y.Begin(), y.End());
  offsetX_ = x[0];
  offsetY_ = y[0];
  scaleX_ = 1./(x[m-1] - x[0]);
  scaleY_ = 1./(y[n-1] - y[0]);
  // copy from StdVector to Vector (no cast available)
  // and scale to [0,1]
  x_.Resize(x.GetSize());
  for(unsigned int i = 0; i < m; ++i)
    x_[i] = (x[i] - offsetX_) * scaleX_;
  y_.Resize(y.GetSize());
  for(unsigned int i = 0; i < n; ++i)
    y_[i] = (y[i] - offsetY_) * scaleY_;

  LOG_DBG2(bicubicappx) << "BCI: x = " << x.ToString();
  LOG_DBG2(bicubicappx) << "BCI: y = " << y.ToString();
  LOG_DBG(bicubicappx) << "BCI: offsetX= " << offsetX_ << " scaleX= " << scaleX_
      << " offsetY= " << offsetY_ << " scaleY= " << scaleY_;
  LOG_DBG2(bicubicappx) << "BCI: x_= " << x_.ToString();
  LOG_DBG2(bicubicappx) << "BCI: y_= " << y_.ToString();
  LOG_DBG3(bicubicappx) << "BCI: data_= " << data_.ToString();
}

BiCubicInterpolate::BiCubicInterpolate(const StdVector<double>& data, const StdVector<double>& a,
    const StdVector<double>& b, const Matrix<double>& coeff, const bool periodic)
  : BiCubicInterpolate(data, a, b, periodic) {

  assert(coeff.GetNumCols() == 16);
  // convert Matrix<double> to StdVector<StdVector<double>>
  coeff_.Resize(coeff.GetNumRows());
  for (unsigned int i=0; i < coeff.GetNumRows(); i++) {
    Vector<double> temp;
    coeff.GetRow(temp, i);
    coeff_[i].Resize(temp.GetSize());
    std::copy(temp.GetPointer(), temp.GetPointer()+temp.GetSize(), coeff_[i].Begin());
  }
}

double BiCubicInterpolate::EvaluateFunc(double x, double y) const {
  double xloc, yloc, dxloc, dyloc;
  unsigned int index = GetLocalValues(x, y, xloc, yloc, dxloc, dyloc);

  double fValue = EvaluatePolynom(index, xloc, yloc) * factor_;

  LOG_DBG2(bicubicappx) << "BCI: Eval bicubic interpolator at point (" << x << ", " << y << "): " << fValue;
  return fValue;
}

Vector<double> BiCubicInterpolate::EvaluatePrime(double x, double y) const {
  Vector<double> fValue(2, 0.0);

  double xloc, yloc, dxloc, dyloc;
  unsigned int index = GetLocalValues(x, y, xloc, yloc, dxloc, dyloc);

  fValue[0] = EvaluatePolynom(index, xloc, yloc, Derivative::X) / dxloc * scaleX_ * factor_;
  fValue[1] = EvaluatePolynom(index, xloc, yloc, Derivative::Y) / dyloc * scaleY_ * factor_;

  LOG_DBG2(bicubicappx) << "BCI::EP df(" << x << ", " << y << "): " << fValue.ToString();
  return fValue;
}

double BiCubicInterpolate::EvaluateDeriv(double x, double y, int dparam) const {
  double fValue = 0.0;

  double xloc, yloc, dxloc, dyloc;
  unsigned int index = GetLocalValues(x, y, xloc, yloc, dxloc, dyloc);

  switch (dparam) {
  case 0:
    fValue = EvaluatePolynom(index, xloc, yloc, Derivative::X) / dxloc * scaleX_ * factor_;
    break;
  case 1:
    fValue = EvaluatePolynom(index, xloc, yloc, Derivative::Y) / dyloc * scaleY_ * factor_;
    break;
  }
  LOG_DBG2(bicubicappx) << "BCI::ED deriv=" << dparam << ", df(" << x << ", " << y << ")= " << fValue;
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

  // Calculation of the coefficients of the interpolation polynomial for all possible intervals
  StdVector<double> dFdx;
  StdVector<double> dFdy;
  StdVector<double> dFdxdy;
  ApproxPartialDeriv(dFdx, dFdy, dFdxdy);

  // for debug, this will write the coefficients to a file
  // these can then be loaded in matlab and used with plot_interp.m to show the interpolating polynomials
  //std::ofstream coeff_file;
  //std::stringstream filename;
  //filename << "bicubicInt_coeff_" << static_cast<const void*>(this) << ".txt";
  //coeff_file.open(filename.str());

  double dx, dy;
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
      dx = x_[i+1] - x_[i];
      dy = y_[j+1] - y_[j];
      for (int q = 0; q < 4; ++q) {
        F[q] *= 1.0;
        Fdx[q] *= dx;
        Fdy[q] *= dy;
        Fdxdy[q] *= dx * dy;
      }

      StdVector<double> coeff(16);
      CalcCoeff(coeff, F, Fdx, Fdy, Fdxdy);
      coeff_[(n-1)*i + j] = coeff;
      LOG_DBG3(bicubicappx) << "CA: coeff_[" << (n-1)*i+j << "]= " << coeff_[(n-1)*i+j].ToString();

      //coeff_file << coeff.ToString(TS_PLAIN, ",") << std::endl;
    }
  }
  //coeff_file.close();
}

void BiCubicInterpolate::ApproxPartialDeriv(StdVector<double>& dFdx, StdVector<double>& dFdy, StdVector<double>& dFdxdy) const {
  // approximation of the derivatives with finite differences at the sample points
  int m = x_.GetSize();
  int n = y_.GetSize();

  dFdx.Resize(numMeas_);
  dFdy.Resize(numMeas_);
  dFdxdy.Resize(numMeas_);

  unsigned int iu, il, ju, jl;
  double dx, dy;

  for(int i = 0; i < m; ++i) {
    for(int j = 0; j < n; ++j) {
      if(periodic_) {
        iu = i+1 < m ? i+1 : 1;
        il = i-1 >= 0 ? i-1 : m-2;
        ju = j+1 < n ? j+1 : 1;
        jl = j-1 >= 0 ? j-1 : n-2;
      } else {
        iu = i+1 < m ? i+1 : i;
        il = i-1 >= 0 ? i-1 : 0;
        ju = j+1 < n ? j+1 : j;
        jl = j-1 >= 0 ? j-1 : 0;
      }
      LOG_DBG3(bicubicappx) << "APD: iu= " << iu << " il= " << il << " ju= " << ju << " jl= " << jl;

      dx = x_[iu] - x_[il];
      if (periodic_ && x_[iu] < x_[il])
        // we cross the boundary of (x[0],x[m-1])
        dx += 1;
      dy = y_[ju] - y_[jl];
      if (periodic_ && y_[ju] < y_[jl])
        // we cross the boundary of (y[0],y[n-1])
        dy += 1;

      // the scaled derivative. for the real one, this has to be multiplied with scaleX_ or scaleY_
      //f(i+1,j)-f(i-1,j) / 2
      dFdx[sub2ind(i,j)] = (data_[sub2ind(iu,j)] - data_[sub2ind(il,j)]) / dx;
      dFdy[sub2ind(i,j)] = (data_[sub2ind(i,ju)] - data_[sub2ind(i,jl)]) / dy;
      //f(i+1,j+1)-f(i+1,j-1)-f(i-1,j+1)+f(i-1,j-1))/4
      dFdxdy[sub2ind(i,j)] = (data_[sub2ind(iu,ju)] - data_[sub2ind(iu,jl)]
                            - data_[sub2ind(il,ju)] + data_[sub2ind(il,jl)]) / (dx * dy);

      LOG_DBG3(bicubicappx) << "APD: i= " << i << ", j= " << j << ", idx= " << sub2ind(i,j)
          << " dFdx= " << dFdx[sub2ind(i,j)] * scaleX_
          << " dFdy= " << dFdy[sub2ind(i,j)] * scaleY_
          << " dFdxdy= " << dFdxdy[sub2ind(i,j)] * scaleX_ * scaleY_;
    }
  }
}


void BiCubicInterpolate::CalcCoeff(StdVector<double>& coeff, const StdVector<double>& F, const StdVector<double>& Fdx, const StdVector<double>& Fdy, const StdVector<double>& Fdxdy) const {
  Vector<double> x(16, 0.0);
  // Solve a linear system in order to get the interpolation coefficients
  for(unsigned int i = 0; i < 4; ++i) {
    x[i] = F[i];
    x[4+i] = Fdx[i];
    x[8+i] = Fdy[i];
    x[12+i] = Fdxdy[i];
  }
  LOG_DBG3(bicubicappx) << "CC: x = " << x.ToString();

  // solve equation system
  for(unsigned int i = 0; i < 16; ++i) {
    coeff[i] = 0.;
    for (unsigned int j = 0; j < 16; ++j) {
      coeff[i] += BCIC[i][j] * x[j];
    }
  }
}

inline unsigned int BiCubicInterpolate::GetLocalValues(double x, double y, double& xloc, double& yloc, double& dxloc, double& dyloc) const {
  // scale to [0,1] like in constructor
  x = (x - offsetX_) * scaleX_;
  y = (y - offsetY_) * scaleY_;
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
  }

  assert(0 <= x && x <= 1);
  assert(0 <= y && y <= 1);

  unsigned int kx, ky;

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
  LOG_DBG3(bicubicappx) << "kx: " << kx << ", ky: " << ky;

  unsigned int index = (y_.GetSize()-1)*kx + ky;

  // find relative x,y in patch
  dxloc = x_[kx + 1] - x_[kx];
  dyloc = y_[ky + 1] - y_[ky];
  xloc = (x - x_[kx]) / dxloc;
  yloc = (y - y_[ky]) / dyloc;

  return index;
}

unsigned int BiCubicInterpolate::sub2ind(unsigned int ii, unsigned int jj) const {
  return jj*x_.GetSize() + ii;
}

}
