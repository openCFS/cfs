// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <fstream>
#include <iostream>
#include <math.h>
#include <algorithm>
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string/classification.hpp>

#include "tools.hh"
#include "MatVec/matrix.hh"
#include "Elements/elements_header.hh"
#include "Domain/elem.hh"
#include "Domain/grid.hh"
#include "DataInOut/WriteInfo.hh"
#include "DataInOut/Logging/cfslog.hh"
#include "General/exception.hh"


DECLARE_LOG(tools)
DEFINE_LOG(tools, "tools")

namespace CoupledField {

  // =========================================
  //   Split a string into a list of strings
  // =========================================
  void SplitStringList( const std::string &list, StdVector<std::string> &strVec,
                        const char delimiter ) {

    UInt lastDelim = 0;
    strVec.Clear();
    UInt i=0;

    // ignore all leading spaces
    while ( i < list.length() && list[i] == ' ' ) {
      i++;
      lastDelim++;
    }

    // get the n-1 entries of the list
    for ( i = 0; i < list.length(); i++ ) {
      if ( list[i] == delimiter ) {
        strVec.Push_back( std::string( list, lastDelim, i-lastDelim ));
        i++;

        // ignore next spaces
        while ( i < list.length() && list[i] == ' ' ) {
          i++;
        }

        lastDelim = i;
      }
    }

    // get the n-th entry
    i = lastDelim;
    while ( i < list.length() && list[i] != ' ' ) {
      i++;
    }

    strVec.Push_back( std::string( list, lastDelim, i-lastDelim ));

  }
  
  Double dist_Mat(const Matrix<Double> &a) {
    Double preSqrt=0;
    UInt i;
    const UInt k=a.GetNumRows();
    // std::cout<<"tools.cc:size of matrix: "<<k<<std::endl;
    for (i=0; i<k; i++)
      preSqrt+= (a[i][0]-a[i][1]) * (a[i][0]-a[i][1]);
    return sqrt(preSqrt);
  }

  /// a-->b
  void calcNormal2Line(Vector<Double> & normal, const Point &a, const Point &b) {
    Double distance=Point::Dist(a,b);
    normal[0]=(b[1]-a[1])/distance;
    normal[1]=(a[0]-b[0])/distance;
  }

  /// a-->b
  void calcNormal2Line_Mat(Vector<Double> & normal, const Matrix<double> &a)
  {
    const double distance=dist_Mat(a);
    // std::cout<<"distance :"<<distance<<std::endl;

    normal[0]=(a[1][1]-a[1][0])/distance;
    normal[1]=(a[0][1]-a[0][0])/distance;
  }


  /// a-->b-->c. no fix orientation.
  void calcNormal2Surface(Vector<Double> &normal, 
                          const Point &a, const Point &b, const Point &c) {
    Point s(a-b);
    Point t(c-b);
    normal[0]=t[1]*s[2]-t[2]*s[1];
    normal[1]=t[2]*s[0]-t[0]*s[2];
    normal[2]=t[0]*s[1]-t[1]*s[0];
    
    Double L2_normal=sqrt(normal[0]*normal[0]+normal[1]*normal[1]+normal[2]*normal[2]);
    normal[0] /= L2_normal;
    normal[1] /= L2_normal;
    normal[2] /= L2_normal;
  }

  /// a-->b-->c. no fix orientation.
  void calcNormal2Surface_Mat(Vector<Double> &normal, const Matrix<Double> &ptCoord){
    
    Point s(ptCoord[0][0], ptCoord[1][0], ptCoord[2][0]);
    Point t(ptCoord[0][2], ptCoord[1][2], ptCoord[2][2]);
    Point b(ptCoord[0][1], ptCoord[1][1], ptCoord[2][1]);
    s -= b;
    t -= b;
    normal[0]=t[1]*s[2]-t[2]*s[1];
    normal[1]=t[2]*s[0]-t[0]*s[2];
    normal[2]=t[0]*s[1]-t[1]*s[0];
    
    Double L2_normal=sqrt(normal[0]*normal[0]+normal[1]*normal[1]+normal[2]*normal[2]);
    normal[0] /= L2_normal;
    normal[1] /= L2_normal;
    normal[2] /= L2_normal;
  }

  UInt defineRefinements(const Double tolElem, const Double tolTotal,
                         const UInt noOfChilds){
    Double tmp = log( tolElem/tolTotal ) / log( (Double) noOfChilds );
    return (UInt)tmp + 1;
  }

  Double CalcArea(Elem * ptE, Grid * ptgrid, const UInt level){
    Double         area = 0;
    BaseFE         * ptelem = ptE->ptElem;
    const StdVector<UInt> & connect = ptE->connect;
    Matrix<Double> ptCoord;
    UInt        nrIntPnts = ptelem->GetNumIntPoints();
    const Vector<Double> & intWeights = ptelem->GetIntWeights();
    UInt        i;
    Double         jacDet;

    ptgrid->GetElemNodesCoord(ptCoord,connect);

    for (i=0; i<nrIntPnts; i++) {
      jacDet = ptelem->CalcJacobianDetAtIp(i+1, ptCoord, NULL);
      area +=jacDet*intWeights[i];
    }

    return area;
  }

  std::string ToValidXML(const std::string& input)
  {
    std::string out = input;
    std::replace_if(out.begin(), out.end(), boost::is_any_of(", |"), '_');
    return out;
  }


  Double NormL2(const Double* data, const UInt size)
  {
    Double result = 0.0;
    for(UInt i = 0; i < size; i++)
      result += data[i] * data[i];

    return std::sqrt(result);
  }


  /** Compares if two doubles are close to each other */
  bool close(Double d1, Double d2)
  {
    return std::abs(d1-d2) < 1e-6;
  }

  /** Compared if two complex are close (if both the real and imaginary part are close) */
  bool close(Complex c1, Complex c2)
  {
    return close(c1.real(), c2.real()) && close(c1.imag(), c2.imag());
  }

  bool IsNoise(Double val)
  {
    return std::abs(val) < 1e-13;
  }
  
  bool IsNoise(Complex val)
  {
    return IsNoise(val.real()) && IsNoise(val.imag());
  }

  bool IsNoise(Integer val)
  {
    return false;
  }

  bool IsNoise(UInt val)
  {
    return false;
  }

  double Average(const double* data, unsigned int size)
  {
    double sum = 0.0;
    for(unsigned int i = 0; i < size; i++)
      sum += *(data + i);

    return sum / size;
  }

  double StandardDeviation(const double* data, unsigned int size)
  {
    // expected value
    double e = Average(data, size);

    // variance
    double v = 0.0;
    for(unsigned int i = 0; i < size; i++)
      v += (*(data + i) - e) * (*(data + i) - e);

    return v / size;

    return sqrt(v);
  }

  void Assign(Matrix<Double>& target, const Matrix<Double>& other, const Double factor)
  {
    const unsigned int orows(other.GetNumRows());
    const unsigned int ocols(other.GetNumCols());
    target.Resize(orows, ocols);
    for(UInt r = 0; r < orows; ++r)
      for(UInt c = 0; c < ocols; ++c)
        target[r][c] = factor * other[r][c];
  }

  void Assign(Matrix<Complex>& target, const Matrix<Complex>& other, const Complex factor)
  {
    const unsigned int orows(other.GetNumRows());
    const unsigned int ocols(other.GetNumCols());
    target.Resize(orows, ocols);
    for(UInt r = 0; r < orows; ++r)
      for(UInt c = 0; c < ocols; ++c)
        target[r][c] = factor * other[r][c];
  }

  void Assign(Matrix<Complex>& target, const Matrix<Double>& other, const Complex factor)
  {
    const unsigned int orows(other.GetNumRows());
    const unsigned int ocols(other.GetNumCols());
    target.Resize(orows, ocols);
    for(UInt r = 0; r < orows; ++r)
      for(UInt c = 0; c < ocols; ++c)
        target[r][c] = factor * other[r][c];
  }

  void Assign(Matrix<Complex>& target, const Matrix<Double>& other, const Double factor)
  {
    const unsigned int orows(other.GetNumRows());
    const unsigned int ocols(other.GetNumCols());
    target.Resize(orows, ocols);
    for(UInt r = 0; r < orows; ++r)
      for(UInt c = 0; c < ocols; ++c)
        target[r][c] = factor * other[r][c];
  }

  void Assign(Vector<Double>& target, const Vector<Double>& other, const Double factor)
  {
    UInt n = other.GetSize();
    target.Resize(n);
    for(UInt i = 0; i < n; i++)
        target[i] = factor * other[i];
  }

  void Assign(Vector<Complex>& target, const Vector<Complex>& other, const Double factor)
  {
    UInt n = other.GetSize();
    target.Resize(n);
    for(UInt i = 0; i < n; i++)
        target[i] = factor * other[i];
  }

  void Assign(Vector<Complex>& target, const Vector<Double>& other, const Double factor)
  {
    UInt n = other.GetSize();
    target.Resize(n);
    for(UInt i = 0; i < n; i++)
        target[i] = factor * other[i];
  }


  int MemoryUsage(bool peak)
  {
    // if the file dows not exist we return 0
    std::ifstream file("/proc/self/status", std::ifstream::in);

    std::string data;
    while(file)
    {
      file >> data;
      if(data == (peak ? "VmPeak:" : "VmSize:"))
      {
        file >> data; // read next value
        try
        {
          return lexical_cast<int>(data);
        }
        catch(bad_lexical_cast &)
        {
          return 0;
        }
      }
    }

    return 0;
  }



  double SmoothMax(double left, double right, double beta)
  {
    assert(beta > 0 || beta == -1.0);
    if(beta == -1.0) return std::max(left, right);

    // the continuous Kreisselmeier and Steinhauser max approximation taken
    // from Sigmund; Morphology-based black and white filters for topology optimization; 2007
    // x = log ( sum(exp(beta * x_i)) / sum 1) / beta

    return std::log(0.5 * (std::exp(left * beta) + std::exp(right * beta))) / beta;
  }

  double SmoothMax(const StdVector<double>& values, double beta)
  {
    assert(beta > 0 || beta == -1.0);
    assert(values.GetSize() > 0);

    double res = 0.0;

    if(beta == -1.0)
    {
      res = values[0];
      for(unsigned int i = 1, n = values.GetSize(); i < n; i++)
        res = std::max(res, values[i]);
    }
    else
    {
      // see  SmoothMax(double left, double right, double beta)
      // x = log ( sum(exp(beta * x_i)) / sum 1) / beta
      double sum = 0.0;
      for(unsigned int i = 0, n = values.GetSize(); i < n; i++)
        sum += std::exp(values[i] * beta);

      res = std::log(sum / (double) values.GetSize()) / beta;
    }

    // LOG_DBG3(tools) << "SmoothMax v=" << values.ToString() << " beta=" << beta << " -> " << res;
    // std::cout << "SmoothMax v=" << values.ToString() << " beta=" << beta << " -> " << res << std::endl;
    return res;
  }

  double SmoothMin(double left, double right, double beta)
  {
    assert(beta > 0 || beta == -1.0);
    assert(right > 0 && left > 0);

    if(beta == -1.0)
      return std::min(left, right);

    // see comment in CalcMaxApproximation()
    return 1.0 - std::log(0.5 * (std::exp((1.0 - left) * beta) + std::exp((1.0 - right) * beta))) / beta;
  }

  double SmoothMin(const StdVector<double>& values, double beta)
  {
    assert(beta > 0 || beta == -1.0);
    assert(values.GetSize() > 0);
    // see  SmoothMax(double left, double right, double beta)

    if(beta == -1.0)
    {
      double t = values[0];
      for(unsigned int i = 1, n = values.GetSize(); i < n; i++)
        t = std::min(t, values[i]);
      return t;
    }

    double sum = 0.0;
    for(unsigned int i = 0, n = values.GetSize(); i < n; i++)
      sum += std::exp((1.0 - values[i]) * beta);

    return 1.0 - std::log(sum / (double) values.GetSize()) / beta;
  }


  double DerivSmoothMax(double left, double right, double beta, int derive)
  {
    assert(derive == -1 || derive == 1); // left or right
    assert(beta > 0);

    double exp_left  = std::exp(left * beta);
    double exp_right = std::exp(right * beta);

    if(derive == -1)
      return exp_left / (exp_left + exp_right);
    else
      return exp_right / (exp_left + exp_right);
  }

  double DerivSmoothMax(const StdVector<double>& values, double beta, unsigned int derive)
  {
    assert(beta > 0);

    double my_exp = -1.0;
    double sum = 0.0;

    if(values.GetSize() == 1)
      return 1.0;

    for(unsigned int i = 0, n = values.GetSize(); i < n; i++)
    {
      double v = std::exp(values[i] * beta);
      sum += v;
      if(derive == i) my_exp = v; // not derive ist not a window based index
    }
    assert(my_exp != -1.0);

    return my_exp / sum;
  }

  double DerivSmoothMin(double left, double right, double beta, int derive)
  {
    assert(derive == -1 || derive == 1); // left or right
    double exp_left  = std::exp((1.0 - left) * beta);
    double exp_right = std::exp((1.0 - right) * beta);

    if(derive == -1)
      return exp_left / (exp_left + exp_right);
    else
      return exp_right / (exp_left + exp_right);
  }

  double DerivSmoothMin(const StdVector<double>& values, double beta, unsigned int derive)
  {
    double my_exp = -1.0;
    double sum = 0.0;

    if(values.GetSize() == 1)
      return 1.0;

    for(unsigned int i = 0, n = values.GetSize(); i < n; i++)
    {
      double v = std::exp((1.0 - values[i]) * beta);
      sum += v;
      if(derive == i) my_exp = v;
    }
    assert(my_exp != -1.0);

    return my_exp / sum;
  }

  double SmoothAbs(double x, double eps)
  {
    assert(eps >= 0);
    return std::sqrt(x*x + eps*eps) - eps;
  }

  double DerivSmoothAbs(double x, double eps)
  {
    assert(eps >= 0);
    assert(abs(x) + eps > 0);
    return x / std::sqrt(x*x + eps*eps);
  }

  VTKStructuredPoints::VTKStructuredPoints(Integer i, Integer j, Integer k, const std::string& scalars, const std::string& vectors)
  {
    assert(i > 0 && j > 0 && k > 0);

    i_max = i;
    j_max = j;
    k_max = k;

    data_.Resize(i);
    for(Integer x = 0; x < i_max; x++)
    {
      data_[x].Resize(j_max);
      for(Integer y = 0; y < j_max; y++)
        data_[x][y].Resize(k_max);
    }

    scalar_label_ = scalars;
    vector_label_ = vectors;
  }

  void VTKStructuredPoints::Set(Integer i, Integer j, Integer k, const Double value)
  {
    data_[i][j][k].first = value;
  }

  void VTKStructuredPoints::Set(Integer i, Integer j, Integer k, const Point& value)
  {
    data_[i][j][k].second = value;
  }

  void VTKStructuredPoints::Set(Integer i, Integer j, Integer k, const Vector<Double>& value)
  {
    assert(value.GetSize() == 2 || value.GetSize() == 3);
    Point& p = data_[i][j][k].second;
    p[0] = value[0];
    p[1] = value[1];
    p[2] = value.GetSize() == 2 ? 0.0 : value[2];
  }


  void VTKStructuredPoints::Write(std::ostream& out) const
  {
    out << "# vtk DataFile Version 2.0\n"
        << "VTK Structured Points from CFS\n"
        << "ASCII\n"
        << "DATASET STRUCTURED_POINTS\n"
        << "DIMENSIONS " << i_max << " " << j_max << " " << k_max << "\n"
        << "ASPECT_RATIO 1 1 1\n"
        << "ORIGIN 0 0 0\n"
        << "POINT_DATA " << (i_max * j_max * k_max) << "\n";

    if(scalar_label_ != "")
    {
      out << "SCALARS " << scalar_label_ << " float\n"
          << "LOOKUP_TABLE default\n";
      for(Integer k = 0; k < k_max; k++)
        for(Integer j = 0; j < j_max; j++)
          for(Integer i = 0; i < i_max; i++)
            out << data_[i][j][k].first << "\n";
    }

    if(vector_label_ != "")
    {
      out << "VECTORS " << vector_label_ << " float\n";
      for(Integer k = 0; k < k_max; k++)
        for(Integer j = 0; j < j_max; j++)
          for(Integer i = 0; i < i_max; i++)
          {
            const Point& vec = data_[i][j][k].second;
            out << vec[0] << " " << vec[1] << " " << vec[2] << "\n";
          }
    }
  }





}// namespace CoupledField
