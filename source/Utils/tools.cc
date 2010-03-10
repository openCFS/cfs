// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <fstream>
#include <iostream>
#include <math.h>
#include <boost/tokenizer.hpp>

#include "tools.hh"
#include "MatVec/matrix.hh"
#include "Elements/elements_header.hh"
#include "Domain/elem.hh"
#include "Domain/grid.hh"
#include "DataInOut/WriteInfo.hh"
#include "General/exception.hh"


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

}// namespace CoupledField
