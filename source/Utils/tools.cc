// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <fstream>
#include <iostream>
#include <math.h>

#include "tools.hh"
#include "Matrix/matrix.hh"
#include "Elements/elements_header.hh"
#include "Domain/elem.hh"
#include "Domain/grid.hh"
#include "DataInOut/WriteInfo.hh"
#include "General/exception.hh"


namespace CoupledField {

  // =======================
  //   Issue Error Message
  // =======================
  void Error( const Char * Text, const Char * const filename,
              const UInt numline ) {
    Info->Error( Text, filename, numline );
  }

  // ====================================================
  //   Issue Error Message (using global string stream)
  // ====================================================
  void Error( const Char *const filename, const UInt numline ) {

    // Obtain error message and clear string stream
    std::string errMsg = "";
    if ( error != NULL ) {
      errMsg = error->str();
      error->str( "" );
    }
    else {
      errMsg = "<Caller did not provide error description>";
    }

    // Delegate work to WriteInfo::Error()
    Info->Error( errMsg.c_str(), filename, numline );
  }

  // =========================
  //   Issue Warning Message
  // =========================
  void Warning( const Char *Text, const Char * const filename,
                const UInt numline ) {
    Info->Warning( Text, filename, numline );
    //std::cerr << "\033[31mWARNING:\033[0m " << Text;
  }

  // ======================================================
  //   Issue Warning Message (using global string stream)
  // ======================================================
  void Warning( const Char *const filename, const UInt numline ) {

    // Obtain error message and clear string stream
    std::string warnMsg = "";
    if ( warning != NULL ) {
      warnMsg = warning->str();
      warning->str( "" );
    }
    else {
      warnMsg = "<Caller did not provide error description>";
    }

    // Delegate work to WriteInfo::Warning()
    Info->Warning( warnMsg.c_str(), filename, numline );
  }

  // =========================================
  //   Split a string into a list of strings
  // =========================================
  void SplitStringList( const std::string &list, StdVector<std::string> &strVec,
                        const Char delimiter ) {

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

  Point & Point::operator+=(const Point& t) {
    UInt i;
    for (i=0; i<3; i++)
      data[i]+=t.data[i];
    return *this;
  }


  Point Point::operator+(const Point& t) {
    UInt i;
    Point ret;
    for (i=0; i<3; i++)
      ret.data[i] = t.data[i] + data[i];
    return ret;
  }


  Point & Point::operator=(const Point& t) {
    UInt i;
    for (i=0; i<3; i++)
      data[i]=t.data[i];
    return *this;
  }

  Point Point::operator-(const Point& t) {
    UInt i;
    for (i=0; i<3; i++)
      data[i]-=t.data[i];
    return *this;
  }

  Point Point::operator*(const Double factor) {
    for (UInt i=0; i<3; i++)
      data[i] *= factor;
    return *this;
  }

  Point& Point::operator*=(const Double factor) {
    for (UInt i=0; i<3; i++)
      data[i] *= factor;
    return *this;
  }


  Point Point::operator/(const Double factor) {
    for (UInt i=0; i<3; i++)
      data[i] /= factor;
    return *this;
  }


  Point& Point::operator/=(const Double factor) {
    for (UInt i=0; i<3; i++)
      data[i] /= factor;
    return *this;
  }

  std::string Point::ToString() const
  {
     std::ostringstream os;
     os << "(" << data[0] << ";" << data[1] << ";" << data[2] << ")";
     return os.str();
  }

  Double dist_Mat(Matrix<Double> a) {
    Double preSqrt=0;
    UInt i;
    UInt k=a.GetSizeRow();
    // std::cout<<"tools.cc:size of matrix: "<<k<<std::endl;
    for (i=0; i<k; i++)
      preSqrt+= (a[i][0]-a[i][1]) * (a[i][0]-a[i][1]);
    return sqrt(preSqrt);
  }

  /// a-->b
  void calcNormal2Line(Vector<Double> & normal,Point a,Point b) {
    Double distance=Point::dist(a,b);
    normal[0]=(b[1]-a[1])/distance;
    normal[1]=(a[0]-b[0])/distance;
  }

  /// a-->b
  void calcNormal2Line_Mat(Vector<Double> & normal,Matrix<Double> a)
  {
    Double distance=dist_Mat(a);
    // std::cout<<"distance :"<<distance<<std::endl;

    normal[0]=(a[1][1]-a[1][0])/distance;
    normal[1]=(a[0][1]-a[0][0])/distance;
  }


  /// a-->b-->c. no fix orientation.
  void calcNormal2Surface(Vector<Double> & normal,Point a,Point b,
                          Point c) {
    Point t,s;
    Double L2_normal;
    s=(a-b);
    t=(c-b);
    normal[0]=t[1]*s[2]-t[2]*s[1];
    normal[1]=t[2]*s[0]-t[0]*s[2];
    normal[2]=t[0]*s[1]-t[1]*s[0];
    L2_normal=sqrt(normal[0]*normal[0]+normal[1]*normal[1]+normal[2]*normal[2]);
    normal[0]=normal[0]/L2_normal;
    normal[1]=normal[1]/L2_normal;
    normal[2]=normal[2]/L2_normal;
  }

  /// a-->b-->c. no fix orientation.
  void calcNormal2Surface_Mat(Vector<Double> & normal,Matrix<Double> ptCoord){
    Point t,s,a,b,c;

    for (UInt i=0;i<3;i++) {
      a[i]=ptCoord[i][0];
      b[i]=ptCoord[i][1];
      c[i]=ptCoord[i][2];
    }

    Double L2_normal;
    s=(a-b);
    t=(c-b);
    normal[0]=t[1]*s[2]-t[2]*s[1];
    normal[1]=t[2]*s[0]-t[0]*s[2];
    normal[2]=t[0]*s[1]-t[1]*s[0];
    L2_normal=sqrt(normal[0]*normal[0]+normal[1]*normal[1]+normal[2]*normal[2]);
    normal[0]=normal[0]/L2_normal;
    normal[1]=normal[1]/L2_normal;
    normal[2]=normal[2]/L2_normal;
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


  // **************** Conversion Function *****************

  Double String2Double( const std::string & val) {


    Double retVal = 0.0;
    char *endp;
    retVal = strtod(val.c_str(), &endp);
    if (! (val.c_str() != endp && *endp == '\0') )
        EXCEPTION("could not convert " << val << " into a double");

    return retVal;
  }

  Integer String2Int( const std::string & val) {
    Integer retVal = 0;
    char *endp;
    retVal = strtol(val.c_str(), &endp, 0);
    if (! (val.c_str() != endp && *endp == '\0') )
        EXCEPTION("could not convert " << val << " into a integer");

    return retVal;
  }

  UInt String2UInt( const std::string & val) {

    UInt retVal=0;
    char *endp;
    retVal = strtoul(val.c_str(), &endp, 0);
    if (! (val.c_str() != endp && *endp == '\0') )
        EXCEPTION("could not convert " << val << " into a unsigned integer");

    return retVal;
  }

  double NormL2(const double* data, unsigned int size)
  {
    double result = 0.0;
    for(unsigned int i = 0; i < size; i++)
      result += *(data + i) * *(data + i);

    return result;
  }


  /** Compares if two doubles are close to each other */
  bool close(Double d1, Double d2)
  {
    return std::abs(d1-d2) < 1e-6;
  }

  /** Compared if two complex are close (if both the real and imaginary part are close) */
  bool close(Complex c1, Complex c2)
  {
    return std::abs(c1.real()-c2.real()) < 1e-6
        && std::abs(c1.imag()-c2.imag()) < 1e-6;
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

  void Assign(Matrix<Double>& target, const Matrix<Double>& other, Double factor)
  {
    target.Resize(other);
    for(UInt r = 0; r < other.GetSizeRow(); r++)
      for(UInt c = 0; c < other.GetSizeCol(); c++)
        target[r][c] = factor * other[r][c];
  }

  void Assign(Matrix<Complex>& target, const Matrix<Complex>& other, Complex factor)
  {
    target.Resize(other);
    for(UInt r = 0; r < other.GetSizeRow(); r++)
      for(UInt c = 0; c < other.GetSizeCol(); c++)
        target[r][c] = factor * other[r][c];
  }

  void Assign(Matrix<Complex>& target, const Matrix<Double>& other, Complex factor)
  {
    target.Resize(other.GetSizeRow(), other.GetSizeCol());
    for(UInt r = 0; r < other.GetSizeRow(); r++)
      for(UInt c = 0; c < other.GetSizeCol(); c++)
        target[r][c] = factor * other[r][c];
  }

  void Assign(Matrix<Complex>& target, const Matrix<Double>& other, Double factor)
  {
    target.Resize(other.GetSizeRow(), other.GetSizeCol());
    for(UInt r = 0; r < other.GetSizeRow(); r++)
      for(UInt c = 0; c < other.GetSizeCol(); c++)
        target[r][c] = factor * other[r][c];
  }




}// namespace CoupledField
