// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <fstream>
#include <iostream>
#include <math.h>
#include <stdio.h> // MemoryUsage
#include <boost/tokenizer.hpp>

#include "tools.hh"
#include "MatVec/matrix.hh"
#include "Elements/elements_header.hh"
#include "Domain/elem.hh"
#include "Domain/grid.hh"
#include "DataInOut/WriteInfo.hh"
#include "General/exception.hh"


namespace CoupledField {

  // =========================
  //   Issue Warning Message
  // =========================
  void Warning( const char *Text, const char * const filename,
                const UInt numline ) {
    Info->Warning( Text, filename, numline );
    //std::cerr << "\033[31mWARNING:\033[0m " << Text;
  }

  // ======================================================
  //   Issue Warning Message (using global string stream)
  // ======================================================
  void Warning( const char *const filename, const UInt numline ) {

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
  
  Point & Point::operator=(const Point& t) {
    for (UInt i=0; i<3; i++)
      data[i]=t.data[i];
    return *this;
  }

  Point & Point::operator+=(const Point& t) {
    UInt i;
    for (i=0; i<3; i++)
      data[i]+=t.data[i];
    return *this;
  }

  Point Point::operator+(const Point &t) {
    return Point(data[0]+t.data[0], data[1]+t.data[1], data[2] + t.data[2]);
  }
  
  Point Point::operator+(const Point &t) const {
    return Point(data[0]+t.data[0], data[1]+t.data[1], data[2] + t.data[2]);
  }

  Point & Point::operator-=(const Point & t) {
    for (UInt i=0; i<3; i++)
      data[i]-=t.data[i];
    return *this;
  }

  Point Point::operator-(const Point &t) {
    return Point(data[0]-t.data[0], data[1]-t.data[1], data[2]-t.data[2]);
  }
  
  Point Point::operator-(const Point &t) const {
    return Point(data[0]-t.data[0], data[1]-t.data[1], data[2]-t.data[2]);
  }

  Point& Point::operator*=(const Double factor) {
    for (UInt i=0; i<3; i++)
      data[i] *= factor;
    return *this;
  }
  
  Point Point::operator*(const Double factor) {
    return Point(data[0]*factor, data[1]*factor, data[2]*factor);
  }
  
  Point Point::operator*(const Double factor) const {
    return Point(data[0]*factor, data[1]*factor, data[2]*factor);
  }
  
  Point& Point::operator/=(const Double factor) {
    for (UInt i=0; i<3; i++)
      data[i] /= factor;
    return *this;
  }
  
  Point Point::operator/(const Double factor) {
    return Point(data[0]/factor, data[1]/factor, data[2]/factor);
  }

  Point Point::operator/(const Double factor) const {
    return Point(data[0]/factor, data[1]/factor, data[2]/factor);
  }

  std::string Point::ToString() const
  {
     std::ostringstream os;
     os << "(" << data[0] << ";" << data[1] << ";" << data[2] << ")";
     return os.str();
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
    Double distance=Point::dist(a,b);
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

  int MemoryUsage()
  {
    //! Process id of this program
    pid_t pid_num = getpid();  // get process id
    std::string pid = lexical_cast<std::string>(pid_num);

    // Get memory consumption
    std::string cmd = "ps -F " + pid + " | grep " + pid;

    // popen is not ANSI-C but shall be on all unix system in the std lib!
    // if this makes problem use defines!

    FILE* pipe = popen(cmd.c_str(), "r");
    if(pipe == NULL) return 0;

    // we have only a single line!
    char line[128]; // the line buffer
    if(fgets(line, 128, pipe) == NULL) return 0;
    pclose(pipe);

    // (ID        PID  PPID  C    SZ   RSS PSR STIME TTY      STAT   TIME CMD)
    // now extract the size
    // fwein     6288  6283  2 135425 111908 1 10:22 ?        Sl     1:57 /usr/lib64/firefox/firefox
    int count = 0;
    char_separator<char> sep(" \t"); // spaces, closing and colon
    tokenizer<char_separator<char> > tok(std::string(line), sep);
    for(tokenizer<char_separator<char> >::iterator beg=tok.begin(); beg!=tok.end(); ++beg)
    {
      // we need the 5th element (nice and red haired)
      // std::cout << "token '" << *beg << "'" << std::endl;
      if(++count < 5) continue;

      return lexical_cast<int>(*beg);
    }

    return 0; // token not found
  }


}// namespace CoupledField
