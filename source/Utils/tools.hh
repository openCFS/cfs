#ifndef TOOLS_2001
#define TOOLS_2001

#include <string>
#include <iostream>
#include <cmath>
#include <boost/lexical_cast.hpp>

#include "General/Environment.hh"

namespace CoupledField {

  template<class TYPE> class Matrix;
  template<class TYPE> class Vector;
  template<class TYPE> class StdVector;
  
  // to prevent using the abs from stdlib.h which is only for int!
  using std::abs;

  // =========================================================================
  //     STRING CONVERSION FUNCTIONS
  // =========================================================================
  //@{



  //! Function for splitting a string into a vector of single entries

  //! This function can be used to split a string into a vector of substrings.
  //! The decomposition into substrings is steered by selecting the delimiter
  //! between the substrings. By default the delimiter is a comma. The
  //! delimiter and all spaces are discarded. Note that the current
  //! implementation cannot deal with spaces as delimiter.
  //! \param list      (input)  single string containing multiple substrings
  //!                           separated by a 'delimiter'
  //! \param strVec    (output) vector of the single substrings
  //! \param delimiter (input)  character used as delimiter
  void SplitStringList( const std::string &list, StdVector<std::string> &strVec,
                        const char delimiter = ',' );

  
  //! Wrap string in braces
  inline std::string Bracket( const std::string& xpr ) {
    return "("+xpr+")";
  }
    //@}

  // =========================================================================
  //  ANGLE CONVERSION
  // =========================================================================
  
   //! Convert grad => rad
   inline Double Grad2Rad( Double rad ) {
     return rad / 180.0 * PI;
   }
   
  //! Convert rad => grad
  inline Double Rad2Grad( Double rad ) {
     return rad / PI * 180.0;
   }
  
  // =========================================================================
  //  COMPLEX CONVERSION
  // =========================================================================
  
  //! Convert (real,imag) => amplitude
  Double RealImagToAmpl( const Complex& in );
  
  //! Convert (real,imag) => phase
  Double RealImagToPhase( const Complex& in );
  
  //! Convert (ampl,phase) => (real,imag)
  Complex AmplPhaseToComplex( Double val, Double phase );
  
  //! Convert (ampl,phase) => real
  Double AmplPhaseToReal( Double val, Double phase );
  
  //! Convert (ampl,phase) => imag
  Double AmplPhaseToImag( Double val, Double phase );
  
  //! Convert (ampl,phase) => real (strings)
  std::string AmplPhaseToReal( const std::string& val, 
                               const std::string& phase );

  //! Convert (ampl,phase) => imag (strings)
  std::string AmplPhaseToImag( const std::string& val, 
                               const std::string& phase );
  
  // ------ vector versions -----
  //! Convert (ampl,phase) => (real,imag) (strings vectors)
  void AmplPhaseToRealImag( const StdVector<std::string>& val, 
                            const StdVector<std::string>& phase,
                            StdVector<std::string>& real,
                            StdVector<std::string>& imag );
    
  // =========================================================================
  //     VARIOUS OTHER METHODS AND CLASSES
  // =========================================================================


  /** Compares if two doubles are close to each other */
  inline bool close(Double d1, Double d2) { return std::abs(d1-d2) < 1e-6; }


  /** separate eps version to be faster! */
  inline bool close(Double d1, Double d2, double eps) { return std::abs(d1-d2) < eps; }

  /** Compared if two complex are close (if both the real and imaginary part are close) */
  inline bool close(Complex c1, Complex c2) { return close(c1.real(), c2.real()) && close(c1.imag(), c2.imag()); }

  inline bool close(Complex c1, Complex c2, double eps) { return close(c1.real(), c2.real(), eps) && close(c1.imag(), c2.imag(), eps); }

  /** identifies numerical noise */
  inline bool IsNoise(Double val) { return std::abs(val) < 1e-13; }

  inline bool IsNoise(Complex val) { return IsNoise(val.real()) && IsNoise(val.imag()); }

  inline bool IsNoise(int val) { return false; }

  inline bool IsNoise(UInt val) { return false; }

  //! power of value
  template<class T>
  T pow(T x, UInt power)
  { T p=x;
    if (!power)
      return 1;

    for (UInt i=2; i<=power; i++)
      p*=x;
    return p;
  }

  /** http://stackoverflow.com/questions/1903954/is-there-a-standard-sign-function-signum-sgn-in-c-c */
  template <typename T> int sgn(T val) {
      return (T(0) < val) - (val < T(0));
  }


  //! calculate distance between two points embedded in matrix

  Double dist_Mat(const Matrix<Double> &a);

  // calculation area or volume of element
  struct Elem;
  class Grid;
  Double CalcArea(Elem * ptElem, Grid * ptgrid);

  // define number of refinement for the element
  UInt defineRefinements(const Double tolElem, const Double tolTotal,
                         const UInt noOfChilds);

//  //! calculate the normal to line with following orientation: a-->b
//  /*!
//    \param normal normal
//    \param a,b pointes
//  */
//  void calcNormal2Line(Vector<Double> & normal, const Point &a, const Point &b);
//
//  //! calculate the normal to line with following orientation: a-->b
//  /*!
//    \param normal normal
//    \param a,b points embedded in Matrix
//  */
//  void calcNormal2Line_Mat(Vector<Double> & normal, const Matrix<Double> &a);
//
//  //! calculate normal to surface element
//  /*!
//    \param normal normal
//    \param a,b,c vertices of element
//  */
//  void calcNormal2Surface(Vector<Double> & normal,
//                          const Point &a, const Point &b, const Point &c);
//
//
//  //! calculate normal to surface element using matrix parameter
//  /*!
//    \param normal normal
//    \param ptCoord matrix containing vertices of surface element
//  */
//  void calcNormal2Surface_Mat(Vector<Double> & normal, const Matrix<Double> &ptCoord);

  /** Assigns the multiple of a matrix to another matrix. The target is resized.
   * This is silly copy & pase code. And is for the non-mixed variants already
   * in Matrix::Assign(). Anybody knows how to do the mixed variant in Matrix? */
  void Assign(Matrix<Double>& target,  const Matrix<Double>&  other, const Double factor);
  void Assign(Matrix<Complex>& target, const Matrix<Complex>& other, const Complex factor);
  void Assign(Matrix<Complex>& target, const Matrix<Double>&  other, const Complex factor);
  void Assign(Matrix<Complex>& target, const Matrix<Double>&  other, const Double factor);

  void Assign(Vector<Double>& target, const Vector<Double>& other, const Double factor);
  void Assign(Vector<Complex>& target, const Vector<Complex>& other, const Double factor);
  void Assign(Vector<Complex>& target, const Vector<Double>& other, const Double factor);

  template<class TYPE, class TYPE2>
  void Add(Matrix<TYPE>& out, const TYPE fac, const Matrix<TYPE2>& other)
  {
   assert(out.GetNumRows() != other.GetNumRows() || out.GetNumCols() != other.GetNumCols());
   for(unsigned int r = 0, rn = out.GetNumRows(); r < rn; r++)
     for(unsigned int c = 0, cn = out.GetNumCols(); c < cn; c++)
       out[r][c] += fac * other[r][c];
  }

  template<class TYPE, class TYPE2>
  void Add(Vector<TYPE>& out, const TYPE2 fac, const Vector<TYPE>& other)
  {
   assert(out.GetSize() != other.GetSize());
   for(unsigned int i = 0, rn = out.GetSize(); i < rn; i++)
       out[i] += fac * other[i];
  }


  /** transforms a complex matrix to its complex conjugate */
  void Conj(Matrix<Complex>& mat);


  /** makes sure the string is a valid xml element and attribute name */
  std::string ToValidXML(const std::string& input);

  /** Calculates the L2 norm of a array. This is for cases where we
   * don't use one of our vectors. E.g. with IPOPT */
  Double NormL2(const Double* data, const UInt size);

  /** Calculate the average of an array */
  double Average(const double* data, unsigned int size);

  /** Calculate the Standard Deviation of an array */
  double StandardDeviation(const double* data, unsigned int size);

  template <class TYPE>
  std::string ToString(const StdVector<Vector<TYPE> >& data, bool new_line = false);

  template <class TYPE>
  std::string ToString(const StdVector<StdVector<TYPE> >& data, bool new_line = false);


  /** converts data arrays to strings such that they can be copy & pasted from log to matlab.
   * Redundant to StdVector::ToString() but there the complex special implementation was not possible */
  template <class TYPE>
  std::string ToString(const TYPE* data, unsigned int size);

  /** Returns the sign of a value 
   * @return 0 if 0 or +/- 1 */ 
  inline int Sign(int a) { return (a == 0) ? 0 : (a < 0 ? -1 : 1); } 
  inline int Sign(Double a) { return (a == 0.0) ? 0 : (a < 0.0 ? -1 : 1); } 

  /** Compares the sign for equality - note that we have three signs! -1, 0, 1 */ 
  inline bool SameSign(Double a, Double b) { return (a < 0.0) == (b < 0.0); } 

  /// prints formatted header including name, version, date
  void PrintCFSHeader(std::ostream & out);
  

  /** Determines the current memory consumption.
   * This is done by calling ps and some post processing from a pipe.
   * Runs clearly only on Unix and is rather expensive
   * @param peak peak memory or current memory
   * @return the memory in KBytes or 0 if there was a problem */
  int MemoryUsage(bool peak);

  /** Calculates the continuous Kreisselmeier and Steinhauser max approxmiation for two values.
   * @param beta -1 is special and makes real max, otherwise beta needs to be > 0 */
  double SmoothMax(double left, double right, double beta);

  /** Calculates the continuous Kreisselmeier and Steinhauser max approximation for arbitrary values. */
  double SmoothMax(const StdVector<double>& values, double beta);

  /** @param deriv -1 for left or 1 for right value to derive for */
  double DerivSmoothMax(double left, double right, double beta, int derive);

  /** @param deriv index within values. */
  double DerivSmoothMax(const StdVector<double>& values, double beta, unsigned int derive);

  /** @see CalcMaxApproximation() */
  double SmoothMin(double left, double right, double beta);

  double SmoothMin(const StdVector<double>& values, double beta);

  /** @see DerivSmoothMax() */
  double DerivSmoothMin(double left, double right, double beta, int derive);

  double DerivSmoothMin(const StdVector<double>& values, double beta, unsigned int derive);

  /** Calculates an approximation of the abs function:A(x) = sqrt(x^2 + eps^2) - eps
   * As used in Poulsen; A new scheme for imposing a minimum length scale in topology optimization; 2003
   * @param eps small, e.g. 10% of x */
  double SmoothAbs(double x, double eps);

  /** derivative of
   * @see CalcAbsApproximation() */
  double DerivSmoothAbs(double x, double eps);

} // end of CoupledField

#endif
