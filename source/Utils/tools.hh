#ifndef TOOLS_2001
#define TOOLS_2001

#include <cmath>
#include <string>
#include <iostream>
#include <boost/lexical_cast.hpp>

#include "General/Environment.hh"
#include "Optimization/EigenInfo.hh"

// C++ is so poor that there even is no really standard for pi - > why not????
#ifndef M_PI
  #include <boost/math/constants/constants.hpp>
  #define M_PI boost::math::constants::pi<double>()
#endif

#include <def_use_cuda.hh>
#include "def_use_openmp.hh"
#ifdef USE_OPENMP
#include <omp.h>
#endif

#include <set> // a gcc7 runner has issues with std::set, try including it late

namespace CoupledField {

  class SingleVector;
  class BaseVector;

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

  
  /** boost based SplitStringList() which consideres almost all whitespaces */
  void SplitStringListWhitespace(const std::string &s, StdVector<std::string> &strVec);

  /** convert a string such that is becomes a valid filename. Might need extensions */
  std::string ConvertToFilename(std::string org);

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
     return rad / 180.0 * M_PI;
   }
   
  //! Convert rad => grad
  inline Double Rad2Grad( Double rad ) {
     return rad / M_PI * 180.0;
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
                               const std::string& phase,
							   bool isInRad = false );

  //! Convert (ampl,phase) => imag (strings)
  std::string AmplPhaseToImag( const std::string& val, 
                               const std::string& phase,
							   bool isInRad = false );

  // ------ vector versions -----
  //! Convert (ampl,phase) => (real,imag) (strings vectors)
  void AmplPhaseToRealImag( const StdVector<std::string>& val, 
                            const StdVector<std::string>& phase,
                            StdVector<std::string>& real,
                            StdVector<std::string>& imag );
    
  // =========================================================================
  //     VARIOUS OTHER METHODS AND CLASSES
  // =========================================================================

  // generate vector with random numbers
  StdVector<double> GenerateRandomVector(size_t size, double first = 0.0, double last = 1.0);

  /** Compares if two doubles are close to each other */
  inline bool close(Double d1, Double d2) { return std::abs(d1-d2) < 1e-6; }

  /** clang complains about std::abs(c1 - c2) if i* is unsigned int but it would be used for templated Matrix::IsSymmetric()
   * @param eps not used but there for compatibelity*/
  inline bool close(int i1, int i2, double eps = 0.0) { return i1 == i2; }
  inline bool close(unsigned int i1, unsigned int i2, double eps = 0.0) { return i1 == i2; }

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

  /** http://stackoverflow.com/questions/1903954/is-there-a-standard-sign-function-signum-sgn-in-c-c */
  template <typename T> int sgn(T val) {
      return (T(0) < val) - (val < T(0));
  }

  // shortcuts for pow with int base - strange this is not in boost?!
  inline double Pow2(double x) { return x * x; }
  inline int Pow2(int x) { return x * x; }
  inline unsigned int Pow2(unsigned int x) { return x * x; }
  inline double Pow3(double x) { return x * x * x; }
  inline int Pow3(int x) { return x * x * x; }
  inline unsigned int Pow3(unsigned int x) { return x * x * x; }


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
   * This is silly copy & paste code. And is for the non-mixed variants already
   * in Matrix::Assign(). Anybody knows how to do the mixed variant in Matrix? */
  void Assign(Matrix<Double>& target,  const Matrix<Double>&  other, const Double factor);
  void Assign(Matrix<Complex>& target, const Matrix<Complex>& other, const Complex factor);
  void Assign(Matrix<Complex>& target, const Matrix<Double>&  other, const Complex factor);
  void Assign(Matrix<Complex>& target, const Matrix<Double>&  other, const Double factor);

  void Assign(Vector<Double>& target, const Vector<Double>& other, const Double factor);
  void Assign(Vector<Complex>& target, const Vector<Complex>& other, const Double factor);
  void Assign(Vector<Complex>& target, const Vector<Double>& other, const Double factor);

  template<class T>
  void Copy(const StdVector<T>& source, Vector<T>& target) {
    target.Resize(source.GetSize());
    if(!source.IsEmpty())
      std::memcpy(target.GetPointer(), source.GetPointer(), source.GetSize() * sizeof(T));
  }

  template<class TYPE, class TYPE2>
  void Add(Matrix<TYPE>& out, const TYPE fac, const Matrix<TYPE2>& other)
  {
   assert(out.GetNumRows() == other.GetNumRows() && out.GetNumCols() == other.GetNumCols());
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

  /** all vectors need to be either complex or real but the same.
   * out += fac1 * vec1 + fac2 * vec2 */
  void Add(BaseVector& out, double fac1, const BaseVector& vec1, double fac2, const BaseVector& vec2);

  /** @see Add() above. If vectors are real, the real part from the scalars is used */
  void Add(BaseVector& out, Complex fac1, const BaseVector& vec1, Complex fac2, const BaseVector& vec2);

  double Inner(const Vector<double>& v1, const Vector<double>& v2);
  Complex Inner(const Vector<Complex>& v1, const Vector<Complex>& v2);
  Complex Inner(const Vector<double>& v1, const Vector<Complex>& v2);
  Complex Inner(const Vector<Complex>& v1, const Vector<double>& v2);

  /** Search for the smallest value within a row
   * @param value set when given
   * @param set the info if given to be used for output
   * @return the 0-based column index */
  unsigned int SearchMinMax(const Matrix<double>& mat, unsigned int row, bool minimum, double* val = NULL, EigenInfo* info = NULL);


  /** transforms a complex matrix to its complex conjugate */
  void Conj(Matrix<Complex>& mat);


  //! Convert a path pattern into a regular expression

  //! Converts a path pattern into a regex by escaping all special regex
  //! characters not used in path patterns.
  std::string PathPatternToRegEx(const std::string & pattern);

  /** makes sure the string is a valid xml element and attribute name */
  std::string ToValidXML(const std::string& input);

  /** Calculates the L2 norm of a array. This is for cases where we
   * don't use one of our vectors. E.g. with IPOPT */
  double NormL2(const Double* data, const UInt size);

  double NormL2(const Double* data, const Double* data2, const UInt size);

  double NormL2(const SingleVector* data, const SingleVector* data2);

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

  /** generic ToString for STL containers (set, vector, ...).
   * A pitty, the containers doen't bring this by themselves :( */
  template <class Cont>
  std::string ToStringCont(const Cont& cont, const std::string& sep = " ")
  {
    std::ostringstream os;
    for(const auto& it : cont)
      os << it << sep;
    return os.str();
  }

  /** We have types for Direction in EnvironmentTypes.hh and ShapeParamElement::Dof
   * but this is an easy int to x,y,z conversion for loops */
  std::string Dof(int d);
  int Dof(const std::string& dim);

  /** Returns the sign of a value 
   * @return 0 if 0 or +/- 1 */ 
  inline int Sign(int a) { return (a == 0) ? 0 : (a < 0 ? -1 : 1); } 
  inline int Sign(Double a) { return (a == 0.0) ? 0 : (a < 0.0 ? -1 : 1); } 

  /** Compares the sign for equality - note that we have three signs! -1, 0, 1 */ 
  inline bool SameSign(Double a, Double b) { return (a < 0.0) == (b < 0.0); } 

  /// prints formatted header including name, version, date
  void PrintCFSHeader(std::ostream & out);
  

  /** Determines the current memory consumption for Linux, macOS and Windows
   * Rather slow on Linux, only current memory on macOS
   * @param peak peak memory, otherwise current memory
   * @return the memory in KBytes or 0 if there was a problem */
  int MemoryUsage(bool peak);

  /** Calculates the continuous Kreisselmeier and Steinhauser max approxmiation for two values.
   * @param beta - -1 is special and makes real max, otherwise beta needs to be > 0
   * @param normalize - applies normalization by factor 0.5 */
  double SmoothMax(double left, double right, double beta, bool normalize = false);

  /** Calculates the continuous Kreisselmeier and Steinhauser max approximation for arbitrary values.
  * @param beta - -1 is special and makes real max, otherwise beta needs to be > 0
  * @param normalize - applies normalization by number of values */
  double SmoothMax(const StdVector<double>& values, double beta, bool normalize = false);

  /** @param deriv -1 for left or 1 for right value to derive for */
  double DerivSmoothMax(double left, double right, double beta, int derive);

  /** @param derive index within values. */
  double DerivSmoothMax(const StdVector<double>& values, double beta, unsigned int derive);

  /** @see SmoothMax(double, double, double, bool) */
  double SmoothMin(double left, double right, double beta, bool normalize = false);

  /** @see SmoothMax(const StdVector<double>&, double, bool) */
  double SmoothMin(const StdVector<double>& values, double beta, bool normalize = false);

  /** @see DerivSmoothMax(double, double, double, int) */
  double DerivSmoothMin(double left, double right, double beta, int derive);

  /** @see DerivSmoothMax(const StdVector<double>&, double, unsigned int) */
  double DerivSmoothMin(const StdVector<double>& values, double beta, unsigned int derive);

  /** Calculates an approximation of the abs function:A(x) = sqrt(x^2 + eps^2) - eps
   * As used in Poulsen; A new scheme for imposing a minimum length scale in topology optimization; 2003
   * @param eps small, e.g. 10% of x */
  double SmoothAbs(double x, double eps);

  /** derivative of
   * @see CalcAbsApproximation() */
  double DerivSmoothAbs(double x, double eps);

  inline unsigned int Product(const StdVector<unsigned int>& vec)
  {
    unsigned int prod = 1;
    for (unsigned int i = 0; i < vec.GetSize(); i++)
      prod *= vec[i];

    return prod;
  }

  /** uses the global domain->GetMathParser() to evaluate an expression */
  double MathParse(const std::string& expr);

  /** returns the weights for numerical quadrature, close Newton-Cotes formula.
   * For equal spacing: 0.0,1.0; 0.0,0.5,1.0; 0,0.3333,0.6666,1.0; ...
   * @return [], [.5,0.5], [1/6,2/6,1/6], ...] for order 0(invalid), 1, 2, ... with -,2,3, ... weights
   * the result has 11 entries there the last, res[10], has 11 entries
   * @see https://de.wikipedia.org/wiki/Newton-Cotes-Formeln */
  StdVector<Vector<double> > GetNewtonCotes();

  /** Multiple subscripts from linear index */
  void Sub2Ind(Vector<unsigned int> size, StdVector<int> sub, unsigned int &ind);

  /** Linear index from multiple subscripts */
  void Ind2Sub(Vector<unsigned int> size, unsigned int ind, StdVector<int> &sub);

  /** Return a linspace similar to numpy.python where the ends are always included
   * s=1, e=4, n=5 -> 1, 1.75, 2.5, 3.25, 4. Order can be reversed! */
  Vector<double> Linspace(double start, double end, int n);

  /** Return a logspace by giving the bases, similar to numpy.logspace
   * e.g. s=1,e=4,n=4: -> 1e1, 1e2, 1e3, 1e4 or s=2, e=-2, n=5 -> 1e2, 1e1, 1e0, 1e-1, 1e-2
   * is actually std::pow(10, Linspace(s,e,n)) */
  Vector<double> LogspaceBase(double start_exponent, double end_exponent, int n);

  // omp_get_thread_num
  inline unsigned int GetThreadNum()
  {
  #ifdef USE_OPENMP
     return omp_get_thread_num();
   #else
     return 0;
   #endif

  }

  inline bool UseOpenMP()
  {
    #ifdef USE_OPENMP
      return true;
    #else
      return false;
    #endif
  }

  inline bool UseCuda()
  {
    #ifdef USE_CUDA
      return true;
    #else
      return false;
    #endif
  }


} // end of CoupledField

#endif
