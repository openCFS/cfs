#ifndef TOOLS_2001
#define TOOLS_2001

/** LightTools.hh is includes almost nothing, here we include slightly more
 * but still try to have not too many dependencies */

#include <cmath>
#include <charconv>
#include <string>
#include <array>
#include <map>
#include <boost/lexical_cast.hpp>

#include <def_use_cuda.hh>
#include "def_use_openmp.hh"
#ifdef USE_OPENMP
  #include <omp.h>
#endif

#include "General/defs.hh"
#include "Utils/LightTools.hh"
#include "Utils/StdVector.hh"

namespace CoupledField 
{
  class SingleVector;
  class BaseVector;
  struct Elem;
  struct EigenInfo;

  template<class TYPE> class Matrix;
  template<class TYPE> class Vector;
  
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

  std::string ToString(const StdVector<Elem*>& data);

  template <class TYPE>
  std::string ToString(const StdVector<Vector<TYPE> >& data, bool new_line = false);

  template <class TYPE>
  std::string ToString(const StdVector<StdVector<TYPE> >& data, bool new_line = false);

  template <class TYPE>
  std::string ToString(const std::vector<TYPE>& data)
  {
    std::ostringstream os;
    os << "[";
    for(unsigned int i = 0; i < data.size(); i++)
      os << data[i] << ((i < data.size()-1) ? ", ": "]");
    return os.str();
  }

  /** converts data arrays to strings such that they can be copy & pasted from log to matlab.
   * Redundant to StdVector::ToString() but there the complex special implementation was not possible */
  template <class TYPE>
  std::string ToString(const TYPE* data, unsigned int size);

  template <class TYPE, size_t N>
  std::string ToString(const std::array<TYPE, N>& data) {
    return ToString(data.data(), data.size());
  }

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
  * @see the SoftMax function (common in AI) and the p-norm are available in FeatureMappingDesign.hh - move when needed
  * @param beta - -1 is special and makes real max, otherwise beta needs to be > 0
  * @param normalize - applies normalization by number of values */
  double SmoothMax(const Vector<double>& values, double beta, bool normalize = false);

  /** @param deriv -1 for left or 1 for right value to derive for */
  double DerivSmoothMax(double left, double right, double beta, int derive);

  /** @param derive index within values. */
  double DerivSmoothMax(const Vector<double>& values, double beta, unsigned int derive);

  /** @see SmoothMax(double, double, double, bool) */
  double SmoothMin(double left, double right, double beta, bool normalize = false);

  /** @see SmoothMax(const StdVector<double>&, double, bool) */
  double SmoothMin(const StdVector<double>& values, double beta, bool normalize = false);
  double SmoothMin(const Vector<double>& values, double beta, bool normalize = false);
  double SmoothMin(const double* values, size_t size, double beta, bool normalize = false);

  /** @see DerivSmoothMax(double, double, double, int) */
  double DerivSmoothMin(double left, double right, double beta, int derive);

  /** @see DerivSmoothMax() */
  double DerivSmoothMin(const Vector<double>& values, double beta, unsigned int derive);
  double DerivSmoothMin(const StdVector<double>& values, double beta, unsigned int derive);
  double DerivSmoothMin(const double* values, size_t size, double beta, unsigned int derive);

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

  /** Convert a map to vector indexed by key
   * Assumes keys are non-negative integers that can be used as vector indices
   * @param inputMap Input map with integer-like keys
   * @return Vector where index corresponds to key value - empty keys have default value (0, nullptr, ...)
   * @throws Exception if any key is negative or keys are not contiguous (gaps are filled with default values) */  
  template<typename KEY_TYPE, typename VALUE_TYPE>
  StdVector<VALUE_TYPE> MapToVector(const std::map<KEY_TYPE, VALUE_TYPE>& inputMap)
  {
    if(inputMap.empty()) 
      return StdVector<VALUE_TYPE>();

    // Check that all keys are >= 0
    for(const auto& pair : inputMap) 
      if(pair.first < 0) 
        throw Exception("Map has negative key " + std::to_string(pair.first));

    // Find the maximum key to determine vector size
    KEY_TYPE maxKey = inputMap.rbegin()->first;
    
    // Create vector with size = maxKey + 1, initialized with default values
    StdVector<VALUE_TYPE> result((unsigned int) (maxKey + 1));
    
    // Fill vector with values at corresponding indices
    for(const auto& pair : inputMap) 
      result[(unsigned int) pair.first] = pair.second;

    return result;
  }

  /** test if a string is of the given value. E.g. int or double without try/catch */
  template<class T>
  bool IsType(const std::string& s) 
  {
    T value;
    auto [ptr, ec] = std::from_chars(s.data(), s.data() + s.size(), value);
    return ec == std::errc() && ptr == s.data() + s.size();
  }

  // omp_get_thread_num
  inline unsigned int GetThreadNum()
  {
  #ifdef USE_OPENMP
     return omp_get_thread_num();
   #else
     return 0;
   #endif

  }

} // end of CoupledField

#endif
