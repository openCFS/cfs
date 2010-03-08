// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef TOOLS_2001
#define TOOLS_2001

#include <string>
#include <iostream>

#include "General/environment.hh"
#include "Utils/Point.hh"

namespace CoupledField {

  template<class TYPE> class Matrix;
  template<class TYPE> class Vector;
  template<class TYPE> class StdVector;
  
  // to prevent using the abs from stdlib.h which is only for int!
  using std::abs;

  // =========================================================================
  //     ERROR / WARNING HANDLING
  // =========================================================================


  //@{
  //! \name Error / Warning Handling

  //! Function for issuing an warning message, will not terminate the program

  //! This function can be used to issue a warning message. It is actually
  //! only a shortcut for calling the Warning() method of the WriteInfo class.
  //! \param Text     Text of the warning message
  //! \param filename This is intended to contain the
  //!                 name of the module/file in which the problem occured. The
  //!                 __FILE__ macro should be inserted in the call. The
  //!                 argument is optional.
  //! \param numline  This is intended to contain the
  //!                 number of the code line of the module/file in which the
  //!                 problem occured. The __LINE__ macro should be inserted in
  //!                 the call. The argument is optional.
  void Warning( const char* Text, const char* const filename = NULL,
                const UInt numline = 0 );


  //! Function for issuing an warning message, will not terminate the program

  //! This function can be used to issue a warning message. This variant of
  //! Warning() obtains the warning message from the global <b>(*warning)</b>
  //! string stream. As in the other Warning() variant the actual work is
  //! delegated to WriteInfo::Warning().
  //! Apropriate usage of this function looks like this
  //! \code
  //! (*warning) << "Index k ='" << k << "' already treated!";
  //! Warning( __FILE__, __LINE__ );
  //! \endcode
  //! \param filename this is intended to contain the
  //!                 name of the module/file in which the problem occured. The
  //!                 __FILE__ macro should be inserted in the call.
  //! \param numline  This is intended to contain the
  //!                 number of the code line of the module/file in which the
  //!                 problem occured. The __LINE__ macro should be inserted in
  //!                 the call.
  void Warning( const char *const filename, const UInt numline );

  //@}

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
  //     VARIOUS OTHER METHODS AND CLASSES
  // =========================================================================

  /** Compares if two doubles are close to each other */
  bool close(Double d1, Double d2);

  /** Compared if two complex are close (if both the real and imaginary part are close) */
  bool close(Complex c1, Complex c2);

  /** identifies numerical noise */
  bool IsNoise(Double val);

  bool IsNoise(Complex val);

  bool IsNoise(int val);

  bool IsNoise(UInt val);

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

  //! calculate distance between two points embedded in matrix

  Double dist_Mat(const Matrix<Double> &a);

  // calculation area or volume of element
  struct Elem;
  class Grid;
  Double CalcArea(Elem * ptElem, Grid * ptgrid);

  // define number of refinement for the element
  UInt defineRefinements(const Double tolElem, const Double tolTotal,
                         const UInt noOfChilds);

  //! calculate the normal to line with following orientation: a-->b
  /*!
    \param normal normal
    \param a,b pointes
  */
  void calcNormal2Line(Vector<Double> & normal, const Point &a, const Point &b);

  //! calculate the normal to line with following orientation: a-->b
  /*!
    \param normal normal
    \param a,b points embedded in Matrix
  */
  void calcNormal2Line_Mat(Vector<Double> & normal, const Matrix<Double> &a);

  //! calculate normal to surface element
  /*!
    \param normal normal
    \param a,b,c vertices of element
  */
  void calcNormal2Surface(Vector<Double> & normal,
                          const Point &a, const Point &b, const Point &c);


  //! calculate normal to surface element using matrix parameter
  /*!
    \param normal normal
    \param ptCoord matrix containing vertices of surface element
  */
  void calcNormal2Surface_Mat(Vector<Double> & normal, const Matrix<Double> &ptCoord);

  /** Assigns the multiple of a matrix to another matrix. The target is resized.
   * This is silly copy & pase code. And is for the non-mixed variants already
   * in Matrix::Assign(). Anybody knows how to do the mixed variant in Matrix? */
  void Assign(Matrix<Double>& target,  const Matrix<Double>&  other, const Double factor);
  void Assign(Matrix<Complex>& target, const Matrix<Complex>& other, const Complex factor);
  void Assign(Matrix<Complex>& target, const Matrix<Double>&  other, const Complex factor);
  void Assign(Matrix<Complex>& target, const Matrix<Double>&  other, const Double factor);


  /// prints formatted header including name, version, date
  void PrintCFSHeader(std::ostream & out);

  /** Calculates the L2 norm of a array. This is for cases where we
   * don't use one of our vectors. E.g. with IPOPT */
  Double NormL2(const Double* data, const UInt size);

  /** Calculate the average of an array */
  double Average(const double* data, unsigned int size);

  /** Calculate the Standard Deviation of an array */
  double StandardDeviation(const double* data, unsigned int size);

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

} // end of CoupledField

#endif
