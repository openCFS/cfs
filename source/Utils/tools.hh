// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef TOOLS_2001
#define TOOLS_2001

#include <string>
#include <iostream>

#include "General/environment.hh"

//! \file tools.hh
//! Defines various methods and classes:
//! - Error / Warning Handling
//! - String conversion functions
//! - Geometric Poin<> class and related helper functions

namespace CoupledField {

  template<class TYPE> class Matrix;
  template<class TYPE> class Vector;
  template<class TYPE> class StdVector;

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
  //! \name String conversion functions

   //! Convert anything to a standard string
  //! This auxilliary method can convert any data type to a standard string,
  //! if the << operator has been overloaded for the respective type.
  template<class T> std::string GenStr( const T &value ) {
    std::ostringstream mystream;
    mystream << value;
    return mystream.str();
  }

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

  //! Converts a string into a double value
  Double String2Double( const std::string & val);

  //! Converts a string into an integer value
  Integer String2Int( const std::string & val);

  //! Converts a string into an unsigned integer value
  UInt String2UInt( const std::string & val);
    //@}

  // =========================================================================
  //     VARIOUS OTHER METHODS AND CLASSES
  // =========================================================================

  /** Compares if two doubles are close to each other */
  bool close(Double d1, Double d2);

  /** Compared if two complex are close (if both the real and imaginary part are close) */
  bool close(Complex c1, Complex c2);

  //! Absolute value of number
  template<class T>
  T abs(T x) { return (x>0 ? x: -x); }

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

  //! class for working with points.
  /*!
    \param dim dimension of the point. 2.. point in 2D, 3.. point in 3D
  */
  class Point {

  public:
    //! constructor
    Point() {
      for(UInt i=0; i<3; i++)
        data[i]=0.0;
    }

    //!destructor
    ~Point(){;}

    /** resets the values */
    void SetZero() {
      for(UInt i=0; i<3; i++)
        data[i]=0.0;
    }

    //!
    Point & operator=(const Point & t);
    //!
    Point & operator+=( const Point & t );
    //!
    Point  operator+(const Point & t);
    //!
    Point  operator-(const Point & t);

    /** scale the point */ 
    Point  operator*(double factor); 

    Point& operator*=(double factor); 

    /** scale the point */ 
    Point  operator/(double factor); 

    /** scale the point */ 
    Point&  operator/=(double factor);   

    //! return coordinate number i
    Double &operator[](UInt i) {
      assert(i < 3);
      return data[i];
    }

    //! return coordinate number i
    Double operator[](UInt i) const {
      assert(i < 3);
      return data[i];
    }

    //! calculate distance between two points
    inline static Double dist(const Point& a, const Point& b) { 
      Double preSqrt = 0.0; 
      for(UInt i=0; i<3; i++) 
        preSqrt+= (a.data[i]-b.data[i]) * (a.data[i]-b.data[i]); 
      return sqrt(preSqrt); 
    } 

    //! calculate distance to another point
    Double dist(const Point& other) const {
      return Point::dist(*this, other);
    }

    /** Lists the content
     * @return the form "(0.3;4.3;0.0)" but no digit control */
    std::string ToString() const;

    Double data[3];
  };


  //! calculate distance between two points embedded in matrix

  Double dist_Mat(Matrix<Double> a);

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
  void calcNormal2Line(Vector<Double> & normal,Point a, Point b);

  //! calculate the normal to line with following orientation: a-->b
  /*!
    \param normal normal
    \param a,b points embedded in Matrix
  */
  void calcNormal2Line_Mat(Vector<Double> & normal, Matrix<Double> a);

  //! calculate normal to surface element
  /*!
    \param normal normal
    \param a,b,c vertices of element
  */
  void calcNormal2Surface(Vector<Double> & normal,Point a,Point b,
                          Point c);


  //! calculate normal to surface element using matrix parameter
  /*!
    \param normal normal
    \param ptCoord matrix containing vertices of surface element
  */
  void calcNormal2Surface_Mat(Vector<Double> & normal,Matrix<Double> ptCoord);

  /** Assigns the multiple of a matrix to another matrix. The target is resized.
   * This is silly copy & pase code. And is for the non-mixed variants already
   * in Matrix::Assign(). Anybody knows how to do the mixed variant in Matrix? */
  void Assign(Matrix<Double>& target,  const Matrix<Double>&  other, Double factor);
  void Assign(Matrix<Complex>& target, const Matrix<Complex>& other, Complex factor);
  void Assign(Matrix<Complex>& target, const Matrix<Double>&  other, Complex factor);
  void Assign(Matrix<Complex>& target, const Matrix<Double>&  other, Double factor);


  /// prints formatted header including name, version, date
  void PrintCFSHeader(std::ostream & out);

  /** Calculates the L2 norm of a array. This is for cases where we
   * don't use one of our vectors. E.g. with IPOPT */
  double NormL2(const double* data, unsigned int size);

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



} // end of CoupledField

#endif
