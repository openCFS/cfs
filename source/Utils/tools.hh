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
  
  //! Function for issuing an error message and terminating program execution.

  //! This function can be used to issue an error message and terminate
  //! execution of the program. In fact this method is only a shortcut for
  //! calling the Error() method of the WriteInfo class.
  //! \param Text     Text of the error message
  //! \param filename This is intended to contain the
  //!                 name of the module/file in which the error occured. The
  //!                 __FILE__ macro should be inserted in the call. The
  //!                 argument is optional.
  //! \param numline  This is intended to contain the
  //!                 number of the code line of the module/file in which the
  //!                 error occured. The __LINE__ macro should be inserted in
  //!                 the call. The argument is optional.
  void Error( const Char *Text, const Char *const filename,
              const UInt numline );


  //! Function for issuing an error message and terminating program execution.

  //! This function can be used to issue an error message and terminate
  //! execution of the program. This variant of Error() obtains the error
  //! message from the global <b>(*error)</b> string stream. As in the
  //! other Error() variant the actual work is delegated to WriteInfo::Error().
  //! Apropriate usage of this function looks like this
  //! \code
  //! (*error) << "Cannot open file '" << fname << "' for reading";
  //! Error( __FILE__, __LINE__ );
  //! \endcode
  //! \param filename This is intended to contain the
  //!                 name of the module/file in which the error occured. The
  //!                 __FILE__ macro should be inserted in the call.
  //! \param numline  This is intended to contain the
  //!                 number of the code line of the module/file in which the
  //!                 error occured. The __LINE__ macro should be inserted in
  //!                 the call.
  void Error( const Char *const filename, const UInt numline );


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
  void Warning( const Char* Text, const Char* const filename = NULL,
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
  void Warning( const Char *const filename, const UInt numline );

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
                        const Char delimiter = ',' );
  
  //! Converts a string into a double value
  Double String2Double( const std::string & val); 
  
  //! Converts a string into an integer value
  Integer String2Int( const std::string & val); 
  
  //! Converts a string into an unsigned integer value
  UInt String2UInt( const std::string & val);
  
  //! Converts a string StdVector into a double StdVector
  void String2Double( StdVector<Double> & retVal, 
                      const StdVector<std::string> & val ); 

  //! Converts a string StdVector into a double vector
  void String2Double( Vector<Double> & retVal, 
                      const StdVector<std::string> & val ); 
  
  //! Converts a double StdVector into a string StdVector
  void Double2String( StdVector<std::string> & retVal, 
                      const StdVector<Double> & val );

  //! Converts a double Vector into a string StdVector
  void Double2String( StdVector<std::string> & retVal, 
                      const Vector<Double> & val );

  //! Converts a string StdVector into an integer StdVector
  void String2Int( StdVector<Integer> & retVal, 
                   const StdVector<std::string> & val );

  //! Converts a string StdVector into an integer Vector
  void String2Int( Vector<Integer> & retVal, 
                   const StdVector<std::string> & val );

  //! Converts an integer StdVector into a string StdVector
  void Int2String( StdVector<std::string> & retVal, 
                   const StdVector<Integer> & val );
  
  //! Converts an integer Vector into a string StdVector
  void Int2String( StdVector<std::string> & retVal, 
                   const Vector<Integer> & val );
  
  //! Converts a string vector into an unsigned integer vector
  void String2UInt( StdVector<UInt> & retVal, 
                    const StdVector<std::string> & val );

  //! Converts a string StdVector into an unsigned integer StdVector
  void String2UInt( Vector<UInt> & retVal, 
                    const StdVector<std::string> & val );
  
  //! Converts an unsigned integer StdVector into a string StdVector
  void UInt2String( StdVector<std::string> & retVal, 
                    const StdVector<UInt> & val );
  
  //! Converts an unsigned integer Vector into a string StdVector
  void UInt2String( StdVector<std::string> & retVal, 
                    const Vector<UInt> & val );

  //@}

  // =========================================================================
  //     VARIOUS OTHER METHODS AND CLASSES
  // =========================================================================
 

  //! Absolute value of number
  template<class T>
  T abs(T x) { return (x>0 ? x: -x); } 

  //! square of value
  template<class T>
  T sqr(T x) { return x*x;}

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
        p[i]=0.0;
    } 

    //!destructor
    ~Point(){;}

    /** resets the values */
    void SetZero() {
      for(UInt i=0; i<3; i++)
        p[i]=0.0;
    } 

    //!
    Point & operator=(const Point & t);
    //!
    Point & operator+=( const Point & t );
    //!
    Point  operator+(const Point & t);
    //!
    Point  operator-(const Point & t); 

    //! return coordinate number i
    Double &operator[](UInt i) {
      assert(i < 3);
      return p[i];
    } 

    //! return coordinate number i
    Double operator[](UInt i) const {
      assert(i < 3);
      return p[i];
    } 
    
    //! calculate distance between two points
    static Double dist(const Point& a, const Point& b);

    //! calculate distance to another point
    Double dist(const Point& other) const;
    
    /** Lists the content 
     * @return the form "(0.3;4.3;0.0)" but no digit control */ 
    std::string ToString() const;

  private:
    Double p[3];
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

  /// prints formatted header including name, version, date
  void PrintCFSHeader(std::ostream & out);

  Double Sin(const Double x);
  Double Cos(const Double x);

  /** Calculates the L2 norm of a array. This is for cases where we
   * don't use one of our vectors. E.g. with IPOPT */ 
  double NormL2(const double* data, unsigned int size);
  
  /** Calculate the average of an array */
  double Average(const double* data, unsigned int size);
  
  /** Calculate the Standard Deviation of an array */
  double StandardDeviation(const double* data, unsigned int size);


  /// prints formatted header including name, version, date
  void PrintCFSHeader(std::ostream & out);

  
  
} // end of CoupledField

#endif
