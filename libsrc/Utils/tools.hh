#ifndef TOOLS_2001
#define TOOLS_2001

#include <string>
#include <iostream>

#include "General/environment.hh"

namespace CoupledField {

  template<class TYPE> class Matrix; 
  template<class TYPE> class Vector;
  template<class TYPE> class StdVector;

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
  void SplitStringList( std::string list, StdVector<std::string> &strVec,
                        Char delimiter = ',' );


  //! Reader of Fnc from conf-file
  /*!
    \param namefnc name of the reading function
  */
  pfn1var FncReader(const std::string namefnc);
  
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
  template<UInt dim>
  class Point {

  public:
    //! constructor
    Point(){
      for(UInt i=0; i<dim; i++)
        p[i]=0.0;
    }

    //!destructor
    ~Point(){;}

    //!
    Point & operator=(const Point & t);
    //!
    Point & operator+(const Point & t);
    //!
    Point & operator-(const Point & t); 

    //! return coordinate number i
    Double &operator[](UInt i){return p[i];} 

    //! return coordinate number i
    Double operator[](UInt i) const {return p[i];} 

  private:
    Double p[dim];
  };

#if defined(__GNUC__) 
  template class Point<2>;
  template class Point<3>;
#endif

  //! calculate distance between two points
  template<UInt dim>
  Double dist(Point<dim> a, Point<dim> b);

  //! calculate distance between two points embedded in matrix

  Double dist_Mat(Matrix<Double> a);

  //! print point in ofstream
  /*!
    \param point point
    \param out pointer to a ofstream
  */
  template<UInt dim>
  void PrintPoint(Point<dim> point, std::ostream * out); 

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
  void calcNormal2Line(Vector<Double> & normal,Point<2> a, Point<2> b);

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
  void calcNormal2Surface(Vector<Double> & normal,Point<3> a,Point<3> b,
                          Point<3> c);


  //! calculate normal to surface element using matrix parameter
  /*!
    \param normal normal
    \param ptCoord matrix containing vertices of surface element
  */
  void calcNormal2Surface_Mat(Vector<Double> & normal,Matrix<Double> ptCoord);


  char * c_string(const std::string & s);

  /// prints formatted header including name, version, date
  void PrintCFSHeader(std::ostream & out);

  Double Sin(const Double x);
  Double Cos(const Double x);


  /// prints formatted header including name, version, date
  void PrintCFSHeader(std::ostream & out);

} // end of CoupledField

#endif
