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
  //! \param Text     text of the error message
  //! \param filename another text string. This is intended to contain the
  //!                 name of the module/file in which the error occured. The
  //!                 __FILE__ macro should be inserted in the call. The
  //!                 argument is optional.
  //! \param numline  another text string. This is intended to contain the
  //!                 number of the code line of the module/file in which the
  //!                 error occured. The __LINE__ macro should be inserted in
  //!                 the call. The argument is optional.
  void Error( const Char * Text, const Char * const filename = NULL,
	      const Integer numline = 0 );


  //! Function for issuing an warning message, will not terminate the program

  //! This function can be used to issue an error message and terminate
  //! execution of the program. In fact this method is only a shortcut for
  //! calling the Warning() method of the WriteInfo class.
  //! \param Text     text of the warning message
  //! \param filename another text string. This is intended to contain the
  //!                 name of the module/file in which the problem occured. The
  //!                 __FILE__ macro should be inserted in the call. The
  //!                 argument is optional.
  //! \param numline  another text string. This is intended to contain the
  //!                 number of the code line of the module/file in which the
  //!                 problem occured. The __LINE__ macro should be inserted in
  //!                 the call. The argument is optional.
  void Warning( const Char * Text, const Char * const filename = NULL,
		const Integer numline = 0 );


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
  T pow(T x, Integer power) 
  { T p=x;
  if (!power)
    return 1;
 
  for (Integer i=2; i<=power; i++)
    p*=x;
  return p;
  }

  //! class for working with points. 
  /*!
    \param dim dimension of the point. 2.. point in 2D, 3.. point in 3D
  */
  template<Integer dim>
  class Point {

  public:
    //! constructor
    Point(){;}
    //!deconstructor
    ~Point(){;}

    //!
    Point & operator=(const Point & t);
    //!
    Point & operator+(const Point & t);
    //!
    Point & operator-(const Point & t); 

    //! return coordinate number i
    Double &operator[](Integer i){return p[i];} 
  private:
    Double p[dim];
  };

#ifdef __GNUC__
  template class Point<2>;
  template class Point<3>;
#endif

  //! calculate distance between two points
  template<Integer dim>
  Double dist(Point<dim> a, Point<dim> b);

  //! calculate distance between two points embedded in matrix

  Double dist_Mat(Matrix<Double> a);

  //! print point in ofstream
  /*!
    \param point point
    \param out pointer to a ofstream
  */
  template<Integer dim>
  void PrintPoint(Point<dim> point, std::ostream * out); 

  // calculation area or volume of element
  struct Elem;
  class Grid;
  Double CalcArea(Elem * ptElem, Grid * ptgrid, const Integer level);

  // define number of refinement for the element
  Integer defineRefinements(const Double tolElem, const Double tolTotal,
			    const Integer noOfChilds);

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
