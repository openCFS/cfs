#ifndef TOOLS_2001
#define TOOLS_2001

#include <string>
#include <iostream>

#include "General/environment.hh"

namespace CoupledField
{
  template<class TYPE> class Matrix; 
  template<class TYPE> class Vector;
  template<class TYPE> class StdVector;

  //! Function for working with errors. The program stops after the error-message. 
  /*!
    \param Text text of the error message
    \param filename name of the file with the error. it can be omitted.
    \param numline number of line with the error. it can be omitted.
  */
  void Error(const Char * Text, const Char * const filename=NULL,
               const Integer numline=0);

  //! Function for output warnings. The execution of the program is continued after the message.
  /*!
    \param Text text of the warning message
    \param filename name of the file with the error. it can be omitted.
    \param numline number of line with the error. it can be omitted.
  */
  void Warning(const Char * Text, const Char * const filename=NULL,
               const Integer numline=0);

  //! Reader of Fnc from conf-file
  /*!
    \param namefnc name of the reading function
  */
 pfn1var FncReader(const std::string namefnc);



  //! Function for splitting up a comma-separated
  //! string into a StdVector of single entries
  //! \param list (input) single string containing multiple strings,
  //! separated by 'delimiter'
  //! \param strVec (output) vector of the single entries
  //! \param delimiter (input) character used as delimiter
  
  void SplitStringList(std::string list, 
		       StdVector<std::string> & strVec,
		       Char delimiter = ',');
  
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
class Point
{
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
void calcNormal2Line_Mat(Vector<Double> & normal,Matrix<Double> a);

  //! calculate normal to surface element
  /*!
    \param normal normal
    \param a,b,c vertices of element
  */
void calcNormal2Surface(Vector<Double> & normal,Point<3> a,Point<3> b, Point<3> c);


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
