#ifndef TOOLS_2001
#define TOOLS_2001

#include "environment.hh"
#include <string>
#include <vector>

namespace CoupledField
{

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

 //! Absolute value of number
template<class T>
T abs(T x) { return (x>0 ? x: -x); } 

  //! square of value
template<class T>
T sqr(T x) { return x*x;}

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

  //! return coordinate number i
  Double &operator[](Integer i){return p[i];} 
private:
  Double p[dim];
};

#ifdef __GNUC__
template class Point<2>;
template class Point<3>;
#endif

  //! calculate distance between two pointes
template<Integer dim>
Double dist(Point<dim> a, Point<dim> b);

  //! print point in ofstream
  /*!
    \param point point
    \param out pointer to a ofstream
  */
template<Integer dim>
void PrintPoint(Point<dim> point, std::ostream * out); 

//! calculate the normal to line with following orientation: a-->b
  /*!
    \param normal normal
    \param a,b pointes
  */
void calcNormal2Line(std::vector<Double> & normal,Point<2> a, Point<2> b);

  //! scalar multiplication of 2 vectors
Double ScalarMult(std::vector<Double> a, std::vector<Double> b);

} // end of CoupledField

#endif
