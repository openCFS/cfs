#ifndef TOOLS_2001
#define TOOLS_2001

#include "environment.hh"
#include <string>
#include <vector>

namespace CoupledField
{

  //! Function for working with errors 
  void Error(const Char * Text, const Char * const filename=NULL,
               const Integer numline=0);

  //! Function for output warnings 
  void Warning(const Char * Text, const Char * const filename=NULL,
               const Integer numline=0);

  //! Reader of Fnc from conf-file
 pfn1var FncReader(const std::string namefnc);

 //! Absolute value of number
template<class T>
T abs(T x) { return (x>0 ? x: -x); } 

  //! square of value
template<class T>
T sqr(T x) { return x*x;}

template<Integer dim>
class Point
{
public:
  Point(){;}
  ~Point(){;}
  
  Point & operator=(const Point & t);
  Point & operator+(const Point & t);

  Double &operator[](Integer i){return p[i];} 
private:
  Double p[dim];
};

#ifdef __GNUC__
template class Point<2>;
template class Point<3>;
#endif

template<Integer dim>
Double dist(Point<dim> a, Point<dim> b);

template<Integer dim>
void PrintPoint(Point<dim> point, std::ostream * out); 

/// a-->b
void calcNormal2Line(std::vector<Double> & normal,Point<2> a,Point<2> b);

Double ScalarMult(std::vector<Double> a, std::vector<Double> b);

} // end of CoupledField

#endif
