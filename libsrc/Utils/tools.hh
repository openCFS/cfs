#ifndef TOOLS_2001
#define TOOLS_2001

#include "environment.hh"

namespace CoupledField
{

  //! Function for working with errors 
  void Error(const Char * Text, const Char * const filename=NULL,
               const Integer numline=0);

 //! Absolute value of number
template<class T>
T abs(T x) { return (x>0 ? x: -x); } 

//! Point in 2D
class Point2D
{
public:

  Double x,y;

  Point2D & operator=(const Point2D & t);
  Point2D & operator+(const Point2D & t);
  Boolean is2D(){return TRUE;}

};

//! Point in 3D
class Point3D
{
public:

  Double x,y,z;
  Point3D  & operator=(const Point3D & t);
  Point3D  & operator+(const Point3D & t);
  Boolean is2D(){return FALSE;}
  
};

//! Print Point3D (overloading)
 void PrintPoint(const Point3D point, std::ostream * out) ;

//! Print Point2D (overloading)
 void PrintPoint(const Point2D point, std::ostream * out) ;

} // end of CoupledField

#endif
