#include <string>
#include <fstream.h>
#include <iostream.h>

#include "general_head.hh"
#include "tools.hh"

namespace CoupledField {

/// Print Point3D
void PrintPoint(const Point3D point, ostream * out) 
{
 (*out)<< "   " << point.x << "   " <<  point.y << "   " <<  point.z ;
}

/// Print Point2D
void PrintPoint(const Point2D point, ostream * out) 
{
 (*out)<< "   " << point.x << "   " << point.y ;
}

void Error(const Char * const Text, const Char * const filename,
                      const Integer numline)
{
 cerr << "ERROR: " << Text;
 if (filename) { cerr <<"( " << filename <<" ";
                 if (numline) cerr << numline;
                 cerr << ")";}
 cerr << endl;
 exit(-1);
}

template Double abs(Double);

Point3D & Point3D::operator+(const Point3D & t)
{
  x+=t.x;
  y+=t.y;
  z+=t.z;
  return *this;
}
 
Point3D & Point3D::operator=(const Point3D & t)
{
  x=t.x;
  y=t.y;
  z=t.z;
  return *this;
}
 
 
Point2D & Point2D::operator+(const Point2D & t)
{
  x+=t.x;
  y+=t.y;
  return *this;
}
 
Point2D & Point2D::operator=(const Point2D & t)
{
  x=t.x;
  y=t.y;
  return *this;
}

} // close namespace CoupledField
