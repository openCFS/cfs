#include <fstream>
#include <iostream>

#include "tools.hh"

namespace CoupledField {


template<Integer dim>
void PrintPoint(Point<dim> point, std::ostream * out) 
 {
   Integer i;
   for(i=0; i<dim; i++) 
  (*out)<< "   " << point[i];
 }

#ifdef __GNUC__
template void PrintPoint(Point<2>, std::ostream *);
template void PrintPoint(Point<3>, std::ostream *);
#endif
 
void Error(const Char * const Text, const Char * const filename,
                      const Integer numline)
{
 std::cerr << "\033[31mERROR:\033[0m " << Text;
 if (filename) { std::cerr <<"( " << filename <<" ";
                 if (numline) std::cerr << numline;
                 std::cerr << ")";}
 std::cerr << std::endl;
 exit(-1);
}

void Warning(const Char * const Text, const Char * const filename,
                      const Integer numline)
{
 std::cerr << "\033[31mWARNING:\033[0m " << Text;
 if (filename) { std::cerr <<"( " << filename <<" ";
                 if (numline) std::cerr << numline;
                 std::cerr << ")";}
 std::cerr << std::endl;
}

template<Integer dim>
Point<dim> & Point<dim>::operator+(const Point<dim>&t)
{
  Integer i;
  for (i=0; i<dim; i++)
    p[i]+=t.p[i];	
  return *this;
}

template<Integer dim>
Point<dim> & Point<dim>::operator=(const Point<dim>&t)
{
  Integer i;
  for (i=0; i<dim; i++)
    p[i]=t.p[i];
  return *this;
}

Double Sin(const Double x)
{
  return sin(x);
}

Double Cos(const Double x)
{
  return cos(x);
}

pfn1var FncReader(const std::string namefnc)
{
  if (namefnc == "sin") return Sin;
  else {
    if (namefnc == "cos") return Cos;
    else Error("This function is not predefined in CFS++");
  }
}

template<Integer dim>
Double dist(Point<dim> a,Point<dim> b)
{
  Double preSqrt=0;
  Integer i;
  for (i=0; i<dim; i++)
    preSqrt+=sqr(a[i]-b[i]);
 return sqrt(preSqrt);
}

#ifdef __GNUC__
template Double dist(Point<2>,Point<2>);
template Double dist(Point<3>,Point<3>);
#endif


/// a-->b
void calcNormal2Line(std::vector<Double> & normal,Point<2> a,Point<2> b)
{
  Double distance=dist(a,b);
  normal[0]=(b[1]-a[1])/distance;
  normal[1]=(a[0]-b[0])/distance;
}

Double ScalarMult(std::vector<Double> a, std::vector<Double> b)
{
  Integer i;
  Double res=0;
  if (a.size()!=b.size()) 
    Error("Dimensions in scalar multiplivcation of 2 vectors are not comparable",__FILE__,__LINE__);

  for (i=0; i<a.size(); i++)
    res+=a[i]*b[i];

  return res;
}

} // close namespace CoupledField
