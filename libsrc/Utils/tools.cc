#include <fstream>
#include <iostream>
#include <math.h>

#include "tools.hh"
#include <Matrix/matrix.hh>
#include <Elements/elements_header.hh>
#include <Domain/elem.hh>
#include <Domain/grid.hh>

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
 
void Error(const Char * Text, const Char * const filename,
                      const Integer numline)
{
 std::cerr << "\033[31mERROR:\033[0m " << Text;

#ifdef DEBUG
 if (filename) { std::cerr <<"( " << filename <<" ";
                 if (numline) std::cerr << numline;
                 std::cerr << ")";}
#endif

 std::cerr << std::endl;
 exit(-1);
}

void Warning(const Char *Text, const Char * const filename,
                      const Integer numline)
{
 std::cerr << "\033[31mWARNING:\033[0m " << Text;

#ifdef DEBUG
 if (filename) { std::cerr <<"( " << filename <<" ";
                 if (numline) std::cerr << numline;
                 std::cerr << ")";}
#endif

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

template<Integer dim>
Point<dim> & Point<dim>::operator-(const Point<dim>&t)
{
  Integer i;
  for (i=0; i<dim; i++)
    p[i]-=t.p[i];	
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


Double dist_Mat(Matrix<Double> a)
{
  Double preSqrt=0;
  Integer i;
  Integer k=a.getSize();
  // std::cout<<"tools.cc:size of matrix: "<<k<<std::endl;
  for (i=0; i<k; i++)
    preSqrt+=sqr(a[i][0]-a[i][1]);
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

/// a-->b
void calcNormal2Line_Mat(std::vector<Double> & normal,Matrix<Double> a)
{
  Double distance=dist_Mat(a);
  // std::cout<<"distance :"<<distance<<std::endl;
  
  normal[0]=(a[1][1]-a[1][0])/distance;
  normal[1]=(a[0][1]-a[0][0])/distance;
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

/// a-->b-->c. no fix orientation.
void calcNormal2Surface(std::vector<Double> & normal,Point<3> a,Point<3> b, Point<3> c)
{
  Point<3> t,s;
  Double L2_normal; 
  s=(a-b);
  t=(c-b);
  normal[0]=t[1]*s[2]-t[2]*s[1];
  normal[1]=t[2]*s[0]-t[0]*s[2];  
  normal[2]=t[0]*s[1]-t[1]*s[0];
  L2_normal=sqrt(sqr(normal[0])+sqr(normal[1])+sqr(normal[2]));
  normal[0]=normal[0]/L2_normal;
  normal[1]=normal[1]/L2_normal;
  normal[2]=normal[2]/L2_normal;  
}

/// a-->b-->c. no fix orientation.
void calcNormal2Surface_Mat(std::vector<Double> & normal,Matrix<Double> ptCoord)
{
  Point<3> t,s,a,b,c;

  for (Integer i=0;i<3;i++)
    {
      a[i]=ptCoord[i][0];
      b[i]=ptCoord[i][1];      
      c[i]=ptCoord[i][2];
    }
  
  Double L2_normal; 
  s=(a-b);
  t=(c-b);
  normal[0]=t[1]*s[2]-t[2]*s[1];
  normal[1]=t[2]*s[0]-t[0]*s[2];  
  normal[2]=t[0]*s[1]-t[1]*s[0];
  L2_normal=sqrt(sqr(normal[0])+sqr(normal[1])+sqr(normal[2]));
  normal[0]=normal[0]/L2_normal;
  normal[1]=normal[1]/L2_normal;
  normal[2]=normal[2]/L2_normal;  
}

char * c_string(const std::string & s)
{
   char * p = new char[s.length()+1];
   s.copy(p, std::string::npos);
   p[s.length()]=0;
   return p;
}


void SetSubVector(std::vector<Double>& mainVec, std::vector<Double>& subVec, Integer position)
{
#ifdef TRACE
  (*trace) << "entering SetSubVector" << std::endl;
#endif

  if (position + subVec.size() > mainVec.size()) 
    Error("Too less space for subvector!",__FILE__,__LINE__);

  if (!(mainVec.size() && subVec.size()) )
    Error("Vectors not initialized!",__FILE__,__LINE__);

  for(Integer i=0; i<subVec.size(); i++)
    mainVec[i+position] = subVec[i];
}

Integer defineRefinements(const Double tolElem, const Double tolTotal, const Integer noOfChilds)
{
  Double tmp=log(tolElem/tolTotal)/log(noOfChilds);
  return (Integer)tmp + 1;
}

Double CalcArea(Elem * ptE, Grid * ptgrid, const Integer level)
{
  Double         area = 0;
  BaseFE         * ptelem = ptE->ptElem;
  const Vector<Integer> & connect = ptE->connect;
  Matrix<Double> ptCoord;
  Integer        nrIntPnts = ptelem->GetNumIntPoints();
  const std::vector<Double> & intWeights = ptelem->GetIntWeights();  
  Integer        i;
  Double         jacDet;

  ptgrid->GetCoordNodesElemMat(connect,ptCoord,level);

  for (i=0; i<nrIntPnts; i++)
    {
      jacDet = ptelem->CalcJacobianDetAtIp(i+1, ptCoord);
      area +=jacDet*intWeights[i];
    } 

  return area;
}

}// namespace CoupledField
