#include <iostream>
#include <fstream>

#include "triangle1.hh"

namespace CoupledField
{
                   
Triangle1 :: Triangle1(enum IntegrationType aintegtype) : GeTriangle(aintegtype)
{
#ifdef TRACE
  (*trace) << "entering Triangle1::Triangle1" << std::endl;
#endif
   ElemType  = TRIANGLE1;   
   Init();
}
  
Triangle1 :: ~Triangle1()
{
#ifdef TRACE
  (*trace) << "entering Triangle1::~Triangle1" << std::endl;
#endif

  ;
}

void Triangle1 :: Init()
{
#ifdef TRACE
  (*trace) << "entering Triangle1::Init" << std::endl;
#endif

  Dim = 2;
  NumNodes = 3;
  NumEdges = 3;
  NumFaces = 1;

  SetIntPoints();
  SetTransformFncAtIntPoints();
  SetDerTransformFncAtIntPoints();
  IsSet=TRUE;
}

Vector<Double> & Triangle1::GetShFncAtIP(const Integer iShFnc)
{
  switch(iShFnc)
    {
    case 1:
      return TransFncAtIP1;
    case 2:
      return TransFncAtIP2;
    case 3:
      return TransFncAtIP3;
    default:
    Error("Shape function does not exist with this number", __FILE__,__LINE__);
    }
}

void  Triangle1::GetGradientShFnc(Vector<Double> & grad, const Integer i, const Integer ip)
{
  grad.Resize(2);
 
  switch(i)
{
 case 1:
    grad[0]=DxTransFncAtIP1[ip]; grad[1]=DyTransFncAtIP1[ip];
    break;
 case 2:
    grad[0]=DxTransFncAtIP2[ip]; grad[1]=DyTransFncAtIP2[ip];
    break;
 case 3:
    grad[0]=DxTransFncAtIP3[ip]; grad[1]=DyTransFncAtIP3[ip];
    break;
}
}

} // end of namespace
