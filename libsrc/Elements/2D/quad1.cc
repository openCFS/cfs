#include <iostream>
#include <fstream>

#include "quad1.hh"

namespace CoupledField
{

Quad1 :: Quad1() : Rectangle()
{
#ifdef TRACE
  (*trace) << "entering Quad1::Quad1" << std::endl;
#endif
  Init();
}
  
Quad1 :: ~Quad1()
{
#ifdef TRACE
  (*trace) << "entering Quad1::~Quad1" << std::endl;
#endif

  ;
}

void Quad1 :: Init()
{
#ifdef TRACE
  (*trace) << "entering Quad1::Init" << std::endl;
#endif

  Dim      = 2;
  NumNodes = 4;
  NumEdges = 4;
  NumFaces = 1;

  SetIntPoints();
  SetTransformFncAtIntPoints(); 
  SetDerTransformFncAtIntPoints();
  IsSet=TRUE;
}


void  Quad1::GetGradientShFnc(Vector<Double> & grad, const Integer i, const Integer ip)
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
 case 4:
    grad[0]=DxTransFncAtIP4[ip]; grad[1]=DyTransFncAtIP4[ip];
    break;
 default:
    Error("Number of shape function is wrong", __FILE__, __LINE__);
} 

}

Vector<Double> & Quad1::GetShFncAtIP(const Integer iShFnc) 
{
  switch(iShFnc)
    { 
    case 1:
      return TransFncAtIP1;
      break;
    case 2:
      return TransFncAtIP2;
      break;
    case 3:
      return TransFncAtIP3;
      break;
    case 4:
      return TransFncAtIP4;
      break;
    default:
      Error("Number of shape function is wrong", __FILE__, __LINE__);
    }
}
}
