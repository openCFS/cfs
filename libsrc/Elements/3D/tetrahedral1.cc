#include <iostream>
#include <fstream>

#include "tetrahedral1.hh"

namespace CoupledField
{

Tetrahedral1::Tetrahedral1():GeTetrahedral()
{ 
#ifdef TRACE
  (*trace) << "entering Tetrahedral1::Tetrahedral1" << std::endl;
#endif
   ElemType = TETRAHEDRAL1;
   Init();
}

Tetrahedral1::~Tetrahedral1()
{
#ifdef TRACE
  (*trace) << "entering Tetrahedral1::~Tetrahedral1" << std::endl;
#endif
 ; 
}

void Tetrahedral1::Init()
{
#ifdef TRACE
  (*trace) << "entering Tetrahedral1::Init" << std::endl;
#endif
  
  Dim = 3;
  NumNodes = 4;
  NumEdges = 6;
  NumFaces = 4;

  SetIntPoints();
  SetTransformFncAtIntPoints();
  SetDerTransformFncAtIntPoints();
  IsSet=TRUE;
}

Vector<Double> & Tetrahedral1::GetShFncAtIP(const Integer iShFnc)
{
  switch(iShFnc)
    {
    case 1:
      return TransFncAtIP1;
    case 2:
      return TransFncAtIP2;
    case 3:
      return TransFncAtIP3;
    case 4:
      return TransFncAtIP4;
    default:
    Error("Shape function does not exist with this number", __FILE__,__LINE__);
    }
}

void  Tetrahedral1::GetGradientShFnc(Vector<Double> & grad, const Integer i, const Integer ip)
{
  grad.Resize(3);

  switch(i)
{
 case 1:
    grad[0]=DxTransFncAtIP1[ip]; grad[1]=DyTransFncAtIP1[ip];  grad[2]=DzTransFncAtIP1[ip];
    break;
 case 2:
    grad[0]=DxTransFncAtIP2[ip]; grad[1]=DyTransFncAtIP2[ip];  grad[2]=DzTransFncAtIP2[ip];
    break;
 case 3:
    grad[0]=DxTransFncAtIP3[ip]; grad[1]=DyTransFncAtIP3[ip]; grad[2]=DzTransFncAtIP3[ip];
    break;
 case 4:
    grad[0]=DxTransFncAtIP4[ip]; grad[1]=DyTransFncAtIP4[ip]; grad[2]=DzTransFncAtIP4[ip];
    break;
 default:
   Error("Wrong number of shape function",__FILE__,__LINE__);
}
}
} // end of namespace

