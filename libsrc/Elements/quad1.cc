#include <iostream.h>
#include <fstream.h>
#include <math.h>

#include <general_head.hh>
#include <utils_head.hh>
#include "baseelem.hh"
#include "rectangle.hh"
#include "quad1.hh"

namespace CoupledField
{

Quad1 :: Quad1(ShortInt aintegtype) : Rectangle(aintegtype)
{
#ifdef TRACE
  (*trace) << "entering Quad1::Quad1" << endl;
#endif
  ElemType  = QUADRILATERAL1;  
  Init();
}
  
Quad1 :: ~Quad1()
{
#ifdef TRACE
  (*trace) << "entering Quad1::~Quad1" << endl;
#endif

  ;
}

void Quad1 :: Init()
{
#ifdef TRACE
  (*trace) << "entering Quad1::Init" << endl;
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
 case 2:
    grad[0]=DxTransFncAtIP2[ip]; grad[1]=DyTransFncAtIP2[ip]; 
 case 3:
    grad[0]=DxTransFncAtIP3[ip]; grad[1]=DyTransFncAtIP3[ip]; 
 case 4:
    grad[0]=DxTransFncAtIP4[ip]; grad[1]=DyTransFncAtIP4[ip];
} 

}

Vector<Double> & Quad1::GetShFncAtIP(const Integer iShFnc) 
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
    }
}
}
