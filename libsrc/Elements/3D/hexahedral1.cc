#include <iostream>
#include <fstream>

#include "hexahedral1.hh"

namespace CoupledField
{

Hexahedral1::Hexahedral1():GeHexahedral()
{ 
#ifdef TRACE
  (*trace) << "entering Hexahedral1::Hexahedral1" << std::endl;
#endif
   Init();
}

Hexahedral1::~Hexahedral1()
{
#ifdef TRACE
  (*trace) << "entering Hexahedral1::~Hexahedral1" << std::endl;
#endif
 ; 
}

void Hexahedral1::Init()
{
#ifdef TRACE
  (*trace) << "entering Hexahedral1::Init" << std::endl;
#endif
  
  Dim = 3;
  NumNodes = 8;
  NumEdges = 12;
  NumFaces = 6;

  SetIntPoints();
  SetTransformFncAtIntPoints();
  SetDerTransformFncAtIntPoints();
  IsSet=TRUE;
}

Vector<Double> & Hexahedral1::GetShFncAtIP(const Integer iShFnc)
{
  switch(iShFnc)
    {
    case 1:
      return TransFnc1AtIP;
    case 2:
      return TransFnc2AtIP;
    case 3:
      return TransFnc3AtIP;
    case 4:
      return TransFnc4AtIP;
    case 5:
      return TransFnc5AtIP;
    case 6:
      return TransFnc6AtIP;
    case 7:
      return TransFnc7AtIP;
    case 8:
      return TransFnc8AtIP;
    default:
    Error("Shape function does not exist with this number", __FILE__,__LINE__);
    }
}

void  Hexahedral1::GetGradientShFnc(Vector<Double> & grad, const Integer i, const Integer ip)
{
  grad.Resize(3);

  switch(i)
    {
    case 1:
      grad[0]=DXiTransFnc1AtIP[ip]; grad[1]=DEtTransFnc1AtIP[ip];  grad[2]=DZeTransFnc1AtIP[ip];
      break;
    case 2:
      grad[0]=DXiTransFnc2AtIP[ip]; grad[1]=DEtTransFnc2AtIP[ip];  grad[2]=DZeTransFnc2AtIP[ip];
      break;
    case 3:
      grad[0]=DXiTransFnc3AtIP[ip]; grad[1]=DEtTransFnc3AtIP[ip]; grad[2]=DZeTransFnc3AtIP[ip];
      break;
    case 4:
      grad[0]=DXiTransFnc4AtIP[ip]; grad[1]=DEtTransFnc4AtIP[ip]; grad[2]=DZeTransFnc4AtIP[ip];
      break;
    case 5:
      grad[0]=DXiTransFnc5AtIP[ip]; grad[1]=DEtTransFnc5AtIP[ip];  grad[2]=DZeTransFnc5AtIP[ip];
      break;
    case 6:
      grad[0]=DXiTransFnc6AtIP[ip]; grad[1]=DEtTransFnc6AtIP[ip];  grad[2]=DZeTransFnc6AtIP[ip];
      break;
    case 7:
      grad[0]=DXiTransFnc7AtIP[ip]; grad[1]=DEtTransFnc7AtIP[ip]; grad[2]=DZeTransFnc7AtIP[ip];
      break;
    case 8:
      grad[0]=DXiTransFnc8AtIP[ip]; grad[1]=DEtTransFnc8AtIP[ip]; grad[2]=DZeTransFnc8AtIP[ip];
      break;
    default:
      Error("Wrong number of shape function",__FILE__,__LINE__);
    }
}

void  Hexahedral1::GetGradientShFncAtCenterOfElem(Vector<Double> & grad, const Integer i)
{
  grad.Resize(3);

  switch(i)
    {
    case 1:
      grad[0]=-0.125; grad[1]=-0.125;  grad[2]=-0.125;
      break;
    case 2:
      grad[0]=0.125; grad[1]=-0.125;  grad[2]=-0.125;
      break;
    case 3:
      grad[0]=0.125; grad[1]=0.125; grad[2]=-0.125;
      break;
    case 4:
      grad[0]=-0.125; grad[1]=0.125; grad[2]=-0.125;
      break;
    case 5:
      grad[0]=-0.125; grad[1]=-0.125;  grad[2]=0.125;
      break;
    case 6:
      grad[0]=0.125; grad[1]=-0.125;  grad[2]=0.125;
      break;
    case 7:
      grad[0]=0.125; grad[1]=0.125; grad[2]=0.125;
      break;
    case 8:
      grad[0]=-0.125; grad[1]=0.125; grad[2]=0.125;
      break;
    default:
      Error("Wrong number of shape function",__FILE__,__LINE__);
    }
}
} // end of namespace

