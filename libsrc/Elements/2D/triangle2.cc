#include <iostream>
#include <fstream>

#include "triangle2.hh"

namespace CoupledField
{
                   
Triangle2 :: Triangle2() : GeTriangle()
{
#ifdef TRACE
  (*trace) << "entering Triangle1::Triangle2" << std::endl;
#endif
   Init();
}
  
Triangle2 :: ~Triangle2()
{
#ifdef TRACE
  (*trace) << "entering Triangle2::~Triangle2" << std::endl;
#endif

  ;
}

void Triangle2 :: Init()
{
#ifdef TRACE
  (*trace) << "entering Triangle2::Init" << std::endl;
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

void Triangle2 :: SetShapeFncAtIntPoints()
{
#ifdef TRACE
  (*trace) << "entering Triangle2::SetShapeFncAtIntPoints" << std::endl;
#endif
  Integer i;
  ShapeFncAtIP1.Resize(NumIntPoints);
  ShapeFncAtIP2.Resize(NumIntPoints);
  ShapeFncAtIP3.Resize(NumIntPoints);
  ShapeFncAtIP4.Resize(NumIntPoints);
  ShapeFncAtIP5.Resize(NumIntPoints);
  ShapeFncAtIP6.Resize(NumIntPoints);
 
  for (i=0; i < NumIntPoints; i++)
    {
      ShapeFncAtIP1[i]=ShapeFnc1(IntPoints[i][0],IntPoints[i][1]);
      ShapeFncAtIP2[i]=ShapeFnc2(IntPoints[i][0],IntPoints[i][1]);
      ShapeFncAtIP3[i]=ShapeFnc3(IntPoints[i][0],IntPoints[i][1]);
      ShapeFncAtIP4[i]=ShapeFnc4(IntPoints[i][0],IntPoints[i][1]);
      ShapeFncAtIP5[i]=ShapeFnc5(IntPoints[i][0],IntPoints[i][1]);
      ShapeFncAtIP6[i]=ShapeFnc6(IntPoints[i][0],IntPoints[i][1]);
    }
}

void Triangle2::SetDerShapeFncAtIntPoints()
{
#ifdef TRACE
  (*trace) << "entering Triangle2::SetDerShapeFncAtIntPoints" << std::endl;
#endif
  Integer i;

  DxShapeFncAtIP1.Resize(NumIntPoints);
  DxShapeFncAtIP2.Resize(NumIntPoints);
  DxShapeFncAtIP3.Resize(NumIntPoints);
  DxShapeFncAtIP4.Resize(NumIntPoints);
  DxShapeFncAtIP5.Resize(NumIntPoints);
  DxShapeFncAtIP6.Resize(NumIntPoints);
 
  for (i=0; i < NumIntPoints; i++)
    {
      DxShapeFncAtIP1[i]=DxShapeFnc1(IntPoints[i][0],IntPoints[i][1]);
      DxShapeFncAtIP2[i]=DxShapeFnc2(IntPoints[i][0],IntPoints[i][1]);
      DxShapeFncAtIP3[i]=DxShapeFnc3(IntPoints[i][0],IntPoints[i][1]);
      DxShapeFncAtIP4[i]=DxShapeFnc4(IntPoints[i][0],IntPoints[i][1]);
      DxShapeFncAtIP5[i]=DxShapeFnc5(IntPoints[i][0],IntPoints[i][1]);
      DxShapeFncAtIP6[i]=DxShapeFnc6(IntPoints[i][0],IntPoints[i][1]);
    }

  DyShapeFncAtIP1.Resize(NumIntPoints);
  DyShapeFncAtIP2.Resize(NumIntPoints);
  DyShapeFncAtIP3.Resize(NumIntPoints);
  DyShapeFncAtIP4.Resize(NumIntPoints);
  DyShapeFncAtIP5.Resize(NumIntPoints);
  DyShapeFncAtIP6.Resize(NumIntPoints);
 
  for (i=0; i < NumIntPoints; i++)
    {
      DyShapeFncAtIP1[i]=DyShapeFnc1(IntPoints[i][0],IntPoints[i][1]);
      DyShapeFncAtIP2[i]=DyShapeFnc2(IntPoints[i][0],IntPoints[i][1]);
      DyShapeFncAtIP3[i]=DyShapeFnc3(IntPoints[i][0],IntPoints[i][1]);
      DyShapeFncAtIP4[i]=DyShapeFnc4(IntPoints[i][0],IntPoints[i][1]);
      DyShapeFncAtIP5[i]=DyShapeFnc5(IntPoints[i][0],IntPoints[i][1]);
      DyShapeFncAtIP6[i]=DyShapeFnc6(IntPoints[i][0],IntPoints[i][1]);
    }
 
}

Vector<Double> & Triangle2::GetShFncAtIP(const Integer iShFnc)
{
 SetShapeFncAtIntPoints();

  switch(iShFnc)
    {
    case 1:
      return ShapeFncAtIP1;
    case 2:
      return ShapeFncAtIP2;
    case 3:
      return ShapeFncAtIP3;
    case 4:
      return ShapeFncAtIP4;
    case 5:
      return ShapeFncAtIP5;
    case 6:
      return ShapeFncAtIP6;  
    default:
    Error("Shape function does not exist with this number", __FILE__,__LINE__);
    }
}

void  Triangle2::GetGradientShFnc(Vector<Double> & grad, const Integer i, const Integer ip)
{
  SetDerShapeFncAtIntPoints();
 
  grad.Resize(2);

  switch(i)
{
 case 1:
    grad[0]=DxShapeFncAtIP1[ip]; grad[1]=DyShapeFncAtIP1[ip];
    break;
 case 2:
    grad[0]=DxShapeFncAtIP2[ip]; grad[1]=DyShapeFncAtIP2[ip];
    break;
 case 3:
    grad[0]=DxShapeFncAtIP3[ip]; grad[1]=DyShapeFncAtIP3[ip];
    break;
 case 4:
    grad[0]=DxShapeFncAtIP4[ip]; grad[1]=DyShapeFncAtIP4[ip];
    break;
 case 5:
   grad[0]=DxShapeFncAtIP5[ip]; grad[1]=DyShapeFncAtIP5[ip];
    break;
  case 6:
    grad[0]=DxShapeFncAtIP6[ip]; grad[1]=DyShapeFncAtIP6[ip];
   break; 
  default:
   Error("Wrong number of shape function",__FILE__,__LINE__);
}
}


void  Triangle2::GetGradientShFncAtCenter(Vector<Double> & grad, const Integer ish)
{
  grad.Resize(2);

  Error("Not implemented",__FILE__,__LINE__);
 //  grad[0]=DxShapeFncAtCenter[ish-1];
//   grad[1]=DyShapeFncAtCenter[ish-1];
}

} // end of namespace
