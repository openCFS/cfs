#include <iostream>
#include <fstream>

#include "quad2.hh"

namespace CoupledField
{

Quad2 :: Quad2(enum IntegrationType aintegtype) : Rectangle(aintegtype)
{
#ifdef TRACE
  (*trace) << "entering Quad2::Quad2" << std::endl;
#endif
  ElemType  = QUADRILATERAL2;  
  Init();
}
  
Quad2 :: ~Quad2()
{
#ifdef TRACE
  (*trace) << "entering Quad2::~Quad2" << std::endl;
#endif

  ;
}

void Quad2 :: Init()
{
#ifdef TRACE
  (*trace) << "entering Quad2::Init" << std::endl;
#endif

  Dim      = 2;
  NumNodes = 9;
  NumEdges = 4;
  NumFaces = 1;

  SetIntPoints();
  SetTransformFncAtIntPoints(); 
  SetDerTransformFncAtIntPoints();
  IsSet=TRUE;
}

void Quad2 :: SetShapeFncAtIntPoints()
{
#ifdef TRACE
  (*trace) << "entering Rectangle::SetShapeFncAtIntPoints" << std::endl;
#endif
  Integer i;
  ShapeFncAtIP1.Resize(NumIntPoints);
  ShapeFncAtIP2.Resize(NumIntPoints);
  ShapeFncAtIP3.Resize(NumIntPoints);
  ShapeFncAtIP4.Resize(NumIntPoints);
  ShapeFncAtIP5.Resize(NumIntPoints);
  ShapeFncAtIP6.Resize(NumIntPoints);
  ShapeFncAtIP7.Resize(NumIntPoints);
  ShapeFncAtIP8.Resize(NumIntPoints);
  ShapeFncAtIP9.Resize(NumIntPoints);
 
  for (i=0; i < NumIntPoints; i++)
    {
      ShapeFncAtIP1[i]=ShapeFnc1(IntPoints[i][0],IntPoints[i][1]);
      ShapeFncAtIP2[i]=ShapeFnc2(IntPoints[i][0],IntPoints[i][1]);
      ShapeFncAtIP3[i]=ShapeFnc3(IntPoints[i][0],IntPoints[i][1]);
      ShapeFncAtIP4[i]=ShapeFnc4(IntPoints[i][0],IntPoints[i][1]);
      ShapeFncAtIP5[i]=ShapeFnc5(IntPoints[i][0],IntPoints[i][1]);
      ShapeFncAtIP6[i]=ShapeFnc6(IntPoints[i][0],IntPoints[i][1]);
      ShapeFncAtIP7[i]=ShapeFnc7(IntPoints[i][0],IntPoints[i][1]);
      ShapeFncAtIP8[i]=ShapeFnc8(IntPoints[i][0],IntPoints[i][1]);
      ShapeFncAtIP9[i]=ShapeFnc9(IntPoints[i][0],IntPoints[i][1]);
    }
}

void Quad2 :: SetDerShapeFncAtIntPoints()
{
#ifdef TRACE
  (*trace) << "entering Quad2::SetDerShapeFncAtIntPoints" << std::endl;
#endif
  Integer i;

  DxShapeFncAtIP1.Resize(NumIntPoints);
  DxShapeFncAtIP2.Resize(NumIntPoints);
  DxShapeFncAtIP3.Resize(NumIntPoints);
  DxShapeFncAtIP4.Resize(NumIntPoints);
  DxShapeFncAtIP5.Resize(NumIntPoints);
  DxShapeFncAtIP6.Resize(NumIntPoints);
  DxShapeFncAtIP7.Resize(NumIntPoints);
  DxShapeFncAtIP8.Resize(NumIntPoints);
  DxShapeFncAtIP9.Resize(NumIntPoints);
 
  for (i=0; i < NumIntPoints; i++)
    {
      DxShapeFncAtIP1[i]=DxShapeFnc1(IntPoints[i][0],IntPoints[i][1]);
      DxShapeFncAtIP2[i]=DxShapeFnc2(IntPoints[i][0],IntPoints[i][1]);
      DxShapeFncAtIP3[i]=DxShapeFnc3(IntPoints[i][0],IntPoints[i][1]);
      DxShapeFncAtIP4[i]=DxShapeFnc4(IntPoints[i][0],IntPoints[i][1]);
      DxShapeFncAtIP5[i]=DxShapeFnc5(IntPoints[i][0],IntPoints[i][1]);
      DxShapeFncAtIP6[i]=DxShapeFnc6(IntPoints[i][0],IntPoints[i][1]);
      DxShapeFncAtIP7[i]=DxShapeFnc7(IntPoints[i][0],IntPoints[i][1]);
      DxShapeFncAtIP8[i]=DxShapeFnc8(IntPoints[i][0],IntPoints[i][1]);
      DxShapeFncAtIP9[i]=DxShapeFnc9(IntPoints[i][0],IntPoints[i][1]);
    }

  DyShapeFncAtIP1.Resize(NumIntPoints);
  DyShapeFncAtIP2.Resize(NumIntPoints);
  DyShapeFncAtIP3.Resize(NumIntPoints);
  DyShapeFncAtIP4.Resize(NumIntPoints);
  DyShapeFncAtIP5.Resize(NumIntPoints);
  DyShapeFncAtIP6.Resize(NumIntPoints);
  DyShapeFncAtIP7.Resize(NumIntPoints);
  DyShapeFncAtIP8.Resize(NumIntPoints);
  DyShapeFncAtIP9.Resize(NumIntPoints);
 
  for (i=0; i < NumIntPoints; i++)
    {
      DyShapeFncAtIP1[i]=DyShapeFnc1(IntPoints[i][0],IntPoints[i][1]);
      DyShapeFncAtIP2[i]=DyShapeFnc2(IntPoints[i][0],IntPoints[i][1]);
      DyShapeFncAtIP3[i]=DyShapeFnc3(IntPoints[i][0],IntPoints[i][1]);
      DyShapeFncAtIP4[i]=DyShapeFnc4(IntPoints[i][0],IntPoints[i][1]);
      DyShapeFncAtIP5[i]=DyShapeFnc5(IntPoints[i][0],IntPoints[i][1]);
      DyShapeFncAtIP6[i]=DyShapeFnc6(IntPoints[i][0],IntPoints[i][1]);
      DyShapeFncAtIP7[i]=DyShapeFnc7(IntPoints[i][0],IntPoints[i][1]);
      DyShapeFncAtIP8[i]=DyShapeFnc8(IntPoints[i][0],IntPoints[i][1]);
      DyShapeFncAtIP9[i]=DyShapeFnc9(IntPoints[i][0],IntPoints[i][1]);
    }
 
}

void  Quad2::GetGradientShFnc(Vector<Double> & grad, const Integer i, const Integer ip)
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
 case 7:
    grad[0]=DxShapeFncAtIP7[ip]; grad[1]=DyShapeFncAtIP7[ip];
    break;
 case 8:
    grad[0]=DxShapeFncAtIP8[ip]; grad[1]=DyShapeFncAtIP8[ip];
    break;
 case 9:
    grad[0]=DxShapeFncAtIP9[ip]; grad[1]=DyShapeFncAtIP9[ip];
    break;
} 

}

Vector<Double> & Quad2::GetShFncAtIP(const Integer iShFnc) 
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
    case 7:
      return ShapeFncAtIP7;
    case 8:
      return ShapeFncAtIP8;
    case 9:
      return ShapeFncAtIP9;
    }
}
}
