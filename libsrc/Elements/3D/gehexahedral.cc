#include <iostream>
#include <fstream>
#include <string>

#include "conffile.hh"
#include "gehexahedral.hh"

namespace CoupledField
{

GeHexahedral::GeHexahedral()
{

  std::string integtype;
  conf->get("hexahedra",integtype,"IntegRules");

  IntegType=String2EnumIntegrationType(integtype.c_str());

}

GeHexahedral::~GeHexahedral()
{
  if (IntWeights) delete IntWeights;
}

void GeHexahedral::SetIntPoints()
{
#ifdef TRACE
  (*trace) << "entering GeHexahedral::SetIntPoints" << std::endl;
#endif
 
  switch(IntegType) 
    {
    case GaussOrder2:
	  
      NumIntPoints=8;
      DegreeInteg=3;
      
      IntPoints.Resize(NumIntPoints, Dim);
      IntWeights=0;

      IntPoints[0][0] = -0.57735026919;
      IntPoints[0][1] = -0.57735026919;
      IntPoints[0][2] = -0.57735026919;
      
      IntPoints[1][0] = 0.57735026919;
      IntPoints[1][1] = -0.57735026919;
      IntPoints[1][2] = -0.57735026919;

      IntPoints[2][0] = 0.57735026919;
      IntPoints[2][1] = 0.57735026919;
      IntPoints[2][2] = -0.57735026919;

      IntPoints[3][0] = -0.57735026919;
      IntPoints[3][1] = 0.57735026919;
      IntPoints[3][2] = -0.57735026919;

      IntPoints[4][0] = -0.57735026919;
      IntPoints[4][1] = -0.57735026919;
      IntPoints[4][2] = 0.57735026919;
      
      IntPoints[5][0] = 0.57735026919;
      IntPoints[5][1] = -0.57735026919;
      IntPoints[5][2] = 0.57735026919;

      IntPoints[6][0] = 0.57735026919;
      IntPoints[6][1] = 0.57735026919;
      IntPoints[6][2] = 0.57735026919;

      IntPoints[7][0] = -0.57735026919;
      IntPoints[7][1] = 0.57735026919;
      IntPoints[7][2] = 0.57735026919;

      break;

   default:
      Error("Integration type is not implemented",__FILE__,__LINE__);
   }
}

void GeHexahedral::SetTransformFncAtIntPoints()
{
#ifdef TRACE
  (*trace) << "entering GeHexahedral::SetTransformFncAtIntPoints" << std::endl;
#endif

 Integer i;
 TransFnc1AtIP.Resize(NumIntPoints);
 TransFnc2AtIP.Resize(NumIntPoints);
 TransFnc3AtIP.Resize(NumIntPoints);
 TransFnc4AtIP.Resize(NumIntPoints);
 TransFnc5AtIP.Resize(NumIntPoints);
 TransFnc6AtIP.Resize(NumIntPoints);
 TransFnc7AtIP.Resize(NumIntPoints);
 TransFnc8AtIP.Resize(NumIntPoints);


  for (i=0; i < NumIntPoints; i++)
    {
      TransFnc1AtIP[i]=TransFnc1(IntPoints[i][0],IntPoints[i][1],IntPoints[i][2]);
      TransFnc2AtIP[i]=TransFnc2(IntPoints[i][0],IntPoints[i][1],IntPoints[i][2]);
      TransFnc3AtIP[i]=TransFnc3(IntPoints[i][0],IntPoints[i][1],IntPoints[i][2]);
      TransFnc4AtIP[i]=TransFnc4(IntPoints[i][0],IntPoints[i][1],IntPoints[i][2]);
      TransFnc5AtIP[i]=TransFnc5(IntPoints[i][0],IntPoints[i][1],IntPoints[i][2]);
      TransFnc6AtIP[i]=TransFnc6(IntPoints[i][0],IntPoints[i][1],IntPoints[i][2]);
      TransFnc7AtIP[i]=TransFnc7(IntPoints[i][0],IntPoints[i][1],IntPoints[i][2]);
      TransFnc8AtIP[i]=TransFnc8(IntPoints[i][0],IntPoints[i][1],IntPoints[i][2]);
    }
}

void GeHexahedral::SetDerTransformFncAtIntPoints()
{
#ifdef TRACE
  (*trace) << "entering GeHexahedral::SetDerTransformFncAtIntPoints" << std::endl;
#endif

  Integer i;

  DXiTransFnc1AtIP.Resize(NumIntPoints);
  DXiTransFnc2AtIP.Resize(NumIntPoints);
  DXiTransFnc3AtIP.Resize(NumIntPoints);
  DXiTransFnc4AtIP.Resize(NumIntPoints);
  DXiTransFnc5AtIP.Resize(NumIntPoints);
  DXiTransFnc6AtIP.Resize(NumIntPoints);
  DXiTransFnc7AtIP.Resize(NumIntPoints);
  DXiTransFnc8AtIP.Resize(NumIntPoints);

  for (i=0; i < NumIntPoints; i++)
    {
      DXiTransFnc1AtIP[i]=TransFnc1dxi(IntPoints[i][0],IntPoints[i][1],IntPoints[i][2]);
      DXiTransFnc2AtIP[i]=TransFnc2dxi(IntPoints[i][0],IntPoints[i][1],IntPoints[i][2]);
      DXiTransFnc3AtIP[i]=TransFnc3dxi(IntPoints[i][0],IntPoints[i][1],IntPoints[i][2]);
      DXiTransFnc4AtIP[i]=TransFnc4dxi(IntPoints[i][0],IntPoints[i][1],IntPoints[i][2]);
      DXiTransFnc5AtIP[i]=TransFnc5dxi(IntPoints[i][0],IntPoints[i][1],IntPoints[i][2]);
      DXiTransFnc6AtIP[i]=TransFnc6dxi(IntPoints[i][0],IntPoints[i][1],IntPoints[i][2]);
      DXiTransFnc7AtIP[i]=TransFnc7dxi(IntPoints[i][0],IntPoints[i][1],IntPoints[i][2]);
      DXiTransFnc8AtIP[i]=TransFnc8dxi(IntPoints[i][0],IntPoints[i][1],IntPoints[i][2]);
    }

  DEtTransFnc1AtIP.Resize(NumIntPoints);
  DEtTransFnc2AtIP.Resize(NumIntPoints);
  DEtTransFnc3AtIP.Resize(NumIntPoints);
  DEtTransFnc4AtIP.Resize(NumIntPoints);
  DEtTransFnc5AtIP.Resize(NumIntPoints);
  DEtTransFnc6AtIP.Resize(NumIntPoints);
  DEtTransFnc7AtIP.Resize(NumIntPoints);
  DEtTransFnc8AtIP.Resize(NumIntPoints);


  for (i=0; i < NumIntPoints; i++)
    {
      DEtTransFnc1AtIP[i]=TransFnc1deta(IntPoints[i][0],IntPoints[i][1],IntPoints[i][2]);
      DEtTransFnc2AtIP[i]=TransFnc2deta(IntPoints[i][0],IntPoints[i][1],IntPoints[i][2]);
      DEtTransFnc3AtIP[i]=TransFnc3deta(IntPoints[i][0],IntPoints[i][1],IntPoints[i][2]);
      DEtTransFnc4AtIP[i]=TransFnc4deta(IntPoints[i][0],IntPoints[i][1],IntPoints[i][2]);
      DEtTransFnc5AtIP[i]=TransFnc5deta(IntPoints[i][0],IntPoints[i][1],IntPoints[i][2]);
      DEtTransFnc6AtIP[i]=TransFnc6deta(IntPoints[i][0],IntPoints[i][1],IntPoints[i][2]);
      DEtTransFnc7AtIP[i]=TransFnc7deta(IntPoints[i][0],IntPoints[i][1],IntPoints[i][2]);
      DEtTransFnc8AtIP[i]=TransFnc8deta(IntPoints[i][0],IntPoints[i][1],IntPoints[i][2]);
    }

  DZeTransFnc1AtIP.Resize(NumIntPoints);
  DZeTransFnc2AtIP.Resize(NumIntPoints);
  DZeTransFnc3AtIP.Resize(NumIntPoints);
  DZeTransFnc4AtIP.Resize(NumIntPoints);
  DZeTransFnc5AtIP.Resize(NumIntPoints);
  DZeTransFnc6AtIP.Resize(NumIntPoints);
  DZeTransFnc7AtIP.Resize(NumIntPoints);
  DZeTransFnc8AtIP.Resize(NumIntPoints);

    for (i=0; i < NumIntPoints; i++)
    {
      DZeTransFnc1AtIP[i]=TransFnc1dzeta(IntPoints[i][0],IntPoints[i][1],IntPoints[i][2]);
      DZeTransFnc2AtIP[i]=TransFnc2dzeta(IntPoints[i][0],IntPoints[i][1],IntPoints[i][2]);
      DZeTransFnc3AtIP[i]=TransFnc3dzeta(IntPoints[i][0],IntPoints[i][1],IntPoints[i][2]);
      DZeTransFnc4AtIP[i]=TransFnc4dzeta(IntPoints[i][0],IntPoints[i][1],IntPoints[i][2]);
      DZeTransFnc5AtIP[i]=TransFnc5dzeta(IntPoints[i][0],IntPoints[i][1],IntPoints[i][2]);
      DZeTransFnc6AtIP[i]=TransFnc6dzeta(IntPoints[i][0],IntPoints[i][1],IntPoints[i][2]);
      DZeTransFnc7AtIP[i]=TransFnc7dzeta(IntPoints[i][0],IntPoints[i][1],IntPoints[i][2]);
      DZeTransFnc8AtIP[i]=TransFnc8dzeta(IntPoints[i][0],IntPoints[i][1],IntPoints[i][2]);
    }
}

void GeHexahedral::CalcJacobian(Jacobian<3> & J, const Integer ip,
                     Point<3> * ptCoord, const Boolean NeedJinv)
{

  if (!IsSet)
 {
 SetTransformFncAtIntPoints();
 SetDerTransformFncAtIntPoints();
 IsSet=TRUE;
 }

 Double aux=0;

J.J[0][0] = DXiTransFnc1AtIP[ip]*ptCoord[0][0] + DXiTransFnc2AtIP[ip]*ptCoord[1][0]
          + DXiTransFnc3AtIP[ip]*ptCoord[2][0] + DXiTransFnc4AtIP[ip]*ptCoord[3][0]
          + DXiTransFnc5AtIP[ip]*ptCoord[4][0] + DXiTransFnc6AtIP[ip]*ptCoord[5][0]
          + DXiTransFnc7AtIP[ip]*ptCoord[6][0] + DXiTransFnc8AtIP[ip]*ptCoord[7][0];

J.J[0][1] = DEtTransFnc1AtIP[ip]*ptCoord[0][0] + DEtTransFnc2AtIP[ip]*ptCoord[1][0]
          + DEtTransFnc3AtIP[ip]*ptCoord[2][0] + DEtTransFnc4AtIP[ip]*ptCoord[3][0]
          + DEtTransFnc5AtIP[ip]*ptCoord[4][0] + DEtTransFnc6AtIP[ip]*ptCoord[5][0]
          + DEtTransFnc7AtIP[ip]*ptCoord[6][0] + DEtTransFnc8AtIP[ip]*ptCoord[7][0]; 

J.J[0][2] = DZeTransFnc1AtIP[ip]*ptCoord[0][0] + DZeTransFnc2AtIP[ip]*ptCoord[1][0]
          + DZeTransFnc3AtIP[ip]*ptCoord[2][0] + DZeTransFnc4AtIP[ip]*ptCoord[3][0]
          + DZeTransFnc5AtIP[ip]*ptCoord[4][0] + DZeTransFnc6AtIP[ip]*ptCoord[5][0]
          + DZeTransFnc7AtIP[ip]*ptCoord[6][0] + DZeTransFnc8AtIP[ip]*ptCoord[7][0];

J.J[1][0] = DXiTransFnc1AtIP[ip]*ptCoord[0][1] + DXiTransFnc2AtIP[ip]*ptCoord[1][1]
          + DXiTransFnc3AtIP[ip]*ptCoord[2][1] + DXiTransFnc4AtIP[ip]*ptCoord[3][1]
          + DXiTransFnc5AtIP[ip]*ptCoord[4][1] + DXiTransFnc6AtIP[ip]*ptCoord[5][1]
          + DXiTransFnc7AtIP[ip]*ptCoord[6][1] + DXiTransFnc8AtIP[ip]*ptCoord[7][1];

J.J[1][1] = DEtTransFnc1AtIP[ip]*ptCoord[0][1] + DEtTransFnc2AtIP[ip]*ptCoord[1][1]
          + DEtTransFnc3AtIP[ip]*ptCoord[2][1] + DEtTransFnc4AtIP[ip]*ptCoord[3][1]
          + DEtTransFnc5AtIP[ip]*ptCoord[4][1] + DEtTransFnc6AtIP[ip]*ptCoord[5][1]
          + DEtTransFnc7AtIP[ip]*ptCoord[6][1] + DEtTransFnc8AtIP[ip]*ptCoord[7][1];

J.J[1][2] = DZeTransFnc1AtIP[ip]*ptCoord[0][1] + DZeTransFnc2AtIP[ip]*ptCoord[1][1]
          + DZeTransFnc3AtIP[ip]*ptCoord[2][1] + DZeTransFnc4AtIP[ip]*ptCoord[3][1]
          + DZeTransFnc5AtIP[ip]*ptCoord[4][1] + DZeTransFnc6AtIP[ip]*ptCoord[5][1]
          + DZeTransFnc7AtIP[ip]*ptCoord[6][1] + DZeTransFnc8AtIP[ip]*ptCoord[7][1];

J.J[2][0] = DXiTransFnc1AtIP[ip]*ptCoord[0][2] + DXiTransFnc2AtIP[ip]*ptCoord[1][2]
          + DXiTransFnc3AtIP[ip]*ptCoord[2][2] + DXiTransFnc4AtIP[ip]*ptCoord[3][2]
          + DXiTransFnc5AtIP[ip]*ptCoord[4][2] + DXiTransFnc6AtIP[ip]*ptCoord[5][2]
          + DXiTransFnc7AtIP[ip]*ptCoord[6][2] + DXiTransFnc8AtIP[ip]*ptCoord[7][2];

J.J[2][1] = DEtTransFnc1AtIP[ip]*ptCoord[0][2] + DEtTransFnc2AtIP[ip]*ptCoord[1][2]
          + DEtTransFnc3AtIP[ip]*ptCoord[2][2] + DEtTransFnc4AtIP[ip]*ptCoord[3][2]
          + DEtTransFnc5AtIP[ip]*ptCoord[4][2] + DEtTransFnc6AtIP[ip]*ptCoord[5][2]
          + DEtTransFnc7AtIP[ip]*ptCoord[6][2] + DEtTransFnc8AtIP[ip]*ptCoord[7][2];

J.J[2][2] = DZeTransFnc1AtIP[ip]*ptCoord[0][2] + DZeTransFnc2AtIP[ip]*ptCoord[1][2]
          + DZeTransFnc3AtIP[ip]*ptCoord[2][2] + DZeTransFnc4AtIP[ip]*ptCoord[3][2]
          + DZeTransFnc5AtIP[ip]*ptCoord[4][2] + DZeTransFnc6AtIP[ip]*ptCoord[5][2]
          + DZeTransFnc7AtIP[ip]*ptCoord[6][2] + DZeTransFnc8AtIP[ip]*ptCoord[7][2];

J.detJ = J.J[0][0]*J.J[1][1]*J.J[2][2]-J.J[0][0]*J.J[1][2]*J.J[2][1]-J.J[0][1]*J.J[1][0]*J.J[2][2]+J.J[0][1]*J.J[1][2]*J.J[2][0]+J.J[0][2]*J.J[1][0]*J.J[2][1]-J.J[0][2]*J.J[1][1]*J.J[2][0];

if (NeedJinv)
{
 aux=1.0/J.detJ;

 J.Jinv[0][0] = J.J[1][1]*J.J[2][2]-J.J[1][2]*J.J[2][1];

 J.Jinv[0][1] = J.J[0][2]*J.J[2][1]- J.J[0][1]*J.J[2][2];

 J.Jinv[0][2] = J.J[0][1]*J.J[1][2]- J.J[1][1]*J.J[0][2];

 J.Jinv[1][0] = J.J[1][2]*J.J[2][0]-J.J[1][0]*J.J[2][2];

 J.Jinv[1][1] = J.J[0][0]*J.J[2][2]-J.J[0][2]*J.J[2][0];

 J.Jinv[1][2] = J.J[0][2]*J.J[1][0]-J.J[0][0]*J.J[1][2];

 J.Jinv[2][0] = J.J[1][0]*J.J[2][1]-J.J[1][1]*J.J[2][0];

 J.Jinv[2][1] = J.J[0][1]*J.J[2][0]-J.J[0][0]*J.J[2][1];

 J.Jinv[2][2] = J.J[0][0]*J.J[1][1]-J.J[0][1]*J.J[1][0];

 J.Jinv*=aux;
}

}

void GeHexahedral::CalcJacobianAtCenterOfElem(Jacobian<3> & J,
                     Point<3> * ptCoord, const Boolean NeedJinv)
{

//   if (!IsSet)
//     {
//       SetTransformFncAtIntPoints();
//       SetDerTransformFncAtIntPoints();
//       IsSet=TRUE;
//     }

 Double aux=0;

J.J[0][0] = 0.125*(-ptCoord[0][0] + ptCoord[1][0]
          + ptCoord[2][0] - ptCoord[3][0]
          - ptCoord[4][0] + ptCoord[5][0]
          + ptCoord[6][0] - ptCoord[7][0]);

J.J[0][1] = 0.125*(-ptCoord[0][0] - ptCoord[1][0]
          + ptCoord[2][0] + ptCoord[3][0]
          - ptCoord[4][0] - ptCoord[5][0]
          + ptCoord[6][0] + ptCoord[7][0]); 

J.J[0][2] = 0.125*(-ptCoord[0][0] - ptCoord[1][0]
          - ptCoord[2][0] - ptCoord[3][0]
          + ptCoord[4][0] + ptCoord[5][0]
          + ptCoord[6][0] + ptCoord[7][0]);

J.J[1][0] = 0.125*(-ptCoord[0][1] + ptCoord[1][1]
          + ptCoord[2][1] - ptCoord[3][1]
          - ptCoord[4][1] + ptCoord[5][1]
          + ptCoord[6][1] - ptCoord[7][1]);

J.J[1][1] = 0.125*(-ptCoord[0][1] - ptCoord[1][1]
          + ptCoord[2][1] + ptCoord[3][1]
          - ptCoord[4][1] - ptCoord[5][1]
          + ptCoord[6][1] + ptCoord[7][1]);

J.J[1][2] = 0.125*(-ptCoord[0][1] - ptCoord[1][1]
          - ptCoord[2][1] - ptCoord[3][1]
          + ptCoord[4][1] + ptCoord[5][1]
          + ptCoord[6][1] + ptCoord[7][1]);

J.J[2][0] = 0.125*(-ptCoord[0][2] + ptCoord[1][2]
          + ptCoord[2][2] - ptCoord[3][2]
          - ptCoord[4][2] + ptCoord[5][2]
          + ptCoord[6][2] - ptCoord[7][2]);

J.J[2][1] = 0.125*(-ptCoord[0][2] - ptCoord[1][2]
          + ptCoord[2][2] + ptCoord[3][2]
          - ptCoord[4][2] - ptCoord[5][2]
          + ptCoord[6][2] + ptCoord[7][2]); 

J.J[2][2] = 0.125*(-ptCoord[0][2] - ptCoord[1][2]
          - ptCoord[2][2] - ptCoord[3][2]
          + ptCoord[4][2] + ptCoord[5][2]
          + ptCoord[6][2] + ptCoord[7][2]);

J.detJ = J.J[0][0]*J.J[1][1]*J.J[2][2]-J.J[0][0]*J.J[1][2]*J.J[2][1]-J.J[0][1]*J.J[1][0]*J.J[2][2]+J.J[0][1]*J.J[1][2]*J.J[2][0]+J.J[0][2]*J.J[1][0]*J.J[2][1]-J.J[0][2]*J.J[1][1]*J.J[2][0];

if (NeedJinv)
{
 aux=1.0/J.detJ;

 J.Jinv[0][0] = J.J[1][1]*J.J[2][2]-J.J[1][2]*J.J[2][1];

 J.Jinv[0][1] = J.J[0][2]*J.J[2][1]- J.J[0][1]*J.J[2][2];

 J.Jinv[0][2] = J.J[0][1]*J.J[1][2]- J.J[1][1]*J.J[0][2];

 J.Jinv[1][0] = J.J[1][2]*J.J[2][0]-J.J[1][0]*J.J[2][2];

 J.Jinv[1][1] = J.J[0][0]*J.J[2][2]-J.J[0][2]*J.J[2][0];

 J.Jinv[1][2] = J.J[0][2]*J.J[1][0]-J.J[0][0]*J.J[1][2];

 J.Jinv[2][0] = J.J[1][0]*J.J[2][1]-J.J[1][1]*J.J[2][0];

 J.Jinv[2][1] = J.J[0][1]*J.J[2][0]-J.J[0][0]*J.J[2][1];

 J.Jinv[2][2] = J.J[0][0]*J.J[1][1]-J.J[0][1]*J.J[1][0];

 J.Jinv*=aux;
}

}

void GeHexahedral::CalcJacobian(Jacobian<2> & J, const Integer ip,
                   Point<2> *  ptCoord, const Boolean NeedJinv)
{
 Error("This element is from 3D", __FILE__, __LINE__);
}

///////////////////////////////////////////////////////////////////////////

Vector<Double> *  GeHexahedral::GetDXiShFncAtIP(const Integer iShFnc)
{
  switch(iShFnc)
{ case 1:
     return &  DXiTransFnc1AtIP;
  case 2:
     return &  DXiTransFnc2AtIP;
  case 3:
     return &  DXiTransFnc3AtIP;
  case 4:
     return &  DXiTransFnc4AtIP;
  case 5:
     return &  DXiTransFnc5AtIP;
  case 6:
     return &  DXiTransFnc6AtIP;
  case 7:
     return &  DXiTransFnc7AtIP;
  case 8:
     return &  DXiTransFnc8AtIP;
}
}

Vector<Double> *  GeHexahedral::GetDEtShFncAtIP(const Integer iShFnc)
{
  switch(iShFnc)
{ case 1:
     return &  DEtTransFnc1AtIP;
  case 2:
     return &  DEtTransFnc2AtIP;
  case 3:
     return &  DEtTransFnc3AtIP;
  case 4:
     return &  DEtTransFnc4AtIP;
  case 5:
     return &  DEtTransFnc5AtIP;
  case 6:
     return &  DEtTransFnc6AtIP;
  case 7:
     return &  DEtTransFnc7AtIP;
  case 8:
     return &  DEtTransFnc8AtIP;
}
}

Vector<Double> *  GeHexahedral::GetDZeShFncAtIP(const Integer iShFnc)
{
  switch(iShFnc)
{ case 1:
     return &  DZeTransFnc1AtIP;
  case 2:
     return &  DZeTransFnc2AtIP;
  case 3:
     return &  DZeTransFnc3AtIP;
  case 4:
     return &  DZeTransFnc4AtIP;
  case 5:
     return &  DZeTransFnc5AtIP;
  case 6:
     return &  DZeTransFnc6AtIP;
  case 7:
     return &  DZeTransFnc7AtIP;
  case 8:
     return &  DZeTransFnc8AtIP;
}
}

} // end of namespace

