#include <iostream>
#include <fstream>
#include <string>

#include "conffile.hh"
#include "getetrahedral.hh"

namespace CoupledField
{

GeTetrahedral::GeTetrahedral()
{

  std::string integtype;
  conf->get("tetrahedra",integtype,"IntegRules");

  IntegType=String2EnumIntegrationType(integtype.c_str());

}

GeTetrahedral::~GeTetrahedral()
{
  if (IntWeights) delete IntWeights;
}

void GeTetrahedral::SetIntPoints()
{
#ifdef TRACE
  (*trace) << "entering GeTetrahedral::SetIntPoints" << std::endl;
#endif

  switch(IntegType)
   {
   case GaussOrder2:
    /// Thomas Hughes "The finite element method", p. 174
    
      NumIntPoints=1;
      DegreeInteg=2;
      IntPoints.Resize(NumIntPoints, Dim);
      IntWeights=NULL;

      IntPoints[0][0]=0.25;
      IntPoints[0][1]=0.25;
      IntPoints[0][2]=0.25;

      if (InfoPrint)
       (*infofile) << " For numerical integration procedures we use Gaussian Quadrature with 1 node, degree of precision is 2 " << std::endl;

      break;

   case GaussOrder3:
  /// Thomas Hughes "The finite element method", p. 174   

      NumIntPoints=4;
      DegreeInteg=3;
      IntPoints.Resize(NumIntPoints, Dim);
      IntWeights=new Vector<Double>(NumIntPoints);

      IntPoints[0][0]=0.5854102;
      IntPoints[0][1]=0.1381966;
      IntPoints[0][2]=0.1381966;

      IntPoints[1][1]=0.5854102;
      IntPoints[1][2]=0.1381966;
      IntPoints[1][0]=0.1381966;

      IntPoints[2][2]=0.5854102;
      IntPoints[2][0]=0.1381966;
      IntPoints[2][1]=0.1381966;

      IntPoints[3][0]=0.1381966;
      IntPoints[3][1]=0.1381966;
      IntPoints[3][2]=0.1381966;

      (*IntWeights)[0]=0.25;
      (*IntWeights)[1]=0.25;
      (*IntWeights)[2]=0.25;
      (*IntWeights)[3]=0.25;

      if (InfoPrint)
       (*infofile) << " For numerical integration procedures we use Gaussian Quadrature with 4 nodes, degree of precision is 3 " << std::endl;

      break;

      case GaussOrder4:
  /// Thomas Hughes "The finite element method", p. 174

      NumIntPoints=5;
      DegreeInteg=4;
      IntPoints.Resize(NumIntPoints, Dim);
      IntWeights=new Vector<Double>(NumIntPoints);

      IntPoints[0][0]=0.25;
      IntPoints[0][1]=0.25;
      IntPoints[0][2]=0.25;

      IntPoints[1][0]=0.3333333;
      IntPoints[1][1]=0.1666667; 
      IntPoints[1][2]=0.1666667; 

      IntPoints[2][1]=0.3333333;
      IntPoints[2][2]=0.1666667;
      IntPoints[2][0]=0.1666667;

      IntPoints[3][2]=0.3333333;
      IntPoints[3][0]=0.1666667;
      IntPoints[3][1]=0.1666667;

      IntPoints[4][0]=0.1666667;
      IntPoints[4][1]=0.1666667;
      IntPoints[4][2]=0.1666667;

      (*IntWeights)[0]=-0.8; 
      (*IntWeights)[1]=0.45;
      (*IntWeights)[2]=0.45;
      (*IntWeights)[3]=0.45;
      (*IntWeights)[4]=0.45;

      if (InfoPrint)
      (*infofile) << " For numerical integration procedures we use Gaussian Quadrature with 5 nodes, degree of precision is 4 " << std::endl;

      break;

      case GaussOrder5:
  /// A.H.Stroud "Approximate Calculation of Multiple Integrals"
  /// Printice Hall 1971

       NumIntPoints=15;
       DegreeInteg=5;
       IntPoints.Resize(NumIntPoints, Dim);
       IntWeights=new Vector<Double>(NumIntPoints);

   IntPoints[0][0]=0.25;  IntPoints[0][1]=0.25;  IntPoints[0][2]=0.25;

   IntPoints[1][0]=0.09197107805272303;  IntPoints[1][1]=0.09197107805272303; IntPoints[1][2]=0.09197107805272303;

  IntPoints[2][0]=0.72408676584183096; IntPoints[2][1]=0.09197107805272303; IntPoints[2][2]=0.09197107805272303;

  IntPoints[3][0]=0.09197107805272303; IntPoints[3][1]=0.72408676584183096; IntPoints[3][2]=0.09197107805272303;   

  IntPoints[4][0]=0.09197107805272303; IntPoints[4][1]=0.09197107805272303; IntPoints[4][2]=0.72408676584183096;

  IntPoints[5][0]=0.44364916731037080; IntPoints[5][1]=0.05635083268962915; IntPoints[5][2]=0.05635083268962915;

  IntPoints[6][0]=0.05635083268962915; IntPoints[6][1]=0.44364916731037080; IntPoints[6][2]=0.05635083268962915;

  IntPoints[7][0]=0.05635083268962915; IntPoints[7][1]=0.05635083268962915; IntPoints[7][2]=0.44364916731037080;

  IntPoints[8][0]=0.05635083268962915; IntPoints[8][1]=0.44364916731037080; IntPoints[8][2]=0.44364916731037080;

  IntPoints[9][0]=0.44364916731037080; IntPoints[9][1]=0.05635083268962915; IntPoints[9][2]=0.44364916731037080;

  IntPoints[10][0]=0.44364916731037080; IntPoints[10][1]=0.44364916731037080; IntPoints[10][2]=0.05635083268962915;

  IntPoints[11][0]=0.31979362782962989; IntPoints[11][1]=0.31979362782962989; IntPoints[11][2]=0.31979362782962989;

  IntPoints[12][0]=0.04061911651111023; IntPoints[12][1]=0.31979362782962989; IntPoints[12][2]=0.31979362782962989;

  IntPoints[13][0]=0.31979362782962989; IntPoints[13][1]=0.04061911651111023; IntPoints[13][2]=0.31979362782962989;

  IntPoints[14][0]=0.31979362782962989; IntPoints[14][1]=0.31979362782962989; IntPoints[14][2]=0.04061911651111023;

 
      (*IntWeights)[0]=0.019753086419753086;
      (*IntWeights)[1]=0.011989513963169772;
      (*IntWeights)[2]=0.011989513963169772;
      (*IntWeights)[3]=0.011989513963169772;
      (*IntWeights)[4]=0.011989513963169772;
      (*IntWeights)[5]=0.008818342151675485;
      (*IntWeights)[6]=0.008818342151675485;
      (*IntWeights)[7]=0.008818342151675485;
      (*IntWeights)[8]=0.008818342151675485;
      (*IntWeights)[9]=0.008818342151675485;
      (*IntWeights)[10]=0.00881834215167548;
      (*IntWeights)[11]=0.011511367871045397;
      (*IntWeights)[12]=0.011511367871045397;
      (*IntWeights)[13]=0.011511367871045397;
      (*IntWeights)[14]=0.011511367871045397;
 
      break;

   default:
      Error("Integration type is not implemented",__FILE__,__LINE__);
   }
}

void GeTetrahedral::SetTransformFncAtIntPoints()
{
#ifdef TRACE
  (*trace) << "entering GeTetrahedral::SetTransformFncAtIntPoints" << std::endl;
#endif

 Integer i;
 TransFncAtIP1.Resize(NumIntPoints);
 TransFncAtIP2.Resize(NumIntPoints);
 TransFncAtIP3.Resize(NumIntPoints);
 TransFncAtIP4.Resize(NumIntPoints);

  for (i=0; i < NumIntPoints; i++)
    {
      TransFncAtIP1[i]=TransFnc1(IntPoints[i][0],IntPoints[i][1],IntPoints[i][2]);
      TransFncAtIP2[i]=TransFnc2(IntPoints[i][0],IntPoints[i][1],IntPoints[i][2]);
      TransFncAtIP3[i]=TransFnc3(IntPoints[i][0],IntPoints[i][1],IntPoints[i][2]);
      TransFncAtIP4[i]=TransFnc4(IntPoints[i][0],IntPoints[i][1],IntPoints[i][2]);
    }
}

void GeTetrahedral::SetDerTransformFncAtIntPoints()
{
#ifdef TRACE
  (*trace) << "entering GeTetrahedral::SetDerTransformFncAtIntPoints" << std::endl;
#endif

  Integer i;

  DxTransFncAtIP1.Resize(NumIntPoints);
  DxTransFncAtIP2.Resize(NumIntPoints);
  DxTransFncAtIP3.Resize(NumIntPoints);
  DxTransFncAtIP4.Resize(NumIntPoints);

  for (i=0; i < NumIntPoints; i++)
    {
      DxTransFncAtIP1[i]=TransFnc1dx(IntPoints[i][0],IntPoints[i][1],IntPoints[i][2]);
      DxTransFncAtIP2[i]=TransFnc2dx(IntPoints[i][0],IntPoints[i][1],IntPoints[i][2]);
      DxTransFncAtIP3[i]=TransFnc3dx(IntPoints[i][0],IntPoints[i][1],IntPoints[i][2]);
      DxTransFncAtIP4[i]=TransFnc4dx(IntPoints[i][0],IntPoints[i][1],IntPoints[i][2]);
    }

  DyTransFncAtIP1.Resize(NumIntPoints);
  DyTransFncAtIP2.Resize(NumIntPoints);
  DyTransFncAtIP3.Resize(NumIntPoints);
  DyTransFncAtIP4.Resize(NumIntPoints);

  for (i=0; i < NumIntPoints; i++)
    {
      DyTransFncAtIP1[i]=TransFnc1dy(IntPoints[i][0],IntPoints[i][1],IntPoints[i][2]);
      DyTransFncAtIP2[i]=TransFnc2dy(IntPoints[i][0],IntPoints[i][1],IntPoints[i][2]);
      DyTransFncAtIP3[i]=TransFnc3dy(IntPoints[i][0],IntPoints[i][1],IntPoints[i][2]);
      DyTransFncAtIP4[i]=TransFnc4dy(IntPoints[i][0],IntPoints[i][1],IntPoints[i][2]);
    }

  DzTransFncAtIP1.Resize(NumIntPoints);
  DzTransFncAtIP2.Resize(NumIntPoints);
  DzTransFncAtIP3.Resize(NumIntPoints);
  DzTransFncAtIP4.Resize(NumIntPoints);

    for (i=0; i < NumIntPoints; i++)
    {
      DzTransFncAtIP1[i]=TransFnc1dz(IntPoints[i][0],IntPoints[i][1],IntPoints[i][2]);
      DzTransFncAtIP2[i]=TransFnc2dz(IntPoints[i][0],IntPoints[i][1],IntPoints[i][2]);
      DzTransFncAtIP3[i]=TransFnc3dz(IntPoints[i][0],IntPoints[i][1],IntPoints[i][2]);
      DzTransFncAtIP4[i]=TransFnc4dz(IntPoints[i][0],IntPoints[i][1],IntPoints[i][2]);
    }
}

void GeTetrahedral::CalcJacobian(Jacobian<Point3D> & J, const Integer ip,
                     const Point3D * const ptCoord, const Boolean NeedJinv)
{

  if (!IsSet)
 {
 SetTransformFncAtIntPoints();
 SetDerTransformFncAtIntPoints();
 IsSet=TRUE;
 }

 Double aux=0;

J.J[0][0] = DxTransFncAtIP1[ip]*ptCoord[0].x + DxTransFncAtIP2[ip]*ptCoord[1].x
          + DxTransFncAtIP3[ip]*ptCoord[2].x + DxTransFncAtIP4[ip]*ptCoord[3].x;

J.J[0][1] = DyTransFncAtIP1[ip]*ptCoord[0].x + DyTransFncAtIP2[ip]*ptCoord[1].x
          + DyTransFncAtIP3[ip]*ptCoord[2].x + DyTransFncAtIP4[ip]*ptCoord[3].x;

J.J[0][2] = DzTransFncAtIP1[ip]*ptCoord[0].x + DzTransFncAtIP2[ip]*ptCoord[1].x
          + DzTransFncAtIP3[ip]*ptCoord[2].x + DzTransFncAtIP4[ip]*ptCoord[3].x;

J.J[1][0] = DxTransFncAtIP1[ip]*ptCoord[0].y + DxTransFncAtIP2[ip]*ptCoord[1].y
          + DxTransFncAtIP3[ip]*ptCoord[2].y + DxTransFncAtIP4[ip]*ptCoord[3].y;

J.J[1][1] = DyTransFncAtIP1[ip]*ptCoord[0].y + DyTransFncAtIP2[ip]*ptCoord[1].y
          + DyTransFncAtIP3[ip]*ptCoord[2].y + DyTransFncAtIP4[ip]*ptCoord[3].y;

J.J[1][2] = DzTransFncAtIP1[ip]*ptCoord[0].y + DzTransFncAtIP2[ip]*ptCoord[1].y
          + DzTransFncAtIP3[ip]*ptCoord[2].y + DzTransFncAtIP4[ip]*ptCoord[3].y;

J.J[2][0] = DxTransFncAtIP1[ip]*ptCoord[0].z + DxTransFncAtIP2[ip]*ptCoord[1].z
          + DxTransFncAtIP3[ip]*ptCoord[2].z + DxTransFncAtIP4[ip]*ptCoord[3].z;

J.J[2][1] = DyTransFncAtIP1[ip]*ptCoord[0].z + DyTransFncAtIP2[ip]*ptCoord[1].z
          + DyTransFncAtIP3[ip]*ptCoord[2].z + DyTransFncAtIP4[ip]*ptCoord[3].z;

J.J[2][2] = DzTransFncAtIP1[ip]*ptCoord[0].z + DzTransFncAtIP2[ip]*ptCoord[1].z
          + DzTransFncAtIP3[ip]*ptCoord[2].z + DzTransFncAtIP4[ip]*ptCoord[3].z;

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

void GeTetrahedral::CalcJacobian(Jacobian<Point2D> & J, const Integer ip,
                     const Point2D * const ptCoord, const Boolean NeedJinv)
{
 Error("This element is from 3D", __FILE__, __LINE__);
}

///////////////////////////////////////////////////////////////////////////

Vector<Double> *  GeTetrahedral::GetDxShFncAtIP(const Integer iShFnc)
{
  switch(iShFnc)
{ case 1:
     return &  DxTransFncAtIP1;
  case 2:
     return &  DxTransFncAtIP2;
  case 3:
     return &  DxTransFncAtIP3;
  case 4:
     return &  DxTransFncAtIP4;
}
}

Vector<Double> *  GeTetrahedral::GetDyShFncAtIP(const Integer iShFnc)
{
  switch(iShFnc)
{ case 1:
     return &  DyTransFncAtIP1;
  case 2:
     return &  DyTransFncAtIP2;
  case 3:
     return &  DyTransFncAtIP3;
  case 4:
     return &  DyTransFncAtIP4;
}
}

Vector<Double> *  GeTetrahedral::GetDzShFncAtIP(const Integer iShFnc)
{
  switch(iShFnc)
{ case 1:
     return &  DzTransFncAtIP1;
  case 2:
     return &  DzTransFncAtIP2;
  case 3:
     return &  DzTransFncAtIP3;
  case 4:
     return &  DzTransFncAtIP4;
}
}

} // end of namespace

