#include <iostream.h>
#include <fstream.h>
#include <math.h>

#include <general_head.hh>
#include <utils_head.hh>
#include "baseelem.hh"
#include "rectangle.hh"

namespace CoupledField
{

Rectangle :: ~Rectangle()
{
#ifdef TRACE
  (*trace) << "entering Rectangle::~Rectangle" << endl;
#endif

  ;
}

void Rectangle:: SetIntPoints()
{
#ifdef TRACE
  (*trace) << "entering Rectangle::SetIntPoints" << endl;
#endif
 
  switch(IntegType) 
    {
    case GaussOrder2:

      NumIntPoints=4;
      DegreeInteg=2;
      IntPoints.Resize(NumIntPoints, Dim);
      IntWeights=0;

      IntPoints[0][0] = -1.0/sqrt(3);
      IntPoints[1][0] =  1.0/sqrt(3);
      IntPoints[2][0] =  1.0/sqrt(3);
      IntPoints[3][0] = -1.0/sqrt(3);
      IntPoints[0][1] = -1.0/sqrt(3);
      IntPoints[1][1] = -1.0/sqrt(3);
      IntPoints[2][1] =  1.0/sqrt(3);
      IntPoints[3][1] =  1.0/sqrt(3);
      if (InfoPrint)
       (*infofile) << " For numerical integration procedures we use Gaussian Quadrature with 4 nodes, degree of precision is 3 " << endl; 
      break;

    case GaussOrder5:

      NumIntPoints=9;
      DegreeInteg=5;
      IntPoints.Resize(NumIntPoints, Dim);
      IntWeights=0;

      IntPoints[0][0] = -1.0/sqrt(3);
      IntPoints[1][0] =  1.0/sqrt(3);
      IntPoints[2][0] =  1.0/sqrt(3);
      IntPoints[3][0] = -1.0/sqrt(3);
      IntPoints[0][1] = -1.0/sqrt(3);
      IntPoints[1][1] = -1.0/sqrt(3);
      IntPoints[2][1] =  1.0/sqrt(3);
      IntPoints[3][1] =  1.0/sqrt(3);
      if (InfoPrint)
       (*infofile) << " For numerical integration procedures we use Gaussian Quadrature with 4 nodes, degree of precision is 3 " << endl;
      break;
 
    default:
      cerr << "Integration type " << IntegType 
	   << " is not implemented\n" << endl; exit(-1);
    }
}

void Rectangle :: SetTransformFncAtIntPoints()
{
#ifdef TRACE
  (*trace) << "entering Rectangle::SetShapeFncAtIntPoints" << endl;
#endif
  Integer i;
  TransFncAtIP1.Resize(NumIntPoints);
  TransFncAtIP2.Resize(NumIntPoints);
  TransFncAtIP3.Resize(NumIntPoints);
  TransFncAtIP4.Resize(NumIntPoints);

  for (i=0; i < NumIntPoints; i++)
    { 
      TransFncAtIP1[i]=TransFnc1(IntPoints[i][0],IntPoints[i][1]);
      TransFncAtIP2[i]=TransFnc2(IntPoints[i][0],IntPoints[i][1]);
      TransFncAtIP3[i]=TransFnc3(IntPoints[i][0],IntPoints[i][1]);
      TransFncAtIP4[i]=TransFnc4(IntPoints[i][0],IntPoints[i][1]);
    }
}

void Rectangle :: SetDerTransformFncAtIntPoints()
{
#ifdef TRACE
  (*trace) << "entering Rectangle::SetDerTransFncAtIntPoints" << endl;
#endif
  Integer i;

  DxTransFncAtIP1.Resize(NumIntPoints);
  DxTransFncAtIP2.Resize(NumIntPoints);
  DxTransFncAtIP3.Resize(NumIntPoints);
  DxTransFncAtIP4.Resize(NumIntPoints);

  for (i=0; i < NumIntPoints; i++)
    { 
      DxTransFncAtIP1[i]=TransFnc1dx(IntPoints[i][0],IntPoints[i][1]);
      DxTransFncAtIP2[i]=TransFnc2dx(IntPoints[i][0],IntPoints[i][1]);
      DxTransFncAtIP3[i]=TransFnc3dx(IntPoints[i][0],IntPoints[i][1]);
      DxTransFncAtIP4[i]=TransFnc4dx(IntPoints[i][0],IntPoints[i][1]);
    }

  DyTransFncAtIP1.Resize(NumIntPoints);
  DyTransFncAtIP2.Resize(NumIntPoints);
  DyTransFncAtIP3.Resize(NumIntPoints);
  DyTransFncAtIP4.Resize(NumIntPoints);

  for (i=0; i < NumIntPoints; i++)
    { 
      DyTransFncAtIP1[i]=TransFnc1dy(IntPoints[i][0],IntPoints[i][1]);
      DyTransFncAtIP2[i]=TransFnc2dy(IntPoints[i][0],IntPoints[i][1]);
      DyTransFncAtIP3[i]=TransFnc3dy(IntPoints[i][0],IntPoints[i][1]);
      DyTransFncAtIP4[i]=TransFnc4dy(IntPoints[i][0],IntPoints[i][1]);
    }
}

void Rectangle::CalcJacobian(Jacobian<Point2D> & J, const Integer ip,
                     const Point2D * const ptCoord, const Boolean NeedJinv)
{

 if (!IsSet)
 { 
 SetTransformFncAtIntPoints();
 SetDerTransformFncAtIntPoints(); 
 IsSet=TRUE;
 }

 Integer i;
 Double aux=0;

 J.J[0][0] = DxTransFncAtIP1[ip]*ptCoord[0].x + DxTransFncAtIP2[ip]*ptCoord[1].x
             + DxTransFncAtIP3[ip]*ptCoord[2].x +  DxTransFncAtIP4[ip]*ptCoord[3].x;
 
 J.J[0][1] = DyTransFncAtIP1[ip]*ptCoord[0].x + DyTransFncAtIP2[ip]*ptCoord[1].x
          + DyTransFncAtIP3[ip]*ptCoord[2].x + DyTransFncAtIP4[ip]*ptCoord[3].x;

 J.J[1][0] = DxTransFncAtIP1[ip]*ptCoord[0].y + DxTransFncAtIP2[ip]*ptCoord[1].y
          + DxTransFncAtIP3[ip]*ptCoord[2].y + DxTransFncAtIP4[ip]*ptCoord[3].y;

 J.J[1][1] = DyTransFncAtIP1[ip]*ptCoord[0].y + DyTransFncAtIP2[ip]*ptCoord[1].y
          + DyTransFncAtIP3[ip]*ptCoord[2].y + DyTransFncAtIP4[ip]*ptCoord[3].y;

 J.detJ = J.J[0][0]*J.J[1][1]-J.J[0][1]*J.J[1][0];

if (NeedJinv)
{
 aux=1.0/J.detJ;

 J.Jinv[0][0] = J.J[1][1];

 J.Jinv[0][1] = - J.Jinv[0][1];

 J.Jinv[1][0] = - J.Jinv[1][0];

 J.Jinv[1][1] = J.J[0][0];

 J.Jinv*=aux;
}
}

void Rectangle::CalcJacobian(Jacobian<Point3D> & J, const Integer ip,
                     const Point3D * const ptCoord, const Boolean NeedJinv)
{
 Error("Not implemented yet", __FILE__, __LINE__);
}

} // end of namespace
