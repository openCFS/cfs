#include <iostream.h>
#include <fstream.h>
#include <math.h>

#include <general_head.hh>
#include <utils_head.hh>
#include "baseelem.hh"
#include "quad1.hh"

namespace CoupledField
{

Quad1 :: Quad1(ShortInt aintegtype) : BaseElem()
{
#ifdef TRACE
  (*trace) << "entering Quad1::Quad1" << endl;
#endif
  ElemType  = QUADRILATERAL1;  
  IntegType = aintegtype;
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
  SetShapeFncAtIntPoints(); 
  SetDShapeFncAtIntPoints();
}


void Quad1:: SetIntPoints()
{
#ifdef TRACE
  (*trace) << "entering Quad1::SetIntPoints" << endl;
#endif
 
  switch(IntegType) 
    {
    case GaussOrder2:
      NumIntPoints=4;
      DegreeInteg=2;
      IntPoints=new  Matrix<Double>(NumIntPoints, Dim);
      IntWeights=0;
      (*IntPoints)[0][0] = -1.0/sqrt(3);
      (*IntPoints)[1][0] =  1.0/sqrt(3);
      (*IntPoints)[2][0] =  1.0/sqrt(3);
      (*IntPoints)[3][0] = -1.0/sqrt(3);
      (*IntPoints)[0][1] = -1.0/sqrt(3);
      (*IntPoints)[1][1] = -1.0/sqrt(3);
      (*IntPoints)[2][1] =  1.0/sqrt(3);
      (*IntPoints)[3][1] =  1.0/sqrt(3);
      if (InfoPrint)
       (*infofile) << " For numerical integration procedures we use Gaussian Quadrature with 4 nodes, degree of precision is 3 " << endl; 
      break;

    case GaussOrder5:
      NumIntPoints=9;
      DegreeInteg=5;
      IntPoints=new  Matrix<Double>(NumIntPoints, Dim);
      IntWeights=9;
      (*IntPoints)[0][0] = -1.0/sqrt(3);
      (*IntPoints)[1][0] =  1.0/sqrt(3);
      (*IntPoints)[2][0] =  1.0/sqrt(3);
      (*IntPoints)[3][0] = -1.0/sqrt(3);
      (*IntPoints)[0][1] = -1.0/sqrt(3);
      (*IntPoints)[1][1] = -1.0/sqrt(3);
      (*IntPoints)[2][1] =  1.0/sqrt(3);
      (*IntPoints)[3][1] =  1.0/sqrt(3);
      if (InfoPrint)
       (*infofile) << " For numerical integration procedures we use Gaussian Quadrature with 4 nodes, degree of precision is 3 " << endl;
      break;
 
    default:
      cerr << "Integration type " << IntegType 
	   << " is not implemented\n" << endl; exit(-1);
    }
 
}

void Quad1 :: SetShapeFncAtIntPoints()
{
#ifdef TRACE
  (*trace) << "entering Quad1::SetShapeFncAtIntPoints" << endl;
#endif
  Integer i;
  ShFncAtIP1.Resize(NumIntPoints);
  ShFncAtIP2.Resize(NumIntPoints);
  ShFncAtIP3.Resize(NumIntPoints);
  ShFncAtIP4.Resize(NumIntPoints);

  for (i=0; i < NumIntPoints; i++)
    { 
      ShFncAtIP1[i]=ShapeFnc1((*IntPoints)[i][0],(*IntPoints)[i][1]);
      ShFncAtIP2[i]=ShapeFnc2((*IntPoints)[i][0],(*IntPoints)[i][1]);
      ShFncAtIP3[i]=ShapeFnc3((*IntPoints)[i][0],(*IntPoints)[i][1]);
      ShFncAtIP4[i]=ShapeFnc4((*IntPoints)[i][0],(*IntPoints)[i][1]);
    }
}

void Quad1 :: SetDShapeFncAtIntPoints()
{
#ifdef TRACE
  (*trace) << "entering Quad1::SetDShapeFncAtIntPoints" << endl;
#endif
  Integer i;

  DxShFncAtIP1.Resize(NumIntPoints);
  DxShFncAtIP2.Resize(NumIntPoints);
  DxShFncAtIP3.Resize(NumIntPoints);
  DxShFncAtIP4.Resize(NumIntPoints);

  for (i=0; i < NumIntPoints; i++)
    { 
      DxShFncAtIP1[i]=ShapeFnc1dx((*IntPoints)[i][0],(*IntPoints)[i][1]);
      DxShFncAtIP2[i]=ShapeFnc2dx((*IntPoints)[i][0],(*IntPoints)[i][1]);
      DxShFncAtIP3[i]=ShapeFnc3dx((*IntPoints)[i][0],(*IntPoints)[i][1]);
      DxShFncAtIP4[i]=ShapeFnc4dx((*IntPoints)[i][0],(*IntPoints)[i][1]);
    }

  DyShFncAtIP1.Resize(NumIntPoints);
  DyShFncAtIP2.Resize(NumIntPoints);
  DyShFncAtIP3.Resize(NumIntPoints);
  DyShFncAtIP4.Resize(NumIntPoints);

  for (i=0; i < NumIntPoints; i++)
    { 
      DyShFncAtIP1[i]=ShapeFnc1dy((*IntPoints)[i][0],(*IntPoints)[i][1]);
      DyShFncAtIP2[i]=ShapeFnc2dy((*IntPoints)[i][0],(*IntPoints)[i][1]);
      DyShFncAtIP3[i]=ShapeFnc3dy((*IntPoints)[i][0],(*IntPoints)[i][1]);
      DyShFncAtIP4[i]=ShapeFnc4dy((*IntPoints)[i][0],(*IntPoints)[i][1]);
    }
}

void Quad1::CalcJacobian(Jacobian<Point2D> & J, const Integer ip,
                     const Point2D * const ptCoord, const Boolean NeedJinv)
{
 Integer i;
 Double aux=0;

 J.J[0][0] = DxShFncAtIP1[ip]*ptCoord[0].x + DxShFncAtIP2[ip]*ptCoord[1].x
             + DxShFncAtIP3[ip]*ptCoord[2].x +  DxShFncAtIP4[ip]*ptCoord[3].x;
 
 J.J[0][1] = DyShFncAtIP1[ip]*ptCoord[0].x + DyShFncAtIP2[ip]*ptCoord[1].x
          + DyShFncAtIP3[ip]*ptCoord[2].x + DyShFncAtIP4[ip]*ptCoord[3].x;

 J.J[1][0] = DxShFncAtIP1[ip]*ptCoord[0].y + DxShFncAtIP2[ip]*ptCoord[1].y
          + DxShFncAtIP3[ip]*ptCoord[2].y + DxShFncAtIP4[ip]*ptCoord[3].y;

 J.J[1][1] = DyShFncAtIP1[ip]*ptCoord[0].y + DyShFncAtIP2[ip]*ptCoord[1].y
          + DyShFncAtIP3[ip]*ptCoord[2].y + DyShFncAtIP4[ip]*ptCoord[3].y;

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

void Quad1::CalcJacobian(Jacobian<Point3D> & J, const Integer ip,
                     const Point3D * const ptCoord, const Boolean NeedJinv)
{
 Error("Not implemented yet", __FILE__, __LINE__);
}

void  Quad1::GetGradientShFnc(Vector<Double> & grad, const Integer i, const Integer ip)
{
  grad.Resize(2);

  switch(i)
{
 case 1:
    grad[0]=DxShFncAtIP1[ip]; grad[1]=DyShFncAtIP1[ip];
 case 2:
    grad[0]=DxShFncAtIP2[ip]; grad[1]=DyShFncAtIP2[ip]; 
 case 3:
    grad[0]=DxShFncAtIP3[ip]; grad[1]=DyShFncAtIP3[ip]; 
 case 4:
    grad[0]=DxShFncAtIP4[ip]; grad[1]=DyShFncAtIP4[ip];
} 

}

Vector<Double> & Quad1::GetShFncAtIP(const Integer iShFnc) 
{
  switch(iShFnc)
    { 
    case 1:
      return ShFncAtIP1;
    case 2:
      return ShFncAtIP2;
    case 3:
      return ShFncAtIP3;
    case 4:
      return ShFncAtIP4;
    }
}
}
