#include <iostream>
#include <fstream>

#include "getriangle.hh"

namespace CoupledField
{
                   
GeTriangle :: ~GeTriangle()
{
#ifdef TRACE
  (*trace) << "entering GeTriangle::~GeTriangle" << std::endl;
#endif

  ;
}

void GeTriangle :: SetIntPoints()
{
#ifdef TRACE
  (*trace) << "entering GeTriangle::SetIntPoints" << std::endl;
#endif

 switch(IntegType)
    {
    case GaussOrder3:

      NumIntPoints=4;
      DegreeInteg=3;
      IntPoints.Resize(NumIntPoints, Dim);
      IntWeights=new Vector<Double>(NumIntPoints);
      IntPoints[0][0] = 1.0/3.0;
      IntPoints[1][0] = 3.0/5.0;
      IntPoints[2][0] = 1.0/5.0;
      IntPoints[3][0] = 1.0/5.0;
      IntPoints[0][1] = 1.0/3.0;
      IntPoints[1][1] = 1.0/5.0;
      IntPoints[2][1] = 3.0/5.0;
      IntPoints[3][1] = 1.0/5.0;
      
      (*IntWeights)[0]=-0.5625;
      (*IntWeights)[1]=0.520833333333333;
      (*IntWeights)[2]=0.520833333333333;
      (*IntWeights)[3]=0.520833333333333;

      if (InfoPrint)
       (*infofile) << " For numerical integration procedures we use Gaussian Quadrature with 4 nodes, degree of precision is 3 " << std::endl;
      break;

    case GaussOrder2:

      NumIntPoints=3;
      DegreeInteg=2;
      IntPoints.Resize(NumIntPoints, Dim);
      IntWeights=new Vector<Double>(NumIntPoints);

      IntPoints[0][0] = 0.166666666666667;
      IntPoints[1][0] = 0.666666666666667; 
      IntPoints[2][0] = 0.166666666666667;
      IntPoints[0][1] = 0.166666666666667;
      IntPoints[1][1] = 0.166666666666667;
      IntPoints[2][1] = 0.666666666666667;

      (*IntWeights)[0]= 0.166666666666667 ;
      (*IntWeights)[1]= 0.166666666666667 ;
      (*IntWeights)[2]= 0.166666666666667 ;

      if (InfoPrint)
    (*infofile) << " For numerical integration procedures we use Gaussian Quadrature with 3 nodes, degree of precision is 2 " << std::endl;
      break;

   case GaussOrder4:

// D. A. Dunavant: High Degree Efficient Symmetrical ...
// Int. J. Numer. Methods in Eng.  Vol 21  S. 1129-1148   1985
// (c) Kaskade

      NumIntPoints=6;
      DegreeInteg=4;
      IntPoints.Resize(NumIntPoints, Dim);
      IntWeights=new Vector<Double>(NumIntPoints);
 
      IntPoints[0][0] = 4.45948490915965e-01;
      IntPoints[1][0] = 1.08103018168070e-01 ;
      IntPoints[2][0] = 4.45948490915965e-01;
      IntPoints[3][0] = 9.15762135097710e-02; 
      IntPoints[4][0] = 8.16847572980459e-01; 
      IntPoints[5][0] = 9.15762135097710e-02;  

      IntPoints[0][1] = 4.45948490915965e-01;
      IntPoints[1][1] = 4.45948490915965e-01;
      IntPoints[2][1] = 1.08103018168070e-01;
      IntPoints[3][1] = 9.15762135097710e-02; 
      IntPoints[4][1] = 9.15762135097710e-02; 
      IntPoints[5][1] = 8.16847572980459e-01;  

      (*IntWeights)[0]= 1.116907948390055e-01*2;
      (*IntWeights)[1]= 1.116907948390055e-01*2;
      (*IntWeights)[2]= 1.116907948390055e-01*2;
      (*IntWeights)[3]= 5.497587182766100e-02*2;
      (*IntWeights)[4]= 5.497587182766100e-02*2; 
      (*IntWeights)[5]= 5.497587182766100e-02*2; 
 
      if (InfoPrint)
    (*infofile) << " For numerical integration procedures we use Gaussian Quadrature with 6 nodes, degree of precision is 4 " << std::endl;
      break;

   case GaussOrder5:
 
// D. A. Dunavant: High Degree Efficient Symmetrical ...
// Int. J. Numer. Methods in Eng.  Vol 21  S. 1129-1148   1985
// (c) Kaskade
 
      NumIntPoints=7;
      DegreeInteg=5;
      IntPoints.Resize(NumIntPoints, Dim);
      IntWeights=new Vector<Double>(NumIntPoints);
 
      IntPoints[0][0] = 3.333333333333333e-01;
      IntPoints[1][0] = 4.701420641051151e-01;
      IntPoints[2][0] = 5.971587178976981e-02;
      IntPoints[3][0] = 4.701420641051151e-01;
      IntPoints[4][0] = 1.012865073234563e-01;
      IntPoints[5][0] = 7.974269853530872e-01;
      IntPoints[6][0] = 1.012865073234563e-01;      
 
      IntPoints[0][1] = 3.333333333333333e-01;
      IntPoints[1][1] = 4.701420641051151e-01;
      IntPoints[2][1] = 4.701420641051151e-01;
      IntPoints[3][1] = 5.971587178976981e-02;
      IntPoints[4][1] = 1.012865073234563e-01;
      IntPoints[5][1] = 1.012865073234563e-01;
      IntPoints[5][1] = 7.974269853530872e-01;
 
      (*IntWeights)[0]= 0.2250300003;
      (*IntWeights)[1]= 0.132394152788506;
      (*IntWeights)[2]= 0.132394152788506;
      (*IntWeights)[3]= 0.132394152788506 ;
      (*IntWeights)[4]= 0.125939180544827;
      (*IntWeights)[5]=  0.125939180544827;
      (*IntWeights)[6]=  0.125939180544827; 
 
      if (InfoPrint)
    (*infofile) << " For numerical integration procedures we use Gaussian Quadrature with 7 nodes, degree of precision is 5 " << std::endl;
      break;
 
    default:
      std::cerr << "Integration type " << IntegType
           << " is not implemented \n" << std::endl; exit(-1);
    }
}

void GeTriangle :: SetTransformFncAtIntPoints()
{
#ifdef TRACE
  (*trace) << "entering GeTriangle::SetShapeFncAtIntPoints" << std::endl;
#endif
 Integer i;
 TransFncAtIP1.Resize(NumIntPoints);
 TransFncAtIP2.Resize(NumIntPoints);
 TransFncAtIP3.Resize(NumIntPoints);

  for (i=0; i < NumIntPoints; i++)
    {
      TransFncAtIP1[i]=TransFnc1(IntPoints[i][0],IntPoints[i][1]);
      TransFncAtIP2[i]=TransFnc2(IntPoints[i][0],IntPoints[i][1]);
      TransFncAtIP3[i]=TransFnc3(IntPoints[i][0],IntPoints[i][1]);
    }
}

void GeTriangle :: SetDerTransformFncAtIntPoints()
{
#ifdef TRACE
  (*trace) << "entering GeTriangle::SetDShapeFnc" << std::endl;
#endif
 Integer i;

  DxTransFncAtIP1.Resize(NumIntPoints);
  DxTransFncAtIP2.Resize(NumIntPoints);
  DxTransFncAtIP3.Resize(NumIntPoints);

  for (i=0; i < NumIntPoints; i++)
    {
      DxTransFncAtIP1[i]=TransFnc1dx(IntPoints[i][0],IntPoints[i][1]);
      DxTransFncAtIP2[i]=TransFnc2dx(IntPoints[i][0],IntPoints[i][1]);
      DxTransFncAtIP3[i]=TransFnc3dx(IntPoints[i][0],IntPoints[i][1]);
    }

  DyTransFncAtIP1.Resize(NumIntPoints);
  DyTransFncAtIP2.Resize(NumIntPoints);
  DyTransFncAtIP3.Resize(NumIntPoints);
 
  for (i=0; i < NumIntPoints; i++)
    {
      DyTransFncAtIP1[i]=TransFnc1dy(IntPoints[i][0],IntPoints[i][1]);
      DyTransFncAtIP2[i]=TransFnc2dy(IntPoints[i][0],IntPoints[i][1]);
      DyTransFncAtIP3[i]=TransFnc3dy(IntPoints[i][0],IntPoints[i][1]);
    } 
}

void GeTriangle::CalcJacobian(Jacobian<Point2D> & J, const Integer ip,
                     const Point2D * const ptCoord, const Boolean NeedJinv)
{

  if (!IsSet)
 {
 SetTransformFncAtIntPoints();
 SetDerTransformFncAtIntPoints();
 IsSet=TRUE;
 }

 Double aux=0;

J.J[0][0] = DxTransFncAtIP1[ip]*ptCoord[0].x + DxTransFncAtIP2[ip]*ptCoord[1].x
             + DxTransFncAtIP3[ip]*ptCoord[2].x ;

J.J[0][1] = DyTransFncAtIP1[ip]*ptCoord[0].x + DyTransFncAtIP2[ip]*ptCoord[1].x
          + DyTransFncAtIP3[ip]*ptCoord[2].x ;
 
J.J[1][0] = DxTransFncAtIP1[ip]*ptCoord[0].y + DxTransFncAtIP2[ip]*ptCoord[1].y
          + DxTransFncAtIP3[ip]*ptCoord[2].y ;
 
J.J[1][1] = DyTransFncAtIP1[ip]*ptCoord[0].y + DyTransFncAtIP2[ip]*ptCoord[1].y
          + DyTransFncAtIP3[ip]*ptCoord[2].y ;
 
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

void GeTriangle::CalcJacobian(Jacobian<Point3D> & J, const Integer ip,
                     const Point3D * const ptCoord, const Boolean NeedJinv)
{
 Error("Not implemented yet", __FILE__, __LINE__);
}

} // end of namespace
