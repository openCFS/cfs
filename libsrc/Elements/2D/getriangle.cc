#include <iostream>
#include <fstream>
#include <string>

#include <DataInOut/conffile.hh>
#include "getriangle.hh"

namespace CoupledField
{

GeTriangle::GeTriangle()
{
#ifdef TRACE
  (*trace) << "entering GeTriangle::GeTriangle" << std::endl;
#endif
 
  std::string integtype="GaussOrder2";
  conf->ifget("triangle",integtype,"IntegRules");

  IntegType=String2EnumIntegrationType(integtype.c_str());

	isSetAtCenter_=FALSE;
}

GeTriangle :: ~GeTriangle()
{
#ifdef TRACE
  (*trace) << "entering GeTriangle::~GeTriangle" << std::endl;
#endif

  if (IntWeights) delete IntWeights;
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
      IntPoints[6][1] = 7.974269853530872e-01;
 
      (*IntWeights)[0]= 0.2250300003;
      (*IntWeights)[1]= 0.132394152788506;
      (*IntWeights)[2]= 0.132394152788506;
      (*IntWeights)[3]= 0.132394152788506 ;
      (*IntWeights)[4]= 0.125939180544827;
      (*IntWeights)[5]=  0.125939180544827;
      (*IntWeights)[6]=  0.125939180544827; 

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

void GeTriangle::CalcJacobian(Jacobian<2> & J, const Integer ip,
                     Point<2> * ptCoord, const Boolean NeedJinv)
{

  if (!IsSet)
 {
 SetTransformFncAtIntPoints();
 SetDerTransformFncAtIntPoints();
 IsSet=TRUE;
 }

 Double aux=0;

J.J[0][0] = DxTransFncAtIP1[ip]*ptCoord[0][0] + DxTransFncAtIP2[ip]*ptCoord[1][0]
             + DxTransFncAtIP3[ip]*ptCoord[2][0] ;

J.J[0][1] = DyTransFncAtIP1[ip]*ptCoord[0][0] + DyTransFncAtIP2[ip]*ptCoord[1][0]
          + DyTransFncAtIP3[ip]*ptCoord[2][0] ;
 
J.J[1][0] = DxTransFncAtIP1[ip]*ptCoord[0][1] + DxTransFncAtIP2[ip]*ptCoord[1][1]
          + DxTransFncAtIP3[ip]*ptCoord[2][1] ;
 
J.J[1][1] = DyTransFncAtIP1[ip]*ptCoord[0][1] + DyTransFncAtIP2[ip]*ptCoord[1][1]
          + DyTransFncAtIP3[ip]*ptCoord[2][1] ;
 
J.detJ = J.J[0][0]*J.J[1][1]-J.J[0][1]*J.J[1][0];


 // calculation of inverse matrix
if (NeedJinv)
{
 aux=1.0/J.detJ;
 
 J.Jinv[0][0] = J.J[1][1];
 
 J.Jinv[0][1] = - J.J[0][1];
 
 J.Jinv[1][0] = - J.J[1][0];
 
 J.Jinv[1][1] = J.J[0][0];
 
 J.Jinv*=aux;
}

}


void GeTriangle::CalcJacobian(Jacobian<3> & J, const Integer ip,
                     Point<3> * ptCoord, const Boolean NeedJinv)
{
 Error("This element is from 2D", __FILE__, __LINE__);
}

void GeTriangle::CalcJacobianAtCenter(Jacobian<2> & J,
              Point<2> * ptCoord, const Boolean NeedJinv)
{

  if (!isSetAtCenter_)
    { 
      SetTransformFncAtCenter();
      SetDerTransformFncAtCenter(); 
      isSetAtCenter_=TRUE;
    }

  Double aux=0;

  // Init
  J.J[0][0]  = 0;
  J.J[0][1]  = 0;
  J.J[1][0]  = 0;
  J.J[1][1]  = 0;

  Integer ish;
  for (ish=0; ish < 3; ish++) {
    J.J[0][0]  += DxTransFncAtCenter[ish]*ptCoord[ish][0];
    J.J[0][1]  += DyTransFncAtCenter[ish]*ptCoord[ish][0];
    J.J[1][0]  += DxTransFncAtCenter[ish]*ptCoord[ish][1];
    J.J[1][1]  += DyTransFncAtCenter[ish]*ptCoord[ish][1];
  }

  
  J.detJ = J.J[0][0]*J.J[1][1]-J.J[0][1]*J.J[1][0];

  if (NeedJinv)
    {
      aux=1.0/J.detJ;

      J.Jinv[0][0] = J.J[1][1];

      J.Jinv[0][1] = - J.J[0][1];

      J.Jinv[1][0] = - J.J[1][0];

      J.Jinv[1][1] = J.J[0][0];

      J.Jinv*=aux;
    }
}

void GeTriangle :: SetTransformFncAtCenter()
{
#ifdef TRACE
  (*trace) << "entering GeTriangle::SetTransformFncAtCenter" << std::endl;
#endif
  
  TransFncAtCenter[0]=TransFnc1(0,0);
  TransFncAtCenter[1]=TransFnc2(0,0);
  TransFncAtCenter[2]=TransFnc3(0,0);
  
}

void GeTriangle :: SetDerTransformFncAtCenter()
{
#ifdef TRACE
  (*trace) << "entering GeTriangle::SetDerTransFncAtCenter" << std::endl;
#endif
 
  DxTransFncAtCenter[0]=TransFnc1dx(0,0);
  DxTransFncAtCenter[1]=TransFnc2dx(0,0);
  DxTransFncAtCenter[2]=TransFnc3dx(0,0);

      
  DyTransFncAtCenter[0]=TransFnc1dy(0,0);
  DyTransFncAtCenter[1]=TransFnc2dy(0,0);
  DyTransFncAtCenter[2]=TransFnc3dy(0,0);

    
}

} // end of namespace
