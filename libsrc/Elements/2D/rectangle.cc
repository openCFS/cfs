#include <iostream>
#include <fstream>
#include <string>

#include <DataInOut/conffile.hh>
#include "rectangle.hh"

namespace CoupledField
{

Rectangle::Rectangle()
{
#ifdef TRACE
  (*trace) << "entering Rectangle::Rectangle" << std::endl;
#endif

  std::string integtype;
  conf->get("rectangle",integtype,"IntegRules");

  IntegType=String2EnumIntegrationType(integtype.c_str());

  IsSet=FALSE;
  isSetAtCenter_=FALSE;
}

Rectangle :: ~Rectangle()
{
#ifdef TRACE
  (*trace) << "entering Rectangle::~Rectangle" << std::endl;
#endif

  if (IntWeights) delete IntWeights;
}

void Rectangle:: SetIntPoints()
{
#ifdef TRACE
  (*trace) << "entering Rectangle::SetIntPoints" << std::endl;
#endif
 
  switch(IntegType) 
    {
    case GaussOrder2:

      NumIntPoints=4;
      DegreeInteg=2;
      IntPoints.Resize(NumIntPoints, Dim);
      IntWeights=0;

      IntPoints[0][0] = -0.57735026919;
      IntPoints[1][0] =  0.57735026919;
      IntPoints[2][0] =  0.57735026919;
      IntPoints[3][0] = -0.57735026919;
      IntPoints[0][1] = -0.57735026919;
      IntPoints[1][1] = -0.57735026919;
      IntPoints[2][1] =  0.57735026919;
      IntPoints[3][1] =  0.57735026919;

      break;

    case GaussOrder5:

      NumIntPoints=9;
      DegreeInteg=5;
      IntPoints.Resize(NumIntPoints, Dim);
      IntWeights=new Vector<Double>(NumIntPoints);

      IntPoints[0][0] = -0.774596669241483;
      IntPoints[1][0] =  0.0;
      IntPoints[2][0] =  0.774596669241483;
      IntPoints[3][0] = -0.774596669241483;
      IntPoints[4][0] = 0.0;
      IntPoints[5][0] =  0.774596669241483;
      IntPoints[6][0] = -0.774596669241483;
      IntPoints[7][0] = 0.0;
      IntPoints[8][0] =  0.774596669241483; 

      IntPoints[0][1] = -0.774596669241483;
      IntPoints[1][1] = -0.774596669241483;
      IntPoints[2][1] = -0.774596669241483;
      IntPoints[3][1] = 0.0;
      IntPoints[4][1] = 0.0;
      IntPoints[5][1] = 0.0;
      IntPoints[6][1] =  0.774596669241483;
      IntPoints[7][1] =  0.774596669241483;
      IntPoints[8][1] =  0.774596669241483;

      (*IntWeights)[0]= 0.308642;
      (*IntWeights)[1]= 0.493827;
      (*IntWeights)[2]= 0.308642;
      (*IntWeights)[3]= 0.493827;
      (*IntWeights)[4]= 0.790123; 
      (*IntWeights)[5]= 0.493827; 
      (*IntWeights)[6]= 0.308642; 
      (*IntWeights)[7]= 0.493827; 
      (*IntWeights)[8]= 0.308642; 

      break;

    case GaussOrder7:

      Error("Type of integration Gauss with order 7 is incorrect", __FILE__, __LINE__);
      NumIntPoints=16;
      DegreeInteg=7;
      IntPoints.Resize(NumIntPoints, Dim);
      IntWeights=new Vector<Double>(NumIntPoints);
 
      IntPoints[0][0] = -0.861136311594053;
      IntPoints[1][0] =  -0.339981043584856;
      IntPoints[2][0] =  0.339981043584856;
      IntPoints[3][0] = 0.861136311594053;
      IntPoints[4][0] = -0.861136311594053;
      IntPoints[5][0] = -0.339981043584856;
      IntPoints[6][0] =   0.339981043584856;
      IntPoints[7][0] =  0.861136311594053;
      IntPoints[8][0] =  -0.861136311594053;
      IntPoints[9][0] =  -0.339981043584856;
      IntPoints[10][0] =  0.339981043584856;
      IntPoints[11][0] =  0.861136311594053;
      IntPoints[12][0] =  -0.861136311594053;
      IntPoints[13][0] =  -0.339981043584856;
      IntPoints[14][0] =  0.339981043584856;
      IntPoints[15][0] =  0.861136311594053;


      IntPoints[0][1] =  -0.861136311594053;
      IntPoints[1][1] =   -0.861136311594053;
      IntPoints[2][1] =  -0.861136311594053;
      IntPoints[3][1] =  -0.861136311594053;
      IntPoints[4][1] =  -0.339981043584856;
      IntPoints[5][1] =  -0.339981043584856;
      IntPoints[6][1] =  -0.339981043584856;
      IntPoints[7][1] =  -0.339981043584856;
      IntPoints[8][1] =  0.339981043584856;
      IntPoints[9][1] =  0.339981043584856;
      IntPoints[10][1] =  0.339981043584856;
      IntPoints[11][1] =  0.339981043584856;
      IntPoints[12][1] =  0.861136311594053;
      IntPoints[13][1] =  0.861136311594053;
      IntPoints[14][1] =  0.861136311594053;
      IntPoints[15][1] =  0.861136311594053;

      (*IntWeights)[0]= 0.121003;
      (*IntWeights)[1]= 0.226852;
      (*IntWeights)[2]= - 0.226852;
      (*IntWeights)[3]= - 0.121003;
      (*IntWeights)[4]= 0.226852;
      (*IntWeights)[5]= 0.425293;
      (*IntWeights)[6]= -0.425293;
      (*IntWeights)[7]= -0.226852;
      (*IntWeights)[8]= -0.226852;
      (*IntWeights)[9]= -0.425293;
      (*IntWeights)[10]= 0.425293;
      (*IntWeights)[11]= 0.226852;
      (*IntWeights)[12]= -0.121003;
      (*IntWeights)[13]= -0.226852;
      (*IntWeights)[14]= 0.226852;
      (*IntWeights)[15]= 0.121003;

      break;
 
    default:
      std::cerr << "Integration type " << IntegType 
	   << " is not implemented\n" << std::endl; exit(-1);
    }
}

void Rectangle :: SetTransformFncAtIntPoints()
{
#ifdef TRACE
  (*trace) << "entering Rectangle::SetTransformFncAtIntPoints" << std::endl;
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
  (*trace) << "entering Rectangle::SetDerTransFncAtIntPoints" << std::endl;
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

void Rectangle :: SetTransformFncAtCenter()
{
#ifdef TRACE
  (*trace) << "entering Rectangle::SetTransformFncAtCenter" << std::endl;
#endif
  
  TransFncAtCenter[0]=TransFnc1(0,0);
  TransFncAtCenter[1]=TransFnc2(0,0);
  TransFncAtCenter[2]=TransFnc3(0,0);
  TransFncAtCenter[3]=TransFnc4(0,0);
  
}

void Rectangle :: SetDerTransformFncAtCenter()
{
#ifdef TRACE
  (*trace) << "entering Rectangle::SetDerTransFncAtCenter" << std::endl;
#endif
 
  DxTransFncAtCenter[0]=TransFnc1dx(0,0);
  DxTransFncAtCenter[1]=TransFnc2dx(0,0);
  DxTransFncAtCenter[2]=TransFnc3dx(0,0);
  DxTransFncAtCenter[3]=TransFnc4dx(0,0);
      
  DyTransFncAtCenter[0]=TransFnc1dy(0,0);
  DyTransFncAtCenter[1]=TransFnc2dy(0,0);
  DyTransFncAtCenter[2]=TransFnc3dy(0,0);
  DyTransFncAtCenter[3]=TransFnc4dy(0,0);
    
}

void Rectangle::CalcJacobian(Jacobian<2> & J, const Integer ip,
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
             + DxTransFncAtIP3[ip]*ptCoord[2][0] +  DxTransFncAtIP4[ip]*ptCoord[3][0];
 
 J.J[0][1] = DyTransFncAtIP1[ip]*ptCoord[0][0] + DyTransFncAtIP2[ip]*ptCoord[1][0]
          + DyTransFncAtIP3[ip]*ptCoord[2][0] + DyTransFncAtIP4[ip]*ptCoord[3][0];

 J.J[1][0] = DxTransFncAtIP1[ip]*ptCoord[0][1] + DxTransFncAtIP2[ip]*ptCoord[1][1]
          + DxTransFncAtIP3[ip]*ptCoord[2][1] + DxTransFncAtIP4[ip]*ptCoord[3][1];

 J.J[1][1] = DyTransFncAtIP1[ip]*ptCoord[0][1] + DyTransFncAtIP2[ip]*ptCoord[1][1]
          + DyTransFncAtIP3[ip]*ptCoord[2][1] + DyTransFncAtIP4[ip]*ptCoord[3][1];

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

void Rectangle::CalcJacobianAtCenter(Jacobian<2> & J,
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
  for (ish=0; ish < 4; ish++) {
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

void Rectangle::CalcJacobian(Jacobian<3> & J, const Integer ip,
                     Point<3> * ptCoord, const Boolean NeedJinv)
{
 Error("Not implemented yet", __FILE__, __LINE__);
}

} // end of namespace
