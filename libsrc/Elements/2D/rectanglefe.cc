#include <iostream>
#include <fstream>
#include <string>

#include <DataInOut/conffile.hh>
#include <Elements/basefe.hh>
#include "rectanglefe.hh"

namespace CoupledField
{

RectangleFE::RectangleFE()
{
#ifdef TRACE
  (*trace) << "entering RectangleFE::RectangleFE" << std::endl;
#endif
  
  Dim = 2;
  NumEdges = 4;
  NumFaces = 1;


  std::string integtype="GaussOrder2";
  conf->ifget("rectangle",integtype,"IntegRules");

  IntegType=String2EnumIntegrationType(integtype.c_str());

}

RectangleFE :: ~RectangleFE()
{
#ifdef TRACE
  (*trace) << "entering RectangleFE::~RectangleFE" << std::endl;
#endif
}


void RectangleFE:: SetIntPoints()
{
#ifdef TRACE
  (*trace) << "entering RectangleFE::SetIntPoints" << std::endl;
#endif
 
  switch(IntegType) 
    {
    case GaussOrder2:

      NumIntPoints=4;
      DegreeInteg=2;
      
      if ( !IntPoints)
	IntPoints = new Vector<Double>[NumIntPoints];

      for(Integer i=0; i<NumIntPoints; i++)
	IntPoints[i].Resize(Dim);

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

      if( !IntPoints)
	IntPoints = new Vector<Double>[NumIntPoints];
      
      for(Integer i=0; i<NumIntPoints; i++)
	IntPoints[i].Resize(Dim);

      IntWeights.Resize(NumIntPoints);

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

      IntWeights[0]= 0.308642;
      IntWeights[1]= 0.493827;
      IntWeights[2]= 0.308642;
      IntWeights[3]= 0.493827;
      IntWeights[4]= 0.790123; 
      IntWeights[5]= 0.493827; 
      IntWeights[6]= 0.308642; 
      IntWeights[7]= 0.493827; 
      IntWeights[8]= 0.308642; 

      break;

    case GaussOrder7:

      Error("Type of integration Gauss with order 7 is incorrect", __FILE__, __LINE__);
      NumIntPoints=16;
      DegreeInteg=7;
      if ( !IntPoints) 
	IntPoints = new Vector<Double>[NumIntPoints];

      for(Integer i=0; i<NumIntPoints; i++)
	IntPoints[i].Resize(Dim);
      IntWeights.Resize(NumIntPoints);
      
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

      IntWeights[0]= 0.121003;
      IntWeights[1]= 0.226852;
      IntWeights[2]= - 0.226852;
      IntWeights[3]= - 0.121003;
      IntWeights[4]= 0.226852;
      IntWeights[5]= 0.425293;
      IntWeights[6]= -0.425293;
      IntWeights[7]= -0.226852;
      IntWeights[8]= -0.226852;
      IntWeights[9]= -0.425293;
      IntWeights[10]= 0.425293;
      IntWeights[11]= 0.226852;
      IntWeights[12]= -0.121003;
      IntWeights[13]= -0.226852;
      IntWeights[14]= 0.226852;
      IntWeights[15]= 0.121003;

      break;
 
    default:
      std::cerr << "Integration type " << IntegType 
	   << " is not implemented\n" << std::endl; 
      exit (-1);
    }
}

void RectangleFE :: SetShapeFncAtIp()
{
#ifdef TRACE
  (*trace) << "entering RectangleFE::SetShapeFncAtIp" << std::endl;
#endif
  
  if (!ShFncAtIp)
    ShFncAtIp = new Vector<Double>[NumIntPoints];

  for( Integer i=0; i<NumIntPoints; i++ )
    CalcShapeFnc( ShFncAtIp[i], IntPoints[i]);

}
  
void RectangleFE :: SetShapeFncDerivAtIp()
{
#ifdef TRACE
  (*trace) << "entering RectangleFE::SetShapeFncDerivAtIp" << std::endl;
#endif

  if( !ShFncDerivAtIp)
    ShFncDerivAtIp = new Matrix<Double>[NumIntPoints];

  for( Integer i=0; i<NumIntPoints; i++ )
    CalcLocalDerivShapeFnc( ShFncDerivAtIp[i], IntPoints[i]);

}


Double RectangleFE :: CalcJacobianDet(const Vector<Double> & LCoord,
				      const Matrix<Double> & CornerCoords)
{
#ifdef TRACE
  (*trace) << "entering RectangleFE::CalcJacobianDet" << std::endl;
#endif

  Matrix<Double> J;

  CalcJacobian( J, LCoord, CornerCoords );
  return J[0][0]*J[1][1]-J[0][1]*J[1][0];

}

Double RectangleFE :: CalcJacobianDetAtIp(const Integer ip, 
					  const Matrix<Double> & CornerCoords)
{
#ifdef TRACE
  (*trace) << "entering RectangleFE::CalcJacobianDetAtIp" << std::endl;
#endif

  Matrix<Double> J;

  CalcJacobianAtIp( J, ip, CornerCoords);
  return J[0][0]*J[1][1]-J[0][1]*J[1][0]; 
}


void RectangleFE :: CalcJacobian(Matrix<Double> & J, 
				 const Vector<Double> & LCoord, 
				 const Matrix<Double> & CornerCoords)
{
#ifdef TRACE
  (*trace) << "entering RectangleFE::CalcJacobian" << std::endl;
#endif
  J.Resize(2,2);

  Matrix<Double> LDeriv;

  CalcLocalDerivShapeFnc(LDeriv, LCoord);
  J = CornerCoords * LDeriv;
}



void RectangleFE :: CalcJacobianAtIp(Matrix<Double> & J, 
				     const Integer ip, 
				     const Matrix<Double> & CornerCoords)
{
#ifdef TRACE
  (*trace) << "entering RectangleFE::CalcJacobianAtIp" << std::endl;
#endif

  J.Resize(2,2);

  J = CornerCoords * ShFncDerivAtIp[ip-1];
}



  
void RectangleFE :: CalcInvJacobian(Matrix<Double> & JInv,
				    const Vector<Double> & LCoord,
				    const Matrix<Double> & CornerCoords)
{
#ifdef TRACE
  (*trace) << "entering RectangleFE::CalcInvJacobian" << std::endl;
#endif
  
  JInv.Resize(2,2);

  Double detJ, aux;
  Matrix<Double> J, LDeriv;

  J.Resize(2,2);

  CalcLocalDerivShapeFnc(LDeriv, LCoord);
  J = CornerCoords * LDeriv;

  detJ = J[0][0]*J[1][1]-J[0][1]*J[1][0];

  aux=1.0/detJ;

  JInv[0][0] = J[1][1];

  JInv[0][1] = - J[0][1];

  JInv[1][0] = - J[1][0];

  JInv[1][1] = J[0][0];

  JInv*=aux;

}

void RectangleFE :: CalcInvJacobianAtIp(Matrix<Double> & JInv,
					const Integer ip,
					const Matrix<Double> & CornerCoords)
{
#ifdef TRACE
  (*trace) << "entering RectangleFE::CalcInvJacobianAtIp" << std::endl;
#endif
  JInv.Resize(2,2);

  Double detJ, aux;
  Matrix<Double> J;

  J.Resize(2,2);

  J = CornerCoords * ShFncDerivAtIp[ip-1];

  detJ = J[0][0]*J[1][1]-J[0][1]*J[1][0];

  aux=1.0/detJ;

  JInv[0][0] = J[1][1];

  JInv[0][1] = - J[0][1];

  JInv[1][0] = - J[1][0];

  JInv[1][1] = J[0][0];

  JInv*=aux;


}



} // end of namespace
