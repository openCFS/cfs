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
  
  Dim_ = 2;
  NumEdges_ = 4;
  NumFaces_ = 1;


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

      NumIntPoints_=4;
      DegreeInteg_=2;
      
      if ( !IntPoints_)
	IntPoints_ = new std::vector<Double>[NumIntPoints_];

      for(Integer i=0; i<NumIntPoints_; i++)
	IntPoints_[i].resize(Dim_);

      IntPoints_[0][0] = -0.57735026919;
      IntPoints_[1][0] =  0.57735026919;
      IntPoints_[2][0] =  0.57735026919;
      IntPoints_[3][0] = -0.57735026919;
      IntPoints_[0][1] = -0.57735026919;
      IntPoints_[1][1] = -0.57735026919;
      IntPoints_[2][1] =  0.57735026919;
      IntPoints_[3][1] =  0.57735026919;

      break;

    case GaussOrder5:

      NumIntPoints_=9;
      DegreeInteg_=5;

      if( !IntPoints_)
	IntPoints_ = new std::vector<Double>[NumIntPoints_];
      
      for(Integer i=0; i<NumIntPoints_; i++)
	IntPoints_[i].resize(Dim_);

      IntWeights_.resize(NumIntPoints_);

      IntPoints_[0][0] = -0.774596669241483;
      IntPoints_[1][0] =  0.0;
      IntPoints_[2][0] =  0.774596669241483;
      IntPoints_[3][0] = -0.774596669241483;
      IntPoints_[4][0] = 0.0;
      IntPoints_[5][0] =  0.774596669241483;
      IntPoints_[6][0] = -0.774596669241483;
      IntPoints_[7][0] = 0.0;
      IntPoints_[8][0] =  0.774596669241483; 

      IntPoints_[0][1] = -0.774596669241483;
      IntPoints_[1][1] = -0.774596669241483;
      IntPoints_[2][1] = -0.774596669241483;
      IntPoints_[3][1] = 0.0;
      IntPoints_[4][1] = 0.0;
      IntPoints_[5][1] = 0.0;
      IntPoints_[6][1] =  0.774596669241483;
      IntPoints_[7][1] =  0.774596669241483;
      IntPoints_[8][1] =  0.774596669241483;

      IntWeights_[0]= 0.308642;
      IntWeights_[1]= 0.493827;
      IntWeights_[2]= 0.308642;
      IntWeights_[3]= 0.493827;
      IntWeights_[4]= 0.790123; 
      IntWeights_[5]= 0.493827; 
      IntWeights_[6]= 0.308642; 
      IntWeights_[7]= 0.493827; 
      IntWeights_[8]= 0.308642; 

      break;

    case GaussOrder7:

      Error("Type of integration Gauss with order 7 is incorrect", __FILE__, __LINE__);
      NumIntPoints_=16;
      DegreeInteg_=7;
      if ( !IntPoints_) 
	IntPoints_ = new std::vector<Double>[NumIntPoints_];

      for(Integer i=0; i<NumIntPoints_; i++)
	IntPoints_[i].resize(Dim_);
      IntWeights_.resize(NumIntPoints_);
      
      IntPoints_[0][0] = -0.861136311594053;
      IntPoints_[1][0] =  -0.339981043584856;
      IntPoints_[2][0] =  0.339981043584856;
      IntPoints_[3][0] = 0.861136311594053;
      IntPoints_[4][0] = -0.861136311594053;
      IntPoints_[5][0] = -0.339981043584856;
      IntPoints_[6][0] =   0.339981043584856;
      IntPoints_[7][0] =  0.861136311594053;
      IntPoints_[8][0] =  -0.861136311594053;
      IntPoints_[9][0] =  -0.339981043584856;
      IntPoints_[10][0] =  0.339981043584856;
      IntPoints_[11][0] =  0.861136311594053;
      IntPoints_[12][0] =  -0.861136311594053;
      IntPoints_[13][0] =  -0.339981043584856;
      IntPoints_[14][0] =  0.339981043584856;
      IntPoints_[15][0] =  0.861136311594053;


      IntPoints_[0][1] =  -0.861136311594053;
      IntPoints_[1][1] =   -0.861136311594053;
      IntPoints_[2][1] =  -0.861136311594053;
      IntPoints_[3][1] =  -0.861136311594053;
      IntPoints_[4][1] =  -0.339981043584856;
      IntPoints_[5][1] =  -0.339981043584856;
      IntPoints_[6][1] =  -0.339981043584856;
      IntPoints_[7][1] =  -0.339981043584856;
      IntPoints_[8][1] =  0.339981043584856;
      IntPoints_[9][1] =  0.339981043584856;
      IntPoints_[10][1] =  0.339981043584856;
      IntPoints_[11][1] =  0.339981043584856;
      IntPoints_[12][1] =  0.861136311594053;
      IntPoints_[13][1] =  0.861136311594053;
      IntPoints_[14][1] =  0.861136311594053;
      IntPoints_[15][1] =  0.861136311594053;

      IntWeights_[0]= 0.121003;
      IntWeights_[1]= 0.226852;
      IntWeights_[2]= - 0.226852;
      IntWeights_[3]= - 0.121003;
      IntWeights_[4]= 0.226852;
      IntWeights_[5]= 0.425293;
      IntWeights_[6]= -0.425293;
      IntWeights_[7]= -0.226852;
      IntWeights_[8]= -0.226852;
      IntWeights_[9]= -0.425293;
      IntWeights_[10]= 0.425293;
      IntWeights_[11]= 0.226852;
      IntWeights_[12]= -0.121003;
      IntWeights_[13]= -0.226852;
      IntWeights_[14]= 0.226852;
      IntWeights_[15]= 0.121003;

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
  
  if (!ShFncAtIp_)
    ShFncAtIp_ = new std::vector<Double>[NumIntPoints_];

  for( Integer i=0; i<NumIntPoints_; i++ )
    CalcShapeFnc( ShFncAtIp_[i], IntPoints_[i]);

}
  
void RectangleFE :: SetShapeFncDerivAtIp()
{
#ifdef TRACE
  (*trace) << "entering RectangleFE::SetShapeFncDerivAtIp" << std::endl;
#endif

  if( !ShFncDerivAtIp_)
    ShFncDerivAtIp_ = new Matrix<Double>[NumIntPoints_];

  for( Integer i=0; i<NumIntPoints_; i++ )
    CalcLocalDerivShapeFnc( ShFncDerivAtIp_[i], IntPoints_[i]);

}


Double RectangleFE :: CalcJacobianDet(const std::vector<Double> & LCoord,
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
				 const std::vector<Double> & LCoord, 
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

  J = CornerCoords * ShFncDerivAtIp_[ip-1];
}



  
void RectangleFE :: CalcInvJacobian(Matrix<Double> & JInv,
				    const std::vector<Double> & LCoord,
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

  J = CornerCoords * ShFncDerivAtIp_[ip-1];

  detJ = J[0][0]*J[1][1]-J[0][1]*J[1][0];

  aux=1.0/detJ;

  JInv[0][0] = J[1][1];

  JInv[0][1] = - J[0][1];

  JInv[1][0] = - J[1][0];

  JInv[1][1] = J[0][0];

  JInv*=aux;


}



} // end of namespace
