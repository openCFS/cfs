#include <iostream>
#include <fstream>
#include <string>

#include "DataInOut/ParamHandling/ConfFile.hh"
#include "DataInOut/ParamHandling/BaseParamHandler.hh"
#include "DataInOut/WriteInfo.hh"
#include "General/environment.hh"
#include "linefe.hh"

namespace CoupledField
{


LineFE :: LineFE()
{
  ENTER_FCN( "LineFE::LineFE" );

  Dim_ = 1;
  NumEdges_ = 1;
  NumFaces_ = 1;
  numChilds_ = 1;
  
#ifndef XMLPARAMS
  std::string integtype="GaussOrder2";
  std::string IntRule;
  if (conf->ifget("IntegRules", IntRule)==TRUE)
      conf->ifget("line",integtype,"IntegRules");
#else
  std::string integtype;
  params->Get( "type", integtype, "integRules", "line" );
#endif

  IntegType=String2EnumIntegrationType(integtype.c_str());
}


LineFE :: ~LineFE()
{
  ENTER_FCN( "LineFE::~LineFE" );
}


void LineFE::SetIntPoints()
{
  ENTER_IFCN( "LineFE::SetIntPoints" );

  switch(IntegType)
    {
    case GaussOrder1:

      NumIntPoints_ = 1;
      DegreeInteg_  = 1;

      if ( !IntPoints_)
	IntPoints_ = new Vector<Double>[NumIntPoints_];
      
      IntWeights_.Resize(NumIntPoints_);
      
      for(Integer i=0; i<NumIntPoints_; i++)
	{	  
	  IntWeights_[i]=1;
	  IntPoints_[i].Resize(Dim_);
	}
      
      IntPoints_[0][0] = 0;
      break;


    case GaussOrder2:

      NumIntPoints_=2;
      DegreeInteg_=2;

      if ( !IntPoints_)
	IntPoints_ = new Vector<Double>[NumIntPoints_];
      
      IntWeights_.Resize(NumIntPoints_);
      
      for(Integer i=0; i<NumIntPoints_; i++)
	{	  
	  IntWeights_[i]=1;
	  IntPoints_[i].Resize(Dim_);
	}
      
      IntPoints_[0][0] = -0.57735026919;
      IntPoints_[1][0] =  0.57735026919;
      break;

    default:
      std::cerr << "Integration type " << IntegType
	   << " is not implemented\n" << std::endl; exit(-1);
    }
}


void LineFE::SetShapeFncAtIp()
{
  ENTER_FCN( "LineFE::SetShapeFncAtIp" );

  if (!ShFncAtIp_)
    ShFncAtIp_ = new Vector<Double>[NumIntPoints_];
  
  for( Integer i=0; i<NumIntPoints_; i++ )
    CalcShapeFnc( ShFncAtIp_[i], IntPoints_[i]);
  
}


void LineFE::SetShapeFncDerivAtIp()
{
  ENTER_FCN( "LineFE::SetShapeFncDerivAtIp" );

  if( !ShFncDerivAtIp_)
    ShFncDerivAtIp_ = new Matrix<Double>[NumIntPoints_];
  
  for( Integer i=0; i<NumIntPoints_; i++ )
    CalcLocalDerivShapeFnc( ShFncDerivAtIp_[i], IntPoints_[i]);
}


Double LineFE::CalcJacobianDet(const Vector<Double> & LCoord,
				 const Matrix<Double> & CornerCoords)
{
  ENTER_FCN( "LineFE::CalcJacobianDet" );

  Matrix<Double> J;
  
  CalcJacobian( J, LCoord, CornerCoords );
  return J[0][0];
  
}


Double LineFE::CalcJacobianDetAtIp(const Integer ip, 
				   const Matrix<Double> & CornerCoords)
{
  ENTER_FCN( "LineFE::CalcJacobianDetAtIp" );

//   Matrix<Double> J;

//   CalcJacobianAtIp( J, ip, CornerCoords);
//   return J[0][0];

  Double length;
  if (CornerCoords.GetSizeRow()==2) 
    {
      //see kaltenbacher, p.23, eq.(2.122)
      Matrix<Double> J;
      CalcJacobianAtIp( J, ip, CornerCoords);
      length = sqrt(J[0][0]*J[0][0] + J[1][0]*J[1][0]);
    }
  else 
    {
      // Length/2 is simply the jacDet for a line
      length=dist_Mat(CornerCoords)/2;
    }

  return length;

  
}

void LineFE:: CalcJacobian(Matrix<Double> & J, 
		  const Vector<Double> & LCoord, 
			  const Matrix<Double> & CornerCoords)
{
  ENTER_FCN( "LineFE::CalcJacobian" );

  J.Resize(1,1);

  Matrix<Double> LDeriv;
  
  CalcLocalDerivShapeFnc(LDeriv, LCoord);
  J = CornerCoords * LDeriv;
}


void LineFE::CalcJacobianAtIp(Matrix<Double> & J, 
		      const Integer ip, 
		      const Matrix<Double> & CornerCoords)
{
  ENTER_FCN( "LineFE::CalcJacobianAtIp" );

  if (CornerCoords.GetSizeRow()==2) {
    // Surface element in 2D 
    J.Resize(CornerCoords.GetSizeRow(),Dim_);
  } 
  else
    J.Resize(1,1);

  J = CornerCoords * ShFncDerivAtIp_[ip-1];
}


void LineFE::CalcInvJacobian(Matrix<Double> & JInv,
			     const Vector<Double> & LCoord,
			     const Matrix<Double> & CornerCoords)
{
  ENTER_FCN( "LineFE::CalcInvJacobian" );
  
  JInv.Resize(1,1);

  Matrix<Double> J, LDeriv;

  J.Resize(1,1);

  CalcLocalDerivShapeFnc(LDeriv, LCoord);
  J = CornerCoords * LDeriv;

  JInv[0][0] = 1 / J[0][0];
}

void LineFE::CalcInvJacobianAtIp(Matrix<Double> & JInv,
				 const Integer ip,
				 const Matrix<Double> & CornerCoords)
{
  ENTER_FCN( "LineFE::CalcInvJacobianAtIp" );
  
  JInv.Resize(1,1);

  Matrix<Double> J;

  J.Resize(1,1);

  J = CornerCoords * ShFncDerivAtIp_[ip-1];

  JInv[0][0] = 1 / J[0][0];
}

} // end of namespace
