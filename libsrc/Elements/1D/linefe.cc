#include <iostream>
#include <fstream>
#include <string>

#include <DataInOut/conffile.hh>
#include <Elements/1D/linefe.hh>

namespace CoupledField
{


LineFE :: LineFE()
{
#ifdef TRACE
  (*trace) << "entering LineFE::LineFE" << std::endl;
#endif

  Dim_ = 1;
  NumEdges_ = 1;
  NumFaces_ = 1;
  
  std::string integtype="GaussOrder2";
  std::string IntRule;
  if (conf->ifget("IntegRules", IntRule)==TRUE)
      conf->ifget("line",integtype,"IntegRules");

  IntegType=String2EnumIntegrationType(integtype.c_str());
}


LineFE :: ~LineFE()
{
#ifdef TRACE
  (*trace) << "entering LineFE::~LineFE" << std::endl;
#endif
}


void LineFE::SetIntPoints()
{
#ifdef TRACE
  (*trace) << "entering LineFE::SetIntPoints" << std::endl;
#endif

  switch(IntegType)
    {
    case GaussOrder1:

      NumIntPoints_ = 1;
      DegreeInteg_  = 1;

      if ( !IntPoints_)
	IntPoints_ = new std::vector<Double>[NumIntPoints_];
      
      IntWeights_.resize(NumIntPoints_);
      
      for(Integer i=0; i<NumIntPoints_; i++)
	{	  
	  IntWeights_[i]=1;
	  IntPoints_[i].resize(Dim_);
	}
      
      IntPoints_[0][0] = 0;
      break;


    case GaussOrder2:

      NumIntPoints_=2;
      DegreeInteg_=2;

      if ( !IntPoints_)
	IntPoints_ = new std::vector<Double>[NumIntPoints_];
      
      IntWeights_.resize(NumIntPoints_);
      
      for(Integer i=0; i<NumIntPoints_; i++)
	{	  
	  IntWeights_[i]=1;
	  IntPoints_[i].resize(Dim_);
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
#ifdef TRACE
  (*trace) << "entering LineFE::SetShapeFncAtIp" << std::endl;
#endif

  if (!ShFncAtIp_)
    ShFncAtIp_ = new std::vector<Double>[NumIntPoints_];
  
  for( Integer i=0; i<NumIntPoints_; i++ )
    CalcShapeFnc( ShFncAtIp_[i], IntPoints_[i]);
  
}


void LineFE::SetShapeFncDerivAtIp()
{
#ifdef TRACE
  (*trace) << "entering LineFE::SetShapeFncDerivAtIp" << std::endl;
#endif

  if( !ShFncDerivAtIp_)
    ShFncDerivAtIp_ = new Matrix<Double>[NumIntPoints_];
  
  for( Integer i=0; i<NumIntPoints_; i++ )
    CalcLocalDerivShapeFnc( ShFncDerivAtIp_[i], IntPoints_[i]);
}


Double LineFE::CalcJacobianDet(const std::vector<Double> & LCoord,
				 const Matrix<Double> & CornerCoords)
{
#ifdef TRACE
  (*trace) << "entering LineFE::CalcJacobianDet" << std::endl;
#endif

  Matrix<Double> J;
  
  CalcJacobian( J, LCoord, CornerCoords );
  return J[0][0];
  
}


Double LineFE::CalcJacobianDetAtIp(const Integer ip, 
				   const Matrix<Double> & CornerCoords)
{
#ifdef TRACE
  (*trace) << "entering LineFE::CalcJacobianDetAtIp" << std::endl;
#endif

//   Matrix<Double> J;

//   CalcJacobianAtIp( J, ip, CornerCoords);
//   return J[0][0];

// Length/2 is simply the jacDet for a line
  Double length;
  length=dist_Mat(CornerCoords);
  return length/2;
}

void LineFE:: CalcJacobian(Matrix<Double> & J, 
		  const std::vector<Double> & LCoord, 
			  const Matrix<Double> & CornerCoords)
{
#ifdef TRACE
  (*trace) << "entering LineFE::CalcJacobian" << std::endl;
#endif

  J.Resize(1,1);

  Matrix<Double> LDeriv;
  
  CalcLocalDerivShapeFnc(LDeriv, LCoord);
  J = CornerCoords * LDeriv;
}


void LineFE::CalcJacobianAtIp(Matrix<Double> & J, 
		      const Integer ip, 
		      const Matrix<Double> & CornerCoords)
{
#ifdef TRACE
  (*trace) << "entering LineFE::CalcJacobianAtIp" << std::endl;
#endif

  J.Resize(1,1);

  J = CornerCoords * ShFncDerivAtIp_[ip-1];
}


void LineFE::CalcInvJacobian(Matrix<Double> & JInv,
			     const std::vector<Double> & LCoord,
			     const Matrix<Double> & CornerCoords)
{
#ifdef TRACE
  (*trace) << "entering LineFE::CalcInvJacobian" << std::endl;
#endif
  
  JInv.Resize(1,1);

  Double detJ, aux;
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
#ifdef TRACE
  (*trace) << "entering LineFE::CalcInvJacobianAtIp" << std::endl;
#endif
  
  JInv.Resize(1,1);

  Double detJ, aux;
  Matrix<Double> J;

  J.Resize(1,1);

  J = CornerCoords * ShFncDerivAtIp_[ip-1];

  JInv[0][0] = 1 / J[0][0];
}

} // end of namespace
