#include <iostream>
#include <fstream>
#include <string>

#include <DataInOut/conffile.hh>
#include <Elements/1D/gelinefe.hh>

namespace CoupledField
{


GeLineFE :: GeLineFE()
{
#ifdef TRACE
  (*trace) << "entering GeLineFE::GeLineFE" << std::endl;
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


GeLineFE :: ~GeLineFE()
{
#ifdef TRACE
  (*trace) << "entering GeLineFE::~GeLineFE" << std::endl;
#endif
}


void GeLineFE::SetIntPoints()
{
#ifdef TRACE
  (*trace) << "entering GeLineFE::SetIntPoints" << std::endl;
#endif

  switch(IntegType)
    {
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


void GeLineFE::SetShapeFncAtIp()
{
#ifdef TRACE
  (*trace) << "entering GeLineFE::SetShapeFncAtIp" << std::endl;
#endif

  if (!ShFncAtIp_)
    ShFncAtIp_ = new std::vector<Double>[NumIntPoints_];
  
  for( Integer i=0; i<NumIntPoints_; i++ )
    CalcShapeFnc( ShFncAtIp_[i], IntPoints_[i]);
  
}


void GeLineFE::SetShapeFncDerivAtIp()
{
#ifdef TRACE
  (*trace) << "entering GeLineFE::SetShapeFncDerivAtIp" << std::endl;
#endif

  if( !ShFncDerivAtIp_)
    ShFncDerivAtIp_ = new Matrix<Double>[NumIntPoints_];
  
  for( Integer i=0; i<NumIntPoints_; i++ )
    CalcLocalDerivShapeFnc( ShFncDerivAtIp_[i], IntPoints_[i]);
}


Double GeLineFE::CalcJacobianDet(const std::vector<Double> & LCoord,
				 const Matrix<Double> & CornerCoords)
{
#ifdef TRACE
  (*trace) << "entering GeLineFE::CalcJacobianDet" << std::endl;
#endif

  Matrix<Double> J;
  
  CalcJacobian( J, LCoord, CornerCoords );
  return J[0][0];
  
}


Double GeLineFE::CalcJacobianDetAtIp(const Integer ip, 
				   const Matrix<Double> & CornerCoords)
{
#ifdef TRACE
  (*trace) << "entering GeLineFE::CalcJacobianDetAtIp" << std::endl;
#endif

  Matrix<Double> J;

  CalcJacobianAtIp( J, ip, CornerCoords);
  return J[0][0];

}

void GeLineFE:: CalcJacobian(Matrix<Double> & J, 
		  const std::vector<Double> & LCoord, 
			  const Matrix<Double> & CornerCoords)
{
#ifdef TRACE
  (*trace) << "entering GeLineFE::CalcJacobian" << std::endl;
#endif

  J.Resize(1,1);

  Matrix<Double> LDeriv;
  
  CalcLocalDerivShapeFnc(LDeriv, LCoord);
  J = CornerCoords * LDeriv;
}


void GeLineFE::CalcJacobianAtIp(Matrix<Double> & J, 
		      const Integer ip, 
		      const Matrix<Double> & CornerCoords)
{
#ifdef TRACE
  (*trace) << "entering GeLineFE::CalcJacobianAtIp" << std::endl;
#endif

  J.Resize(1,1);

  J = CornerCoords * ShFncDerivAtIp_[ip-1];
}


void GeLineFE::CalcInvJacobian(Matrix<Double> & JInv,
			     const std::vector<Double> & LCoord,
			     const Matrix<Double> & CornerCoords)
{
#ifdef TRACE
  (*trace) << "entering GeLineFE::CalcInvJacobian" << std::endl;
#endif
  
  JInv.Resize(1,1);

  Double detJ, aux;
  Matrix<Double> J, LDeriv;

  J.Resize(1,1);

  CalcLocalDerivShapeFnc(LDeriv, LCoord);
  J = CornerCoords * LDeriv;

  JInv[0][0] = 1 / J[0][0];
}

void GeLineFE::CalcInvJacobianAtIp(Matrix<Double> & JInv,
				 const Integer ip,
				 const Matrix<Double> & CornerCoords)
{
#ifdef TRACE
  (*trace) << "entering GeLineFE::CalcInvJacobianAtIp" << std::endl;
#endif
  
  JInv.Resize(1,1);

  Double detJ, aux;
  Matrix<Double> J;

  J.Resize(1,1);

  J = CornerCoords * ShFncDerivAtIp_[ip-1];

  JInv[0][0] = 1 / J[0][0];
}

} // end of namespace
