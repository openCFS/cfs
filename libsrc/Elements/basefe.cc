#include <iostream>
#include <fstream>
#include <string>

#include "basefe.hh"
 
namespace CoupledField
{

 BaseFE :: BaseFE()
{
#ifdef TRACE
  (*trace) << "entering BaseFE::BaseFE" << std::endl;
#endif
}
 
BaseFE :: ~BaseFE()
{
#ifdef TRACE
  (*trace) << "entering BaseFE::~BaseFE" << std::endl;
#endif
 
  if( ShFncAtIp_ ) delete[] ShFncAtIp_;
  if( ShFncDerivAtIp_ ) delete[] ShFncDerivAtIp_;
  if( IntPoints_ ) delete[] IntPoints_;
}

void BaseFE :: GetShFnc(std::vector<double> & S, 
			const std::vector<Double> & LCoord)
{
#ifdef TRACE
  (*trace) << "entering BaseFE::GetShFnc" << std::endl;
#endif

  CalcShapeFnc(S, LCoord);

}

void BaseFE :: GetShFncAtIp(std::vector<Double> & S, 
			    const Integer ip)
{
#ifdef TRACE
  (*trace) << "entering BaseFE::GetShFncAtIp" << std::endl;
#endif

  S.resize(NumNodes_);
  
  S = ShFncAtIp_[ip-1];
  
}

void BaseFE :: GetGlobDerivShFnc(Matrix<Double> & Deriv, 
				 const std::vector<Double> & LCoord,
				 const Matrix<Double> & CornerCoords)
{
#ifdef TRACE
  (*trace) << "entering BaseFE::GetGlobDerivShFnc" << std::endl;
#endif


  Deriv.Resize(Dim_,Dim_);
  Matrix<Double> LDeriv, JInv;

  CalcLocalDerivShapeFnc(LDeriv, LCoord);
  CalcInvJacobian(JInv, LCoord, CornerCoords);

  Deriv = LDeriv * JInv;
}

void BaseFE :: GetGlobDerivShFncAtIp(Matrix<Double> & Deriv, 
				     const Integer ip,
				     const Matrix<Double> & CornerCoords)
{
#ifdef TRACE
  (*trace) << "entering BaseFE::GetGlobDerivShFncAtIp" << std::endl;
#endif

  Deriv.Resize(Dim_,Dim_);
  Matrix<Double> JInv;

  CalcInvJacobianAtIp(JInv, ip, CornerCoords);

  Deriv = ShFncDerivAtIp_[ip-1] * JInv;
}



enum IntegrationType BaseFE::String2EnumIntegrationType(const Char * inttype)
{
#ifdef TRACE
  (*trace) << "entering BaseFE::String2EnumIntegrationType" << std::endl;
#endif

 enum IntegrationType result=null;
 if (!std::strcmp(inttype,"GaussOrder2")) result=GaussOrder2;
 if (!std::strcmp(inttype,"GaussOrder3")) result=GaussOrder3;
 if (!std::strcmp(inttype,"GaussOrder4")) result=GaussOrder4;
 if (!std::strcmp(inttype,"GaussOrder5")) result=GaussOrder5;
 if (!std::strcmp(inttype,"GaussOrder7")) result=GaussOrder7;

 if (result==null) Error("Check your config file. Your integration type is wrong",__FILE__,__LINE__);

 return result;
}

}
