#include <iostream>
#include <fstream>
#include <string>

#include "basefe.hh"
#include <Utils/tools.hh> 
 
namespace CoupledField
{

 BaseFE :: BaseFE()
{
#ifdef TRACE
  (*trace) << "entering BaseFE::BaseFE" << std::endl;
#endif
  
  // initializing dynamic objects
  ShFncAtIp_      = NULL;
  ShFncDerivAtIp_ = NULL; 
  IntPoints_      = NULL; 
  
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


  //  Deriv.Resize(NumNodes_,Dim_);
  Matrix<Double> LDeriv, JInv;

  CalcLocalDerivShapeFnc(LDeriv, LCoord);

  CalcInvJacobian(JInv, LCoord, CornerCoords);
 
  Deriv = LDeriv * JInv;
}

void BaseFE :: GetGlobDerivShFncAtIp(Matrix<Double> & Deriv, 
				     const Integer ip,
				     const Matrix<Double> & CornerCoords,
				     Double & jacDet)
{
#ifdef TRACE
  (*trace) << "entering BaseFE::GetGlobDerivShFncAtIp" << std::endl;
#endif

  //  Deriv.Resize(NumNodes_,Dim_);
  Matrix<Double> JInv;

  CalcInvJacobianAtIp(JInv, ip, CornerCoords);

  Deriv = ShFncDerivAtIp_[ip-1] * JInv;

  // det(A) = 1 / det(A^(-1))
  jacDet = 1.0 / JInv.Det();

}


void BaseFE :: GetGlobDerivShFncAtIp(Matrix<Double> & Deriv, 
				     const Integer ip,
				     const Matrix<Double> & CornerCoords)
{
#ifdef TRACE
  (*trace) << "entering BaseFE::GetGlobDerivShFncAtIp" << std::endl;
#endif

  //  Deriv.Resize(NumNodes_,Dim_);
  Matrix<Double> JInv;

  CalcInvJacobianAtIp(JInv, ip, CornerCoords);

  Deriv = ShFncDerivAtIp_[ip-1] * JInv;
}



void BaseFE :: CalcJacobian(Matrix<Double> & J, 
			    const std::vector<Double> & LCoord, 
			    const Matrix<Double> & CornerCoords)
{
#ifdef TRACE
  (*trace) << "entering BaseFE::CalcJacobian" << std::endl;
#endif
  //  J.Resize(Dim_,Dim_);

  Matrix<Double> LDeriv;

  CalcLocalDerivShapeFnc(LDeriv, LCoord);
  J = CornerCoords * LDeriv;
}


void BaseFE :: CalcJacobianAtIp(Matrix<Double> & J, 
				const Integer ip, 
				const Matrix<Double> & CornerCoords)
{
#ifdef TRACE
  (*trace) << "entering BaseFE::CalcJacobianAtIp" << std::endl;
#endif

  //  J.Resize(Dim_,Dim_);

  J = CornerCoords * ShFncDerivAtIp_[ip-1];
}


Double BaseFE :: CalcJacobianDet(const std::vector<Double> & LCoord,
				      const Matrix<Double> & CornerCoords)
{
#ifdef TRACE
  (*trace) << "entering BaseFE::CalcJacobianDet" << std::endl;
#endif

  Matrix<Double> J;

  CalcJacobian( J, LCoord, CornerCoords );
  return J.Det();
}





Double BaseFE :: CalcJacobianDetAtIp(const Integer ip, 
				     const Matrix<Double> & CornerCoords)
{
#ifdef TRACE
  (*trace) << "entering BaseFE::CalcJacobianDetAtIp" << std::endl;
#endif

  Matrix<Double> J;

  CalcJacobianAtIp( J, ip, CornerCoords);

  return J.Det();
}


void BaseFE :: CalcInvJacobian(Matrix<Double> & JInv,
			       const std::vector<Double> & LCoord,
			       const Matrix<Double> & CornerCoords)
{
#ifdef TRACE
  (*trace) << "entering BaseFE::CalcInvJacobian" << std::endl;
#endif
  
  Matrix<Double> J, LDeriv;

  //  J.Resize(Dim_,Dim_);

  CalcLocalDerivShapeFnc(LDeriv, LCoord);

  J = CornerCoords * LDeriv;

  J.Invert(JInv);
}



 
void BaseFE :: CalcInvJacobianAtIp(Matrix<Double> & JInv,
				   const Integer ip,
				   const Matrix<Double> & CornerCoords)
{
#ifdef TRACE
  (*trace) << "entering BaseFE::CalcInvJacobianAtIp" << std::endl;
#endif
  JInv.Resize(Dim_,Dim_);

  Double detJ, aux;
  Matrix<Double> J;

  //  J.Resize(Dim_,Dim_);

  J = CornerCoords * ShFncDerivAtIp_[ip-1];

  J.Invert(JInv);
}


void BaseFE :: SetShapeFncAtIp()
{
#ifdef TRACE
  (*trace) << "entering BaseFE::SetShapeFncAtIp" << std::endl;
#endif
  
  if (!ShFncAtIp_)
    ShFncAtIp_ = new std::vector<Double>[NumIntPoints_];


  for( Integer i=0; i<NumIntPoints_; i++ )      
      CalcShapeFnc( ShFncAtIp_[i], IntPoints_[i]);
}
  
void BaseFE :: SetShapeFncDerivAtIp()
{
#ifdef TRACE
  (*trace) << "entering BaseFE::SetShapeFncDerivAtIp" << std::endl;
#endif

  if( !ShFncDerivAtIp_)
    ShFncDerivAtIp_ = new Matrix<Double>[NumIntPoints_];

  for( Integer i=0; i<NumIntPoints_; i++ )
    CalcLocalDerivShapeFnc( ShFncDerivAtIp_[i], IntPoints_[i]);

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



// calculate global derivates of edge shape functions in integration point ip
void BaseFE::GetEdgeGlobDerivShFncAtIp(std::vector< Matrix<Double>* > & deriv, 
			       const Integer ip,
			       const Matrix<Double> & cornerCoords)
{
#ifdef TRACE
  (*trace) << "entering BaseFE::GetEdgeGlobDerivShFncAtIp" << std::endl;
#endif

  // vector of coordinates of the desired integration point
  std::vector<Double> lCoord;
  lCoord.resize(Dim_);

  for (Integer i=0; i<Dim_; i++)
    lCoord[i] = IntPoints_[ip-1][i];
  
  GetEdgeGlobalDerivShapeFnc(deriv, lCoord, cornerCoords);
}

void BaseFE::CalcEdgeShapeFncAtIp(Matrix<Double> & shape, 
				  const Integer ip,
				  const Matrix<Double> & cornerCoords)
{
  std::vector<Double> lCoord;
  lCoord.resize(Dim_);

  for (Integer i=0; i<Dim_; i++)
    lCoord[i] = IntPoints_[ip-1][i];
  
  CalcEdgeShapeFnc(shape, lCoord);
}


} // end namespace CoupledField
