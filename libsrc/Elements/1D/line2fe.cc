#include <iostream>
#include <fstream>

#include "line2fe.hh"

namespace CoupledField
{


Line2FE :: Line2FE() : LineFE()
{
#ifdef TRACE
  (*trace) << "entering Line2FE::Line2FE" << std::endl;
#endif

  Init();
}
  
Line2FE :: ~Line2FE()
{
#ifdef TRACE
  (*trace) << "entering Line2FE::~Line2FE" << std::endl;
#endif
}

void Line2FE :: Init()
{
#ifdef TRACE
  (*trace) << "entering Line2FE::Init" << std::endl;
#endif  
  NumNodes_ = 3;
  SetIntPoints();
  SetCornerCoords();
  SetShapeFncAtIp();
  SetShapeFncDerivAtIp();  
}

void Line2FE :: SetCornerCoords()
{
#ifdef TRACE
  (*trace) << "entering Line2FE::SetCornerCoords" << std::endl;
#endif

  LCornerCoords_.Resize(Dim_,NumNodes_);
  
  LCornerCoords_[0][0] =   1;
  LCornerCoords_[0][1] =   0;
  LCornerCoords_[0][2] =  -1;

}

void Line2FE :: CalcShapeFnc(std::vector<Double> & Shape, 
			     const std::vector<Double> & LCoord)
{
#ifdef TRACE
  (*trace) << "entering Line2FE::CalcShapeFnc" << std::endl;
#endif

  Shape.resize(NumNodes_);

  Shape[0] = 0.5*LCoord[0]*(LCoord[0]-1);

  Shape[1] = 1.0 - LCoord[0]*LCoord[0];
  
  Shape[2] = 0.5*LCoord[0]*(LCoord[0]+1);

}


void Line2FE :: CalcLocalDerivShapeFnc(Matrix<Double> & LDeriv, 
				       const std::vector<Double> & LCoord)
{
#ifdef TRACE
  (*trace) << "entering Line2FE::CalcLocalDerivShapeFnc" << std::endl;
#endif

  LDeriv.Resize(NumNodes_,Dim_);

  LDeriv[0][0] = 0.5*(2*LCoord[0] - 1);
  LDeriv[1][0] = -2.0*LCoord[0];
  LDeriv[2][0] = 0.5*(2*LCoord[0] + 1); 

}


} // end of namespace
