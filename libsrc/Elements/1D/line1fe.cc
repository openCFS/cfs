#include <iostream>
#include <fstream>

#include "line1fe.hh"

namespace CoupledField
{


Line1FE :: Line1FE() : GeLineFE()
{
#ifdef TRACE
  (*trace) << "entering Line1FE::Line1FE" << std::endl;
#endif

  Init();
}
  
Line1FE :: ~Line1FE()
{
#ifdef TRACE
  (*trace) << "entering Line1FE::~Line1FE" << std::endl;
#endif
}

void Line1FE :: Init()
{
#ifdef TRACE
  (*trace) << "entering Line1FE::Init" << std::endl;
#endif  
  NumNodes_ = 2;
  SetIntPoints();
  SetCornerCoords();
  SetShapeFncAtIp();
  SetShapeFncDerivAtIp();  
}

void Line1FE :: SetCornerCoords()
{
#ifdef TRACE
  (*trace) << "entering Line1FE::SetCornerCoords" << std::endl;
#endif

  LCornerCoords_.Resize(Dim_,NumNodes_);
  
  LCornerCoords_[0][0] =   1;
  LCornerCoords_[0][1] =  -1;

}

void Line1FE :: CalcShapeFnc(std::vector<Double> & Shape, 
			     const std::vector<Double> & LCoord)
{
#ifdef TRACE
  (*trace) << "entering Line1FE::CalcShapeFnc" << std::endl;
#endif

  Shape.resize(NumNodes_);
  
  for( Integer i=0; i<NumNodes_; i++)
    Shape[i] = 0.5*(1+LCornerCoords_[0][i] * LCoord[0]);

}


void Line1FE :: CalcLocalDerivShapeFnc(Matrix<Double> & LDeriv, 
				       const std::vector<Double> & LCoord)
{
#ifdef TRACE
  (*trace) << "entering Line1FE::CalcLocalDerivShapeFnc" << std::endl;
#endif

  LDeriv.Resize(NumNodes_,Dim_);

  for( Integer i=0; i<NumNodes_; i++)
    {
      LDeriv[i][0] = 0.5*LCornerCoords_[0][i];

    }
  
}


} // end of namespace
