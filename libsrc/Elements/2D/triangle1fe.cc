#include <iostream>
#include <fstream>

#include "triangle1fe.hh"

namespace CoupledField
{

Triangle1FE :: Triangle1FE() : TriangleFE()
{
#ifdef TRACE
  (*trace) << "entering Triangle1FE::Triangle1FE" << std::endl;
#endif

  Init();
}
  
Triangle1FE :: ~Triangle1FE()
{
#ifdef TRACE
  (*trace) << "entering Triangle1FE::~Triangle1FE" << std::endl;
#endif
}

void Triangle1FE :: Init()
{
#ifdef TRACE
  (*trace) << "entering Triangle1FE::Init" << std::endl;
#endif  
  NumNodes_ = 3;
  SetIntPoints();
  SetCornerCoords();
  SetShapeFncAtIp();
  SetShapeFncDerivAtIp();  
}

void Triangle1FE :: SetCornerCoords()
{
#ifdef TRACE
  (*trace) << "entering Triangle1FE::SetCornerCoords" << std::endl;
#endif

  LCornerCoords_.Resize(Dim_,NumNodes_);
  
  LCornerCoords_[0][0] = 0;
  LCornerCoords_[1][0] = 0;
  LCornerCoords_[0][1] = 1;
  LCornerCoords_[1][1] = 0;
  LCornerCoords_[0][2] = 0;
  LCornerCoords_[1][2] = 1;

}

void Triangle1FE :: CalcShapeFnc(std::vector<Double> & Shape, 
			     const std::vector<Double> & LCoord)
{
#ifdef TRACE
  (*trace) << "entering Triangle1FE::CalcShapeFnc" << std::endl;
#endif

  Shape.resize(NumNodes_);

  Shape[0] = 1.0 - LCoord[0] - LCoord[1];

  if (Shape[0] < 0)
    Error("Local coordinates are not inside tetrahedral element!",__FILE__,__LINE__);

  for( Integer i=1; i<NumNodes_; i++)
    Shape[i] = LCoord[i-1];

#ifdef DEBUG
  //  (*debug) << "LCoord \n " << LCoord << std::endl;
  //  (*debug) << "Shape \n " << Shape << std::endl;
#endif

}


void Triangle1FE :: CalcLocalDerivShapeFnc(Matrix<Double> & LDeriv, 
				       const std::vector<Double> & LCoord)
{
#ifdef TRACE
  (*trace) << "entering Triangle1FE::CalcLocalDerivShapeFnc" << std::endl;
#endif

  LDeriv.Resize(NumNodes_,Dim_);

  LDeriv[0][0] = -1;
  LDeriv[0][1] = -1;
  LDeriv[1][0] =  1;
  LDeriv[1][1] =  0;
  LDeriv[2][0] =  0;
  LDeriv[2][1] =  1; 
}


} // end of namespace


