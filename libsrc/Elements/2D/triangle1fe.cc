#include <iostream>
#include <fstream>

#include "triangle1fe.hh"

namespace CoupledField
{

Triangle1FE :: Triangle1FE() : TriangleFE()
{
  ENTER_FCN( "Triangle1FE::Triangle1FE" );

  Init();
}
  
Triangle1FE :: ~Triangle1FE()
{
  ENTER_FCN( "Triangle1FE::~Triangle1FE" );
}

void Triangle1FE :: Init()
{
  ENTER_FCN( "Triangle1FE::Init" );

  NumNodes_ = 3;
  SetIntPoints();
  SetCornerCoords();
  SetShapeFncAtIp();
  SetShapeFncDerivAtIp();  
}

void Triangle1FE :: SetCornerCoords()
{
  ENTER_FCN( "Triangle1FE::SetCornerCoords" );

  LCornerCoords_.Resize(Dim_,NumNodes_);
  
  LCornerCoords_[0][0] = 0;
  LCornerCoords_[1][0] = 0;
  LCornerCoords_[0][1] = 1;
  LCornerCoords_[1][1] = 0;
  LCornerCoords_[0][2] = 0;
  LCornerCoords_[1][2] = 1;

}

void Triangle1FE :: CalcShapeFnc(Vector<Double> & Shape, 
				 const Vector<Double> & LCoord)
{
  ENTER_FCN( "Triangle1FE::CalcShapeFnc" );

  Shape.Resize(NumNodes_);

  Shape[0] = 1.0 - LCoord[0] - LCoord[1];

  if (Shape[0] < 0)
    Error("Local coordinates are not inside triangular element!",__FILE__,__LINE__);

  for( Integer i=1; i<NumNodes_; i++)
    Shape[i] = LCoord[i-1];

#ifdef DEBUG
  //  (*debug) << "LCoord \n " << LCoord << std::endl;
  //  (*debug) << "Shape \n " << Shape << std::endl;
#endif

}


void Triangle1FE :: CalcLocalDerivShapeFnc(Matrix<Double> & LDeriv, 
				       const Vector<Double> & LCoord)
{
  ENTER_FCN( "Triangle1FE::CalcLocalDerivShapeFnc" );

  LDeriv.Resize(NumNodes_,Dim_);

  LDeriv[0][0] = -1;
  LDeriv[0][1] = -1;
  LDeriv[1][0] =  1;
  LDeriv[1][1] =  0;
  LDeriv[2][0] =  0;
  LDeriv[2][1] =  1; 
}


} // end of namespace


