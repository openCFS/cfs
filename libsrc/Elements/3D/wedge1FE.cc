#include <iostream>
#include <fstream>

#include <General/environment.hh>
#include "wedge1FE.hh"

namespace CoupledField
{

Wedge1FE::Wedge1FE():WedgeFE()
{ 
  ENTER_FCN( "Wedge1FE::Wedge1FE" );
  Init();
}


Wedge1FE::~Wedge1FE()
{
  ENTER_FCN( "Wedge1FE::~Wedge1FE" );
}

void Wedge1FE::Init()
{
  ENTER_IFCN( "Wedge1FE::Init" );

  NumNodes_ = 6;
  NumEdges_ = 9;

  // first set integration points and corner coords ...
  SetIntPoints();
  SetCornerCoords();

  // ... then calc shape function values at integration points
  SetShapeFncAtIp();
  SetShapeFncDerivAtIp();
  // SetEdgeVertices();
}


void Wedge1FE::SetCornerCoords()
{
  ENTER_IFCN( "Wedge1FE::SetCornerCoords" );

  LCornerCoords_.Resize(Dim_,NumNodes_);
  
  LCornerCoords_[0][0] =   1.0;
  LCornerCoords_[1][0] =   0.0;
  LCornerCoords_[2][0] =  -1.0;

  LCornerCoords_[0][1] =  -0.5;
  LCornerCoords_[1][1] =  -0.5*1.732050808;
  LCornerCoords_[2][1] =  -1.0;

  LCornerCoords_[0][2] =  -0.5;
  LCornerCoords_[1][2] =   0.5*1.732050808;
  LCornerCoords_[2][2] =  -1.0;

  LCornerCoords_[0][3] =   1.0;
  LCornerCoords_[1][3] =   0.0;
  LCornerCoords_[2][3] =   1.0;

  LCornerCoords_[0][4] =  -0.5;
  LCornerCoords_[1][4] =  -0.5*1.732050808;
  LCornerCoords_[2][4] =   1.0;

  LCornerCoords_[0][5] =  -0.5;
  LCornerCoords_[1][5] =   0.5*1.732050808;
  LCornerCoords_[2][5] =   1.0;


}


void Wedge1FE :: CalcShapeFnc(Vector<Double> & Shape, 
			     const Vector<Double> & LCoord)
{
  ENTER_IFCN( "Wedge1FE::CalcShapeFnc" );

  Shape.Resize(NumNodes_);

  //"Wedge Elements"
  // From Finite Element Library. www.mathsoft.cse.clrc.ac.uk/felib3/Docs/html/Level-0/wdg6/wdg6.html.
  // corner nodes

  Shape[0] = 0.1666666667 * (1 + 2*LCoord[0]) * (1-LCoord[2]);
  Shape[1] = 0.1666666667 * (1 - LCoord[0] - 1.732050808 * LCoord[1]) * (1-LCoord[2]);
  Shape[2] = 0.1666666667 * (1 - LCoord[0] + 1.732050808 * LCoord[1]) * (1-LCoord[2]);
  Shape[3] = 0.1666666667 * (1 + 2*LCoord[0]) * (1+LCoord[2]);
  Shape[4] = 0.1666666667 * (1 - LCoord[0] - 1.732050808 * LCoord[1]) * (1+LCoord[2]);
  Shape[5] = 0.1666666667 * (1 - LCoord[0] + 1.732050808 * LCoord[1]) * (1+LCoord[2]);

}


void Wedge1FE :: CalcLocalDerivShapeFnc(Matrix<Double> & LDeriv, 
				       const Vector<Double> & LCoord)
{
  ENTER_IFCN( "Wedge1FE::CalcLocalDerivShapeFnc" );

  LDeriv.Resize(NumNodes_,Dim_);

  LDeriv.Init();
  
  LDeriv[0][0] = 0.3333333333 * (1 - LCoord[2]);
  LDeriv[0][1] = 0.0;
  LDeriv[0][2] = 0.1666666667 * (-1 - 2*LCoord[0]);
  
  LDeriv[1][0] = 0.1666666667 * (-1 + LCoord[2]);
  LDeriv[1][1] = 1.7320508080 / 6.0 * (-1 + LCoord[2]);
  LDeriv[1][2] = 0.1666666667 * (-1 + LCoord[0] + 1.7320508080*LCoord[1]);
  
  LDeriv[2][0] = 0.1666666667 * (-1 + LCoord[2]);
  LDeriv[2][1] = 1.7320508080 / 6.0 * (1 - LCoord[2]);
  LDeriv[2][2] = 0.1666666667 * (-1 + LCoord[0] - 1.7320508080*LCoord[1]);
  
  LDeriv[3][0] = 0.3333333333 * (1 + LCoord[2]);
  LDeriv[3][1] = 0.0;
  LDeriv[3][2] = 0.1666666667 * (1 + 2*LCoord[0]);
  
  LDeriv[4][0] = 0.1666666667 * (-1 - LCoord[2]);
  LDeriv[4][1] = 1.7320508080 / 6.0 * (-1 - LCoord[2]);
  LDeriv[4][2] = 0.1666666667 * (1 - LCoord[0] - 1.7320508080*LCoord[1]);
  
  LDeriv[5][0] = 0.1666666667 * (-1 - LCoord[2]);
  LDeriv[5][1] = 1.7320508080 / 6.0 * (1 + LCoord[2]);
  LDeriv[5][2] = 0.1666666667 * (1 - LCoord[0] + 1.7320508080*LCoord[1]);
}
  



} // end of namespace

