#include <iostream>
#include <fstream>

#include <General/environment.hh>
#include "tetra1FE.hh"

namespace CoupledField
{

Tetra1FE::Tetra1FE():TetraFE()
{ 
#ifdef TRACE
  (*trace) << "entering Tetra1FE::Tetra1FE" << std::endl;
#endif
   Init();
}


Tetra1FE::~Tetra1FE()
{
#ifdef TRACE
  (*trace) << "entering Tetra1FE::~Tetra1FE" << std::endl;
#endif
 ; 
}

void Tetra1FE::Init()
{
#ifdef TRACE
  (*trace) << "entering Tetra1FE::Init" << std::endl;
#endif
  
  Dim_ = 3;
  NumNodes_ = 4;

  // first set integration points and corner coords ...
  SetIntPoints();
  SetCornerCoords();

  // then calc shape function values at integration points
  SetShapeFncAtIp();
  SetShapeFncDerivAtIp();
}


void Tetra1FE::SetCornerCoords()
{
#ifdef TRACE
  (*trace) << "entering Tetra1FE::SetCornerCoords" << std::endl;
#endif

  LCornerCoords_.Resize(Dim_,NumNodes_);
  
  LCornerCoords_[0][0] =  0;
  LCornerCoords_[1][0] =  0;
  LCornerCoords_[2][0] =  0;

  LCornerCoords_[0][1] =  1;
  LCornerCoords_[1][1] =  0;
  LCornerCoords_[2][1] =  0;

  LCornerCoords_[0][2] =  0;
  LCornerCoords_[1][2] =  1;
  LCornerCoords_[2][2] =  0;

  LCornerCoords_[0][3] =  0;
  LCornerCoords_[1][3] =  0;
  LCornerCoords_[2][3] =  1;
}


std::ostream& operator<< (std::ostream & outStr, std::vector<Double> xOut)
{
  for (Integer i=0; i<xOut.size(); i++)
    outStr <<  " " << xOut[i];
  return outStr;
}


void Tetra1FE :: CalcShapeFnc(std::vector<Double> & Shape, 
			     const std::vector<Double> & LCoord)
{
#ifdef TRACE
  (*trace) << "entering Tetra1FE::CalcShapeFnc" << std::endl;
#endif

  Shape.resize(NumNodes_);

  // see Hughes p. 170
  
  Shape[0] = 1.0 - LCoord[0] - LCoord[1] - LCoord[2];

  if (Shape[0] < 0)
    Error("Local coordinates are not inside tetrahedral element!",__FILE__,__LINE__);

  for( Integer i=1; i<NumNodes_; i++)
    Shape[i] = LCoord[i-1];

#ifdef DEBUG
  *debug << "LCoord \n " << LCoord << std::endl;
  *debug << "Shape \n " << Shape << std::endl;
#endif
}


void Tetra1FE :: CalcLocalDerivShapeFnc(Matrix<Double> & LDeriv, 
				       const std::vector<Double> & LCoord)
{
#ifdef TRACE
  (*trace) << "entering Tetra1FE::CalcLocalDerivShapeFnc" << std::endl;
#endif

  LDeriv.Resize(NumNodes_,Dim_);

  LDeriv.Init();
  
  for( Integer i=0; i<Dim_; i++)
    LDeriv[0][i] = -1.0;

  for( Integer i=1; i < NumNodes_; i++)
    LDeriv[i][i-1] = 1.0;

#ifdef DEBUG
  (*debug) << "LDeriv \n " << LDeriv << std::endl;
#endif

}
  



} // end of namespace

