#include <iostream>
#include <fstream>

#include "triangle2fe.hh"

namespace CoupledField
{

Triangle2FE :: Triangle2FE() : TriangleFE()
{
#ifdef TRACE
  (*trace) << "entering Triangle2FE::Triangle2FE" << std::endl;
#endif

  Init();
}
  
Triangle2FE :: ~Triangle2FE()
{
#ifdef TRACE
  (*trace) << "entering Triangle2FE::~Triangle2FE" << std::endl;
#endif
}

void Triangle2FE :: Init()
{
#ifdef TRACE
  (*trace) << "entering Triangle2FE::Init" << std::endl;
#endif
  NumNodes_ = 6;
  SetIntPoints();
  SetCornerCoords();
  SetShapeFncAtIp();
  SetShapeFncDerivAtIp();  
}

void Triangle2FE :: SetCornerCoords()
{
#ifdef TRACE
  (*trace) << "entering Triangle2FE::SetCornerCoords" << std::endl;
#endif

  LCornerCoords_.Resize(Dim_,NumNodes_);
  
  LCornerCoords_[0][0] = 0;
  LCornerCoords_[1][0] = 0;
  LCornerCoords_[0][1] = 1;
  LCornerCoords_[1][1] = 0;
  LCornerCoords_[0][2] = 0;
  LCornerCoords_[1][2] = 1;

  LCornerCoords_[0][3] = 0.5;
  LCornerCoords_[1][3] = 0;
  LCornerCoords_[0][4] = 0.5;
  LCornerCoords_[1][4] = 0.5;
  LCornerCoords_[0][5] = 0;
  LCornerCoords_[1][5] = 0.5;
}

void Triangle2FE :: CalcShapeFnc(std::vector<Double> & Shape, 
			     const std::vector<Double> & LCoord)
{
#ifdef TRACE
  (*trace) << "entering Triangle2FE::CalcShapeFnc" << std::endl;
#endif
   // From Zienkiewicz, The Finite Element Method. Vol 1, page 128.
 // corner nodes
  Shape.resize(NumNodes_);

  Double t = 1.0 - LCoord[0] - LCoord[1]; // Define the third component of the triangular coord.

  Shape[0] = t * (2*t - 1);
  Shape[1] = LCoord[0]*(2*LCoord[0] - 1);
  Shape[2] = LCoord[1]*(2*LCoord[1] - 1);
  Shape[3] = 4 * LCoord[0] * t;
  Shape[4] = 4 * LCoord[0] * LCoord[1];
  Shape[5] = 4 * LCoord[1] * t;

//   if (Shape[0] < 0)
//     Error("Local coordinates are not inside tetrahedral element!",__FILE__,__LINE__);

#ifdef DEBUG
  //  (*debug) << "LCoord \n " << LCoord << std::endl;
  //  (*debug) << "Shape \n " << Shape << std::endl;
#endif

}


void Triangle2FE :: CalcLocalDerivShapeFnc(Matrix<Double> & LDeriv, 
				       const std::vector<Double> & LCoord)
{
#ifdef TRACE
  (*trace) << "entering Triangle2FE::CalcLocalDerivShapeFnc" << std::endl;
#endif

  LDeriv.Resize(NumNodes_,Dim_);

  LDeriv[0][0] =  4*LCoord[0] + 4*LCoord[1] - 3.0;
  LDeriv[0][1] =  4*LCoord[0] + 4*LCoord[1] - 3.0;

  LDeriv[1][0] =  4*LCoord[0]-1;
  LDeriv[1][1] =  0;

  LDeriv[2][0] =  0;
  LDeriv[2][1] =  4*LCoord[1]-1; 

  LDeriv[3][0] =  4*(1 - 2*LCoord[0] - LCoord[1]);
  LDeriv[3][1] = -4*LCoord[0];

  LDeriv[4][0] =  4* LCoord[1];
  LDeriv[4][1] =  4* LCoord[0];

  LDeriv[5][0] = -4*LCoord[1];
  LDeriv[5][1] =  4*(1 - 2*LCoord[1] - LCoord[0]); 
}


} // end of namespace


