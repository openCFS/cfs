#include <iostream>
#include <fstream>

#include "quad1fe.hh"

namespace CoupledField
{


Quad1FE :: Quad1FE() : RectangleFE()
{
#ifdef TRACE
  (*trace) << "entering Quad1FE::Quad1FE" << std::endl;
#endif

  Init();
}
  
Quad1FE :: ~Quad1FE()
{
#ifdef TRACE
  (*trace) << "entering Quad1FE::~Quad1FE" << std::endl;
#endif
}

void Quad1FE :: Init()
{
#ifdef TRACE
  (*trace) << "entering Quad1FE::Init" << std::endl;
#endif  
  NumNodes_ = 4;
  SetIntPoints();
  SetCornerCoords();
  SetShapeFncAtIp();
  SetShapeFncDerivAtIp();  
}

void Quad1FE :: SetCornerCoords()
{
#ifdef TRACE
  (*trace) << "entering Quad1FE::SetCornerCoords" << std::endl;
#endif

  LCornerCoords_.Resize(Dim_,NumNodes_);
  
  LCornerCoords_[0][0] =  1;
  LCornerCoords_[1][0] =  1;
  LCornerCoords_[0][1] = -1;
  LCornerCoords_[1][1] =  1;
  LCornerCoords_[0][2] = -1;
  LCornerCoords_[1][2] = -1;
  LCornerCoords_[0][3] =  1;
  LCornerCoords_[1][3] = -1;

}

void Quad1FE :: CalcShapeFnc(std::vector<Double> & Shape, 
			     const std::vector<Double> & LCoord)
{
#ifdef TRACE
  (*trace) << "entering Quad1FE::CalcShapeFnc" << std::endl;
#endif

  Shape.resize(NumNodes_);
  
  for( Integer i=0; i<NumNodes_; i++)
    Shape[i] = 0.25 * (1 + LCornerCoords_[0][i] * LCoord[0])
                    * (1 + LCornerCoords_[1][i] * LCoord[1]);

}


void Quad1FE :: CalcLocalDerivShapeFnc(Matrix<Double> & LDeriv, 
				       const std::vector<Double> & LCoord)
{
#ifdef TRACE
  (*trace) << "entering Quad1FE::CalcLocalDerivShapeFnc" << std::endl;
#endif

  LDeriv.Resize(NumNodes_,Dim_);

  for( Integer i=0; i<NumNodes_; i++)
    {
      LDeriv[i][0] = 0.25 * LCornerCoords_[0][i] 
                          * (1 + LCornerCoords_[1][i] * LCoord[1] );
      LDeriv[i][1] = 0.25 * (1 + LCornerCoords_[0][i] * LCoord[0] )
                          * LCornerCoords_[i][1];
    }
  
}


} // end of namespace


