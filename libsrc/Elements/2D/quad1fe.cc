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
  NumNodes = 4;
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

  LCornerCoords.Resize(Dim,NumNodes);
  
  LCornerCoords[0][0] =  1;
  LCornerCoords[1][0] =  1;
  LCornerCoords[0][1] = -1;
  LCornerCoords[1][1] =  1;
  LCornerCoords[0][2] = -1;
  LCornerCoords[1][2] = -1;
  LCornerCoords[0][3] =  1;
  LCornerCoords[1][3] = -1;

}

void Quad1FE :: CalcShapeFnc(Vector<Double> & LShape, 
			     const Vector<Double> & LCoord)
{
#ifdef TRACE
  (*trace) << "entering Quad1FE::CalcShapeFnc" << std::endl;
#endif

  LShape.Resize(NumNodes);
  
  for( Integer i=0; i<NumNodes; i++)
    LShape[i] = 0.25 * (1 + LCornerCoords[0][i] * LCoord[0])
                     * (1 + LCornerCoords[1][i] * LCoord[1]);

}


void Quad1FE :: CalcLocalDerivShapeFnc(Matrix<Double> & LDeriv, 
				       const Vector<Double> & LCoord)
{
#ifdef TRACE
  (*trace) << "entering Quad1FE::CalcLocalDerivShapeFnc" << std::endl;
#endif

  LDeriv.Resize(NumNodes,Dim);

  for( Integer i=0; i<NumNodes; i++)
    {
      LDeriv[i][0] = 0.25 * LCornerCoords[0][i] 
                          * (1 + LCornerCoords[1][i] * LCoord[1] );
      LDeriv[i][1] = 0.25 * (1 + LCornerCoords[0][i] * LCoord[0] )
                          * LCornerCoords[i][1];
    }
  
}


} // end of namespace


