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
  
  LCornerCoords_[0][0] =  -1;
  LCornerCoords_[1][0] =  -1;
  LCornerCoords_[0][1] = 1;
  LCornerCoords_[1][1] =  -1;
  LCornerCoords_[0][2] = 1;
  LCornerCoords_[1][2] = 1;
  LCornerCoords_[0][3] = -1;
  LCornerCoords_[1][3] = 1;


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
                          * LCornerCoords_[1][i];
    }
  
}


Double Quad1FE::CalcMeanStrain(Matrix<Double> &cornerCoords, Array<Double> &displacements)
{
#ifdef TRACE
  (*trace) << "entering Quad1FE::CalcDistortion" << std::endl;
#endif

  Double factor;
  Double eps1, eps2, eps4, eps5, eps11, eps12, eps21, eps22, eps41, eps42, eps51, eps52;
  Double length1, length2, length11, length12, length21, length22;


  // calculate inital size of element
  length11 = abs(cornerCoords[0][0] - cornerCoords[0][1]);
  length12 = abs(cornerCoords[0][3] - cornerCoords[0][2]);
  length1 = (length11+length12) * 0.5;
  
  length21 = abs(cornerCoords[1][0] - cornerCoords[1][3]);
  length22 = abs(cornerCoords[1][1] - cornerCoords[1][2]);
  length2 = (length21+length22) * 0.5;

   // calculate absolute change of size
  eps11 = displacements[0][0] - displacements[0][1];
  eps12 = displacements[0][3] - displacements[0][2];
  eps1 = (eps11+eps12) * 0.5;

  eps21 = displacements[1][0] - displacements[1][3];
  eps22 = displacements[1][1] - displacements[1][2];
  eps2 = (eps21+eps22) * 0.5;
  
  eps41 = displacements[0][2] - displacements[0][1];
  eps42 = displacements[0][3] - displacements[0][0]; 
  eps4 = (eps41+eps42)*0.5;
    
  eps51 = displacements[1][1] - displacements[1][0];
  eps52 = displacements[1][3] - displacements[1][2];
  eps5= (eps51+eps52)*0.5;  

  factor =  0.2 * ((eps1*eps1/(length1*length1))
		 + (eps2*eps2/(length2*length2)) 
		 + (eps5*eps5/(length1*length2))
		 + (eps4*eps4/(length1*length1)));

  return factor;
}

} // end of namespace


